use std::collections::{BTreeMap, HashMap};
use std::fmt::Write;

use super::{Error, Result};

use super::raw::RawHelixerGenome;
use super::{BlockID, Sequence, SequenceID, Species, SpeciesID};

pub struct HelixerIndex {
    species: Vec<Species>,
    species_name_idx: HashMap<String, SpeciesID>,

    sequences: Vec<Sequence>,
    sequence_name_idx: Vec<HashMap<String, SequenceID>>, // Species ID -> Map
    species_sequences: Vec<Vec<SequenceID>>,             // Species ID -> [Sequence_ID]

    block_offsets: Vec<(u64, u64)>, // Block Idx -> Position range
    sequence_blocks: Vec<(Vec<BlockID>, Vec<BlockID>)>, // SequenceID -> ([BlockID], [BlockID])
}

impl HelixerIndex {
    pub fn new(genome: &RawHelixerGenome) -> Result<HelixerIndex> {
        let all_species = genome.get_species()?;
        let all_sequences = genome.get_seqids()?;
        let all_startends = genome.get_start_ends()?;

        Self::build_from_slices(&all_species, &all_sequences, &all_startends)
    }

    fn build_from_slices(
        all_species: &[String],
        all_sequences: &[String],
        all_startends: &[(u64, u64)],
    ) -> Result<HelixerIndex> {
        let mut species = Vec::new();
        let mut species_name_idx = HashMap::new();

        let mut sequences = Vec::new();
        let mut sequence_name_idx = Vec::new();
        let mut species_sequences = Vec::new();

        let mut block_offsets = Vec::with_capacity(all_startends.len());
        let mut sequence_blocks_trees = Vec::new(); // Build initially as BTreeMaps

        let mut maybe_last_species_name = None;
        let mut maybe_last_sequence_name = None;

        let mut species_id = SpeciesID::new(usize::max_value());
        let mut sequence_id = SequenceID::new(usize::max_value());

        for ((species_name, sequence_name), (start, end)) in all_species
            .iter()
            .zip(all_sequences.iter())
            .zip(all_startends.iter())
        {
            if maybe_last_species_name.is_none() || maybe_last_species_name.unwrap() != species_name
            // New species
            {
                species_id = SpeciesID::new(species.len());
                species.push(Species::new(species_name.clone(), species_id));
                species_name_idx.insert(species_name.clone(), species_id);

                sequence_name_idx.push(HashMap::new());
                species_sequences.push(Vec::new());

                maybe_last_species_name = Some(species_name);
                maybe_last_sequence_name = None;
            }

            if maybe_last_sequence_name.is_none()
                || maybe_last_sequence_name.unwrap() != sequence_name
            // New sequence (or species)
            {
                sequence_id = SequenceID::new(sequences.len());
                sequences.push(Sequence::new(
                    sequence_name.clone(),
                    sequence_id,
                    species_id,
                ));

                sequence_name_idx[species_id.0].insert(sequence_name.clone(), sequence_id);

                species_sequences[species_id.0].push(sequence_id);
                sequence_blocks_trees.push((BTreeMap::new(), BTreeMap::new()));

                maybe_last_sequence_name = Some(sequence_name);
            }

            let block_id = BlockID::new(block_offsets.len());
            block_offsets.push((*start, *end));

            if start < end {
                if let Some(prev_id) = sequence_blocks_trees[sequence_id.0]
                    .0
                    .insert(*start, block_id)
                {
                    let mut msg = String::new();
                    write!(
                        &mut msg,
                        "Block Start {} at index {} already occurred at index {}",
                        start, block_id.0, prev_id.0
                    )
                    .unwrap();
                    return Err(Error::DuplicateValue(msg));
                }
            } else if start > end {
                if let Some(prev_id) = sequence_blocks_trees[sequence_id.0]
                    .1
                    .insert(*end, block_id)
                {
                    let mut msg = String::new();
                    write!(
                        &mut msg,
                        "Block End {} at index {} already occurred at index {}",
                        end, block_id.0, prev_id.0
                    )
                    .unwrap();
                    return Err(Error::DuplicateValue(msg));
                }
            } else {
                return Err("Zero length block".into());
            }
        }

        let (sequence_blocks, sequence_lengths) =
            Self::build_sequence_blocks(sequence_blocks_trees, &block_offsets)?;

        for (idx, sequence_length) in sequence_lengths.into_iter().enumerate() {
            sequences[idx].set_length(sequence_length);
        }

        Ok(HelixerIndex {
            species,
            species_name_idx,
            sequences,
            sequence_name_idx,
            species_sequences,
            block_offsets,
            sequence_blocks,
        })
    }

    fn build_sequence_blocks_fwd(
        sequence_block_tree: BTreeMap<u64, BlockID>,
    ) -> Result<Vec<BlockID>> {
        Ok(sequence_block_tree.into_iter().map(|(_k, v)| v).collect())
    }

    fn build_sequence_blocks_rev(
        sequence_block_tree: BTreeMap<u64, BlockID>,
    ) -> Result<Vec<BlockID>> {
        Ok(sequence_block_tree
            .into_iter()
            .rev()
            .map(|(_k, v)| v)
            .collect())
    }

    fn check_sequence_block_contiguity(
        block_indexes: &[BlockID],
        block_offsets: &[(u64, u64)],
    ) -> Result<(u64, u64)> {
        let mut maybe_prev = None;
        for (id, curr) in block_indexes.iter().map(|id| (id, block_offsets[id.0])) {
            if let Some(prev) = maybe_prev {
                if curr.0 != prev {
                    let mut msg = String::new();
                    write!(
                        &mut msg,
                        "Gap between blocks {} to {} at index {}",
                        prev, curr.0, id.0
                    )
                    .unwrap();
                    return Err(Error::InvalidValue(msg));
                }
                maybe_prev = Some(curr.1);
            }
        }

        if block_indexes.len() == 0 {
            Ok((0, 0))
        } else {
            Ok((
                block_offsets[block_indexes.first().unwrap().0].0,
                block_offsets[block_indexes.last().unwrap().0].1,
            ))
        }
    }

    fn build_sequence_blocks(
        sequence_block_trees: Vec<(BTreeMap<u64, BlockID>, BTreeMap<u64, BlockID>)>,
        block_offsets: &[(u64, u64)],
    ) -> Result<(Vec<(Vec<BlockID>, Vec<BlockID>)>, Vec<u64>)> {
        let sequence_blocks_with_lengths = sequence_block_trees
            .into_iter()
            .map(|(fwd_tree, rev_tree)| {
                let fwd_vec = Self::build_sequence_blocks_fwd(fwd_tree)?;
                let rev_vec = Self::build_sequence_blocks_rev(rev_tree)?;

                if fwd_vec.len() != rev_vec.len() {
                    return Err(
                        "Forward / Reverse block count mismatched (perhaps using filtered data?)"
                            .into(),
                    );
                }

                let (fwd_start, fwd_end) =
                    Self::check_sequence_block_contiguity(&fwd_vec, block_offsets)?;
                let (rev_start, rev_end) =
                    Self::check_sequence_block_contiguity(&rev_vec, block_offsets)?;

                if fwd_start != 0 {
                    return Err(
                        "Forward blocks not starting at offset zero (perhaps using filtered data?)"
                            .into(),
                    );
                }

                if rev_end != 0 {
                    return Err(
                        "Reverse blocks not ending at offset zero (perhaps using filtered data?)"
                            .into(),
                    );
                }

                if fwd_end != rev_start {
                    return Err(
                        "Forward / Reverse max offset mismatched (perhaps using filtered data?)"
                            .into(),
                    );
                }

                Ok((fwd_vec, rev_vec, fwd_end))
            })
            .collect::<Result<Vec<(Vec<BlockID>, Vec<BlockID>, u64)>>>()?;

        let mut sequence_blocks = Vec::with_capacity(sequence_blocks_with_lengths.len());
        let mut sequence_lengths = Vec::with_capacity(sequence_blocks_with_lengths.len());

        for (fwd_blocks, rev_blocks, length) in sequence_blocks_with_lengths.into_iter() {
            sequence_blocks.push((fwd_blocks, rev_blocks));
            sequence_lengths.push(length);
        }

        Ok((sequence_blocks, sequence_lengths))
    }

    pub fn get_all_species(&self) -> &[Species] {
        &self.species
    }

    pub fn get_species_by_id(&self, id: SpeciesID) -> &Species {
        &self.species[id.0]
    }

    pub fn get_species_by_name(&self, name: &str) -> Option<&Species> {
        self.species_name_idx
            .get(name)
            .and_then(|id| self.species.get(id.0))
    }

    pub fn get_all_sequences(&self) -> &[Sequence] {
        &self.sequences
    }

    pub fn get_sequence_by_id(&self, id: SequenceID) -> &Sequence {
        &self.sequences[id.0]
    }

    pub fn get_sequence_by_species_id_and_sequence_name(
        &self,
        species_id: SpeciesID,
        sequence_name: &str,
    ) -> Option<&Sequence> {
        self.sequence_name_idx[species_id.0]
            .get(sequence_name)
            .and_then(|id| self.sequences.get(id.0))
    }

    pub fn get_sequences_for_all_species(&self) -> &[Vec<SequenceID>] {
        &self.species_sequences
    }

    pub fn get_sequences_for_species(&self, id: SpeciesID) -> &[SequenceID] {
        &self.species_sequences[id.0]
    }

    pub fn get_all_block_offsets(&self) -> &[(u64, u64)] {
        &self.block_offsets
    }

    pub fn get_all_block_ids(&self) -> &[(Vec<BlockID>, Vec<BlockID>)] {
        &self.sequence_blocks
    }

    pub fn get_block_ids_for_sequence(&self, id: SequenceID) -> &(Vec<BlockID>, Vec<BlockID>) {
        &self.sequence_blocks[id.0]
    }

    pub fn dump(&self) {
        println!("Species Count: {}", self.species.len());

        for (species_idx, species) in self.species.iter().enumerate() {
            println!("Species [{}] is {}", species_idx, species.get_name());

            let sequence_idxs = &self.species_sequences[species_idx];
            println!("    Sequence Count: {}", sequence_idxs.len());

            for (sequence_idx, sequence) in sequence_idxs
                .iter()
                .map(|idx| (idx, &self.sequences[idx.0]))
            {
                println!(
                    "    Sequence [{}] is {}",
                    sequence_idx.0,
                    sequence.get_name()
                );

                let (fwd_blocks, rev_blocks) = &self.sequence_blocks[sequence_idx.0];
                print!("        Fwd Block Count: {} -", fwd_blocks.len());

                for block in fwd_blocks.iter() {
                    let (block_start, block_end) = self.block_offsets[block.0];
                    print!(" {}({}:{})", block.0, block_start, block_end);
                }
                println!();

                print!("        Rev Block Count: {} -", rev_blocks.len());
                for block in rev_blocks.iter() {
                    let (block_start, block_end) = self.block_offsets[block.0];
                    print!(" {}({}:{})", block.0, block_start, block_end);
                }
                println!();
            }
        }
    }
}
