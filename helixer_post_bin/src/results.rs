use std::path::Path;

pub mod conv;
pub mod error;
pub mod index;
pub mod iter;
pub mod raw;

pub use crate::results::error::Error;
pub type Result<T> = std::result::Result<T, Error>;

use self::index::HelixerIndex;
use self::raw::{RawHelixerGenome, RawHelixerPredictions};
use crate::results::conv::{
    Bases, ClassPrediction, ClassReference, PhasePrediction, PhaseReference, Transitions,
};
use crate::results::iter::{BlockedDataset1D, BlockedDataset2D};

pub struct HelixerResults {
    predictions: RawHelixerPredictions,
    genome: RawHelixerGenome,

    index: HelixerIndex,
}

impl HelixerResults {
    pub fn new(predictions_path: &Path, genome_path: &Path) -> Result<HelixerResults> {
        let predictions = RawHelixerPredictions::new(predictions_path)?;
        let (blocks, blocksize) = predictions.get_blocks_and_blocksize()?;
        let genome = RawHelixerGenome::new(genome_path, blocks, blocksize)?;

        let index = HelixerIndex::new(&genome)?;

        Ok(HelixerResults {
            predictions,
            genome,
            index,
        })
    }

    pub fn get_raw_predictions(&self) -> &RawHelixerPredictions {
        &self.predictions
    }

    pub fn get_raw_genome(&self) -> &RawHelixerGenome {
        &self.genome
    }

    pub fn get_index(&self) -> &HelixerIndex {
        &self.index
    }

    // Delegate species/sequence/block lookups to Index

    pub fn get_all_species(&self) -> &[Species] {
        self.index.get_all_species()
    }

    pub fn get_species_by_id(&self, id: SpeciesID) -> &Species {
        self.index.get_species_by_id(id)
    }

    pub fn get_species_by_name(&self, name: &str) -> Option<&Species> {
        self.index.get_species_by_name(name)
    }

    pub fn get_all_sequences(&self) -> &[Sequence] {
        self.index.get_all_sequences()
    }

    pub fn get_sequence_by_id(&self, id: SequenceID) -> &Sequence {
        self.index.get_sequence_by_id(id)
    }

    pub fn get_sequence_by_species_id_and_sequence_name(
        &self,
        species_id: SpeciesID,
        sequence_name: &str,
    ) -> Option<&Sequence> {
        self.index
            .get_sequence_by_species_id_and_sequence_name(species_id, sequence_name)
    }

    pub fn get_sequences_for_all_species(&self) -> &[Vec<SequenceID>] {
        self.index.get_sequences_for_all_species()
    }

    pub fn get_sequences_for_species(&self, id: SpeciesID) -> &[SequenceID] {
        self.index.get_sequences_for_species(id)
    }

    pub fn get_all_block_offsets(&self) -> &[(u64, u64)] {
        self.index.get_all_block_offsets()
    }

    pub fn get_all_block_ids(&self) -> &[(Vec<BlockID>, Vec<BlockID>)] {
        self.index.get_all_block_ids()
    }

    pub fn get_block_ids_for_sequence(&self, id: SequenceID) -> &(Vec<BlockID>, Vec<BlockID>) {
        self.index.get_block_ids_for_sequence(id)
    }

    // Wrapped dataset accessors for large datasets, delegate smaller datasets to standard collection converters

    pub fn get_class_predictions(&self) -> Result<BlockedDataset2D<f32, ClassPrediction>> {
        Ok(BlockedDataset2D::new(
            &self.index,
            self.predictions.get_class_raw()?,
        ))
    }

    pub fn get_phase_predictions(&self) -> Result<BlockedDataset2D<f32, PhasePrediction>> {
        Ok(BlockedDataset2D::new(
            &self.index,
            self.predictions.get_phase_raw()?,
        ))
    }

    pub fn get_x(&self) -> Result<BlockedDataset2D<f32, Bases>> {
        Ok(BlockedDataset2D::new(&self.index, self.genome.get_x_raw()?))
    }

    /*
        /data/err_samples        	Dataset {268164/Inf}			Sida: Mostly TRUE, occasional FALSE
        /data/fully_intergenic_samples Dataset {268164/Inf}			Sida: Mostly FALSE, occasional TRUE
        /data/gene_lengths       	Dataset {268164/Inf, 20000}		Longest spanning gene model
        /data/is_annotated       	Dataset {268164/Inf}			Sida: All TRUE
        /data/sample_weights     	Dataset {268164/Inf, 20000}
        /data/seqids             	Dataset {268164/Inf}
        /data/species            	Dataset {268164/Inf}
        /data/start_ends         	Dataset {268164/Inf, 2}
    */

    pub fn get_err_samples(&self) -> Result<Vec<bool>> {
        self.genome.get_err_samples()
    }

    pub fn get_fully_intergenic_samples(&self) -> Result<Vec<bool>> {
        self.genome.get_fully_intergenic_samples()
    }

    pub fn get_gene_lengths(&self) -> Result<BlockedDataset1D<u32>> {
        Ok(BlockedDataset1D::new(
            &self.index,
            self.genome.get_gene_lengths_raw()?,
        ))
    }

    pub fn get_is_annotated(&self) -> Result<Vec<bool>> {
        self.genome.get_is_annotated()
    }

    pub fn get_sample_weights(&self) -> Result<BlockedDataset1D<i8>> {
        Ok(BlockedDataset1D::new(
            &self.index,
            self.genome.get_sample_weights_raw()?,
        ))
    }

    pub fn get_transitions(&self) -> Result<BlockedDataset2D<i8, Transitions>> {
        Ok(BlockedDataset2D::new(
            &self.index,
            self.genome.get_transitions_raw()?,
        ))
    }

    pub fn get_class_reference(&self) -> Result<Option<BlockedDataset2D<i8, ClassReference>>> {
        match self.genome.get_y_raw()? {
            Some(dataset) => Ok(Some(BlockedDataset2D::new(&self.index, dataset))),
            None => Ok(None),
        }
    }

    pub fn get_class_reference_as_pseudo_predictions(
        &self,
    ) -> Result<Option<BlockedDataset2D<i8, ClassPrediction>>> {
        match self.genome.get_y_raw()? {
            Some(dataset) => Ok(Some(BlockedDataset2D::new(&self.index, dataset))),
            None => Ok(None),
        }
    }

    pub fn get_phase_reference(&self) -> Result<Option<BlockedDataset2D<i8, PhaseReference>>> {
        match self.genome.get_phases_raw()? {
            Some(dataset) => Ok(Some(BlockedDataset2D::new(&self.index, dataset))),
            None => Ok(None),
        }
    }

    pub fn get_phase_reference_as_pseudo_predictions(
        &self,
    ) -> Result<Option<BlockedDataset2D<i8, PhasePrediction>>> {
        match self.genome.get_phases_raw()? {
            Some(dataset) => Ok(Some(BlockedDataset2D::new(&self.index, dataset))),
            None => Ok(None),
        }
    }
}

// Newtype for Species ID

#[derive(Copy, Clone)]
pub struct SpeciesID(usize);

impl SpeciesID {
    fn new(id: usize) -> SpeciesID {
        SpeciesID(id)
    }

    pub fn inner(&self) -> usize {
        self.0
    }
}

pub struct Species {
    name: String,
    id: SpeciesID,
}

impl Species {
    fn new(name: String, id: SpeciesID) -> Species {
        Species { name, id }
    }

    pub fn get_name(&self) -> &str {
        &self.name
    }

    pub fn get_id(&self) -> SpeciesID {
        self.id
    }
}

// Newtype for Sequence ID

#[derive(Copy, Clone)]
pub struct SequenceID(usize);

impl SequenceID {
    fn new(id: usize) -> SequenceID {
        SequenceID(id)
    }

    pub fn inner(&self) -> usize {
        self.0
    }
}

pub struct Sequence {
    name: String,
    length: u64,
    id: SequenceID,
    species_id: SpeciesID,
}

impl Sequence {
    fn new(name: String, id: SequenceID, species_id: SpeciesID) -> Sequence {
        Sequence {
            name,
            length: 0,
            id,
            species_id,
        }
    }

    fn set_length(&mut self, length: u64) {
        self.length = length;
    }

    pub fn get_name(&self) -> &str {
        &self.name
    }

    pub fn get_length(&self) -> u64 {
        self.length
    }

    pub fn get_id(&self) -> SequenceID {
        self.id
    }

    pub fn get_species_id(&self) -> SpeciesID {
        self.species_id
    }
}

// Newtype for Block ID

#[derive(Copy, Clone)]
pub struct BlockID(usize);

impl BlockID {
    fn new(id: usize) -> BlockID {
        BlockID(id)
    }

    pub fn inner(&self) -> usize {
        self.0
    }
}
