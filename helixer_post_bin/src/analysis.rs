use crate::results::{Species, Sequence};
use crate::results::conv::{ClassPrediction, PhasePrediction, ArrayConvInto};
use std::io::Write;
use crate::analysis::rater::{SequenceRating, SequenceRater};
use crate::gff::GffWriter;
use crate::analysis::window::BasePredictionWindowThresholdIterator;
use crate::analysis::hmm::{PredictionHmm, HmmStateRegion};
use crate::analysis::gff_conv::hmm_solution_to_gff;
use crate::analysis::extractor::{BasePredictionExtractor, ComparisonExtractor};

pub mod extractor;
pub mod gff_conv;
pub mod hmm;
pub mod rater;
pub mod window;


pub struct Analyzer<'a, TC: ArrayConvInto<ClassPrediction>, TP: ArrayConvInto<PhasePrediction>>
{
    bp_extractor: BasePredictionExtractor<'a, TC, TP>,
    comp_extractor: ComparisonExtractor<'a>,
    window_size: usize,
    edge_threshold: f32,
    peak_threshold: f32,
    min_coding_length: usize,
}

impl<'a, TC: ArrayConvInto<ClassPrediction>, TP: ArrayConvInto<PhasePrediction>> Analyzer<'a, TC, TP>
{
    pub fn new(bp_extractor: BasePredictionExtractor<'a, TC, TP>, comp_extractor: ComparisonExtractor<'a>,
               window_size: usize, edge_threshold: f32, peak_threshold: f32, min_coding_length: usize) -> Analyzer<'a, TC, TP>
    {
        Analyzer { bp_extractor, comp_extractor, window_size, edge_threshold, peak_threshold, min_coding_length }
    }

    fn process_sequence_1d<W:Write>(&self, species: &Species, seq: &Sequence, rev: bool, bp_iter: BasePredictionWindowThresholdIterator<TC, TP>,
                                    gene_idx: &mut usize, rater: &mut SequenceRater, gff_writer: &mut GffWriter<W>) -> (usize, usize)
    {
        let mut window_count = 0;
        let mut window_length_total = 0;

        for (bp_vec, _total_vec, start_pos, _peak) in bp_iter
        {
            window_count += 1;
            window_length_total += bp_vec.len();

            let end_pos = start_pos + bp_vec.len();

            let hmm = PredictionHmm::new(bp_vec);
            let maybe_solution = hmm.solve();

            if let Some(solution) = maybe_solution
            {
                //solution.dump(start_pos);
                let solution_regions = solution.trace_regions();
                let genes = HmmStateRegion::split_genes(solution_regions);

                for (gene_regions, coding_length) in genes.iter()
                { rater.rate_regions(start_pos, &gene_regions, *coding_length < self.min_coding_length); }

                let gff_records = hmm_solution_to_gff(genes, species.get_name(), seq.get_name(), "HelixerPost",
                                                      rev, start_pos, seq.get_length(), self.min_coding_length, gene_idx);
                gff_writer.write_records(&gff_records).expect("Failed to write to GFF");
            } else { panic!("No solution at {} {} - {}", seq.get_name(), start_pos, end_pos); }
        }

        (window_count, window_length_total)
    }

    pub fn process_sequence<W: Write>(&self, species: &Species, seq: &Sequence,
                                      fwd_rating: &mut SequenceRating, rev_rating: &mut SequenceRating, gff_writer: &mut GffWriter<W>) -> (usize, usize)
    {
        let id = seq.get_id();
        println!("  BP_Extractor for Sequence {} - ID {}", seq.get_name(), id.inner());

        let mut gene_idx = 1;

        let fwd_bp_iter =
            BasePredictionWindowThresholdIterator::new(self.bp_extractor.fwd_iterator(id), self.window_size, self.edge_threshold, self.peak_threshold).unwrap();
        let mut fwd_comp_rater = SequenceRater::new(self.comp_extractor.fwd_iterator(id), seq.get_length() as usize);

        let (fwd_window_count, fwd_window_length_total) = self.process_sequence_1d(species, seq, false, fwd_bp_iter, &mut gene_idx, &mut fwd_comp_rater, gff_writer);
        let fwd_seq_rating = fwd_comp_rater.calculate_stats();
        println!("Forward for Sequence {} - ID {}", seq.get_name(), id.inner());
        fwd_seq_rating.dump();

        fwd_rating.accumulate(&fwd_seq_rating);

        let rev_bp_iter =
            BasePredictionWindowThresholdIterator::new(self.bp_extractor.rev_iterator(id), self.window_size, self.edge_threshold, self.peak_threshold).unwrap();
        let mut rev_comp_rater = SequenceRater::new(self.comp_extractor.rev_iterator(id), seq.get_length() as usize);

        let (rev_window_count, rev_window_length_total) = self.process_sequence_1d(species, seq, true, rev_bp_iter, &mut gene_idx, &mut rev_comp_rater, gff_writer);
        let rev_seq_rating = rev_comp_rater.calculate_stats();
        println!("Reverse for Sequence {} - ID {}", seq.get_name(), id.inner());
        rev_seq_rating.dump();

        rev_rating.accumulate(&rev_seq_rating);

        (fwd_window_count+rev_window_count, fwd_window_length_total+rev_window_length_total)
      }
}
