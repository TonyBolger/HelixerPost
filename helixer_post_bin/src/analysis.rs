use crate::results::{SequenceID, SpeciesID, HelixerResults, Result};
use std::ops::Range;
use crate::results::conv::{ClassPrediction, Bases, PhasePrediction};
use crate::results::iter::{BlockedDataset2D, BlockedDataset2DIter};

pub mod gff_conv;
pub mod hmm;
pub mod window;

pub struct SequenceRegion
{
    species_id: SpeciesID,
    sequence_id: SequenceID,

    rc: bool,
    position: Range<usize>,

}

impl SequenceRegion
{
    pub fn new(species_id: SpeciesID, sequence_id: SequenceID, rc: bool, position: Range<usize>) -> SequenceRegion
    {
        SequenceRegion { species_id, sequence_id, position, rc }
    }

    pub fn get_species(&self) -> SpeciesID { self.species_id }

    pub fn get_sequence(&self) -> SequenceID { self.sequence_id }

    pub fn get_rc(&self) -> bool { self.rc }

    pub fn get_range(&self) -> &Range<usize> { &self.position }
}

/*

    Functionality Split:

    1) Dataset iterator / accessor (here)
    2) Prediction Window (window module)
    3) Window iteration (window module)

 */


// Hard coded as Bases / Predictions for now. Would make more sense to have all base-level datasets optionally available (and seekable)
pub struct BasePredictionExtractor<'a>
{
    helixer_res: &'a HelixerResults,

    bases_blocked_dataset: BlockedDataset2D<'a, f32, Bases>,

    class_pred_blocked_dataset: BlockedDataset2D<'a, f32, ClassPrediction>,
    phase_pred_blocked_dataset: BlockedDataset2D<'a, f32, PhasePrediction>,
}

impl<'a> BasePredictionExtractor<'a>
{
    pub fn new(helixer_res: &'a HelixerResults) -> Result<BasePredictionExtractor<'a>>
    {
        let bases_blocked_dataset = helixer_res.get_x()?;
        let class_pred_blocked_dataset = helixer_res.get_class_predictions()?;
        let phase_pred_blocked_dataset = helixer_res.get_phase_predictions()?;

        Ok( BasePredictionExtractor { helixer_res, bases_blocked_dataset, class_pred_blocked_dataset, phase_pred_blocked_dataset } )
    }

    pub fn fwd_iterator(&'a self, sequence_id: SequenceID) -> BasePredictionIterator<'a>
    {
        let sequence = self.helixer_res.get_index().get_sequence_by_id(sequence_id);

        let base_iter = self.bases_blocked_dataset.fwd_iter(sequence_id);
        let class_pred_iter = self.class_pred_blocked_dataset.fwd_iter(sequence_id);
        let phase_pred_iter = self.phase_pred_blocked_dataset.fwd_iter(sequence_id);

        BasePredictionIterator::new(self, sequence.get_species_id(), sequence_id, false, base_iter, class_pred_iter, phase_pred_iter)
    }

    pub fn rev_iterator(&'a self, sequence_id: SequenceID) -> BasePredictionIterator<'a>
    {
        let sequence = self.helixer_res.get_index().get_sequence_by_id(sequence_id);

        let base_iter = self.bases_blocked_dataset.rev_iter(sequence_id);
        let class_pred_iter = self.class_pred_blocked_dataset.rev_iter(sequence_id);
        let phase_pred_iter = self.phase_pred_blocked_dataset.rev_iter(sequence_id);

        BasePredictionIterator::new(self, sequence.get_species_id(), sequence_id, true, base_iter, class_pred_iter, phase_pred_iter)
    }
}


pub struct BasePredictionIterator<'a>
{
    extractor: &'a BasePredictionExtractor<'a>,

    species_id: SpeciesID,
    sequence_id: SequenceID,
    rc: bool,

    base_iter: BlockedDataset2DIter<'a, f32, Bases>,

    class_pred_iter: BlockedDataset2DIter<'a, f32, ClassPrediction>,
    phase_pred_iter: BlockedDataset2DIter<'a, f32, PhasePrediction>,
}

impl<'a> BasePredictionIterator<'a>
{
    fn new(extractor: &'a BasePredictionExtractor<'a>, species_id: SpeciesID, sequence_id: SequenceID, rc: bool, base_iter: BlockedDataset2DIter<'a, f32, Bases>,
           class_pred_iter: BlockedDataset2DIter<'a, f32, ClassPrediction>, phase_pred_iter: BlockedDataset2DIter<'a, f32, PhasePrediction>) -> BasePredictionIterator<'a>
    {
        BasePredictionIterator { extractor, species_id, sequence_id, rc, base_iter, class_pred_iter, phase_pred_iter }
    }

    pub fn get_extractor(&self) -> &BasePredictionExtractor { self.extractor }

    pub fn get_species_id(&self) -> SpeciesID { self.species_id }

    pub fn get_sequence_id(&self) -> SequenceID { self.sequence_id }

    pub fn get_rc(&self) -> bool { self.rc }
}

impl<'a> Iterator for BasePredictionIterator<'a>
{
    type Item = (Bases, ClassPrediction, PhasePrediction);

    fn next(&mut self) -> Option<Self::Item>
        {
            let bases = self.base_iter.next();
            let class_pred = self.class_pred_iter.next();
            let phase_pred = self.phase_pred_iter.next();

            if bases.is_none() && class_pred.is_none() && phase_pred.is_none()
                { return None }

            if bases.is_none() || class_pred.is_none() || phase_pred.is_none()
                { panic!("Different lengths in base/class_pred/phase_pred iterators") }

            let bases = bases.expect("Unexpected end of base iter");
            let class_pred = class_pred.expect("Unexpected end of class_pred iter");
            let phase_pred = phase_pred.expect("Unexpected end of phase_pred iter");

            return Some((bases, class_pred, phase_pred))
        }
}




