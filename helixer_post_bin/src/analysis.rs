use crate::results::{SequenceID, SpeciesID, HelixerResults, Result};
use std::ops::Range;
use crate::results::conv::{Prediction, Bases};
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
    pred_blocked_dataset: BlockedDataset2D<'a, f32, Prediction>,
}

impl<'a> BasePredictionExtractor<'a>
{
    pub fn new(helixer_res: &'a HelixerResults) -> Result<BasePredictionExtractor<'a>>
    {
        let bases_blocked_dataset = helixer_res.get_x()?;
        let pred_blocked_dataset = helixer_res.get_predictions()?;

        Ok( BasePredictionExtractor { helixer_res, bases_blocked_dataset, pred_blocked_dataset } )
    }

    pub fn fwd_iterator(&'a self, sequence_id: SequenceID) -> BasePredictionIterator<'a>
    {
        let sequence = self.helixer_res.get_index().get_sequence_by_id(sequence_id);

        let base_iter = self.bases_blocked_dataset.fwd_iter(sequence_id);
        let pred_iter = self.pred_blocked_dataset.fwd_iter(sequence_id);

        BasePredictionIterator::new(self, sequence.get_species_id(), sequence_id, false, base_iter, pred_iter)
    }

    pub fn rev_iterator(&'a self, sequence_id: SequenceID) -> BasePredictionIterator<'a>
    {
        let sequence = self.helixer_res.get_index().get_sequence_by_id(sequence_id);

        let base_iter = self.bases_blocked_dataset.rev_iter(sequence_id);
        let pred_iter = self.pred_blocked_dataset.rev_iter(sequence_id);

        BasePredictionIterator::new(self, sequence.get_species_id(), sequence_id, true, base_iter, pred_iter)
    }
}


pub struct BasePredictionIterator<'a>
{
    extractor: &'a BasePredictionExtractor<'a>,

    species_id: SpeciesID,
    sequence_id: SequenceID,
    rc: bool,

    base_iter: BlockedDataset2DIter<'a, f32, Bases>,
    pred_iter: BlockedDataset2DIter<'a, f32, Prediction>,
}

impl<'a> BasePredictionIterator<'a>
{
    fn new(extractor: &'a BasePredictionExtractor<'a>, species_id: SpeciesID, sequence_id: SequenceID, rc: bool,
           base_iter: BlockedDataset2DIter<'a, f32, Bases>, pred_iter: BlockedDataset2DIter<'a, f32, Prediction>) -> BasePredictionIterator<'a>
    {
        BasePredictionIterator { extractor, species_id, sequence_id, rc, base_iter, pred_iter }
    }

    pub fn get_extractor(&self) -> &BasePredictionExtractor { self.extractor }

    pub fn get_species_id(&self) -> SpeciesID { self.species_id }

    pub fn get_sequence_id(&self) -> SequenceID { self.sequence_id }

    pub fn get_rc(&self) -> bool { self.rc }
}

impl<'a> Iterator for BasePredictionIterator<'a>
{
    type Item = (Bases, Prediction);

    fn next(&mut self) -> Option<Self::Item>
        {
            let bases = self.base_iter.next();
            let pred = self.pred_iter.next();

            if bases.is_none() && pred.is_none()
                { return None }

            if bases.is_none() || pred.is_none()
                { panic!("Different lengths in base/pred iterators") }

            let bases = bases.expect("Unexpected end of base iter");
            let pred = pred.expect("Unexpected end of base iter");

            return Some((bases, pred))
        }
}




