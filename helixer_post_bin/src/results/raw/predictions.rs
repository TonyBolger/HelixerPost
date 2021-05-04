
use super::super::{Error, Result};
use hdf5::{Dataset, File};
use std::path::Path;

pub struct RawHelixerPredictions
{
    predictions_file: File
}


const PREDICTIONS_DATASIZE: usize = 4;

impl RawHelixerPredictions
{
    pub fn new(predictions_file_path: &Path) -> Result<RawHelixerPredictions>
    {
        let predictions_file = File::open(predictions_file_path)?;
        Ok(RawHelixerPredictions { predictions_file })
    }

    pub fn get_predictions_raw(&self) -> Result<Dataset>
    {
        let pred_dataset = self.predictions_file.dataset("predictions")?;

        let shape = pred_dataset.shape();
        if shape.len()!=3
            { return Err(Error::MismatchedDimensions(shape.len(), 3)); }

        if shape[2]!= PREDICTIONS_DATASIZE
            { return Err(Error::MismatchedDataSize(shape[2], PREDICTIONS_DATASIZE)); }

        Ok(pred_dataset)
    }

    pub fn get_blocks_and_blocksize(&self) -> Result<(usize, usize)>
    {
        let pred_dataset = self.get_predictions_raw()?;
        let shape = pred_dataset.shape();

        let blocks = shape[0];
        let blocksize = shape[1];

        Ok((blocks, blocksize))
    }

}