use ndarray::ArrayView1;

pub trait ArrayConv<T>
{
    fn conv(array: ArrayView1<'_, T>) -> Self;
}

#[derive(Clone, Copy)]
pub struct Prediction
{
    values: [f32; 4] // Ordering is intergenic, utr, coding, intron
}

impl Prediction
{
    pub fn get(&self) -> &[f32;4] { &self.values }

    pub fn get_intergenic(&self) -> f32 { self.values[0] }

    pub fn get_utr(&self) -> f32 { self.values[1] }

    pub fn get_coding(&self) -> f32 { self.values[2] }

    pub fn get_intron(&self) -> f32 { self.values[3] }

    pub fn get_genic(&self) -> f32 { 1.0 - self.values[0] }
}

impl ArrayConv<f32> for Prediction
{
    fn conv(array: ArrayView1<'_, f32>) -> Self {
        let values: [f32;4] = [ array[0], array[1], array[2], array[3] ];
        Prediction { values }
    }
}

#[derive(Clone, Copy)]
pub struct Bases
{
    values: [f32; 4] // Ordering is C, A, T, G
}

impl Bases
{
    pub fn get(&self) -> &[f32; 4] { &self.values }
}

impl ArrayConv<f32> for Bases
{
    fn conv(array: ArrayView1<'_, f32>) -> Self {
        let values: [f32;4] = [ array[0], array[1], array[2], array[3] ];
        Bases { values }
    }
}

pub struct Transitions
{
    values: [i8; 6] // Ordering is +TR, +CDS, +In, -TR, -CDS, -In
                    // [transcription start site, start codon, donor splice site, transcription stop site, stop codon, acceptor splice site]
}

impl Transitions
{
    pub fn get(&self) -> &[i8; 6] { &self.values }
}

impl ArrayConv<i8> for Transitions
{
    fn conv(array: ArrayView1<'_, i8>) -> Self {
        let values: [i8;6] = [ array[0], array[1], array[2], array[3], array[4], array[5] ];
        Transitions { values }
    }
}



pub struct Annotation
{
    values: [i8; 4] // Ordering is intergenic, utr, coding, intron
}

impl Annotation
{
    pub fn get(&self) -> &[i8; 4] { &self.values }
}

impl ArrayConv<i8> for Annotation
{
    fn conv(array: ArrayView1<'_, i8>) -> Self {
        let values: [i8;4] = [ array[0], array[1], array[2], array[3] ];
        Annotation { values }
    }
}
