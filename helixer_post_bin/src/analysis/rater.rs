

/*
    Precision = TP / ( TP + FP )
    Recall = TP / ( TP + FN )

    F1 = (2 * Precision * Recall) / (Precision + Recall)
    F1 = TP / ( TP + ( FP + FN ) / 2 )
    F1 = (2 * TP) / (2 * TP + FP + FN)
 */



use crate::analysis::ComparisonIterator;
use crate::analysis::hmm::{HmmStateRegion, HmmAnnotationLabel};
use crate::results::conv::{ClassReference, PhaseReference, ClassPrediction, PhasePrediction};

pub struct ConfusionMatrix<const N: usize>
{
    count: [[u64; N]; N]
}


impl<const N: usize> ConfusionMatrix<N>
{
    pub fn new() -> ConfusionMatrix<N>
    {
        let count = [[0; N]; N];

        ConfusionMatrix { count }
    }

    pub fn size(&self) -> usize { N }

    pub fn accumulate(&mut self, other: & Self)
    {
        for p in 0..N
            {
            for r in 0..N
                { self.count[r][p]+=other.count[r][p]; }
            }
    }

    pub fn increment(&mut self, ref_idx: usize, pred_idx: usize)
    {
        self.count[ref_idx][pred_idx]+=1;
    }


    pub fn get_tp(&self, idx: usize) -> u64
    {
        self.count[idx][idx]
    }

    pub fn get_fp(&self, idx: usize) -> u64
    {
        let mut false_pos = 0;
        for r in 0..idx
            { false_pos +=self.count[r][idx]; }

        for r in idx+1..N
            { false_pos +=self.count[r][idx]; }

        false_pos
    }

    pub fn get_fn(&self, idx: usize) -> u64
    {
        let mut false_neg = 0;
        for p in 0..idx
            { false_neg += self.count[idx][p]; }

        for p in idx+1 .. N
            { false_neg += self.count[idx][p]; }

        false_neg
    }

    pub fn dump(&self)
    {
        for i in 0..N
        {
            for j in 0..N
                { print!("{}\t\t", self.count[i][j]); }
            println!();
        }
    }
}




#[derive(Clone, Copy)]
enum Annotation
{
    OutsideWindow,
    Intergenic,
    UTR,
    CodingPhase0,
    CodingPhase1,
    CodingPhase2,
    Intron
}

impl Annotation
{
    pub fn to_str(&self) -> &str
    {
        match self
        {
            Annotation::OutsideWindow => "OutsideWindow",
            Annotation::Intergenic => "Intergenic",
            Annotation::UTR => "UTR",
            Annotation::CodingPhase0 => "CodingPhase0",
            Annotation::CodingPhase1 => "CodingPhase1",
            Annotation::CodingPhase2 => "CodingPhase2",
            Annotation::Intron => "Intron",
        }
    }

    fn get_window_idx(self) -> usize
    {
        match self
        {
            Annotation::OutsideWindow => 0,
            _ => 1
        }
    }

    fn get_class_idx(self) -> usize
    {
        match self
        {
            Annotation::OutsideWindow => 0,
            Annotation::Intergenic => 0,
            Annotation::UTR => 1,
            Annotation::CodingPhase0 => 2,
            Annotation::CodingPhase1 => 2,
            Annotation::CodingPhase2 => 2,
            Annotation::Intron => 3,
        }
    }

    fn get_phase_idx(self) -> usize
    {
        match self
        {
            Annotation::OutsideWindow => 0,
            Annotation::Intergenic => 0,
            Annotation::UTR => 0,
            Annotation::CodingPhase0 => 1,
            Annotation::CodingPhase1 => 2,
            Annotation::CodingPhase2 => 3,
            Annotation::Intron => 0,
        }
    }
}

impl ClassReference
{
    fn get_window_idx(&self) -> usize
    {
        match self.get_max_idx()
            {
            0 => 0,
            _ => 1
            }
    }

    fn get_class_idx(&self) -> usize { self.get_max_idx() }
}

impl PhaseReference
{
    fn get_phase_idx(&self) -> usize { self.get_max_idx() }
}


impl ClassPrediction
{
    fn get_class_idx(&self) -> usize { self.get_max_idx() }
}

impl PhasePrediction
{
    fn get_phase_idx(&self) -> usize { self.get_max_idx() }
}





pub struct SequenceRater<'a>
{
    comp_iterator: ComparisonIterator<'a>,
    annotation: Vec<Annotation>
}

impl<'a> SequenceRater<'a>
{
    pub fn new(comp_iterator: ComparisonIterator<'a>, seq_length: usize) -> SequenceRater<'a>
    {
        let mut annotation = Vec::with_capacity(seq_length);
        annotation.resize(seq_length, Annotation::OutsideWindow);

        SequenceRater { comp_iterator, annotation }
    }

    pub fn rate_regions(&mut self, start_offset: usize, gene_regions: &[HmmStateRegion])
    {
        let mut coding_annotation = Annotation::CodingPhase0;

        for region in gene_regions
            {
            if region.get_annotation_label()!=HmmAnnotationLabel::Coding
                {
                let annotation = match region.get_annotation_label()
                    {
                    HmmAnnotationLabel::Intergenic => Annotation::Intergenic,
                    HmmAnnotationLabel::UTR5 => Annotation::UTR,
                    HmmAnnotationLabel::Intron => Annotation::Intron,
                    HmmAnnotationLabel::UTR3 => Annotation::UTR,
                    _ => panic!("Unexpected Hmm Annotation Label {}", region.get_annotation_label().to_str())
                   };

                for pos in region.get_start_pos()+start_offset..region.get_end_pos()+start_offset
                    { self.annotation[pos]=annotation; }
                }
            else
                {
                for pos in region.get_start_pos()+start_offset..region.get_end_pos()+start_offset
                    {
                    self.annotation[pos]=coding_annotation;

                    coding_annotation = match coding_annotation
                        {
                        Annotation::CodingPhase0 => Annotation::CodingPhase2,
                        Annotation::CodingPhase1 => Annotation::CodingPhase0,
                        Annotation::CodingPhase2 => Annotation::CodingPhase1,
                        _ => panic!("Unexpected Coding Annotation Label {}", coding_annotation.to_str())
                        }
                    }

                }

            }
    }

    pub fn calculate_stats(self) -> SequenceRating
    {
        let mut rating = SequenceRating::new();
        rating.rate(self.comp_iterator, self.annotation);

        rating
    }
}




// Window States = Outside / Inside * Intergenic / Genic
// Class = Intergenic / UTR / Coding / Intron
// Phase = NonCoding / Ph0 / Ph1 / Ph2

pub struct SequenceRating
{
    ref_ml_class_confusion: ConfusionMatrix<4>,
    ref_ml_phase_confusion: ConfusionMatrix<4>,

    ref_hp_class_confusion: ConfusionMatrix<4>,
    ref_hp_phase_confusion: ConfusionMatrix<4>,

    ml_hp_class_confusion: ConfusionMatrix<4>,
    ml_hp_phase_confusion: ConfusionMatrix<4>,

    ref_window_confusion: ConfusionMatrix<2>
}

impl SequenceRating
{
    fn new() -> SequenceRating
    {
        SequenceRating { ref_ml_class_confusion: ConfusionMatrix::<4>::new(), ref_ml_phase_confusion: ConfusionMatrix::<4>::new(),
            ref_hp_class_confusion: ConfusionMatrix::<4>::new(), ref_hp_phase_confusion: ConfusionMatrix::<4>::new(),
            ml_hp_class_confusion: ConfusionMatrix::<4>::new(), ml_hp_phase_confusion: ConfusionMatrix::<4>::new(),
            ref_window_confusion: ConfusionMatrix::<2>::new()}
    }

    fn rate(&mut self, comp_iterator: ComparisonIterator, annotation: Vec<Annotation>)
    {
        for ((class_ref, phase_ref, class_ml, phase_ml), annotation) in
            comp_iterator.zip(annotation.into_iter())
        {
            let ref_class_idx = class_ref.get_class_idx();
            let ref_phase_idx = phase_ref.get_phase_idx();

            let ml_class_idx = class_ml.get_class_idx();
            let ml_phase_idx = phase_ml.get_phase_idx();

            let hp_class_idx = annotation.get_class_idx();
            let hp_phase_idx = annotation.get_phase_idx();

            let ref_window_idx = class_ref.get_window_idx();
            let hp_window_idx = annotation.get_window_idx();

            self.ref_ml_class_confusion.increment(ref_class_idx, ml_class_idx);
            self.ref_ml_phase_confusion.increment(ref_phase_idx, ml_phase_idx);

            self.ref_hp_class_confusion.increment(ref_class_idx, hp_class_idx);
            self.ref_hp_phase_confusion.increment(ref_phase_idx, hp_phase_idx);

            self.ml_hp_class_confusion.increment(ml_class_idx, hp_class_idx);
            self.ml_hp_phase_confusion.increment(ml_phase_idx, hp_phase_idx);

            self.ref_window_confusion.increment(ref_window_idx, hp_window_idx);
        }
    }

    pub fn dump(&self)
    {
        println!("Ref vs HP Window");
        self.ref_window_confusion.dump();
        println!();

        println!("Ref vs ML class");
        self.ref_ml_class_confusion.dump();
        println!("Ref vs HP class");
        self.ref_hp_class_confusion.dump();
        println!("ML vs HP class");
        self.ml_hp_class_confusion.dump();
        println!();

        println!("Ref vs ML phase");
        self.ref_ml_phase_confusion.dump();
        println!("Ref vs HP Phase");
        self.ref_hp_phase_confusion.dump();
        println!("ML vs HP phase");
        self.ml_hp_phase_confusion.dump();
        println!();
    }
}


















#[cfg(test)]
mod tests {
    use crate::analysis::rater::ConfusionMatrix;

    #[test]
    fn test_confusion_empty_matrix()
    {
        let empty_matrix = ConfusionMatrix::<3>::new();
        assert_eq!(empty_matrix.size(), 3);

        for i in 0..3
            {
            assert_eq!(empty_matrix.get_tp(i), 0);
            assert_eq!(empty_matrix.get_fp(i), 0);
            assert_eq!(empty_matrix.get_fn(i), 0);
            }
    }

    #[test]
    fn test_confusion_perfect_matrix()
    {
        let mut matrix = ConfusionMatrix::<3>::new();

        for i in 0..3
            { matrix.increment(i,i); }

        for i in 0..3
        {
            assert_eq!(matrix.get_tp(i), 1);
            assert_eq!(matrix.get_fp(i), 0);
            assert_eq!(matrix.get_fn(i), 0);
        }
    }

    #[test]
    fn test_confusion_fp_matrix()
    {
        for i in 0..3
            {
            let mut matrix = ConfusionMatrix::<3>::new();

            for j in 0..3
                { matrix.increment(j, i); } // Same pred for all refs

            assert_eq!(matrix.get_tp(i), 1);
            assert_eq!(matrix.get_fp(i), 2);
            assert_eq!(matrix.get_fn(i), 0);
            }
    }


    #[test]
    fn test_confusion_fn_matrix()
    {
        for i in 0..3
        {
            let mut matrix = ConfusionMatrix::<3>::new();

            for j in 0..3
                { matrix.increment(i, j); } // Same pred for all refs

            assert_eq!(matrix.get_tp(i), 1);
            assert_eq!(matrix.get_fp(i), 0);
            assert_eq!(matrix.get_fn(i), 2);
        }
    }
}





