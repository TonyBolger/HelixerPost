use crate::results::conv::{Bases, ClassPrediction, PhasePrediction};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Clone, Copy)]
struct ClassPredPenalty
{
    penalty: [f64; 4] // Ordering is intergenic, utr, coding, intron
}

#[allow(dead_code)]
impl ClassPredPenalty
{
    pub fn get_intergenic_penalty(&self) -> f64 { self.penalty[0] }

    pub fn get_utr_penalty(&self) -> f64 { self.penalty[1] }

    pub fn get_coding_penalty(&self) -> f64 { self.penalty[2] }

    pub fn get_intron_penalty(&self) -> f64 { self.penalty[3] }
}

const CLASS_PRED_PROB_FLOOR: f64 = 0.000_000_001; // Prevent infinite penalties

impl From<&ClassPrediction> for ClassPredPenalty
{
    fn from(pred: &ClassPrediction) -> Self {

        let raw_pred = pred.get();

        let mut penalty = [0.0 ; 4];
        for i in 0..4
            {
            let adjusted_pred = raw_pred[i] as f64;

            let adjusted_pred = if adjusted_pred > CLASS_PRED_PROB_FLOOR { adjusted_pred } else { CLASS_PRED_PROB_FLOOR };
            penalty[i] = -f64::log2(adjusted_pred);
            }

        let mut min_penalty = penalty[0];
        for i in 1..4
            { if penalty[i] < min_penalty { min_penalty = penalty[i] } }


        for i in 0..4
            { penalty[i]-=min_penalty; }

        ClassPredPenalty { penalty }
    }
}


#[derive(Clone, Copy)]
struct PhasePredPenalty
{
    penalty: [f64; 4] // Ordering is non-coding, phase 0, phase 1, phase 2
}

#[allow(dead_code)]
impl PhasePredPenalty
{
    pub fn get_non_coding_penalty(&self) -> f64 { self.penalty[0] }

    pub fn get_phase0_penalty(&self) -> f64 { self.penalty[1] }

    pub fn get_phase1_penalty(&self) -> f64 { self.penalty[2] }

    pub fn get_phase2_penalty(&self) -> f64 { self.penalty[3] }
}

//const PHASE_PRED_PROB_FLOOR: f64 = 0.000_000_001; // Prevent infinite penalties
const PHASE_PRED_PROB_FLOOR: f64 = 0.5; // Limit impact of incorrect phase

impl From<&PhasePrediction> for PhasePredPenalty
{
    fn from(pred: &PhasePrediction) -> Self {

        let raw_pred = pred.get();

        let mut penalty = [0.0 ; 4];
        for i in 0..4
        {
            let adjusted_pred = raw_pred[i] as f64;

            let adjusted_pred = if adjusted_pred > PHASE_PRED_PROB_FLOOR { adjusted_pred } else { PHASE_PRED_PROB_FLOOR };
            penalty[i] = -f64::log2(adjusted_pred);
        }

        let mut min_penalty = penalty[0];
        for i in 1..4
        { if penalty[i] < min_penalty { min_penalty = penalty[i] } }

        for i in 0..4
        { penalty[i]-=min_penalty; }

        PhasePredPenalty { penalty }
    }
}


#[derive(Clone, Copy)]
struct PredPenalty
{
    penalty: [f64; 6] // Ordering is intergenic, utr, coding_phase0, coding_phase1, coding_phase2, intron
}

#[allow(dead_code)]
impl PredPenalty
{
    pub fn get_intergenic_penalty(&self) -> f64 { self.penalty[0] }

    pub fn get_utr_penalty(&self) -> f64 { self.penalty[1] }

    pub fn get_coding_phase0_penalty(&self) -> f64 { self.penalty[2] }

    pub fn get_coding_phase1_penalty(&self) -> f64 { self.penalty[3] }

    pub fn get_coding_phase2_penalty(&self) -> f64 { self.penalty[4] }

    pub fn get_intron_penalty(&self) -> f64 { self.penalty[5] }
}

const PHASE_RETAIN: f64 = 1.0;      // Adjust as needed

const PHASE_DILUTE: f64 = 1.0 - PHASE_RETAIN;
const PRED_PROB_FLOOR: f64 =  0.000_000_001; // Prevent infinite penalties

impl From<(&ClassPrediction, &PhasePrediction)> for PredPenalty
{
    fn from((class_pred, phase_pred): (&ClassPrediction, &PhasePrediction)) -> Self {

        let phase0 = phase_pred.get_phase0() as f64;
        let phase1 = phase_pred.get_phase1() as f64;
        let phase2 = phase_pred.get_phase2() as f64;

        // Approach 1: rescale total phase to match coding and blend to dilution target (mean coding or total coding)
        let coding = class_pred.get_coding() as f64;
        let total_coding_phase = phase0+phase1+phase2;
        let (phase0, phase1, phase2) = if total_coding_phase > 0.0  // Prevent div by zero risk
            {
            let phase_scale = coding / total_coding_phase;
            (phase0 * phase_scale, phase1 * phase_scale, phase2 * phase_scale)
            }
        else { ( coding / 3.0, coding / 3.0, coding / 3.0 ) };

        //let dilution_target = coding / 3.0;
        let dilution_target = coding;

        let phase0 = phase0 * PHASE_RETAIN + dilution_target * PHASE_DILUTE;
        let phase1 = phase1 * PHASE_RETAIN + dilution_target * PHASE_DILUTE;
        let phase2 = phase2 * PHASE_RETAIN + dilution_target * PHASE_DILUTE;

        /*
        // Approach 2: rescale total phase to 1, dilute towards 1, then scale by coding
        let total_coding_phase = phase0+phase1+phase2;
        let (phase0, phase1, phase2) = if total_coding_phase > 0.0  // Prevent div by zero risk
            {
            let phase_scale = 1.0 / total_coding_phase;
            (phase0 * phase_scale, phase1 * phase_scale, phase2 * phase_scale)
            }
        else { ( 1.0/3.0, 1.0/3.0, 1.0/3.0 ) };

        let coding = class_pred.get_coding() as f64;
        let phase0 = (phase0 * PHASE_RETAIN + 1.0 * PHASE_DILUTE) * coding;
        let phase1 = (phase1 * PHASE_RETAIN + 1.0 * PHASE_DILUTE) * coding;
        let phase2 = (phase2 * PHASE_RETAIN + 1.0 * PHASE_DILUTE) * coding;
*/

        // Convert to combined prob (not necessarily 1-hot)
        let raw_probs = [class_pred.get_intergenic() as f64, class_pred.get_utr() as f64, phase0, phase1, phase2, class_pred.get_intron() as f64];

        // Calculate penalty from probs
        let mut penalty = [0.0 ; 6];
        for i in 0..6
        {
            let adjusted_pred = raw_probs[i];

            let adjusted_pred = if adjusted_pred > PRED_PROB_FLOOR { adjusted_pred } else { PRED_PROB_FLOOR };
            penalty[i] = -f64::log2(adjusted_pred);
        }

        let mut min_penalty = penalty[0];
        for i in 1..6
            { if penalty[i] < min_penalty { min_penalty = penalty[i] } }

        for i in 0..6
            { penalty[i]-=min_penalty; }

        PredPenalty { penalty }
    }

}


const BASE_PROB_FLOOR: f64 = 0.000_000_001; // Prevent infinite penalties

#[derive(Clone, Copy)]
pub struct BasesPenalty
{
    penalty: [f64; 4] // Ordering is C, A, T, G
}

impl BasesPenalty
{
    pub fn get_c(&self) -> f64
    {
        self.penalty[0]
    }

    pub fn get_a(&self) -> f64
    {
        self.penalty[1]
    }

    pub fn get_t(&self) -> f64 { self.penalty[2] }

    pub fn get_g(&self) -> f64
    {
        self.penalty[3]
    }

    pub fn as_str(&self) -> char
    {
        if self.penalty[0]<BASE_PROB_FLOOR
            { 'C' }
        else if self.penalty[1]<BASE_PROB_FLOOR
            { 'A' }
        else if self.penalty[2]<BASE_PROB_FLOOR
            { 'T' }
        else
            { 'G' }
    }
}

impl From<&Bases> for BasesPenalty
{
    fn from(bases: &Bases) -> Self {

        let raw_bases = bases.get();

        let mut penalty = [0.0 ; 4];
        for i in 0..4
        {
            let adjusted_base = raw_bases[i] as f64;

            let adjusted_base = if adjusted_base > BASE_PROB_FLOOR { adjusted_base } else { BASE_PROB_FLOOR };
            penalty[i] = -f64::log2(adjusted_base);
        }

        let mut min_penalty = penalty[0];
        for i in 1..4
            { if penalty[i] < min_penalty { min_penalty = penalty[i] } }


        for i in 0..4
            { penalty[i]-=min_penalty; }

        BasesPenalty { penalty }
    }

}


struct TransitionContext<'a>
{
    class_pred_pen: &'a [ClassPredPenalty],
    phase_pred_pen: &'a [PhasePredPenalty],
    pred_pen: &'a [PredPenalty],

    base_pen: &'a [BasesPenalty],
    offset: usize
}

#[allow(dead_code)]
impl<'a> TransitionContext<'a>
{
    fn new(class_pred_pen: &'a [ClassPredPenalty], phase_pred_pen: &'a [PhasePredPenalty], pred_pen: &'a [PredPenalty], base_pen: &'a [BasesPenalty], offset: usize) -> TransitionContext<'a>
    {
        TransitionContext { class_pred_pen, phase_pred_pen, pred_pen, base_pen, offset}
    }

    fn get_class_pred(&self, position: usize) -> Option<&ClassPredPenalty>
    {
        self.class_pred_pen.get(self.offset + position)
    }

    fn get_phase_pred(&self, position: usize) -> Option<&PhasePredPenalty>
    {
        self.phase_pred_pen.get(self.offset + position)
    }

    fn get_pred(&self, position: usize) -> Option<&PredPenalty>
    {
        self.pred_pen.get(self.offset + position)
    }

    /*
    fn get_upstream(&self, len: usize) -> Option<&[BasesPenalty]>
    {
        if len <= self.offset
            { Some(&self.base_pen[self.offset-len .. self.offset]) }
        else { None }
    }
*/

    fn get_downstream(&self, len: usize) -> Option<&[BasesPenalty]>
    {
        if self.offset + len <= self.base_pen.len()
        { Some(&self.base_pen[self.offset .. self.offset + len]) }
        else
        { None }
    }

    fn get_ctx(&self, ulen: usize, dlen: usize) -> Option<&[BasesPenalty]>
    {
        if ulen <= self.offset && self.offset + dlen <= self.base_pen.len()
        { Some(&self.base_pen[self.offset-ulen .. self.offset + dlen]) }
        else
        { None }
    }

    fn get_donor_penalty_ag_gt(&self, can_splice: bool) -> Option<f64>
    {
        if let (Some(ds), true) = (self.get_ctx(2, 2), can_splice)
        {
            let pen =
//                ds[0].get_a() +
//                ds[1].get_g() +
                ds[2].get_g() +
                    ds[3].get_t();
            Some(pen * DONOR_WEIGHT + DONOR_FIXED_PENALTY)
        }
        else
        { None }
    }

    fn get_acceptor_penalty_ag_g(&self) -> Option<f64>
    {
        if let Some(us) = self.get_ctx(2,1)
        {
            let pen =
                us[0].get_a() +
                    us[1].get_g();
//                us[2].get_g();
            Some(pen*ACCEPTOR_WEIGHT)
        }
        else
        { None }
    }

}

























#[allow(dead_code)]
#[derive(Clone, Copy, Eq, PartialEq)]
pub enum HmmAnnotationLabel
{
    Intergenic,
    UTR5,
    Start,
    Coding,
    Intron,
    Stop,
    UTR3,
}

impl HmmAnnotationLabel
{
    fn to_str(&self) -> &str
    {
        match self
        {
            HmmAnnotationLabel::Intergenic => "Intergenic",
            HmmAnnotationLabel::UTR5 => "UTR5",
            HmmAnnotationLabel::Start => "Start",
            HmmAnnotationLabel::Coding => "Coding",
            HmmAnnotationLabel::Intron => "Intron",
            HmmAnnotationLabel::Stop => "Stop",
            HmmAnnotationLabel::UTR3 => "UTR3",
        }
    }
}

/*
#[allow(dead_code)]
#[derive(Clone, Copy, Eq, PartialEq)]
enum HmmIntron
{
    Intron(usize)
}


impl HmmIntron
{
    fn len(self) -> usize
    {
        match self
        {
            HmmIntron::Intron(l) => l
        }
    }
}
*/

const HMM_STATES: usize = 33;

#[derive(Clone, Copy, Eq, PartialEq, Ord)]
enum HmmState
{
    Intergenic = 0,
    UTR5 = 1,
    UTR5IntronDSS = 2,
    UTR5Intron = 3,
    Start0 = 4, // After A
    Start0IntronDSS = 5,
    Start0Intron = 6,
    Start1 = 7, // After AT
    Start1IntronDSS = 8,
    Start1Intron = 9,
    Start2 = 10, // After ATG
    Coding0 = 11,
    Coding0IntronDSS = 12,
    Coding0Intron = 13,
    Coding1 = 14,
    Coding1IntronDSS = 15,
    Coding1Intron = 16,
    Coding2 = 17,
    Coding2IntronDSS = 18,
    Coding2Intron = 19,
    Stop0T = 20,
    Stop0TIntronDSS = 21,
    Stop0TIntron = 22,
    Stop1TA = 23,
    Stop1TAIntronDSS = 24,
    Stop1TAIntron = 25,
    Stop1TG = 26,
    Stop1TGIntronDSS = 27,
    Stop1TGIntron = 28,
    Stop2 = 29,
    UTR3 = 30,
    UTR3IntronDSS = 31,
    UTR3Intron = 32
}

impl std::cmp::PartialOrd for HmmState
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering>
    {
        let s = *self as u8;
        let o = *other as u8;

        Some(s.cmp(&o))
    }
}

// Allow/Prevent introns within each state and/or at state transitions

const CAN_SPLICE_UTR5: bool = true;
const CAN_SPLICE_UTR5_START: bool = true;
const CAN_SPLICE_START: bool = true;
const CAN_SPLICE_START_CODING: bool = true;
const CAN_SPLICE_CODING: bool = true;
const CAN_SPLICE_CODING_STOP: bool = true;
const CAN_SPLICE_STOP: bool = true;
const CAN_SPLICE_STOP_UTR3: bool = true;
const CAN_SPLICE_UTR3: bool = true;

const START_WEIGHT: f64 = 1_000.0;

const DONOR_FIXED_PENALTY: f64 = 0.0;
const DONOR_WEIGHT: f64 = 1.0;

const ACCEPTOR_WEIGHT: f64 = 1.0;

const STOP_WEIGHT: f64 = 1_000.0;

//const MIN_INTRON_LENGTH: usize = 50;


pub fn show_config()
{
    println!("HMM Config");
    println!("  Splicing Flags: U:{} US:{} S:{} SC:{} C:{} CS:{} S:{} SU:{} U:{}",
             CAN_SPLICE_UTR5, CAN_SPLICE_UTR5_START, CAN_SPLICE_START,
             CAN_SPLICE_START_CODING, CAN_SPLICE_CODING, CAN_SPLICE_CODING_STOP,
             CAN_SPLICE_STOP, CAN_SPLICE_STOP_UTR3, CAN_SPLICE_UTR3);

    println!("  Splicing - Fixed Penalty: Donor {}, Weights: Donor {}, Acceptor {}", DONOR_FIXED_PENALTY, DONOR_WEIGHT, ACCEPTOR_WEIGHT);
    println!("  Coding - Weights: Start {}, Stop {}", START_WEIGHT, STOP_WEIGHT);
    println!("  Phase Mode: Implementation 1, Dilute to Total, Retention: {}", PHASE_RETAIN);
    println!();
}




impl HmmState
{
    fn get_annotation_label(self) -> HmmAnnotationLabel
    {
        let start_base = HmmAnnotationLabel::Coding; //Start;
        let stop_base = HmmAnnotationLabel::Coding; //Stop;

        match self
        {
            HmmState::Intergenic => HmmAnnotationLabel::Intergenic,

            HmmState::UTR5 => HmmAnnotationLabel::UTR5,
            HmmState::UTR5IntronDSS => HmmAnnotationLabel::Intron,
            HmmState::UTR5Intron => HmmAnnotationLabel::Intron,

            HmmState::Start0 => start_base,
            HmmState::Start0IntronDSS => HmmAnnotationLabel::Intron,
            HmmState::Start0Intron => HmmAnnotationLabel::Intron,
            HmmState::Start1 => start_base,
            HmmState::Start1IntronDSS => HmmAnnotationLabel::Intron,
            HmmState::Start1Intron => HmmAnnotationLabel::Intron,
            HmmState::Start2 => start_base,

            HmmState::Coding0 => HmmAnnotationLabel::Coding,
            HmmState::Coding0IntronDSS => HmmAnnotationLabel::Intron,
            HmmState::Coding0Intron => HmmAnnotationLabel::Intron,
            HmmState::Coding1 => HmmAnnotationLabel::Coding,
            HmmState::Coding1IntronDSS => HmmAnnotationLabel::Intron,
            HmmState::Coding1Intron => HmmAnnotationLabel::Intron,
            HmmState::Coding2 => HmmAnnotationLabel::Coding,
            HmmState::Coding2IntronDSS => HmmAnnotationLabel::Intron,
            HmmState::Coding2Intron => HmmAnnotationLabel::Intron,

            HmmState::Stop0T => stop_base,
            HmmState::Stop0TIntronDSS => HmmAnnotationLabel::Intron,
            HmmState::Stop0TIntron => HmmAnnotationLabel::Intron,
            HmmState::Stop1TA => stop_base,
            HmmState::Stop1TAIntronDSS => HmmAnnotationLabel::Intron,
            HmmState::Stop1TAIntron => HmmAnnotationLabel::Intron,
            HmmState::Stop1TG => stop_base,
            HmmState::Stop1TGIntronDSS => HmmAnnotationLabel::Intron,
            HmmState::Stop1TGIntron => HmmAnnotationLabel::Intron,
            HmmState::Stop2 => stop_base,

            HmmState::UTR3 => HmmAnnotationLabel::UTR3,
            HmmState::UTR3IntronDSS => HmmAnnotationLabel::Intron,
            HmmState::UTR3Intron => HmmAnnotationLabel::Intron,
        }
    }


    fn get_state_penalty(self, class_pred: &ClassPredPenalty, phase_pred: &PhasePredPenalty, pred: &PredPenalty) -> f64
    {
        fn coding_phase0(_class_pred: &ClassPredPenalty, _phase_pred: &PhasePredPenalty, pred: &PredPenalty) -> f64 {
            //class_pred.get_coding_penalty() + phase_pred.get_phase0_penalty()
            pred.get_coding_phase0_penalty()
        }

        fn coding_phase1(_class_pred: &ClassPredPenalty, _phase_pred: &PhasePredPenalty, pred: &PredPenalty) -> f64 {
            //class_pred.get_coding_penalty() + phase_pred.get_phase1_penalty()
            pred.get_coding_phase1_penalty()
        }

        fn coding_phase2(_class_pred: &ClassPredPenalty, _phase_pred: &PhasePredPenalty, pred: &PredPenalty) -> f64 {
            //class_pred.get_coding_penalty() + phase_pred.get_phase2_penalty()
            pred.get_coding_phase2_penalty()
        }

        match self
        {
            HmmState::Intergenic=> class_pred.get_intergenic_penalty(),

            HmmState::UTR5=> class_pred.get_utr_penalty(),
            HmmState::UTR5IntronDSS=> class_pred.get_intron_penalty(),
            HmmState::UTR5Intron=> class_pred.get_intron_penalty(),

            HmmState::Start0 => coding_phase0(class_pred, phase_pred, pred),
            HmmState::Start0IntronDSS => class_pred.get_intron_penalty(),
            HmmState::Start0Intron => class_pred.get_intron_penalty(),

            HmmState::Start1 => coding_phase2(class_pred, phase_pred, pred),
            HmmState::Start1IntronDSS => class_pred.get_intron_penalty(),
            HmmState::Start1Intron => class_pred.get_intron_penalty(),

            HmmState::Start2 => coding_phase1(class_pred, phase_pred, pred),

            HmmState::Coding0 => coding_phase0(class_pred, phase_pred, pred),
            HmmState::Coding0IntronDSS => class_pred.get_intron_penalty(),
            HmmState::Coding0Intron => class_pred.get_intron_penalty(),

            HmmState::Coding1 => coding_phase2(class_pred, phase_pred, pred),
            HmmState::Coding1IntronDSS => class_pred.get_intron_penalty(),
            HmmState::Coding1Intron => class_pred.get_intron_penalty(),

            HmmState::Coding2 => coding_phase1(class_pred, phase_pred, pred),
            HmmState::Coding2IntronDSS => class_pred.get_intron_penalty(),
            HmmState::Coding2Intron => class_pred.get_intron_penalty(),

            HmmState::Stop0T => coding_phase0(class_pred, phase_pred, pred),
            HmmState::Stop0TIntronDSS => class_pred.get_intron_penalty(),
            HmmState::Stop0TIntron => class_pred.get_intron_penalty(),

            HmmState::Stop1TA => coding_phase2(class_pred, phase_pred, pred),
            HmmState::Stop1TAIntronDSS => class_pred.get_intron_penalty(),
            HmmState::Stop1TAIntron => class_pred.get_intron_penalty(),

            HmmState::Stop1TG => coding_phase2(class_pred, phase_pred, pred),
            HmmState::Stop1TGIntronDSS => class_pred.get_intron_penalty(),
            HmmState::Stop1TGIntron => class_pred.get_intron_penalty(),

            HmmState::Stop2 => coding_phase1(class_pred, phase_pred, pred),

            HmmState::UTR3=> class_pred.get_utr_penalty(),
            HmmState::UTR3IntronDSS=> class_pred.get_intron_penalty(),
            HmmState::UTR3Intron=> class_pred.get_intron_penalty(),
        }
    }

    fn get_base_count(self) -> usize
    {
        match self
        {
            HmmState::UTR5IntronDSS => 39,
            HmmState::Start0IntronDSS => 39,
            HmmState::Start1IntronDSS => 39,
            HmmState::Coding0IntronDSS => 39,
            HmmState::Coding1IntronDSS => 39,
            HmmState::Coding2IntronDSS => 39,
            HmmState::Stop0TIntronDSS => 39,
            HmmState::Stop1TAIntronDSS => 39,
            HmmState::Stop1TGIntronDSS => 39,
            HmmState::UTR3IntronDSS => 39,
            _ => 1
        }
    }

    fn populate_successor_states_and_transition_penalties(self, trans_ctx: TransitionContext, successors: &mut Vec<(HmmState, f64)>)
    {
        match self
            {
            HmmState::Intergenic =>
                    {
                    successors.push((HmmState::Intergenic, 0.0));
                    successors.push((HmmState::UTR5, 0.0));
                    },

            HmmState::UTR5 =>
                    {
                    successors.push((HmmState::UTR5, 0.0));

                    if let Some(ds) = trans_ctx.get_downstream(1)
                        { successors.push((HmmState::Start0, ds[0].get_a()*START_WEIGHT)); }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_UTR5)
                        { successors.push((HmmState::UTR5IntronDSS, donor_penalty)); }

                    },

            HmmState::UTR5IntronDSS => { successors.push((HmmState::UTR5Intron, 0.0))},

            HmmState::UTR5Intron =>
                    {
                    successors.push((HmmState::UTR5Intron, 0.0));

                    if let Some(acceptor_penalty) = trans_ctx.get_acceptor_penalty_ag_g()
                        {
                        successors.push((HmmState::UTR5, acceptor_penalty));

                        if let (Some(ds), true) = (trans_ctx.get_downstream(1), CAN_SPLICE_UTR5_START)
                            { successors.push((HmmState::Start0, acceptor_penalty+ds[0].get_a()*START_WEIGHT)); }
                        }
                    },

            HmmState::Start0 =>
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        { successors.push((HmmState::Start1, ds[0].get_t()*START_WEIGHT)); }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_START)
                        { successors.push((HmmState::Start0IntronDSS, donor_penalty)); }
                    },

            HmmState::Start0IntronDSS => { successors.push((HmmState::Start0Intron, 0.0))},

            HmmState::Start0Intron =>
                    {
                    successors.push((HmmState::Start0Intron, 0.0));

                    if let (Some(acceptor_penalty), Some(ds)) = (trans_ctx.get_acceptor_penalty_ag_g(), trans_ctx.get_downstream(1))
                        { successors.push((HmmState::Start1, acceptor_penalty + ds[0].get_t()*START_WEIGHT)); }
                    },

            HmmState::Start1 =>
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        { successors.push((HmmState::Start2, ds[0].get_g()*START_WEIGHT)); }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_START)
                        { successors.push((HmmState::Start1IntronDSS, donor_penalty)); }
                    },


            HmmState::Start1IntronDSS => { successors.push((HmmState::Start1Intron, 0.0))},

            HmmState::Start1Intron =>
                    {
                    successors.push((HmmState::Start1Intron, 0.0));

                    if let (Some(acceptor_penalty), Some(ds)) = (trans_ctx.get_acceptor_penalty_ag_g(), trans_ctx.get_downstream(1))
                        {
                        successors.push((HmmState::Start2, acceptor_penalty + ds[0].get_g()*START_WEIGHT));
                        }
                    },

             HmmState::Start2 =>
                    {
                    successors.push((HmmState::Coding0, 0.0));

                    if let Some(ds) = trans_ctx.get_downstream(1)
                        { successors.push((HmmState::Stop0T, ds[0].get_t()*STOP_WEIGHT)); }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_START_CODING)
                        { successors.push((HmmState::Coding2IntronDSS, donor_penalty)); }
                    },

             HmmState::Coding0 =>
                    {
                    successors.push((HmmState::Coding1,0.0));

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_CODING)
                        { successors.push((HmmState::Coding0IntronDSS, donor_penalty)); }
                    }

            HmmState::Coding0IntronDSS => { successors.push((HmmState::Coding0Intron, 0.0))},

            HmmState::Coding0Intron =>
                    {
                    successors.push((HmmState::Coding0Intron, 0.0));

                    if let Some(acceptor_penalty) = trans_ctx.get_acceptor_penalty_ag_g()
                        { successors.push((HmmState::Coding1, acceptor_penalty)); }
                    },

            HmmState::Coding1 =>
                    {
                    successors.push((HmmState::Coding2, 0.0));

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_CODING)
                        { successors.push((HmmState::Coding1IntronDSS, donor_penalty)); }
                    }

            HmmState::Coding1IntronDSS => { successors.push((HmmState::Coding1Intron, 0.0))},

            HmmState::Coding1Intron =>
                    {
                    successors.push((HmmState::Coding1Intron, 0.0));

                    if let Some(acceptor_penalty) = trans_ctx.get_acceptor_penalty_ag_g()
                        { successors.push((HmmState::Coding2, acceptor_penalty)); }
                    },

            HmmState::Coding2 =>
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        {
                        successors.push((HmmState::Coding0, ds[0].get_a()*STOP_WEIGHT));
                        successors.push((HmmState::Coding0, ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding0, ds[0].get_g()*STOP_WEIGHT));
                        successors.push((HmmState::Stop0T, ds[0].get_t()*STOP_WEIGHT));
                        }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_CODING)
                        { successors.push((HmmState::Coding2IntronDSS, donor_penalty)); }
                    }

            HmmState::Coding2IntronDSS => { successors.push((HmmState::Coding2Intron, 0.0))},

            HmmState::Coding2Intron =>
                    {
                    successors.push((HmmState::Coding2Intron, 0.0));

                    if let Some(acceptor_penalty) = trans_ctx.get_acceptor_penalty_ag_g()
                        {
                        successors.push((HmmState::Coding0, acceptor_penalty));

                        if let (Some(ds), true) = (trans_ctx.get_downstream(1), CAN_SPLICE_CODING_STOP)
                            { successors.push((HmmState::Stop0T, ds[0].get_t()*STOP_WEIGHT+acceptor_penalty)); }
                        }
                    },

            HmmState::Stop0T =>  // Equivalent to Coding0, but potentially a stop codon (Txx)
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        {
                        successors.push((HmmState::Stop1TA, ds[0].get_a()*STOP_WEIGHT));
                        successors.push((HmmState::Stop1TG, ds[0].get_g()*STOP_WEIGHT));

                        successors.push((HmmState::Coding1, ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding1, ds[0].get_t()*STOP_WEIGHT));
                        }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_STOP)
                        { successors.push((HmmState::Stop0TIntronDSS, donor_penalty)); }
                    },

            HmmState::Stop0TIntronDSS => { successors.push((HmmState::Stop0TIntron, 0.0))},

            HmmState::Stop0TIntron =>
                    {
                    successors.push((HmmState::Stop0TIntron, 0.0));

                    if let (Some(acceptor_penalty), Some(ds)) = (trans_ctx.get_acceptor_penalty_ag_g(), trans_ctx.get_downstream(1))
                        {
                        successors.push((HmmState::Stop1TA, acceptor_penalty + ds[0].get_a()*STOP_WEIGHT));
                        successors.push((HmmState::Stop1TG, acceptor_penalty + ds[0].get_g()*STOP_WEIGHT));

                        successors.push((HmmState::Coding1, acceptor_penalty + ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding1, acceptor_penalty + ds[0].get_t()*STOP_WEIGHT));
                        }
                    }

            HmmState::Stop1TA =>  // Equivalent to Coding1, but potentially a stop codon (TAx)
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        {
                        successors.push((HmmState::Stop2, ds[0].get_a()*STOP_WEIGHT));
                        successors.push((HmmState::Stop2, ds[0].get_g()*STOP_WEIGHT));

                        successors.push((HmmState::Coding2, ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding2, ds[0].get_t()*STOP_WEIGHT));
                        }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_STOP)
                        { successors.push((HmmState::Stop1TAIntronDSS, donor_penalty)); }

                    },

            HmmState::Stop1TAIntronDSS => { successors.push((HmmState::Stop1TAIntron, 0.0))},

            HmmState::Stop1TAIntron =>
                    {
                    successors.push((HmmState::Stop1TAIntron, 0.0));

                    if let (Some(acceptor_penalty), Some(ds)) = (trans_ctx.get_acceptor_penalty_ag_g(), trans_ctx.get_downstream(1))
                        {
                        successors.push((HmmState::Stop2, acceptor_penalty + ds[0].get_a()*STOP_WEIGHT));
                        successors.push((HmmState::Stop2, acceptor_penalty + ds[0].get_g()*STOP_WEIGHT));

                        successors.push((HmmState::Coding2, acceptor_penalty + ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding2, acceptor_penalty + ds[0].get_t()*STOP_WEIGHT));
                        }
                    }

            HmmState::Stop1TG => // Equivalent to Coding1, but potentially a stop codon (TGx)
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        {
                        successors.push((HmmState::Stop2, ds[0].get_a()*STOP_WEIGHT));

                        successors.push((HmmState::Coding2, ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding2, ds[0].get_g()*STOP_WEIGHT));
                        successors.push((HmmState::Coding2, ds[0].get_t()*STOP_WEIGHT));
                        }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_STOP)
                        { successors.push((HmmState::Stop1TGIntronDSS, donor_penalty)); }
                    },

            HmmState::Stop1TGIntronDSS => { successors.push((HmmState::Stop1TGIntron, 0.0))},

            HmmState::Stop1TGIntron =>
                    {
                    successors.push((HmmState::Stop1TGIntron, 0.0));

                    if let (Some(acceptor_penalty), Some(ds)) = (trans_ctx.get_acceptor_penalty_ag_g(), trans_ctx.get_downstream(1))
                        {
                        successors.push((HmmState::Stop2, acceptor_penalty + ds[0].get_a()*STOP_WEIGHT));

                        successors.push((HmmState::Coding2, acceptor_penalty + ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding2, acceptor_penalty + ds[0].get_g()*STOP_WEIGHT));
                        successors.push((HmmState::Coding2, acceptor_penalty + ds[0].get_t()*STOP_WEIGHT));
                        }
                    }

            HmmState::Stop2 =>
                    {
                    successors.push((HmmState::UTR3, 0.0));

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_STOP_UTR3)
                        { successors.push((HmmState::UTR3IntronDSS, donor_penalty)); }
                    },

            HmmState::UTR3 =>
                    {
                    successors.push((HmmState::UTR3, 0.0));
                    successors.push((HmmState::Intergenic, 0.0));

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_UTR3)
                        { successors.push((HmmState::UTR3IntronDSS, donor_penalty)); }
                    },

            HmmState::UTR3IntronDSS => { successors.push((HmmState::UTR3Intron, 0.0))},

            HmmState::UTR3Intron =>
                    {
                    successors.push((HmmState::UTR3Intron, 0.0));

                    if let Some(acceptor_penalty) = trans_ctx.get_acceptor_penalty_ag_g()
                        { successors.push((HmmState::UTR3, acceptor_penalty)); }
                    },
            }
    }


}








const PENALTY_SCALE: f64=1_000_000.0;   // Convert to u64 to avoid FP annoyances

#[derive(Copy, Clone, PartialEq, Eq)]
pub struct HmmEval
{
    start_position: usize,   // Number of bases produced before this state
    end_position: usize, // Number of bases produced including this state
    state: HmmState,
    previous_state: HmmState,

    penalty: u64,
}

impl HmmEval
{
    fn new_root() -> HmmEval
    {
        HmmEval { start_position: 0, end_position: 0, state: HmmState::Intergenic, previous_state: HmmState::Intergenic, penalty: 0}
    }

    fn new_successor(start_position: usize, end_position: usize, state: HmmState, previous_state: HmmState, penalty: u64) -> HmmEval
    {
        HmmEval { start_position, end_position, state, previous_state, penalty }
    }
}

impl Ord for HmmEval {
    fn cmp(&self, other: &Self) -> Ordering {
        // Order on penalty (lowest), then position (highest)
        other.penalty.cmp(&self.penalty).
            then_with(|| self.end_position.cmp(&other.end_position))
            //.then_with(|| self.state.cmp(&other.state))
    }
}

impl PartialOrd for HmmEval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}




const MAX_EVALS: u64 = 100_000_000;

pub struct PredictionHmm
{
    class_pred_pen: Vec<ClassPredPenalty>,
    phase_pred_pen: Vec<PhasePredPenalty>,

    pred_pen: Vec<PredPenalty>,
    bases_pen: Vec<BasesPenalty>,

    best_eval: Vec<Option<HmmEval>>,

    eval_heap: BinaryHeap<HmmEval>,
}




impl PredictionHmm
{
    pub fn new(bp_vector: Vec<(Bases, ClassPrediction, PhasePrediction)>) -> PredictionHmm
    {
        let mut class_pred_pen = Vec::with_capacity(bp_vector.len());
        let mut phase_pred_pen = Vec::with_capacity(bp_vector.len());
        let mut pred_pen = Vec::with_capacity(bp_vector.len());

        let mut bases_pen = Vec::with_capacity(bp_vector.len());

        for (bases, class_pred, phase_pred) in bp_vector.iter()
        {
            class_pred_pen.push(class_pred.into());
            phase_pred_pen.push(phase_pred.into());
            pred_pen.push((class_pred, phase_pred).into());

            bases_pen.push(bases.into());
        }

        let total_states = (bp_vector.len()+1) * HMM_STATES;
        let best_eval = vec![None; total_states];

        let eval_heap = BinaryHeap::new();
        PredictionHmm { class_pred_pen, phase_pred_pen, pred_pen, bases_pen, best_eval, eval_heap }
    }

    fn consider_eval(&mut self, eval: HmmEval)
    {
        let idx = eval.end_position * HMM_STATES + (eval.state as usize);

        let maybe_old_eval = &self.best_eval[idx];

        if let Some(old_eval) = maybe_old_eval
            {
            if eval.penalty >= old_eval.penalty
                { return }
            }

        self.best_eval[idx]=Some(eval);
        self.eval_heap.push(eval);
    }

    fn is_eval_current(&self, eval: &HmmEval) -> bool
    {
        let idx = eval.end_position * HMM_STATES + (eval.state as usize);
        let best_eval = &self.best_eval[idx].expect("Eval from heap not in best_eval");

        best_eval == eval
    }

    fn process_eval(&mut self, eval: &HmmEval)
    {
        if !self.is_eval_current(eval)
            { return; }

        let trans_ctx = TransitionContext::new(&self.class_pred_pen, &self.phase_pred_pen, &self.pred_pen, &self.bases_pen, eval.end_position);

        let mut successors = Vec::with_capacity(HMM_STATES);
        eval.state.populate_successor_states_and_transition_penalties(trans_ctx, &mut successors);

        for (next_state, trans_penalty) in successors.into_iter()
            {
            let mut extra_penalty = trans_penalty;

            let start_position = eval.end_position;
            let end_position = start_position + next_state.get_base_count();

            if end_position <= self.class_pred_pen.len()  // Drop 'long' state picked near end
                {
                for pos in start_position..end_position
                    { extra_penalty += next_state.get_state_penalty(&self.class_pred_pen[pos], &self.phase_pred_pen[pos], &self.pred_pen[pos]); }

                let penalty = eval.penalty + ((extra_penalty * PENALTY_SCALE) as u64);

                let next_eval = HmmEval::new_successor(start_position, end_position, next_state, eval.state, penalty);

                self.consider_eval(next_eval);
                }
            }
    }

    pub fn solve(mut self) -> Option<PredictionHmmSolution>
    {
        let initial_eval = HmmEval::new_root();
        self.consider_eval(initial_eval);

        let mut evals = 0;

        while evals < MAX_EVALS
        {
            if let Some(eval) = self.eval_heap.pop()
                {
                if eval.end_position == self.class_pred_pen.len()
                    { return Some(PredictionHmmSolution::new(self, eval)); }

                self.process_eval(&eval);
                }
            else
                { break; }  // Nothing left to do

            evals+=1;
        }

        panic!("MAX_EVALS exceeded - raise limit or window thresholds");
//        None
    }
}

pub struct HmmStateRegion
{
    start_pos: usize,
    end_pos: usize,
    annotation_label: HmmAnnotationLabel
}

impl HmmStateRegion
{
    fn new(start_pos: usize, end_pos: usize, base_state: HmmAnnotationLabel) -> HmmStateRegion
    {
        HmmStateRegion { start_pos, end_pos, annotation_label: base_state }
    }

    pub fn get_start_pos(&self) -> usize { self.start_pos }

    pub fn get_end_pos(&self) -> usize { self.end_pos }

    pub fn get_annotation_label(&self) -> HmmAnnotationLabel { self.annotation_label }

    pub fn len(&self) -> usize { self.end_pos - self.start_pos }
}


pub struct PredictionHmmSolution
{
    hmm: PredictionHmm,
    eval: HmmEval
}

impl PredictionHmmSolution
{
    fn new(hmm: PredictionHmm, eval: HmmEval) -> PredictionHmmSolution
    {
        PredictionHmmSolution { hmm, eval }
    }

    pub fn trace_regions(&self) -> Vec<HmmStateRegion>
    {
        let mut regions = Vec::new();

        let mut eval = &self.eval;
        let mut region_end_pos = eval.end_position;

        // State positions are 1 above base positions - state.end_pos = 0 is the dummy start, state end_pos 1 is the first real assignment

        // End position is exclusive

        while eval.end_position > 0
            {
            if eval.state.get_annotation_label()!=eval.previous_state.get_annotation_label() // If prev state changes annotation label
                {
                regions.push(HmmStateRegion::new(eval.start_position, region_end_pos, eval.state.get_annotation_label()));
                region_end_pos = eval.start_position;
                }

            let prev_position = eval.start_position;
            let idx = prev_position * HMM_STATES + (eval.previous_state as usize);


            eval = self.hmm.best_eval.get(idx).unwrap().as_ref().unwrap();

            }

        if region_end_pos > 1
            { regions.push(HmmStateRegion::new(0, region_end_pos-1, eval.state.get_annotation_label())); }

        regions.reverse();

        regions
    }


    pub fn split_genes(regions: Vec<HmmStateRegion>) -> Vec<Vec<HmmStateRegion>>
    {
        let mut vec_of_vecs = Vec::new();

        let mut current_vec = Vec::new();
        for region in regions
            {
            if region.annotation_label == HmmAnnotationLabel::Intergenic
                {
                if current_vec.len()>0
                    {
                    vec_of_vecs.push(current_vec);
                    current_vec=Vec::new();
                    }
                }
            else
                {
                current_vec.push(region);
                }
            }

        if current_vec.len()>0
            { vec_of_vecs.push(current_vec); }

        vec_of_vecs
    }




    pub fn dump(&self, position: usize)
    {
        println!("Solution Penalty: {} over {} bp starting at {}", self.eval.penalty, self.hmm.class_pred_pen.len(), position);

        let regions = self.trace_regions();

        for region in regions.iter()
            {
//            println!("{} to {} is {}", region.start_pos+position+1, region.end_pos+position, region.state.to_str()); // Biologist coordinates

            let mut seq = String::with_capacity(region.end_pos-region.start_pos);
            for idx in region.start_pos .. region.end_pos
                {
                seq.push(self.hmm.bases_pen[idx].as_str());
                }

            println!("{} to {} aka {} to {} is {} - {}", region.start_pos, region.end_pos, region.start_pos+position, region.end_pos+position, region.annotation_label.to_str(), seq);

            }

//        if position > 30000
//            { panic!("First 30k"); }
        panic!("Show only one region");
    }




}




/*

    base[0]=A
    base[1]=T
    base[2]=G


0   UTR5    - : A
1   Start1  A : T
    Start2  AT : G
    Coding0 ATG : -
    Coding1 ATG - : -
    Coding2 ATG -- : -
    Coding0 ATG --- : T
    StopT   ATG --- T : A
    StopTA  ATG --- TA : A
    Stop3   ATG --- TAA : -
    UTR3    ATG --- TAA -

                           0120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120120
    191 to 345 is Coding - ATGGAGGATCAAGTTGGGTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTTGGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATCAGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCT

                           12012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012
    427 to 708 is Coding - TCCAGTCAAAGTACAAATCGAGAGATGCTATGTGGTACTTCTTCTCTCGTAGAGAAAACAACAAAGGGAATCGACAGAGCAGGACAACGGTTTCTGGTAAATGGAAGCTTACCGGAGAATCTGTTGAGGTCAAGGACCAGTGGGGATTTTGTAGTGAGGGCTTTCGTGGTAAGATTGGTCATAAAAGGGTTTTGGTGTTCCTCGATGGAAGATACCCTGACAAAACCAAATCTGATTGGGTTATCCACGAGTTCCACTACGACCTCTTACCAGAACATCAG

                           012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012
   917 to 1037 is Coding - AGGACATATGTCATCTGCAGACTTGAGTACAAGGGTGATGATGCGGACATTCTATCTGCTTATGCAATAGATCCCACTCCCGCTTTTGTCCCCAATATGACTAGTAGTGCAGGTTCTGTG

                           012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012
  1137 to 1527 is Coding - GTCAACCAATCACGTCAACGAAATTCAGGATCTTACAACACTTACTCTGAGTATGATTCAGCAAATCATGGCCAGCAGTTTAATGAAAACTCTAACATTATGCAGCAGCAACCACTTCAAGGATCATTCAACCCTCTCCTTGAGTATGATTTTGCAAATCACGGCGGTCAGTGGCTGAGTGACTATATCGACCTGCAACAGCAAGTTCCTTACTTGGCACCTTATGAAAATGAGTCGGAGATGATTTGGAAGCATGTGATTGAAGAAAATTTTGAGTTTTTGGTAGATGAAAGGACATCTATGCAACAGCATTACAGTGATCACCGGCCCAAAAAACCTGTGTCTGGGGTTTTGCCTGATGATAGCAGTGATACTGAAACTGGATCAATG

                           012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012
  1605 to 1758 is Coding - ATTTTCGAAGACACTTCGAGCTCCACTGATAGTGTTGGTAGTTCAGATGAACCGGGCCATACTCGTATAGATGATATTCCATCATTGAACATTATTGAGCCTTTGCACAATTATAAGGCACAAGAGCAACCAAAGCAGCAGAGCAAAGAAAAG

                         012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012012
1870 to 2059 is Coding - GTGATAAGTTCGCAGAAAAGCGAATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCTCCATCCACCAACACGGTGAAGCAGAGCTGGATTGTTTTGGAGAATGCACAGTGGAACTATCTCAAGAACATGATCATTGGTGTCTTGTTGTTCATCTCCGTCATTAGTTGGATCATTCTTGTTGGT

2059 to 2062 is Stop - TAA

3760 to 3762 is Start
3763 to 3913 is Coding
3996 to 4276 is Coding
4486 to 4605 is Coding
4706 to 5095 is Coding
5174 to 5326 is Coding
5439 to 5627 is Coding
5628 to 5630 is Stop

Chr1	phytozomev10	CDS	3760	3913
Chr1	phytozomev10	CDS	3996	4276
Chr1	phytozomev10	CDS	4486	4605
Chr1	phytozomev10	CDS	4706	5095
Chr1	phytozomev10	CDS	5174	5326
Chr1	phytozomev10	CDS	5439	5630

23136 to 23349 is UTR5
23350 to 23514 is Intron
23515 to 23524 is UTR5
23525 to 24451 is Coding
24542 to 24655 is Coding
24752 to 24962 is Coding
25041 to 25435 is Coding
25524 to 25743 is Coding
25825 to 25997 is Coding
26081 to 26203 is Coding
26292 to 26452 is Coding
26543 to 26776 is Coding
26862 to 27012 is Coding
27099 to 27281 is Coding
27372 to 27533 is Coding
27618 to 27713 is Coding
27803 to 28431 is Coding
28708 to 28805 is Coding
28890 to 29080 is Coding
29193 to 30065 is Coding
30147 to 30311 is Coding
30410 to 30816 is Coding
30902 to 31079 is Coding
31080 to 31200 is UTR3


                                                                    YYYYYYYYNCAG                      YYYYYYYYNCAG
 1234567890123456789012345678901234567890123456789012345678901234567890123456789 012345678901234567890123456789012 34567
 GTATATATATATATATTATGCTTAGTGTCTTTTTTTTTTTTGTTGAAACTATCTAATCATATTTGGTATATATATGTAG ATTCTTGAAGCCTTGACTGCCGCCTCGTGCCAG GAAAC

                                                                                            YTRAC




Chr1    HelixerPost     CDS     8594    8646    .       -       0       ID=Athaliana_Chr1_006990.1.CDS.1;Parent=Athaliana_Chr1_006990.1
Chr1    HelixerPost     CDS     8417    8464    .       -       2       ID=Athaliana_Chr1_006990.1.CDS.2;Parent=Athaliana_Chr1_006990.1
Chr1    HelixerPost     CDS     8236    8325    .       -       2       ID=Athaliana_Chr1_006990.1.CDS.3;Parent=Athaliana_Chr1_006990.1

Chr1    HelixerPost     CDS     7762    7886    .       -       2       ID=Athaliana_Chr1_006990.1.CDS.4;Parent=Athaliana_Chr1_006990.1

Chr1    HelixerPost     CDS     7564    7649    .       -       1       ID=Athaliana_Chr1_006990.1.CDS.5;Parent=Athaliana_Chr1_006990.1
Chr1    HelixerPost     CDS     7384    7450    .       -       0       ID=Athaliana_Chr1_006990.1.CDS.6;Parent=Athaliana_Chr1_006990.1
Chr1    HelixerPost     CDS     7159    7232    .       -       1       ID=Athaliana_Chr1_006990.1.CDS.7;Parent=Athaliana_Chr1_006990.1

Chr1    HelixerPost     CDS     6428    6655    .       -       0       ID=Athaliana_Chr1_006991.1.CDS.1;Parent=Athaliana_Chr1_006991.1


Chr1    phytozomev10    CDS     8571    8666    .       -       0       ID=AT1G01020.1.TAIR10.CDS.1;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     8417    8464    .       -       0       ID=AT1G01020.1.TAIR10.CDS.2;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     8236    8325    .       -       0       ID=AT1G01020.1.TAIR10.CDS.3;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     7942    7987    .       -       0       ID=AT1G01020.1.TAIR10.CDS.4;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     7762    7835    .       -       2       ID=AT1G01020.1.TAIR10.CDS.5;Parent=AT1G01020.1.TAIR10;pacid=19655142

Chr1    phytozomev10    CDS     7564    7649    .       -       0       ID=AT1G01020.1.TAIR10.CDS.6;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     7384    7450    .       -       1       ID=AT1G01020.1.TAIR10.CDS.7;Parent=AT1G01020.1.TAIR10;pacid=19655142
Chr1    phytozomev10    CDS     7157    7232    .       -       0       ID=AT1G01020.1.TAIR10.CDS.8;Parent=AT1G01020.1.TAIR10;pacid=19655142

Chr1    phytozomev10    CDS     6915    7069    .       -       2       ID=AT1G01020.1.TAIR10.CDS.9;Parent=AT1G01020.1.TAIR10;pacid=19655142











Chr1    HelixerPost     CDS     32547   32670   .       -       0       ID=Athaliana_Chr1_006988.1.CDS.1;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     32431   32477   .       -       1       ID=Athaliana_Chr1_006988.1.CDS.2;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     32282   32347   .       -       0       ID=Athaliana_Chr1_006988.1.CDS.3;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     32088   32195   .       -       0       ID=Athaliana_Chr1_006988.1.CDS.4;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     31933   31998   .       -       0       ID=Athaliana_Chr1_006988.1.CDS.5;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     31693   31813   .       -       0       ID=Athaliana_Chr1_006988.1.CDS.6;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     31521   31602   .       -       1       ID=Athaliana_Chr1_006988.1.CDS.7;Parent=Athaliana_Chr1_006988.1
Chr1    HelixerPost     CDS     31382   31424   .       -       2       ID=Athaliana_Chr1_006988.1.CDS.8;Parent=Athaliana_Chr1_006988.1

Chr1    phytozomev10    CDS     32547   32670   .       -       0       ID=AT1G01050.1.TAIR10.CDS.1;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     32431   32459   .       -       2       ID=AT1G01050.1.TAIR10.CDS.2;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     32282   32347   .       -       0       ID=AT1G01050.1.TAIR10.CDS.3;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     32088   32195   .       -       0       ID=AT1G01050.1.TAIR10.CDS.4;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     31933   31998   .       -       0       ID=AT1G01050.1.TAIR10.CDS.5;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     31693   31813   .       -       0       ID=AT1G01050.1.TAIR10.CDS.6;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     31521   31602   .       -       2       ID=AT1G01050.1.TAIR10.CDS.7;Parent=AT1G01050.1.TAIR10;pacid=19652974
Chr1    phytozomev10    CDS     31382   31424   .       -       1       ID=AT1G01050.1.TAIR10.CDS.8;Parent=AT1G01050.1.TAIR10;pacid=19652974




Chr1    HelixerPost     gene    23136   31200   .       +       .       ID=Athaliana_Chr1_000002
Chr1    HelixerPost     mRNA    23136   31200   .       +       .       ID=Athaliana_Chr1_000002.1;Parent=Athaliana_Chr1_000002

Chr1    HelixerPost     exon    23136   23265   .       +       .       ID=Athaliana_Chr1_000002.1.exon.1;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     five_prime_UTR  23136   23264   .       +       .       ID=Athaliana_Chr1_000002.1.five_prime_UTR.1;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     CDS     23265   23265   .       +       0       ID=Athaliana_Chr1_000002.1.CDS.1;Parent=Athaliana_Chr1_000002.1

Chr1    HelixerPost     exon    23267   23339   .       +       .       ID=Athaliana_Chr1_000002.1.exon.2;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     CDS     23267   23339   .       +       2       ID=Athaliana_Chr1_000002.1.CDS.2;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     exon    23515   24451   .       +       .       ID=Athaliana_Chr1_000002.1.exon.3;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     CDS     23515   24451   .       +       1       ID=Athaliana_Chr1_000002.1.CDS.3;Parent=Athaliana_Chr1_000002.1





Chr1    HelixerPost     mRNA    23136   31200   .       +       .       ID=Athaliana_Chr1_000002.1;Parent=Athaliana_Chr1_000002
Chr1    HelixerPost     exon    23136   23344   .       +       .       ID=Athaliana_Chr1_000002.1.exon.1;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     five_prime_UTR  23136   23344   .       +       .       ID=Athaliana_Chr1_000002.1.five_prime_UTR.1;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     exon    23515   24451   .       +       .       ID=Athaliana_Chr1_000002.1.exon.2;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     five_prime_UTR  23515   23524   .       +       .       ID=Athaliana_Chr1_000002.1.five_prime_UTR.2;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     CDS     23525   24451   .       +       0       ID=Athaliana_Chr1_000002.1.CDS.1;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     exon    24546   24655   .       +       .       ID=Athaliana_Chr1_000002.1.exon.3;Parent=Athaliana_Chr1_000002.1
Chr1    HelixerPost     CDS     24546   24655   .       +       0       ID=Athaliana_Chr1_000002.1.CDS.2;Parent=Athaliana_Chr1_000002.1






 */
