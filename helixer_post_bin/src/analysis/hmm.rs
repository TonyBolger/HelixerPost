use crate::results::conv::{Bases, Prediction};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::num::NonZeroU32;

#[derive(Clone, Copy)]
struct PredictionPenalty
{
    penalty: [f64; 4] // Ordering is intergenic, utr, coding, intron
}

impl PredictionPenalty
{
    pub fn get_intergenic_penalty(&self) -> f64 { self.penalty[0] }

    pub fn get_utr_penalty(&self) -> f64 { self.penalty[1] }

    pub fn get_coding_penalty(&self) -> f64 { self.penalty[2] }

    pub fn get_intron_penalty(&self) -> f64 { self.penalty[3] }
}

const PREDICTION_PROB_FLOOR: f64 = 0.000_000_001; // Prevent infinite penalties

impl From<&Prediction> for PredictionPenalty
{
    fn from(pred: &Prediction) -> Self {

        let raw_pred = pred.get();

        let mut penalty = [0.0 ; 4];
        for i in 0..4
            {
            let adjusted_pred = raw_pred[i] as f64;

            let adjusted_pred = if adjusted_pred > PREDICTION_PROB_FLOOR { adjusted_pred } else { PREDICTION_PROB_FLOOR };
            penalty[i] = -f64::log2(adjusted_pred);
            }

        let mut min_penalty = penalty[0];
        for i in 1..4
            { if penalty[i] < min_penalty { min_penalty = penalty[i] } }


        for i in 0..4
            { penalty[i]-=min_penalty; }

        PredictionPenalty { penalty }
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
    base_pen: &'a [BasesPenalty],
    offset: usize
}

impl<'a> TransitionContext<'a>
{
    fn new(base_pen: &'a [BasesPenalty],  offset: usize)-> TransitionContext<'a>
    {
        TransitionContext { base_pen, offset}
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

#[allow(dead_code)]
#[derive(Clone, Copy, Eq, PartialEq)]
enum HmmIntron
{
    None,
    Intron(NonZeroU32)
}


impl HmmIntron
{
    fn len(self) -> usize
    {
        match self
        {
            HmmIntron::None => 0,
            HmmIntron::Intron(l) => l.get() as usize
        }
    }
}


const HMM_STATES: usize = 23;

#[derive(Clone, Copy, Eq, PartialEq, Ord)]
enum HmmState
{
    Intergenic = 0,
    UTR5 = 1,
    UTR5Intron = 2,
    Start1 = 3, // After A
    Start1Intron = 4,
    Start2 = 5, // After AT
    Start2Intron = 6,
    Start3 = 7,
    Coding0 = 8,
    Coding0Intron = 9,
    Coding1 = 10,
    Coding1Intron = 11,
    Coding2 = 12,
    Coding2Intron = 13,
    StopT = 14,
    StopTIntron = 15,
    StopTA = 16,
    StopTAIntron = 17,
    StopTG = 18,
    StopTGIntron = 19,
    Stop3 = 20,
    UTR3 = 21,
    UTR3Intron = 22
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

const START_WEIGHT: f64 = 10_000.0;

const DONOR_FIXED_PENALTY: f64 = 100.0;
const DONOR_WEIGHT: f64 = 10.0;

const ACCEPTOR_WEIGHT: f64 = 10.0;

const STOP_WEIGHT: f64 = 10_000.0;

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
            HmmState::UTR5Intron => HmmAnnotationLabel::Intron,

            HmmState::Start1 => start_base,
            HmmState::Start1Intron => HmmAnnotationLabel::Intron,
            HmmState::Start2 => start_base,
            HmmState::Start2Intron => HmmAnnotationLabel::Intron,
            HmmState::Start3 => start_base,

            HmmState::Coding0 => HmmAnnotationLabel::Coding,
            HmmState::Coding0Intron => HmmAnnotationLabel::Intron,
            HmmState::Coding1 => HmmAnnotationLabel::Coding,
            HmmState::Coding1Intron => HmmAnnotationLabel::Intron,
            HmmState::Coding2 => HmmAnnotationLabel::Coding,
            HmmState::Coding2Intron => HmmAnnotationLabel::Intron,

            HmmState::StopT => stop_base,
            HmmState::StopTIntron => HmmAnnotationLabel::Intron,
            HmmState::StopTA => stop_base,
            HmmState::StopTAIntron => HmmAnnotationLabel::Intron,
            HmmState::StopTG => stop_base,
            HmmState::StopTGIntron => HmmAnnotationLabel::Intron,
            HmmState::Stop3 => stop_base,

            HmmState::UTR3 => HmmAnnotationLabel::UTR3,
            HmmState::UTR3Intron => HmmAnnotationLabel::Intron,
        }
    }

    fn get_state_penalty(self, pred: &PredictionPenalty) -> f64
    {
        match self
        {
            HmmState::Intergenic=> pred.get_intergenic_penalty(),

            HmmState::UTR5=> pred.get_utr_penalty(),
            HmmState::UTR5Intron=> pred.get_intron_penalty(),

            HmmState::Start1=> pred.get_coding_penalty(),
            HmmState::Start1Intron => pred.get_intron_penalty(),
            HmmState::Start2=> pred.get_coding_penalty(),
            HmmState::Start2Intron => pred.get_intron_penalty(),
            HmmState::Start3=> pred.get_coding_penalty(),

            HmmState::Coding0=> pred.get_coding_penalty(),
            HmmState::Coding0Intron => pred.get_intron_penalty(),
            HmmState::Coding1=> pred.get_coding_penalty(),
            HmmState::Coding1Intron => pred.get_intron_penalty(),
            HmmState::Coding2=> pred.get_coding_penalty(),
            HmmState::Coding2Intron => pred.get_intron_penalty(),

            HmmState::StopT => pred.get_coding_penalty(),
            HmmState::StopTIntron => pred.get_intron_penalty(),
            HmmState::StopTA => pred.get_coding_penalty(),
            HmmState::StopTAIntron => pred.get_intron_penalty(),
            HmmState::StopTG => pred.get_coding_penalty(),
            HmmState::StopTGIntron => pred.get_intron_penalty(),
            HmmState::Stop3 => pred.get_coding_penalty(),

            HmmState::UTR3=> pred.get_utr_penalty(),
            HmmState::UTR3Intron=> pred.get_intron_penalty(),
        }
    }


    /*
    fn psstp_intron_donor(state: HmmPrimaryState, can_flag: bool)
    {

    }
    */

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
                        { successors.push((HmmState::Start1, ds[0].get_a()*START_WEIGHT)); }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_UTR5)
                        { successors.push((HmmState::UTR5Intron, donor_penalty)); }
                    },

            HmmState::UTR5Intron =>
                    {
                    successors.push((HmmState::UTR5Intron, 0.0));

                    if let Some(acceptor_penalty) = trans_ctx.get_acceptor_penalty_ag_g()
                        {
                        successors.push((HmmState::UTR5, acceptor_penalty));

                        if let (Some(ds), true) = (trans_ctx.get_downstream(1), CAN_SPLICE_UTR5_START)
                            { successors.push((HmmState::Start1, acceptor_penalty+ds[0].get_a()*START_WEIGHT)); }
                        }
                    },

            HmmState::Start1 =>
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        { successors.push((HmmState::Start2, ds[0].get_t()*START_WEIGHT)); }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_START)
                        { successors.push((HmmState::Start1Intron, donor_penalty)); }
                    },

            HmmState::Start1Intron =>
                    {
                    successors.push((HmmState::Start1Intron, 0.0));

                    if let (Some(acceptor_penalty), Some(ds)) = (trans_ctx.get_acceptor_penalty_ag_g(), trans_ctx.get_downstream(1))
                        { successors.push((HmmState::Start2, acceptor_penalty + ds[0].get_t()*START_WEIGHT)); }
                    },

            HmmState::Start2 =>
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        { successors.push((HmmState::Start3, ds[0].get_g()*START_WEIGHT)); }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_START)
                        { successors.push((HmmState::Start2Intron, donor_penalty)); }
                    },

            HmmState::Start2Intron =>
                    {
                    successors.push((HmmState::Start2Intron, 0.0));

                    if let (Some(acceptor_penalty), Some(ds)) = (trans_ctx.get_acceptor_penalty_ag_g(), trans_ctx.get_downstream(1))
                        {
                        successors.push((HmmState::Start3, acceptor_penalty + ds[0].get_g()*START_WEIGHT));
                        }
                    },

             HmmState::Start3 =>
                    {
                    successors.push((HmmState::Coding1, 0.0));

                    if let Some(ds) = trans_ctx.get_downstream(1)
                        { successors.push((HmmState::StopT, ds[0].get_t()*STOP_WEIGHT)); }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_START_CODING)
                        { successors.push((HmmState::Coding0Intron, donor_penalty)); }
                    },

             HmmState::Coding0 =>
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        {
                            successors.push((HmmState::Coding1, ds[0].get_a()*STOP_WEIGHT));
                            successors.push((HmmState::Coding1, ds[0].get_c()*STOP_WEIGHT));
                            successors.push((HmmState::Coding1, ds[0].get_g()*STOP_WEIGHT));
                            successors.push((HmmState::StopT, ds[0].get_t()*STOP_WEIGHT));
                        }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_CODING)
                        { successors.push((HmmState::Coding0Intron, donor_penalty)); }
                    }

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
                        { successors.push((HmmState::Coding1Intron, donor_penalty)); }
                    }

            HmmState::Coding1Intron =>
                    {
                    successors.push((HmmState::Coding1Intron, 0.0));

                    if let Some(acceptor_penalty) = trans_ctx.get_acceptor_penalty_ag_g()
                        { successors.push((HmmState::Coding2, acceptor_penalty)); }
                    },

            HmmState::Coding2 =>
                    {
                    successors.push((HmmState::Coding0, 0.0));

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_CODING)
                        { successors.push((HmmState::Coding2Intron, donor_penalty)); }
                    }

            HmmState::Coding2Intron =>
                    {
                    successors.push((HmmState::Coding2Intron, 0.0));

                    if let Some(acceptor_penalty) = trans_ctx.get_acceptor_penalty_ag_g()
                        {
                        successors.push((HmmState::Coding0, acceptor_penalty));

                        if let (Some(ds), true) = (trans_ctx.get_downstream(1), CAN_SPLICE_CODING_STOP)
                            { successors.push((HmmState::StopT, ds[0].get_t()*STOP_WEIGHT+acceptor_penalty)); }
                        }
                    },

            HmmState::StopT =>  // Equivalent to Coding1, but potentially a stop codon (Txx)
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        {
                        successors.push((HmmState::StopTA, ds[0].get_a()*STOP_WEIGHT));
                        successors.push((HmmState::StopTG, ds[0].get_g()*STOP_WEIGHT));

                        successors.push((HmmState::Coding2, ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding2, ds[0].get_t()*STOP_WEIGHT));
                        }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_STOP)
                        { successors.push((HmmState::StopTIntron, donor_penalty)); }
                    },

            HmmState::StopTIntron =>
                    {
                    successors.push((HmmState::StopTIntron, 0.0));

                    if let (Some(acceptor_penalty), Some(ds)) = (trans_ctx.get_acceptor_penalty_ag_g(), trans_ctx.get_downstream(1))
                        {
                        successors.push((HmmState::StopTA, acceptor_penalty + ds[0].get_a()*STOP_WEIGHT));
                        successors.push((HmmState::StopTG, acceptor_penalty + ds[0].get_g()*STOP_WEIGHT));

                        successors.push((HmmState::Coding2, acceptor_penalty + ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding2, acceptor_penalty + ds[0].get_t()*STOP_WEIGHT));
                        }
                    }

            HmmState::StopTA =>  // Equivalent to Coding2, but potentially a stop codon (TAx)
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        {
                        successors.push((HmmState::Stop3, ds[0].get_a()*STOP_WEIGHT));
                        successors.push((HmmState::Stop3, ds[0].get_g()*STOP_WEIGHT));

                        successors.push((HmmState::Coding0, ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding0, ds[0].get_t()*STOP_WEIGHT));
                        }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_STOP)
                        { successors.push((HmmState::StopTAIntron, donor_penalty)); }

                    },

            HmmState::StopTAIntron =>
                    {
                    successors.push((HmmState::StopTAIntron, 0.0));

                    if let (Some(acceptor_penalty), Some(ds)) = (trans_ctx.get_acceptor_penalty_ag_g(), trans_ctx.get_downstream(1))
                        {
                        successors.push((HmmState::Stop3, acceptor_penalty + ds[0].get_a()*STOP_WEIGHT));
                        successors.push((HmmState::Stop3, acceptor_penalty + ds[0].get_g()*STOP_WEIGHT));

                        successors.push((HmmState::Coding0, acceptor_penalty + ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding0, acceptor_penalty + ds[0].get_t()*STOP_WEIGHT));
                        }
                    }

            HmmState::StopTG => // Equivalent to Coding2, but potentially a stop codon (TGx)
                    {
                    if let Some(ds) = trans_ctx.get_downstream(1)
                        {
                        successors.push((HmmState::Stop3, ds[0].get_a()*STOP_WEIGHT));

                        successors.push((HmmState::Coding0, ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding0, ds[0].get_g()*STOP_WEIGHT));
                        successors.push((HmmState::Coding0, ds[0].get_t()*STOP_WEIGHT));
                        }

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_STOP)
                        { successors.push((HmmState::StopTGIntron, donor_penalty)); }
                    },

            HmmState::StopTGIntron =>
                    {
                    successors.push((HmmState::StopTGIntron, 0.0));

                    if let (Some(acceptor_penalty), Some(ds)) = (trans_ctx.get_acceptor_penalty_ag_g(), trans_ctx.get_downstream(1))
                        {
                        successors.push((HmmState::Stop3, acceptor_penalty + ds[0].get_a()*STOP_WEIGHT));

                        successors.push((HmmState::Coding0, acceptor_penalty + ds[0].get_c()*STOP_WEIGHT));
                        successors.push((HmmState::Coding0, acceptor_penalty + ds[0].get_g()*STOP_WEIGHT));
                        successors.push((HmmState::Coding0, acceptor_penalty + ds[0].get_t()*STOP_WEIGHT));
                        }
                    }

            HmmState::Stop3 =>
                    {
                    successors.push((HmmState::UTR3, 0.0));

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_STOP_UTR3)
                        { successors.push((HmmState::UTR3Intron, donor_penalty)); }
                    },

            HmmState::UTR3 =>
                    {
                    successors.push((HmmState::UTR3, 0.0));
                    successors.push((HmmState::Intergenic, 0.0));

                    if let Some(donor_penalty) = trans_ctx.get_donor_penalty_ag_gt(CAN_SPLICE_UTR3)
                        { successors.push((HmmState::UTR3Intron, donor_penalty)); }
                    },

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
    position: usize,   // Number of bases _already_ produced
    intron: HmmIntron, // Intron (if any) between this and prior state (not populated _yet_)
    state: HmmState,
    previous_state: HmmState,

    penalty: u64,
}

impl HmmEval
{
    fn new_root() -> HmmEval
    {
        HmmEval { position: 0, intron: HmmIntron::None, state: HmmState::Intergenic, previous_state: HmmState::Intergenic, penalty: 0}
    }

    fn new_successor(prev: &HmmEval, state: HmmState, extra_penalty: f64) -> HmmEval
    {
        let position = prev.position+1;
        let previous_state = prev.state;
        let penalty = prev.penalty+((extra_penalty * PENALTY_SCALE) as u64);

        HmmEval { position, intron: HmmIntron::None, state, previous_state, penalty }
    }
}

impl Ord for HmmEval {
    fn cmp(&self, other: &Self) -> Ordering {
        // Order on penalty (lowest), then position (highest)
        other.penalty.cmp(&self.penalty).
            then_with(|| self.position.cmp(&other.position))
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
    pred_pen: Vec<PredictionPenalty>,
    bases_pen: Vec<BasesPenalty>,

    best_eval: Vec<Option<HmmEval>>,

    eval_heap: BinaryHeap<HmmEval>,
}


impl PredictionHmm
{
    pub fn new(bp_vector: Vec<(Bases, Prediction)>) -> PredictionHmm
    {
        let mut pred_pen = Vec::with_capacity(bp_vector.len());
        let mut bases_pen = Vec::with_capacity(bp_vector.len());
        for (bases, predictions) in bp_vector.iter()
        {
            pred_pen.push(predictions.into());
            bases_pen.push(bases.into());
        }

        let total_states = (bp_vector.len()+1) * HMM_STATES;
        let best_eval = vec![None; total_states];

        let eval_heap = BinaryHeap::new();
        PredictionHmm { pred_pen, bases_pen, best_eval, eval_heap }
    }

    fn consider_eval(&mut self, eval: HmmEval)
    {
        let idx = eval.position * HMM_STATES + (eval.state as usize);

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
        let idx = eval.position * HMM_STATES + (eval.state as usize);
        let best_eval = &self.best_eval[idx].expect("Eval from heap not in best_eval");

        best_eval == eval
    }

    fn process_eval(&mut self, eval: &HmmEval)
    {
        if !self.is_eval_current(eval)
            { return; }

        let trans_ctx = TransitionContext::new(&self.bases_pen, eval.position);

        let mut successors = Vec::with_capacity(HMM_STATES);
        eval.state.populate_successor_states_and_transition_penalties(trans_ctx, &mut successors);

        for (next_state, trans_penalty) in successors.into_iter()
            {
            let state_penalty = next_state.get_state_penalty(&self.pred_pen[eval.position]);
            let next_eval = HmmEval::new_successor(eval, next_state, trans_penalty + state_penalty);

            self.consider_eval(next_eval);
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
                if eval.position == self.pred_pen.len()
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
        let mut region_end_pos = eval.position;

        // State positions are 1 above base positions - state.pos = 0 is the dummy start, state pos 1 is the first real assignment

        // End position is exclusive

        while eval.position > 0
            {
            if eval.intron!=HmmIntron::None // Break region if current state contains an intron or ...
                {
                regions.push(HmmStateRegion::new(eval.position-1, region_end_pos-1, eval.state.get_annotation_label()));

                let intron_start = eval.position - eval.intron.len();
                let intron_end = eval.position - 1;
                regions.push(HmmStateRegion::new(intron_start-1, intron_end-1, HmmAnnotationLabel::Intron));

                region_end_pos = intron_start;
                }
            else if eval.state.get_annotation_label()!=eval.previous_state.get_annotation_label() // ... if next state changes annotation label
                {
                regions.push(HmmStateRegion::new(eval.position-1, region_end_pos-1, eval.state.get_annotation_label()));
                region_end_pos = eval.position;
                }

            let prev_position = eval.position - 1;
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
        println!("Solution Penalty: {} over {} bp starting at {}", self.eval.penalty, self.hmm.pred_pen.len(), position);

        let regions = self.trace_regions();

        for region in regions.iter()
            {
//            println!("{} to {} is {}", region.start_pos+position+1, region.end_pos+position, region.state.to_str()); // Biologist coordinates

            let mut seq = String::with_capacity(region.end_pos-region.start_pos);
            for idx in region.start_pos .. region.end_pos
                {
                seq.push(self.hmm.bases_pen[idx].as_str());
                }

            println!("{} to {} is {} - {}", region.start_pos+position+1, region.end_pos+position, region.annotation_label.to_str(), seq);

            }

        if position > 30000
            { panic!("First 30k"); }
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













 */
