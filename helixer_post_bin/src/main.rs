
use std::process::exit;
use helixer_post_bin::results::{HelixerResults, Sequence, Species};
use helixer_post_bin::analysis::BasePredictionExtractor;
use helixer_post_bin::analysis::window::BasePredictionWindowThresholdIterator;
use helixer_post_bin::analysis::hmm::PredictionHmm;
use helixer_post_bin::gff::GffWriter;
use std::fs::File;
use std::io::{BufWriter, Write};

/*
fn dump_seq(helixer_res: &HelixerResults, seq: &Sequence)
{
    let id = seq.get_id();
    println!("  Data for Sequence {} - ID {}", seq.get_name(), id.inner());

    let bases_blocked_dataset = helixer_res.get_x().unwrap();
    let bases_fwd = bases_blocked_dataset.fwd_iter(id).count();
    let bases_rev = bases_blocked_dataset.rev_iter(id).count();
    println!("      Got {} / {} bases", bases_fwd, bases_rev);

    let gl_blocked_dataset = helixer_res.get_gene_lengths().unwrap();
    let gl_fwd = gl_blocked_dataset.fwd_iter(id).count();
    let gl_rev = gl_blocked_dataset.fwd_iter(id).count();
    println!("      Got {} / {} gene lengths", gl_fwd, gl_rev);

    let sw_blocked_dataset = helixer_res.get_sample_weights().unwrap();
    let sw_fwd = sw_blocked_dataset.fwd_iter(id).count();
    let sw_rev = sw_blocked_dataset.rev_iter(id).count();
    println!("      Got {} / {} sample weights", sw_fwd, sw_rev);

    let anno_blocked_dataset = helixer_res.get_y().unwrap();
    let anno_fwd = anno_blocked_dataset.fwd_iter(id).count();
    let anno_rev = anno_blocked_dataset.rev_iter(id).count();
    println!("      Got {} / {} annos", anno_fwd, anno_rev);

    let trans_blocked_dataset = helixer_res.get_transitions().unwrap();
    let trans_fwd = trans_blocked_dataset.fwd_iter(id).count();
    let trans_rev = trans_blocked_dataset.rev_iter(id).count();
    println!("      Got {} / {} trans", trans_fwd, trans_rev);

    let pred_blocked_dataset = helixer_res.get_predictions().unwrap();
    let pred_fwd = pred_blocked_dataset.fwd_iter(id).count();
    let pred_rev = pred_blocked_dataset.fwd_iter(id).count();
    println!("      Got {} / {} preds", pred_fwd, pred_rev);
}
*/

fn dump_extractor<W: Write>(extractor: &BasePredictionExtractor,
                            species: &Species, seq: &Sequence, window_size: usize, edge_threshold: f32, peak_threshold: f32, gff_writer: &mut GffWriter<W>) -> (usize, usize)
{
    let id = seq.get_id();
    println!("  BP_Extractor for Sequence {} - ID {}", seq.get_name(), id.inner());

    let mut window_count = 0;
    let mut window_length_total = 0;

    let fwd_iter =
        BasePredictionWindowThresholdIterator::new(extractor.fwd_iterator(id), window_size, edge_threshold, peak_threshold).unwrap();

    let mut gene_idx = 1;
    for (bp_vec, _total_vec, start_pos, _peak) in fwd_iter
        {
        window_count+=1;
        window_length_total+=bp_vec.len();

        let end_pos = start_pos + bp_vec.len();

        let hmm = PredictionHmm::new(bp_vec);
        let maybe_solution = hmm.solve();

        if let Some(solution) = maybe_solution
            {
            let gff_records = solution.as_gff(species.get_name(), seq.get_name(), "HelixerPost",
                                              false, start_pos, seq.get_length(), &mut gene_idx);
            gff_writer.write_records(&gff_records).expect("Failed to write to GFF");
            }
        else
            { panic!("No solution at {} {} - {}", seq.get_name(), start_pos, end_pos); }
        }

    let rev_iter =
        BasePredictionWindowThresholdIterator::new(extractor.rev_iterator(id), window_size, edge_threshold, peak_threshold).unwrap();

    for (bp_vec, _total_vec, start_pos, _peak) in rev_iter
        {
        window_count+=1;
        window_length_total+=bp_vec.len();

        let end_pos = start_pos + bp_vec.len();

        let hmm = PredictionHmm::new(bp_vec);
        let maybe_solution = hmm.solve();

        if let Some(solution) = maybe_solution
            {
            let gff_records = solution.as_gff(species.get_name(), seq.get_name(), "HelixerPost",
                                                  true, start_pos, seq.get_length(), &mut gene_idx);
            gff_writer.write_records(&gff_records).expect("Failed to write to GFF");
            }
         else
            { panic!("No solution at {} {} - {}", seq.get_name(), start_pos, end_pos); }

        }

    (window_count, window_length_total)
}


fn main()
{
    let arg_vec = std::env::args().collect::<Vec<_>>(); // Arg iterator into vector

    if arg_vec.len() != 7
    {
        println!("HelixerPost <genome.h5> <predictions.h5> <windowSize> <edgeThresh> <peakThresh> <gff>");
        exit(1);
    }

    let genome_path = arg_vec[1].as_str();
    let predictions_path = arg_vec[2].as_str();
    let window_size = arg_vec[3].parse().unwrap();
    let edge_threshold = arg_vec[4].parse().unwrap();
    let peak_threshold = arg_vec[5].parse().unwrap();
    let gff_filename = &arg_vec[6];

    let helixer_res = HelixerResults::new(predictions_path.as_ref(), genome_path.as_ref()).unwrap();

    let is_err_sample = helixer_res.get_raw_genome().get_err_samples().unwrap();
    let is_fully_ig = helixer_res.get_raw_genome().get_fully_intergenic_samples().unwrap();
    let is_anno = helixer_res.get_raw_genome().get_is_annotated().unwrap();

    println!("Is err_sample - T {} vs F {} ", is_err_sample.iter().filter(|x| **x).count(), is_err_sample.iter().filter(|x| !**x).count());
    println!("Is fully_ig - T {} vs F {} ", is_fully_ig.iter().filter(|x| **x).count(), is_fully_ig.iter().filter(|x| !**x).count());
    println!("Is anno - T {} vs F {} ", is_anno.iter().filter(|x| **x).count(), is_anno.iter().filter(|x| !**x).count());

    let extractor = BasePredictionExtractor::new(&helixer_res).unwrap();

    let mut total_count = 0;
    let mut total_length = 0;

    let gff_file = File::create(gff_filename).unwrap();
    let mut gff_writer = GffWriter::new(BufWriter::new(gff_file));

    for species in helixer_res.get_all_species()
    {
        let id = species.get_id();
        println!("Sequences for Species {} - {}", species.get_name(), id.inner() );
        for seq_id in helixer_res.get_sequences_for_species(id)
        {
            let seq = helixer_res.get_sequence_by_id(*seq_id);

//            dump_seq(&helixer_res, seq);
            let (count, length) = dump_extractor(&extractor, species, seq, window_size, edge_threshold, peak_threshold, &mut gff_writer);

            total_count+=count;
            total_length+=length;
        }
    }

    println!("Total: {}bp in {} sequences", total_length, total_count);


}