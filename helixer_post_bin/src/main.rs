
use std::process::exit;
use helixer_post_bin::results::{HelixerResults, Sequence, Species};
use helixer_post_bin::analysis::{BasePredictionExtractor, ComparisonExtractor};
use helixer_post_bin::analysis::window::BasePredictionWindowThresholdIterator;
use helixer_post_bin::analysis::hmm::{PredictionHmm, show_config, HmmStateRegion};
use helixer_post_bin::gff::GffWriter;
use std::fs::File;
use std::io::{BufWriter, Write};
use helixer_post_bin::analysis::gff_conv::hmm_solution_to_gff;
use helixer_post_bin::analysis::rater::{SequenceRater, SequenceRating};


fn process_sequence<W: Write>(bp_extractor: &BasePredictionExtractor, comp_extractor: &ComparisonExtractor,
                              species: &Species, seq: &Sequence, window_size: usize, edge_threshold: f32, peak_threshold: f32, min_coding_length: usize,
                              fwd_rating: &mut SequenceRating, rev_rating: &mut SequenceRating, gff_writer: &mut GffWriter<W>) -> (usize, usize)
{
    let id = seq.get_id();
    println!("  BP_Extractor for Sequence {} - ID {}", seq.get_name(), id.inner());

    let mut window_count = 0;
    let mut window_length_total = 0;

    let fwd_bp_iter =
        BasePredictionWindowThresholdIterator::new(bp_extractor.fwd_iterator(id), window_size, edge_threshold, peak_threshold).unwrap();
    let mut fwd_comp_rater = SequenceRater::new(comp_extractor.fwd_iterator(id), seq.get_length() as usize);

    let mut gene_idx = 1;
    for (bp_vec, _total_vec, start_pos, _peak) in fwd_bp_iter
        {
        window_count+=1;
        window_length_total+=bp_vec.len();

        let end_pos = start_pos + bp_vec.len();

        let hmm = PredictionHmm::new(bp_vec);
        let maybe_solution = hmm.solve();

        if let Some(solution) = maybe_solution
            {
            //solution.dump(start_pos);
            let solution_regions = solution.trace_regions();
            let genes = HmmStateRegion::split_genes(solution_regions);

            for (gene_regions, coding_length) in genes.iter()
                { fwd_comp_rater.rate_regions(start_pos, &gene_regions, *coding_length < min_coding_length); }

            let gff_records = hmm_solution_to_gff(genes, species.get_name(), seq.get_name(), "HelixerPost",
                                              false, start_pos, seq.get_length(), min_coding_length, &mut gene_idx);
            gff_writer.write_records(&gff_records).expect("Failed to write to GFF");
            }
        else
            { panic!("No solution at {} {} - {}", seq.get_name(), start_pos, end_pos); }
        }

    let fwd_seq_rating = fwd_comp_rater.calculate_stats();
    println!("Forward for Sequence {} - ID {}", seq.get_name(), id.inner());
    fwd_seq_rating.dump();

    fwd_rating.accumulate(&fwd_seq_rating);

    let rev_bp_iter =
        BasePredictionWindowThresholdIterator::new(bp_extractor.rev_iterator(id), window_size, edge_threshold, peak_threshold).unwrap();
    let mut rev_comp_rater = SequenceRater::new(comp_extractor.rev_iterator(id), seq.get_length() as usize);

    for (bp_vec, _total_vec, start_pos, _peak) in rev_bp_iter
        {
        window_count+=1;
        window_length_total+=bp_vec.len();

        let end_pos = start_pos + bp_vec.len();
//        rev_comp_rater.rate_window(start_pos, end_pos);

        let hmm = PredictionHmm::new(bp_vec);
        let maybe_solution = hmm.solve();

        if let Some(solution) = maybe_solution
            {
            let solution_regions = solution.trace_regions();
            let genes = HmmStateRegion::split_genes(solution_regions);

            for (gene_regions,coding_length) in genes.iter()
                { rev_comp_rater.rate_regions(start_pos, &gene_regions, *coding_length < min_coding_length); }


            let gff_records = hmm_solution_to_gff(genes, species.get_name(), seq.get_name(), "HelixerPost",
                                                  true, start_pos, seq.get_length(), min_coding_length, &mut gene_idx);
            gff_writer.write_records(&gff_records).expect("Failed to write to GFF");
            }
         else
            { panic!("No solution at {} {} - {}", seq.get_name(), start_pos, end_pos); }

        }

    let rev_seq_rating = rev_comp_rater.calculate_stats();
    println!("Reverse for Sequence {} - ID {}", seq.get_name(), id.inner());
    rev_seq_rating.dump();

    rev_rating.accumulate(&rev_seq_rating);

    (window_count, window_length_total)
}


fn main()
{
    let arg_vec = std::env::args().collect::<Vec<_>>(); // Arg iterator into vector

    if arg_vec.len() != 8
    {
        println!("HelixerPost <genome.h5> <predictions.h5> <windowSize> <edgeThresh> <peakThresh> <minCodingLength> <gff>");
        exit(1);
    }

    let genome_path = arg_vec[1].as_str();
    let predictions_path = arg_vec[2].as_str();
    let window_size = arg_vec[3].parse().unwrap();
    let edge_threshold = arg_vec[4].parse().unwrap();
    let peak_threshold = arg_vec[5].parse().unwrap();
    let min_coding_length = arg_vec[6].parse().unwrap();
    let gff_filename = &arg_vec[7];

    let helixer_res = HelixerResults::new(predictions_path.as_ref(), genome_path.as_ref()).expect("Failed to open input files");

    let bp_extractor = BasePredictionExtractor::new(&helixer_res).expect("Failed to open Base / ClassPrediction / PhasePrediction Datasets");
    let comp_extractor = ComparisonExtractor::new(&helixer_res).expect("Failed to open ClassReference / PhaseReference / ClassPrediction / PhasePrediction Datasets");

    let mut total_count = 0;
    let mut total_length = 0;

    show_config();

    let gff_file = File::create(gff_filename).unwrap();
    let mut gff_writer = GffWriter::new(BufWriter::new(gff_file));

    for species in helixer_res.get_all_species()
    {
        let mut fwd_species_rating = SequenceRating::new();
        let mut rev_species_rating = SequenceRating::new();

        let id = species.get_id();
        println!("Sequences for Species {} - {}", species.get_name(), id.inner() );
        for seq_id in helixer_res.get_sequences_for_species(id)
        {
            let seq = helixer_res.get_sequence_by_id(*seq_id);

//            dump_seq(&helixer_res, seq);
            let (count, length) = process_sequence(&bp_extractor, &comp_extractor, species, seq,
                                                   window_size, edge_threshold, peak_threshold, min_coding_length,
                                                   &mut fwd_species_rating, &mut rev_species_rating, &mut gff_writer);

            total_count+=count;
            total_length+=length;
        }

        println!("Forward for Species {} - {}", species.get_name(), id.inner() );
        fwd_species_rating.dump();

        println!("Reverse for Species {} - {}", species.get_name(), id.inner() );
        rev_species_rating.dump();

        let mut species_rating = SequenceRating::new();
        species_rating.accumulate(&fwd_species_rating);
        species_rating.accumulate(&rev_species_rating);

        println!("Total for Species {} - {}", species.get_name(), id.inner() );
        species_rating.dump();
    }

    println!("Total: {}bp across {} windows", total_length, total_count);


}