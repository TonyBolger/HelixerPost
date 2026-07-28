use helixer_post_bin::analysis::extractor::{BasePredictionExtractor, ComparisonExtractor};
use helixer_post_bin::analysis::hmm::show_hmm_config;
use helixer_post_bin::analysis::rater::SequenceRating;
use helixer_post_bin::analysis::Analyzer;
use helixer_post_bin::gff::GffWriter;
use helixer_post_bin::results::raw::RawHelixerPredictions;
use helixer_post_bin::results::HelixerResults;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::exit;
use std::thread;

#[derive(Clone, Copy)]
struct Parameters {
    window_size: usize,
    edge_threshold: f32,
    peak_threshold: f32,
    min_coding_length: usize,
}

struct SequenceOutput {
    sequence_idx: usize,
    gff: Vec<u8>,
    count: usize,
    length: usize,
    fwd_rating: SequenceRating,
    rev_rating: SequenceRating,
}

fn process_sequence(
    genome_path: &Path,
    predictions_path: &Path,
    sequence_idx: usize,
    params: Parameters,
) -> Result<SequenceOutput, String> {
    let helixer_res = HelixerResults::new(predictions_path, genome_path)
        .map_err(|err| format!("failed to open input files: {:?}", err))?;
    let bp_extractor = BasePredictionExtractor::new_from_prediction(&helixer_res)
        .map_err(|err| format!("failed to open prediction datasets: {:?}", err))?;
    let comp_extractor = ComparisonExtractor::new(&helixer_res)
        .map_err(|err| format!("failed to open comparison datasets: {:?}", err))?;
    let analyzer = Analyzer::new(
        bp_extractor,
        comp_extractor,
        params.window_size,
        params.edge_threshold,
        params.peak_threshold,
        params.min_coding_length,
    );
    let seq = helixer_res
        .get_all_sequences()
        .get(sequence_idx)
        .ok_or_else(|| format!("sequence index {} is out of range", sequence_idx))?;
    let species = helixer_res.get_species_by_id(seq.get_species_id());
    let mut fwd_rating = SequenceRating::new();
    let mut rev_rating = SequenceRating::new();
    let mut gff_writer = GffWriter::new(BufWriter::new(Vec::new()));
    gff_writer
        .write_region_header(seq.get_name(), seq.get_length())
        .map_err(|err| format!("failed to write sequence header: {}", err))?;
    let (count, length) = analyzer.process_sequence(
        species,
        seq,
        &mut fwd_rating,
        &mut rev_rating,
        &mut gff_writer,
    );
    let gff = gff_writer
        .into_inner()
        .map_err(|err| format!("failed to finalize sequence output: {}", err))?;
    Ok(SequenceOutput { sequence_idx, gff, count, length, fwd_rating, rev_rating })
}

fn sequence_groups(sequence_count: usize, worker_count: usize) -> Vec<Vec<usize>> {
    let groups = worker_count.min(sequence_count);
    (0..groups)
        .map(|worker| {
            let start = worker * sequence_count / groups;
            let end = (worker + 1) * sequence_count / groups;
            (start..end).collect()
        })
        .collect()
}

fn parse_args() -> Result<(PathBuf, PathBuf, Parameters, PathBuf, usize), String> {
    let args = std::env::args().collect::<Vec<_>>();
    if args.len() != 8 && args.len() != 10 {
        return Err("HelixerPost <genome.h5> <predictions.h5> <windowSize> <edgeThresh> <peakThresh> <minCodingLength> <gff> [--workers N]".into());
    }
    let workers = if args.len() == 10 {
        if args[8] != "--workers" { return Err("optional arguments must be --workers N".into()); }
        args[9].parse::<usize>().map_err(|_| "worker count must be a positive integer".to_string())?
    } else { 1 };
    if workers == 0 { return Err("worker count must be at least one".into()); }
    Ok((
        PathBuf::from(&args[1]), PathBuf::from(&args[2]),
        Parameters {
            window_size: args[3].parse().map_err(|_| "invalid windowSize".to_string())?,
            edge_threshold: args[4].parse().map_err(|_| "invalid edgeThresh".to_string())?,
            peak_threshold: args[5].parse().map_err(|_| "invalid peakThresh".to_string())?,
            min_coding_length: args[6].parse().map_err(|_| "invalid minCodingLength".to_string())?,
        },
        PathBuf::from(&args[7]), workers,
    ))
}

fn main() {
    let (genome_path, predictions_path, params, gff_filename, workers) = match parse_args() {
        Ok(value) => value,
        Err(message) => { eprintln!("{}", message); exit(1); }
    };

    let helixer_res = HelixerResults::new(&predictions_path, &genome_path)
        .expect("Failed to open input files");

    show_hmm_config();

    // There should only ever be one species for the gff output
    assert_eq!(
        helixer_res.get_all_species().len(),
        1,
        "Error: Multiple Species are not allowed for GFF output."
    );
    let rhg = RawHelixerPredictions::new(Path::new(&predictions_path))
        .expect("Error: Something went wrong while accessing the predictions.");
    let model_md5sum = rhg.get_model_md5sum().ok();
    let species_name = helixer_res.get_all_species().first().map(|x| x.get_name());
    // The per-sequence analyzers are worker-local, but reference availability
    // is a property of the shared prediction datasets and is needed for the
    // serial aggregate rating report below.
    let has_ref = ComparisonExtractor::new(&helixer_res)
        .expect("Failed to open comparison datasets")
        .has_ref();

    let mut final_gff = GffWriter::new(BufWriter::new(Vec::new()));
    final_gff
        .write_global_header(species_name, model_md5sum)
        .expect(&*format!(
            "Error: Could not write header to file {}.",
            gff_filename.display()
        ));
    let sequence_count = helixer_res.get_all_sequences().len();
    let mut outputs = Vec::new();
    for group in sequence_groups(sequence_count, workers) {
        let genome_path = genome_path.clone();
        let predictions_path = predictions_path.clone();
        outputs.push(thread::spawn(move || {
            group.into_iter().map(|sequence_idx| {
                process_sequence(&genome_path, &predictions_path, sequence_idx, params)
            }).collect::<Result<Vec<_>, _>>()
        }));
    }
    let mut sequence_outputs = Vec::with_capacity(sequence_count);
    for output in outputs {
        match output.join() {
            Ok(Ok(mut output)) => sequence_outputs.append(&mut output),
            Ok(Err(message)) => { eprintln!("worker failed: {}", message); exit(1); }
            Err(_) => { eprintln!("worker panicked; final output was not written"); exit(1); }
        }
    }
    sequence_outputs.sort_by_key(|output| output.sequence_idx);
    let mut total_count = 0;
    let mut total_length = 0;
    for species in helixer_res.get_all_species() {
        let mut fwd_species_rating = SequenceRating::new();
        let mut rev_species_rating = SequenceRating::new();

        let id = species.get_id();
        println!(
            "Sequences for Species {} - {}",
            species.get_name(),
            id.inner()
        );
        for seq_id in helixer_res.get_sequences_for_species(id) {
            let output = &sequence_outputs[seq_id.inner()];
            final_gff.write_bytes(&output.gff).expect("Failed to buffer GFF output");
            total_count += output.count;
            total_length += output.length;
            fwd_species_rating.accumulate(&output.fwd_rating);
            rev_species_rating.accumulate(&output.rev_rating);
        }

        println!(
            "Forward for Species {} - {}",
            species.get_name(),
            id.inner()
        );
        fwd_species_rating.dump(has_ref);

        println!(
            "Reverse for Species {} - {}",
            species.get_name(),
            id.inner()
        );
        rev_species_rating.dump(has_ref);

        let mut species_rating = SequenceRating::new();
        species_rating.accumulate(&fwd_species_rating);
        species_rating.accumulate(&rev_species_rating);

        println!("Total for Species {} - {}", species.get_name(), id.inner());
        species_rating.dump(has_ref);
    }

    let final_bytes = final_gff.into_inner().expect("Failed to finalize GFF");
    let temporary = gff_filename.with_extension(format!("{}.tmp", gff_filename.extension().and_then(|x| x.to_str()).unwrap_or("gff")));
    {
        let mut file = File::create(&temporary).expect("Failed to create temporary GFF");
        file.write_all(&final_bytes).expect("Failed to write temporary GFF");
        file.flush().expect("Failed to flush temporary GFF");
    }
    std::fs::rename(&temporary, &gff_filename).expect("Failed to commit final GFF");

    println!("Total: {}bp across {} windows", total_length, total_count);
}
