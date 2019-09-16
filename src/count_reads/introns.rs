use super::chunked_genome::{Chunk, ChunkedGenome};
use crate::bam_ext::{open_bam, BamRecordExtensions};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::collections::{HashMap};
use crate::BamError;

type IntronsOnOneChromosome = HashMap<(u32, u32), u32>;
pub type IntronResult = HashMap<String, IntronsOnOneChromosome>;

fn combine_intron_results(mut a: IntronResult, b: IntronResult) -> IntronResult {
    for (key, source) in b.iter() {
        let target = &mut a
            .entry(key.to_string())
            .or_insert(HashMap::<(u32, u32), u32>::new());
        for (coordinates, counts) in source.iter() {
            *target.entry(*coordinates).or_insert(0) += counts;
        }
    }
    a
}

fn count_introns(
    mut bam: bam::IndexedReader,
    tid: u32,
    start: u32,
    stop: u32,
) -> Result<IntronsOnOneChromosome, BamError> {
    let mut result = HashMap::new();
    bam.fetch(tid, start, stop)?;
    let mut read: bam::Record = bam::Record::new();
    while let Ok(_) = bam.read(&mut read) {
        // do not count multiple blocks matching in one gene multiple times
        if ((read.pos() as u32) < start) || ((read.pos() as u32) >= stop) {
            continue;
        }
        let mut intron_count = 0;
        for (start, stop) in read.introns().iter() {
            if stop < start {
                panic!("stop < start")
            }
            *result.entry((*start, *stop)).or_insert(0 as u32) += 1;
            intron_count += 1;
        }
        *result
            .entry((std::u32::MAX, std::u32::MAX - intron_count))
            .or_insert(0 as u32) += 1;
    }
    Ok(result)
}

/// find all introns from a bam
/// result is a {Reference: {(start, stop): count}}
///
pub fn py_count_introns(
    filename: &str,
    index_filename: Option<&str>,
) -> Result<IntronResult, BamError> {
    //check whether the bam file can be openend
    //and we need it for the chunking
    let bam = open_bam(filename, index_filename)?;

    //perform the counting
    let cg = ChunkedGenome::new_without_tree(bam); // can't get the ParallelBridge to work with our lifetimes.
    let it: Vec<Chunk> = cg.iter().collect();
    let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
    pool.install(|| {
        let result = it
            .into_par_iter()
            .map(|chunk| {
                let bam = open_bam(filename, index_filename).unwrap();

                let introns = count_introns(bam, chunk.tid, chunk.start, chunk.stop);
                let mut result: IntronResult = HashMap::new();
                match introns {
                    Ok(introns) => {
                        result.insert(chunk.chr.to_string(), introns);
                    }
                    Err(_) => (),
                };
                result
            })
            .reduce(|| IntronResult::new(), combine_intron_results);
        //.fold(HashMap::<String, u32>::new(), add_hashmaps);
        Ok(result)
    })
}
