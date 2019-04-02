use failure::Error;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::collections::HashMap;
use crate::bam_ext::open_bam;

fn add_hashmaps(mut a: HashMap<u32, u64>, b: HashMap<u32, u64>) -> HashMap<u32, u64> {
    for (k, v) in b.iter() {
        let x = a.entry(*k).or_insert(0);
        *x += v;
    }
    a
}

//used as a
#[derive(Hash, PartialEq, Eq)]
struct ReadAtPos {
    is_reverse: bool,
    //after a long decision, we've decided that
    //it's best to count based on the seq
    //to detect pcr artificats
    // raw_cigar: Vec<u32>
    seq: Vec<u8>,
}

impl std::convert::From<&bam::Record> for ReadAtPos {
    fn from(read: &bam::Record) -> ReadAtPos {
        ReadAtPos {
            is_reverse: read.is_reverse(),
            //raw_cigar: read.raw_cigar().to_vec()
            seq: read.seq().as_bytes(),
        }
    }
}

/// Calculate how often for each k you observe
/// exactly k identically aligned reads (ie. 'repetitions')
/// where identically is defined by read position, strand
/// and cigar string.
///
/// Returns a HashMap {k_reads: number_of_occurences}
/// This is useful to estimate library complexity
pub fn py_calculate_duplicate_distribution(
    filename: &str,
    index_filename: Option<&str>,
) -> Result<HashMap<u32, u64>, Error> {
    let bam = open_bam(filename, index_filename)?;

    let it = 0..bam.header().target_count();
    let result = it
        .into_par_iter()
        .map(|tid| {
            let mut bam2 = open_bam(filename, index_filename).unwrap();
            
            bam2.fetch(tid, 0, bam2.header().target_len(tid).unwrap())
                .unwrap();
            let mut counts: HashMap<u32, u64> = HashMap::new();
            let mut read: bam::Record = bam::Record::new();
            let mut reads_here: HashMap<ReadAtPos, u32> = HashMap::new();
            let mut last_pos = -1;
            while let Ok(_) = bam2.read(&mut read) {
                //*counts.entry(9999999).or_insert(0) += 1; // record total
                if read.pos() != last_pos {
                    if last_pos != -1 {
                        for (_k, v) in reads_here.iter() {
                            *counts.entry(*v).or_insert(0) += 1;
                        }
                        reads_here.clear();
                    }
                    last_pos = read.pos();
                }
                let key: ReadAtPos = ReadAtPos::from(&read); //(read.is_reverse(), read.raw_cigar().to_vec());
                *reads_here.entry(key).or_insert(0) += 1;
            }
            for (_k, v) in reads_here.iter() {
                *counts.entry(*v).or_insert(0) += 1;
            }
            counts
        })
        .reduce(HashMap::<u32, u64>::new, add_hashmaps);
    Ok(result)
}
