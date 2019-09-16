use super::chunked_genome::{Chunk, ChunkedGenome};
use super::{add_hashmaps, OurTree};
use crate::bam_ext::{open_bam, BamRecordExtensions};
use crate::BamError;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::collections::{HashMap, HashSet};

fn add_dual_hashmaps(
    a: (HashMap<String, u32>, HashMap<String, u32>),
    b: (HashMap<String, u32>, HashMap<String, u32>),
) -> (HashMap<String, u32>, HashMap<String, u32>) {
    (add_hashmaps(a.0, b.0), add_hashmaps(a.1, b.1))
}

///count_reads_in_region_unstranded
///
///counts the unstranded reads in a region,
///matching the to the tree entries.
///
///if each_read_counts_once is set, each read can count only for one gene - the first
///one where a block hits
///otherwise reads having blocks that can not be assigned to just one gene count for
///all they're hitting
fn count_reads_in_region_unstranded(
    mut bam: bam::IndexedReader,
    tree: &OurTree,
    tid: u32,
    start: u32,
    stop: u32,
    gene_count: u32,
    each_read_counts_once: bool,
) -> Result<(Vec<u32>, u32), BamError> {
    let mut result = vec![0; gene_count as usize];
    let mut multimapper_dedup: HashMap<u32, HashSet<Vec<u8>>> = HashMap::new();
    let mut gene_nos_seen = HashSet::<u32>::new();
    let mut outside_count = 0;
    let mut read: bam::Record = bam::Record::new();
    bam.fetch(tid, start, stop)?;
    while let Ok(_) = bam.read(&mut read) {
        // do not count multiple blocks matching in one gene multiple times
        gene_nos_seen.clear();
        let mut hit = false;
        let mut skipped = false;
        if ((read.pos() as u32) < start) || ((read.pos() as u32) >= stop) {
            skipped = true;
        }
        if !skipped {
            let blocks = read.blocks();
            for iv in blocks.iter() {
                if (iv.1 < start) || iv.0 >= stop || ((iv.0 < start) && (iv.1 >= start)) {
                    // if this block is outside of the region
                    // don't count it at all.
                    // if it is on a block boundary
                    // only count it for the left side.
                    // which is ok, since we place the blocks to the right
                    // of our intervals.
                    continue;
                }
                for r in tree.find(iv.0..iv.1) {
                    hit = true;
                    let entry = r.data();
                    let gene_no = (*entry).0;
                    let nh = read.aux(b"NH");
                    let nh = nh.map_or(1, |aux| aux.integer());
                    if nh == 1 {
                        gene_nos_seen.insert(gene_no);
                    } else {
                        let hs = multimapper_dedup
                            .entry(gene_no)
                            .or_insert_with(HashSet::new);
                        hs.insert(read.qname().to_vec());
                    }
                    /*if gene_ids[gene_no as usize] == "FBgn0037275" {
                    println!(
                        "{}, {}, {}",
                        start,
                        stop,
                        std::str::from_utf8(read.qname()).unwrap()
                    );
                    }*/
                    if each_read_counts_once {
                        break; // enable this (part 1 of 2) for each read hitting only once
                    }
                }
                if each_read_counts_once {
                    if hit {
                        //enable this (part 2 of 2) for each read hitting only once
                        break;
                    }
                }
            }
        }
        if !hit && !skipped {
            outside_count += 1;
        }
        for gene_no in gene_nos_seen.iter() {
            result[*gene_no as usize] += 1;
        }
    }
    for (gene_no, hs) in multimapper_dedup.iter() {
        result[*gene_no as usize] += hs.len() as u32;
    }
    Ok((result, outside_count))
}
///
/// count_reads_in_region_stranded
//
/// counts the stranded reads in a region,
/// matching them to the tree entries.
/// returns two vectors to be translated from gene_no
/// to gene_id: matching, reverse matching
fn count_reads_in_region_stranded(
    mut bam: bam::IndexedReader,
    tree: &OurTree,
    tid: u32,
    start: u32,
    stop: u32,
    gene_count: u32,
    each_read_counts_once: bool,
) -> Result<(Vec<u32>, Vec<u32>, u32), BamError> {
    let mut result_forward = vec![0; gene_count as usize];
    let mut result_reverse = vec![0; gene_count as usize];
    let mut multimapper_dedup_forward: HashMap<u32, HashSet<Vec<u8>>> = HashMap::new();
    let mut multimapper_dedup_reverse: HashMap<u32, HashSet<Vec<u8>>> = HashMap::new();
    let mut gene_nos_seen_forward = HashSet::<u32>::new();
    let mut gene_nos_seen_reverse = HashSet::<u32>::new();
    let mut outside_count = 0;
    let mut read: bam::Record = bam::Record::new();
    bam.fetch(tid, start, stop)?;
    while let Ok(_) = bam.read(&mut read) {
        // do not count multiple blocks matching in one gene multiple times
        gene_nos_seen_forward.clear();
        gene_nos_seen_reverse.clear();
        let mut hit = false;
        let mut skipped = false;
        if ((read.pos() as u32) < start) || ((read.pos() as u32) >= stop) {
            skipped = true;
        }
        if !skipped {
            let blocks = read.blocks();
            for iv in blocks.iter() {
                for r in tree.find(iv.0..iv.1) {
                    hit = true;
                    let entry = r.data();
                    let gene_no = (*entry).0;
                    let strand = (*entry).1; // this is 1 or -1
                    let nh = read.aux(b"NH");
                    let nh = nh.map_or(1, |aux| aux.integer());
                    if ((strand == 1) && !read.is_reverse()) || ((strand != 1) && read.is_reverse())
                    {
                        // read is in correct orientation
                        if nh == 1 {
                            gene_nos_seen_forward.insert(gene_no);
                        } else {
                            let hs = multimapper_dedup_forward
                                .entry(gene_no)
                                .or_insert_with(HashSet::new);
                            hs.insert(read.qname().to_vec());
                        }
                    } else {
                        // read is in inverse orientation
                        if nh == 1 {
                            gene_nos_seen_reverse.insert(gene_no);
                        } else {
                            let hs = multimapper_dedup_reverse
                                .entry(gene_no)
                                .or_insert_with(HashSet::new);
                            hs.insert(read.qname().to_vec());
                        }
                    }
                    if each_read_counts_once {
                        break; // enable this (part 1 of 2) for each read hitting only once
                    }
                }
                if each_read_counts_once {
                    if hit {
                        //enable this (part 2 of 2) for each read hitting only once
                        break;
                    }
                }
            }
        }
        if !hit && !skipped {
            outside_count += 1;
        }
        for gene_no in gene_nos_seen_forward.iter() {
            result_forward[*gene_no as usize] += 1;
        }
        for gene_no in gene_nos_seen_reverse.iter() {
            result_reverse[*gene_no as usize] += 1;
        }
    }
    for (gene_no, hs) in multimapper_dedup_forward.iter() {
        result_forward[*gene_no as usize] += hs.len() as u32;
    }
    for (gene_no, hs) in multimapper_dedup_reverse.iter() {
        result_reverse[*gene_no as usize] += hs.len() as u32;
    }
    Ok((result_forward, result_reverse, outside_count))
}
/// python wrapper for py_count_reads_unstranded
pub fn py_count_reads_unstranded(
    filename: &str,
    index_filename: Option<&str>,
    trees: HashMap<String, (OurTree, Vec<String>)>,
    gene_trees: HashMap<String, (OurTree, Vec<String>)>,
    each_read_counts_once: bool,
) -> Result<HashMap<String, u32>, BamError> {
    //check whether the bam file can be openend
    //and we need it for the chunking
    let bam = open_bam(filename, index_filename)?;

    //perform the counting
    let cg = ChunkedGenome::new(gene_trees, bam); // can't get the ParallelBridge to work with our lifetimes.
    let it: Vec<Chunk> = cg.iter().collect();
    let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
    pool.install(|| {
        let result = it
            .into_par_iter()
            .map(|chunk| {
                let bam = open_bam(filename, index_filename).unwrap();
                let (tree, gene_ids) = trees.get(&chunk.chr).unwrap();

                let counts = count_reads_in_region_unstranded(
                    bam,
                    //&chunk.tree,
                    tree,
                    chunk.tid,
                    chunk.start,
                    chunk.stop,
                    gene_ids.len() as u32,
                    each_read_counts_once,
                );
                let mut total = 0;
                let mut outside = 0;
                let mut result: HashMap<String, u32> = match counts {
                    Ok(counts) => {
                        let mut res = HashMap::new();
                        for (gene_no, cnt) in counts.0.iter().enumerate() {
                            let gene_id = &gene_ids[gene_no];
                            res.insert(gene_id.to_string(), *cnt);
                            total += cnt;
                        }
                        outside += counts.1;
                        res
                    }
                    _ => HashMap::new(),
                };
                result.insert("_total".to_string(), total);
                result.insert("_outside".to_string(), outside);
                result.insert(format!("_{}", chunk.chr), total);
                result
            })
            .reduce(HashMap::<String, u32>::new, add_hashmaps);
        //.fold(HashMap::<String, u32>::new(), add_hashmaps);
        Ok(result)
    })
}

fn to_hashmap(counts: Vec<u32>, gene_ids: &Vec<String>, chr: &str) -> HashMap<String, u32> {
    let mut total = 0;
    let mut result = HashMap::new();
    for (gene_no, cnt) in counts.iter().enumerate() {
        let gene_id = &gene_ids[gene_no];
        *result.entry(gene_id.to_string()).or_insert(0) += *cnt;
        total += cnt;
    }
    result.insert("_total".to_string(), total);
    result.insert(format!("_{}", chr), total);
    result
}

/// python wrapper for py_count_reads_stranded
pub fn py_count_reads_stranded(
    filename: &str,
    index_filename: Option<&str>,
    trees: HashMap<String, (OurTree, Vec<String>)>,
    gene_trees: HashMap<String, (OurTree, Vec<String>)>,
    each_read_counts_once: bool,
) -> Result<(HashMap<String, u32>, HashMap<String, u32>), BamError> {
    //check whether the bam file can be openend
    //and we need it for the chunking
    let bam = open_bam(filename, index_filename)?;

    //perform the counting
    let cg = ChunkedGenome::new(gene_trees, bam); // can't get the ParallelBridge to work with our lifetimes.
    let it: Vec<Chunk> = cg.iter().collect();
    let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
    pool.install(|| {
        let result = it
            .into_par_iter()
            .map(|chunk| {
                let bam = open_bam(filename, index_filename).unwrap();
                let (tree, gene_ids) = trees.get(&chunk.chr).unwrap();

                let both_counts = count_reads_in_region_stranded(
                    bam,
                    tree,
                    chunk.tid,
                    chunk.start,
                    chunk.stop,
                    gene_ids.len() as u32,
                    each_read_counts_once,
                );
                let both_counts = both_counts.unwrap_or_else(|_| (Vec::new(), Vec::new(), 0));

                let mut result = (
                    to_hashmap(both_counts.0, gene_ids, &chunk.chr),
                    to_hashmap(both_counts.1, gene_ids, &chunk.chr),
                );
                result.0.insert("_outside".to_string(), both_counts.2);
                result
            })
            .reduce(
                || (HashMap::<String, u32>::new(), HashMap::<String, u32>::new()),
                add_dual_hashmaps,
            );
        //.fold(HashMap::<String, u32>::new(), add_hashmaps);
        Ok(result)
    })
}
