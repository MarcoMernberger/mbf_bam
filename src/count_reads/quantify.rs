use super::chunked_genome::{Chunk, ChunkedGenome};
use super::OurTree;
use crate::bam_ext::{open_bam, BamRecordExtensions};
use crate::BamError;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;

///Quantify reads into gene_stable_id -> [(position, no-of-reads-with-this.pos)]
///
///
/// Please note one single alignment only counts once for a gene
/// but the same read being aligned multiple times will be counted multiple
/// times.
///
/// this is the high level distributo over chunks of chromosomes and aggregate
/// part
pub fn py_quantify_gene_reads(
    filename: &str,
    index_filename: Option<&str>,
    trees: HashMap<String, (OurTree, Vec<String>)>,
    gene_trees: HashMap<String, (OurTree, Vec<String>)>,
) -> Result<
    (
        HashMap<String, Vec<(u32, u32)>>,
        HashMap<String, Vec<(u32, u32)>>,
    ),
    BamError,
> {
    //check whether the bam file can be openend
    //and we need it for the chunking
    let bam = open_bam(filename, index_filename)?;

    //perform the counting
    let cg = ChunkedGenome::new(gene_trees, bam); // can't get the ParallelBridge to work with our lifetimes.
    let it: Vec<Chunk> = cg.iter().collect();
    let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
    let result = pool.install(|| {
        let result = it
            .into_par_iter()
            .map(|chunk| {
                let bam = open_bam(filename, index_filename).unwrap();
                let (tree, gene_ids) = trees.get(&chunk.chr).unwrap();

                let both_counts = quantify_gene_reads(
                    bam,
                    tree,
                    chunk.tid,
                    chunk.start,
                    chunk.stop,
                    gene_ids.len() as u32,
                );
                let both_counts = both_counts.unwrap_or_else(|_| (Vec::new(), Vec::new()));

                (
                    to_hashmap(both_counts.0, gene_ids),
                    to_hashmap(both_counts.1, gene_ids),
                )
            })
            .reduce(
                || {
                    (
                        HashMap::<String, HashMap<u32, u32>>::new(),
                        HashMap::<String, HashMap<u32, u32>>::new(),
                    )
                },
                add_dual_hashmaps,
            );

        result
    });
    Ok((convert_hashmap(result.0), convert_hashmap(result.1)))
}

fn convert_hashmap(input: HashMap<String, HashMap<u32, u32>>) -> HashMap<String, Vec<(u32, u32)>> {
    input
        .into_iter()
        .map(|(key, value)| (key, value.into_iter().collect()))
        .collect()
}

fn vec_hashmap<K: Hash + Eq, V>(count: usize) -> Vec<HashMap<K, V>> {
    let mut r = Vec::new();
    for _ in 0..count {
        r.push(HashMap::new());
    }
    r
}

/// quantify_gene_reads
//
/// collect two vecs of Hashmaps(position->count)
/// first one is 'congruent with gene direction',
/// second one is 'incongruent with gene direction'

fn quantify_gene_reads(
    mut bam: bam::IndexedReader,
    tree: &OurTree,
    tid: u32,
    chunk_start: u32,
    chunk_stop: u32,
    gene_count: u32,
) -> Result<(Vec<HashMap<u32, u32>>, Vec<HashMap<u32, u32>>), BamError> {
    let mut result_forward: Vec<HashMap<u32, u32>> = vec_hashmap(gene_count as usize);
    let mut result_reverse: Vec<HashMap<u32, u32>> = vec_hashmap(gene_count as usize);

    let mut read: bam::Record = bam::Record::new();
    bam.fetch(tid, chunk_start, chunk_stop)?;
    let mut seen = HashSet::new();
    while let Ok(_) = bam.read(&mut read) {
        if ((read.pos() as u32) < chunk_start) || ((read.pos() as u32) >= chunk_stop) {
            continue;
        } else {
            seen.clear();
            let blocks = read.blocks();
            for iv in blocks.iter() {
                for r in tree.find(iv.0..iv.1) {
                    let entry = r.data();
                    let gene_no = (*entry).0;
                    // do not count multiple blocks matching in one gene multiple times
                    if seen.contains(&gene_no) {
                        continue;
                    }
                    seen.insert(gene_no);
                    let strand = (*entry).1; // this is 1 or -1
                    let congruent = ((strand == 1) && !read.is_reverse())
                        || ((strand != 1) && read.is_reverse());
                    if congruent {
                        let counter = result_forward[gene_no as usize]
                            .entry(read.pos() as u32)
                            .or_insert(0);
                        *counter += 1
                    } else {
                        let counter = result_reverse[gene_no as usize]
                            .entry(read.pos() as u32)
                            .or_insert(0);
                        *counter += 1
                    }
                }
            }
        }
    }
    Ok((result_forward, result_reverse))
}

fn to_hashmap(
    counts: Vec<HashMap<u32, u32>>,
    gene_ids: &Vec<String>,
) -> HashMap<String, HashMap<u32, u32>> {
    counts
        .into_iter()
        .enumerate()
        .map(|(gene_no, cnt)| ((&gene_ids[gene_no]).to_string(), cnt))
        .collect()
    /*
    let mut result = HashMap::new();
    for (gene_no, cnt) in counts.iter().enumerate() {
        let gene_id = &gene_ids[gene_no];
        let key = gene_id.to_string();
        if result.contains_key(key) {
            panic!("Were not disjoint as you thought");
        }
        result.insert(key, cnt);
    }
    result
    */
}
fn add_dual_hashmaps(
    a: (
        HashMap<String, HashMap<u32, u32>>,
        HashMap<String, HashMap<u32, u32>>,
    ),
    b: (
        HashMap<String, HashMap<u32, u32>>,
        HashMap<String, HashMap<u32, u32>>,
    ),
) -> (
    HashMap<String, HashMap<u32, u32>>,
    HashMap<String, HashMap<u32, u32>>,
) {
    (add_hashmaps(a.0, b.0), add_hashmaps(a.1, b.1))
}

fn add_hashmaps(
    mut a: HashMap<String, HashMap<u32, u32>>,
    b: HashMap<String, HashMap<u32, u32>>,
) -> HashMap<String, HashMap<u32, u32>> {
    for (k, v) in b.into_iter() {
        let entry = a.entry(k).or_insert_with(|| HashMap::new());
        for (key, value) in v.into_iter() {
            let o = entry.insert(key, value);
            if o.is_some() {
                panic!("still not disjoint");
            }
        }
    }
    a
}
