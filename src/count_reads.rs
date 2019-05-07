use bio::data_structures::interval_tree::IntervalTree;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyList, PyTuple};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::prelude::*;
/// Read counting functions
///
/// These are smart, fast read counters
/// that avoid the more common pitfalls of htseq data
/// in terms of handling reads matching multiple genes
/// or being multi mapped
///
use std::collections::{HashMap, HashSet};
use std::str;

use crate::bam_ext::{open_bam, BamRecordExtensions};
use crate::BamError;

fn add_hashmaps(mut a: HashMap<String, u32>, b: HashMap<String, u32>) -> HashMap<String, u32> {
    for (k, v) in b.iter() {
        let x = a.entry(k.to_string()).or_insert(0);
        *x += v;
    }
    a
}
fn add_dual_hashmaps(
    a: (HashMap<String, u32>, HashMap<String, u32>),
    b: (HashMap<String, u32>, HashMap<String, u32>),
) -> (HashMap<String, u32>, HashMap<String, u32>) {
    (add_hashmaps(a.0, b.0), add_hashmaps(a.1, b.1))
}

/// OurTree stores an interval tree
/// a gene_no (ie. an index into a vector of gene_ids)
/// and the strand (+1/ -1, 0)
pub type OurTree = IntervalTree<u32, (u32, i8)>;

/// build_tree converts an list of python intervals
/// into an OurTree and a Vec of gene_ids
///
/// python intervals are a list of tuples
/// (gene_stable_id, strand, start, stop)
/// - each reference sequence has it's own list.
pub fn build_tree(iv_obj: &PyAny) -> Result<(OurTree, Vec<String>), PyErr> {
    let iv_list: &PyList = iv_obj.extract()?;
    let mut tree = IntervalTree::new();
    let mut gene_ids = Vec::new();
    for (gene_no, iv_entry_obj) in iv_list.iter().enumerate() {
        let iv_tuple: &PyTuple = iv_entry_obj.extract()?;
        let lgene_id: String = iv_tuple.get_item(0).extract()?;
        gene_ids.push(lgene_id);
        let lstrand: i8 = iv_tuple.get_item(1).extract()?;
        let lstart: &PyList = iv_tuple.get_item(2).extract()?;
        let lend: &PyList = iv_tuple.get_item(3).extract()?;
        for (ls, le) in lstart.iter().zip(lend.iter()) {
            let ls: u32 = ls.extract()?;
            let le: u32 = le.extract()?;
            tree.insert(ls..le, (gene_no as u32, lstrand))
        }
    }
    Ok((tree, gene_ids))
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
/// counts the unstranded reads in a region,
/// matching the to the tree entries.
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

struct ChunkedGenome {
    trees: Option<HashMap<String, (OurTree, Vec<String>)>>,
    bam: bam::IndexedReader,
    chromosomes: Vec<String>,
}

impl ChunkedGenome {
    ///create a new chunked genome for iteration
    ///if you pass in a tree, it is guranteed that the splits happen
    ///between entries of the tree, not inside.
    fn new(
        trees: HashMap<String, (OurTree, Vec<String>)>,
        bam: bam::IndexedReader,
    ) -> ChunkedGenome {
        let chrs_in_tree_and_bam = trees
            .keys()
            .map(|x| x.clone())
            .filter(|x| bam.header().tid(x.as_bytes()).is_some())
            .collect();
        ChunkedGenome {
            chromosomes: chrs_in_tree_and_bam,
            trees: Some(trees),
            bam,
        }
    }
    fn new_without_tree(bam: bam::IndexedReader) -> ChunkedGenome {
        ChunkedGenome {
            trees: None,
            chromosomes: bam
                .header()
                .target_names()
                .iter()
                .map(|x| str::from_utf8(x).unwrap().to_string())
                .collect(),
            bam,
        }
    }
    fn iter(&self) -> ChunkedGenomeIterator {
        ChunkedGenomeIterator {
            cg: &self,
            it: self.chromosomes.iter(),
            last_start: 0,
            last_tid: 0,
            last_chr_length: 0,
            last_chr: "".to_string(),
        }
    }
}

struct ChunkedGenomeIterator<'a> {
    cg: &'a ChunkedGenome,
    it: std::slice::Iter<'a, String>,
    last_start: u32,
    last_chr: String,
    last_tid: u32,
    last_chr_length: u32,
}
#[derive(Debug)]
struct Chunk {
    chr: String,
    tid: u32,
    start: u32,
    stop: u32,
}

impl<'a> Iterator for ChunkedGenomeIterator<'a> {
    type Item = Chunk;
    fn next(&mut self) -> Option<Chunk> {
        let chunk_size = 1_000_000;
        if self.last_start >= self.last_chr_length {
            let next_chr = match self.it.next() {
                Some(x) => x,
                None => return None,
            };
            let tid = self.cg.bam.header().tid(next_chr.as_bytes()).unwrap();
            let chr_length = self.cg.bam.header().target_len(tid).unwrap();
            self.last_tid = tid;
            self.last_chr_length = chr_length;
            self.last_chr = next_chr.to_string();
            self.last_start = 0;
        }

        let mut stop = self.last_start + chunk_size;
        if self.cg.trees.is_some() {
            let (next_tree, _next_gene_ids) =
                self.cg.trees.as_ref().unwrap().get(&self.last_chr).unwrap();
            loop {
                ////this has been adjusted not to cut genes in half
                //cut gene in half?
                //option 0 for that is to pass in the gene intervals as well
                //just for constructing the chunks
                let overlapping = next_tree.find(stop..stop + 1).next();
                match overlapping {
                    None => break,
                    Some(entry) => {
                        let iv = entry.interval();
                        if iv.end + 1 < stop {
                            panic!("WHAT?");
                        }
                        stop = iv.end + 1;
                    }
                }
            }
        }
        let c = Chunk {
            chr: self.last_chr.clone(),
            tid: self.last_tid,
            start: self.last_start,
            stop,
        };
        self.last_start = stop;
        Some(c)
    }
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
