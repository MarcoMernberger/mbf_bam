#![feature(nll)]
extern crate pyo3;
extern crate rust_htslib;
#[macro_use]
extern crate failure;
extern crate bio;

use crate::bam_ext::BamRecordExtensions;
use bio::data_structures::interval_tree::IntervalTree;
use failure::Error;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyObjectRef, PyTuple};
use pyo3::wrap_pyfunction;
use pyo3::{exceptions, PyErr, PyResult};
use rayon::iter::ParallelBridge;
use rayon::prelude::ParallelIterator;
use rayon::prelude::*;

use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::collections::{HashMap, HashSet};

mod bam_ext;
mod count_reads;
mod duplicate_distribution;

fn add_hashmaps(mut a: HashMap<String, u32>, b: HashMap<String, u32>) -> HashMap<String, u32> {
    for (k, v) in b.iter() {
        let x = a.entry(k.to_string()).or_insert(0);
        *x += v;
    }
    a
}

#[derive(Debug, Fail)]
enum BamError {
    #[fail(display = "unknown error: {}", msg)]
    UnknownError { msg: String },
}

impl std::convert::Into<PyErr> for BamError {
    fn into(self: BamError) -> PyErr {
        match self {
            BamError::UnknownError { msg } => exceptions::ValueError::py_err(msg),
        }
    }
}

impl std::convert::From<PyErr> for BamError {
    fn from(error: PyErr) -> BamError {
        BamError::UnknownError {
            msg: format!("Python error {:?}", error),
        }
    }
}

impl std::convert::From<rust_htslib::bam::ReaderPathError> for BamError {
    fn from(_error: rust_htslib::bam::ReaderPathError) -> BamError {
        BamError::UnknownError {
            msg: "Could not read bam file".to_string(),
        }
    }
}
impl std::convert::From<rust_htslib::bam::IndexedReaderPathError> for BamError {
    fn from(_error: rust_htslib::bam::IndexedReaderPathError) -> BamError {
        BamError::UnknownError {
            msg: "Could not read bam file".to_string(),
        }
    }
}

impl std::convert::From<rust_htslib::bam::FetchError> for BamError {
    fn from(_error: rust_htslib::bam::FetchError) -> BamError {
        BamError::UnknownError {
            msg: "Fetch error".to_string(),
        }
    }
}

#[pyfunction]
/// calculate_duplicate_distribution(filename, (index_filename) /)
/// --
/// python wrapper for py_calculate_duplicate_distribution
pub fn calculate_duplicate_distribution(
    filename: &str,
    index_filename: Option<&str>,
) -> PyResult<HashMap<u32, u64>> {
    match duplicate_distribution::py_calculate_duplicate_distribution(filename, index_filename) {
        Ok(x) => Ok(x),
        Err(x) => Err(exceptions::ValueError::py_err(format!("{}", x))),
    }
}

type OurTree = IntervalTree<u32, (u32, i8)>;

fn build_tree(iv_obj: &PyObjectRef) -> Result<(OurTree, Vec<String>), PyErr> {
    let iv_list: &PyList = iv_obj.extract()?;
    let mut tree = IntervalTree::new();
    let mut gene_ids = Vec::new();
    let mut last_gene_no = 0;
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
        last_gene_no = gene_no;
    }
    Ok((tree, gene_ids))
}

fn count_reads(
    mut bam: bam::IndexedReader,
    tree: &OurTree,
    tid: u32,
    start: u32,
    stop: u32,
    gene_count: u32,
    gene_ids: &Vec<String>
) -> Result<Vec<u32>, BamError> {
    let mut result = vec![0; gene_count as usize];
    let mut multimapper_dedup: HashMap<u32, HashSet<Vec<u8>>> = HashMap::new();
    let mut gene_nos_seen = HashSet::<u32>::new();
    let mut read: bam::Record = bam::Record::new();
    bam.fetch(tid, start, stop)?;
    while let Ok(_) = bam.read(&mut read) {
        let blocks = read.blocks();
        // do not count multiple blocks matching in one gene multiple times
        gene_nos_seen.clear();
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
                let entry = r.data();
                let gene_no = (*entry).0;
                let nh = read.aux(b"NH");
                let nh = nh.map_or(1, |aux| aux.integer());
                if nh == 1 {
                    gene_nos_seen.insert(gene_no);
                } else {
                    let hs = multimapper_dedup
                        .entry(gene_no)
                        .or_insert_with(|| HashSet::new());
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
            }
        }
        for gene_no in gene_nos_seen.iter() {
            result[*gene_no as usize] += 1;
        }
    }
    for (gene_no, hs) in multimapper_dedup.iter() {
        result[*gene_no as usize] += hs.len() as u32;
    }
    Ok(result)
}

struct ChunkedGenome {
    trees: HashMap<String, (OurTree, Vec<String>)>,
    bam: bam::IndexedReader,
}

impl ChunkedGenome {
    fn new(
        trees: HashMap<String, (OurTree, Vec<String>)>,
        bam: bam::IndexedReader,
    ) -> ChunkedGenome {
        ChunkedGenome { trees, bam }
    }

    fn iter(&self) -> ChunkedGenomeIterator {
        ChunkedGenomeIterator {
            cg: &self,
            it: self.trees.iter(),
            last_start: 0,
            last_tid: 0,
            last_chr_length: 0,
            last_chr: "".to_string(),
        }
    }
}

struct ChunkedGenomeIterator<'a> {
    cg: &'a ChunkedGenome,
    it: std::collections::hash_map::Iter<'a, String, (OurTree, Vec<String>)>,
    last_start: u32,
    last_chr: String,
    last_tid: u32,
    last_chr_length: u32,
}
struct Chunk<'a> {
    chr: String,
    tid: u32,
    tree: &'a OurTree,
    gene_ids: &'a Vec<String>,
    start: u32,
    stop: u32,
}

impl<'a> Iterator for ChunkedGenomeIterator<'a> {
    type Item = Chunk<'a>;
    fn next(&mut self) -> Option<Chunk<'a>> {
        let chunk_size = 1000000;
        if self.last_start >= self.last_chr_length {
            let (next_chr, (next_tree, next_gene_ids)) = match self.it.next() {
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

        let (next_tree, next_gene_ids) = self.cg.trees.get(&self.last_chr).unwrap();
        let mut stop = self.last_start + chunk_size;
        while true {
            //TODO: migh have to extend this for exon based counting to not
            //cut gene in half?
            let overlapping = next_tree.find(stop..stop + 1).next();
            match (overlapping) {
                None => break,
                Some(entry) => {
                    let iv = entry.interval();
                    stop = iv.end + 1;
                }
            }
        }
        let c = Chunk {
            chr: self.last_chr.clone(),
            tid: self.last_tid,
            tree: next_tree,
            gene_ids: next_gene_ids,
            start: self.last_start,
            stop: stop,
        };
        self.last_start = stop;
        Some(c)
    }
}

/// python wrapper for py_count_reads_unstranded
#[pyfunction]
pub fn count_reads_unstranded(
    filename: &str,
    index_filename: Option<&str>,
    intervals: &PyDict,
) -> PyResult<HashMap<String, u32>> {
    //check whether the bam file can be opend
    let bam = match index_filename {
        Some(ifn) => bam::IndexedReader::from_path_and_index(filename, ifn),
        _ => bam::IndexedReader::from_path(filename),
    };
    let bam = match bam {
        Ok(x) => x,
        Err(e) => {
            return Err(BamError::UnknownError {
                msg: format!("Could not read bam: {}", e),
            }
            .into());
        }
    };
    //convert the intervals into our interval trees
    let trees: Result<HashMap<String, (OurTree, Vec<String>)>, BamError> = intervals
        .iter()
        .map(|(chr, iv_obj)| {
            let chr_str: String = chr.extract()?;
            let (tree, gene_list) = build_tree(iv_obj)?;
            Ok((chr_str, (tree, gene_list)))
        })
        .collect();
    let trees = match trees {
        Ok(trees) => trees,
        Err(x) => return Err(x.into()),
    };
    //perform the counting
    let cg = ChunkedGenome::new(trees, bam);
    let it: Vec<Chunk> = cg.iter().collect();
    let result = it
        .into_par_iter()
        .map(|chunk| {
            let bam = match index_filename {
                Some(ifn) => bam::IndexedReader::from_path_and_index(filename, ifn).unwrap(),
                _ => bam::IndexedReader::from_path(filename).unwrap(),
            };

            let counts = count_reads(
                bam,
                &chunk.tree,
                chunk.tid,
                chunk.start,
                chunk.stop,
                chunk.gene_ids.len() as u32,
                &chunk.gene_ids,
            );
            let mut total = 0;
            let mut result: HashMap<String, u32> = match counts {
                Ok(counts) => {
                    let mut res = HashMap::new();
                    for (gene_no, cnt) in counts.iter().enumerate() {
                        let gene_id = &chunk.gene_ids[gene_no];
                        res.insert(gene_id.to_string(), *cnt);
                        total += cnt;
                    }
                    res
                }
                _ => HashMap::new(),
            };
            result.insert("_total".to_string(), total);
            result.insert(format!("_{}", chunk.chr), total);
            result
        })
        .reduce(|| HashMap::<String, u32>::new(), add_hashmaps);
    //.fold(HashMap::<String, u32>::new(), add_hashmaps);
    Ok(result)
}

/// This module is a python module implemented in Rust.
#[pymodule]
fn mbf_bam(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(calculate_duplicate_distribution))?;
    m.add_wrapped(wrap_pyfunction!(count_reads_unstranded))?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
//tests are in the callers until we can actually specify that we need mbf_align (and it's sample
//data) for the testing.
