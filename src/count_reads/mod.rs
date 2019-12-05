mod chunked_genome;
mod counters;
mod introns;
mod quantify;
pub mod by_barcode;

pub use counters::{py_count_reads_stranded, py_count_reads_unstranded};
pub use introns::{py_count_introns, IntronResult};
pub use quantify::py_quantify_gene_reads;

use bio::data_structures::interval_tree::IntervalTree;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyList, PyTuple};
use std::collections::HashMap;

// OurTree stores an interval tree
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

fn add_hashmaps(mut a: HashMap<String, u32>, b: HashMap<String, u32>) -> HashMap<String, u32> {
    for (k, v) in b.iter() {
        let x = a.entry(k.to_string()).or_insert(0);
        *x += v;
    }
    a
}
