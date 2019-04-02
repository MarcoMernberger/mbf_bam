#![feature(nll)]
extern crate pyo3;
extern crate rust_htslib;
#[macro_use]
extern crate failure;
extern crate bio;

//use failure::Error;
use pyo3::prelude::*;
use pyo3::types::{PyDict};
use pyo3::wrap_pyfunction;
use pyo3::{exceptions, PyErr, PyResult};
use std::collections::HashMap;

mod bam_ext;
mod count_reads;
mod duplicate_distribution;

#[derive(Debug, Fail)]
pub enum BamError {
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
/// python wrapper for py_count_reads_unstranded
#[pyfunction]
pub fn count_reads_unstranded(
    filename: &str,
    index_filename: Option<&str>,
    intervals: &PyDict,
) -> PyResult<HashMap<String, u32>> {
    //check whether the bam file can be opend
    //convert the intervals into our interval trees
    let trees: Result<HashMap<String, (count_reads::OurTree, Vec<String>)>, BamError> = intervals
        .iter()
        .map(|(chr, iv_obj)| {
            let chr_str: String = chr.extract()?;
            let (tree, gene_list) = count_reads::build_tree(iv_obj)?;
            Ok((chr_str, (tree, gene_list)))
        })
        .collect();
    let trees = match trees {
        Ok(trees) => trees,
        Err(x) => return Err(x.into()),
    };

    match count_reads::py_count_reads_unstranded(filename, index_filename, trees) {
        Ok(x) => return Ok(x),
        Err(y) => return Err(y.into()),
    };
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
