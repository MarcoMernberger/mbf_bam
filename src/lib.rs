#![feature(nll)]
extern crate pyo3;
extern crate rust_htslib;
#[macro_use]
extern crate failure;
extern crate bio;

//use failure::Error;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::wrap_pyfunction;
use pyo3::{exceptions, PyErr, PyResult};
use std::collections::HashMap;

mod bam_ext;
mod bam_manipulation;
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

impl std::convert::From<rust_htslib::bam::WriterPathError> for BamError {
    fn from(_error: rust_htslib::bam::WriterPathError) -> BamError {
        BamError::UnknownError {
            msg: "Could not read file".to_string(),
        }
    }
}

impl std::convert::From<rust_htslib::bam::WriteError> for BamError {
    fn from(_error: rust_htslib::bam::WriteError) -> BamError {
        BamError::UnknownError {
            msg: "WriteError".to_string(),
        }
    }
}

impl std::convert::From<rust_htslib::bam::index::IndexBuildError> for BamError {
    fn from(_error: rust_htslib::bam::index::IndexBuildError) -> BamError {
        BamError::UnknownError {
            msg: "WriteError".to_string(),
        }
    }
}

impl std::convert::From<std::io::Error> for BamError {
    fn from(error: std::io::Error) -> BamError {
        BamError::UnknownError {
            msg: format!("IO Error: {:?}", error).to_string(),
        }
    }
}
impl std::convert::From<rust_htslib::bam::AuxWriteError> for BamError {
    fn from(_error: rust_htslib::bam::AuxWriteError) -> BamError {
        BamError::UnknownError {
            msg: "AuxWriteError".to_string(),
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
// /convert the intervals into our interval trees
fn py_intervals_to_trees(
    intervals: &PyDict,
) -> PyResult<HashMap<String, (count_reads::OurTree, Vec<String>)>> {
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
    Ok(trees)
}

/// python wrapper for py_count_reads_unstranded
#[pyfunction]
pub fn count_reads_unstranded(
    filename: &str,
    index_filename: Option<&str>,
    intervals: &PyDict,
    gene_intervals: &PyDict,
    each_read_counts_once: Option<bool>,
) -> PyResult<HashMap<String, u32>> {
    let trees = py_intervals_to_trees(intervals)?;
    let gene_trees = py_intervals_to_trees(gene_intervals)?;
    match count_reads::py_count_reads_unstranded(
        filename,
        index_filename,
        trees,
        gene_trees,
        each_read_counts_once.unwrap_or(false),
    ) {
        Ok(x) => Ok(x),
        Err(y) => Err(y.into()),
    }
}

/// python wrapper for py_count_reads_stranded
#[pyfunction]
pub fn count_reads_stranded(
    filename: &str,
    index_filename: Option<&str>,
    intervals: &PyDict,
    gene_intervals: &PyDict,
    each_read_counts_once: Option<bool>,
) -> PyResult<(HashMap<String, u32>, HashMap<String, u32>)> {
    let trees = py_intervals_to_trees(intervals)?;
    let gene_trees = py_intervals_to_trees(gene_intervals)?;
    let res = match count_reads::py_count_reads_stranded(
        filename,
        index_filename,
        trees,
        gene_trees,
        each_read_counts_once.unwrap_or(false),
    ) {
        Ok(x) => x,
        Err(y) => return Err(y.into()),
    };
    Ok(res)
}
// python wrapper for py_count_reads_unstranded
#[pyfunction]
pub fn count_reads_primary_only_right_strand_only_by_barcode(
    filename: &str,
    index_filename: Option<&str>,
    intervals: &PyDict,
    gene_intervals: &PyDict,
    umi_strategy: String,
) -> PyResult<(
    Vec<String>,
    Vec<String>,
    Vec<(u32, u32, u32)>
    )>{
    let trees = py_intervals_to_trees(intervals)?;
    let gene_trees = py_intervals_to_trees(gene_intervals)?;
    let umi_strategy = match umi_strategy.as_ref() {
        "straight" => count_reads::by_barcode::UmiStrategy::Straight,
        _ => return Err(BamError::UnknownError{msg: "invalid umi_strategy".to_string()}.into()),
    };
    match count_reads::by_barcode::py_count_reads_primary_only_right_strand_only_by_barcode(
        filename,
        index_filename,
        trees,
        gene_trees,
        umi_strategy,
    ) {
        Ok(x) => Ok(x),
        Err(y) => Err(y.into()),
    }
}


/// python wrapper for py_count_introns
#[pyfunction]
pub fn count_introns(
    filename: &str,
    index_filename: Option<&str>,
) -> PyResult<count_reads::IntronResult> {
    let res = match count_reads::py_count_introns(filename, index_filename) {
        Ok(x) => x,
        Err(y) => return Err(y.into()),
    };
    Ok(res)
}
///
/// python wrapper for py_substract_bam
#[pyfunction]
pub fn subtract_bam(
    output_filename: &str,
    minuend_filename: &str,
    subtrahend_filename: &str,
) -> PyResult<()> {
    let res = match bam_manipulation::py_substract_bam(
        output_filename,
        minuend_filename,
        subtrahend_filename,
    ) {
        Ok(x) => x,
        Err(y) => return Err(y.into()),
    };
    Ok(res)
}

/// python wrapper for py_quantify_gene_reads
#[pyfunction]
pub fn quantify_gene_reads(
    filename: &str,
    index_filename: Option<&str>,
    intervals: &PyDict,
    gene_intervals: &PyDict,
) -> PyResult<(
    HashMap<String, Vec<(u32, u32)>>,
    HashMap<String, Vec<(u32, u32)>>,
)> {
    let trees = py_intervals_to_trees(intervals)?;
    let gene_trees = py_intervals_to_trees(gene_intervals)?;
    let res = match count_reads::py_quantify_gene_reads(filename, index_filename, trees, gene_trees)
    {
        Ok(x) => x,
        Err(y) => return Err(y.into()),
    };
    Ok(res)
}

/// python wrapper for py_annotate_barcodes_from_fastq
#[pyfunction]
pub fn annotate_barcodes_from_fastq(
    output_filename: &str,
    input_filename: &str,
    fastq2_filenames: Vec<&str>,
    barcodes: Vec<(String, usize, usize)>,
) -> PyResult<()> {
    match bam_manipulation::py_annotate_barcodes_from_fastq(
        output_filename,
        input_filename,
        fastq2_filenames,
        barcodes,
    ) {
        Ok(x) => Ok(x),
        Err(y) => return Err(y.into()),
    }
}

/// This module is a python module implemented in Rust.
#[pymodule]
fn mbf_bam(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(calculate_duplicate_distribution))?;
    m.add_wrapped(wrap_pyfunction!(count_reads_unstranded))?;
    m.add_wrapped(wrap_pyfunction!(count_reads_stranded))?;
    m.add_wrapped(wrap_pyfunction!(count_reads_primary_only_right_strand_only_by_barcode))?;
    m.add_wrapped(wrap_pyfunction!(count_introns))?;
    m.add_wrapped(wrap_pyfunction!(subtract_bam))?;
    m.add_wrapped(wrap_pyfunction!(quantify_gene_reads))?;
    m.add_wrapped(wrap_pyfunction!(annotate_barcodes_from_fastq))?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
//tests are in the callers until we can actually specify that we need mbf_align (and it's sample
//data) for the testing.
