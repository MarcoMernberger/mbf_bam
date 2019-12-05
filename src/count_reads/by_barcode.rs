///
/// count reads, but split them by (multiple) barcodes
/// ie. singe cell stuff
 use super::chunked_genome::{Chunk, ChunkedGenome};
use super::{OurTree};
use crate::bam_ext::{open_bam, BamRecordExtensions};
use crate::BamError;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::collections::{HashMap, HashSet};
/// python wrapper count_reads_primary_only_right_strand_only_by_barcode
///
pub enum UmiStrategy {
    Straight,
//    Hamming(u8),
}
pub fn py_count_reads_primary_only_right_strand_only_by_barcode(
    filename: &str,
    index_filename: Option<&str>,
    trees: HashMap<String, (OurTree, Vec<String>)>,
    gene_trees: HashMap<String, (OurTree, Vec<String>)>,
    umi_strategy: UmiStrategy
) -> Result<(
        Vec<String>, //barcodes
        Vec<String>, //genes
        Vec<(u32, u32, u32)> //barcode_id, gene_id, count

        ), BamError> {
    //check whether the bam file can be openend
    //and we need it for the chunking
    let bam = open_bam(filename, index_filename)?;

    //perform the counting
    let cg = ChunkedGenome::new(gene_trees, bam); // can't get the ParallelBridge to work with our lifetimes.
    let it: Vec<Chunk> = cg.iter().collect();
    let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
    let with_barcodes = pool.install(|| {
        let result = it
            .into_par_iter()
            .map(|chunk| {
                let bam = open_bam(filename, index_filename).unwrap();
                let (tree, gene_ids) = trees.get(&chunk.chr).unwrap();

            let mut result: Vec<(String, String, u32)> = Vec::new();
                match count_reads_primary_only_right_strand_only_by_barcode(
                    bam,
                    //&chunk.tree,
                    tree,
                    chunk.tid,
                    chunk.start,
                    chunk.stop,
                    gene_ids.len() as u32,
                    &umi_strategy,
                ) {
                    Ok(counts) => {
                        for (gene_id, barcode, count) in counts.into_iter() {
                            let barcode = std::str::from_utf8(&barcode).unwrap().to_string();
                            result.push((
                                    gene_ids[gene_id as usize].clone(),
                                    barcode,
                                    count));
                        }
                    }, 
                    _ => {}
                }
                result
            })
            .reduce(Vec::<(String, String, u32)>::new, |mut a, b | { a.extend(b); a});
        //.fold(HashMap::<String, u32>::new(), add_hashmaps);
        result
    });

    let mut genes = HashMap::new();
    let mut barcodes = HashMap::new();
    let mut counts: Vec<(u32, u32, u32)> = Vec::new();
    for (gene_name, barcode, count) in with_barcodes.into_iter() {
        let l = genes.len() as u32;
        let gene_id = *(genes.entry(gene_name).or_insert(l));
        let l = barcodes.len() as u32;
        let barcode_id = *(barcodes.entry(barcode).or_insert(l));

        counts.push((gene_id, barcode_id, count));
    }
    let genes = to_lookup_array(genes);
    let barcodes = to_lookup_array(barcodes);
    Ok((genes, barcodes, counts))
}

fn to_lookup_array(mut input: HashMap<String, u32>) -> Vec<String> {
    let mut flat: Vec<(String, u32)> = input.drain().collect();
    flat.sort_by(|a, b| (a.1.cmp(&b.1)));
    flat.into_iter().map(|(k,_) | k).collect()
}

/// count_reads_primary_only_right_strand_only_by_barcode
///
/// counts the correct strand reads in genes,
/// umi deduped,
/// and split by (cell) barcodes in the XC tag.
fn count_reads_primary_only_right_strand_only_by_barcode
(

    mut bam: bam::IndexedReader,
    tree: &OurTree,
    tid: u32,
    start: u32,
    stop: u32,
    _gene_count: u32,
    umi_strategy: &UmiStrategy
) -> Result<Vec<(
            u32, //gene id
            Vec<u8>, //barcode
            u32, //count
             )>, BamError> 
    {
    let mut read: bam::Record = bam::Record::new();
    bam.fetch(tid, start, stop)?;
    let mut positions: HashMap<
        (u32, Vec<u8>, i32, bool), HashSet<Vec<u8>>> = HashMap::new();;
    while let Ok(_) = bam.read(&mut read) {
        let mut skipped = false;
        if ((read.pos() as u32) < start) || ((read.pos() as u32) >= stop) {
            skipped = true;
        }
        if read.is_secondary() { continue }
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
                    let entry = r.data();
                    let gene_no = (*entry).0;
                    let strand = (*entry).1; // this is 1 or -1
                    //if we are on the right strand...
                    if ((strand == 1) && !read.is_reverse()) || ((strand != 1) && read.is_reverse()) {
                        let barcode = read.aux(b"XC").ok_or_else(|| BamError::UnknownError{msg: "missing XC tag".to_string()})?;
                        let barcode = match barcode {
                            bam::record::Aux::String(m) => m.to_vec(), 
                            _ => return Err(BamError::UnknownError{msg: "XC did not contain string".to_string()}),
                        };

                        let umi = read.aux(b"XM").ok_or_else(|| BamError::UnknownError{msg: "missing XC tag".to_string()})?;
                        let umi = match umi {
                            bam::record::Aux::String(m) => m.to_vec(), 
                            _ => return Err(BamError::UnknownError{msg: "XM did not contain string".to_string()}),
                        };
                        let real_position = read.pos(); //TODO

                        positions.entry(
                            (gene_no, barcode, real_position, read.is_reverse())).or_insert_with(|| HashSet::new()).insert(
                                umi);
                    }
                }
               }
            }
        }
    
    let mut by_gene_barcode: HashMap<(u32, Vec<u8>), u32> = HashMap::new();
    for ((gene_no, barcode, _position, _reverse), umis) in positions.into_iter() {
        let e = by_gene_barcode.entry((gene_no, barcode)).or_insert(0);
        match umi_strategy {
            UmiStrategy::Straight => {
                *e += umis.len() as u32;
            }
        }
        
    }
    let mut result = Vec::new();
    for ((gene_no, barcode), count) in by_gene_barcode.into_iter() {
        result.push((gene_no, barcode, count));
    }
    Ok(result)
}

