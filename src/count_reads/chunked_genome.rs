use super::OurTree;
use rust_htslib::bam;
use rust_htslib::prelude::*;
use std::collections::HashMap;
use std::str;

pub struct ChunkedGenome {
    trees: Option<HashMap<String, (OurTree, Vec<String>)>>,
    bam: bam::IndexedReader,
    chromosomes: Vec<String>,
}

impl ChunkedGenome {
    ///create a new chunked genome for iteration
    ///if you pass in a tree, it is guranteed that the splits happen
    ///between entries of the tree, not inside.
    pub fn new(
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
    pub fn new_without_tree(bam: bam::IndexedReader) -> ChunkedGenome {
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
    pub fn iter(&self) -> ChunkedGenomeIterator {
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

pub struct ChunkedGenomeIterator<'a> {
    cg: &'a ChunkedGenome,
    it: std::slice::Iter<'a, String>,
    last_start: u32,
    last_chr: String,
    last_tid: u32,
    last_chr_length: u32,
}
#[derive(Debug)]
pub struct Chunk {
    pub chr: String,
    pub tid: u32,
    pub start: u32,
    pub stop: u32,
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
