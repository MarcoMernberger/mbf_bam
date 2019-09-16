
///
pub fn py_quantify_gene_reads(
    filename: &str,
    index_filename: Option<&str>,
    trees: HashMap<String, (OurTree, Vec<String>)>,
    gene_trees: HashMap<String, (OurTree, Vec<String>)>,
) -> Result<(HashMap<String, (u32, u32), HashMap<String, (u32, u32)), BamError> {
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

