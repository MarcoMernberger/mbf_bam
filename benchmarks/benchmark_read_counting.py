from pathlib import Path
import pysam
import pickle
import time
import pypipegraph as ppg
import mbf_genomes
import os
from old_read_counters import (
    GeneUnstrandedPython,
    GeneStrandedPython,
    ExonStrandedPython,
    ExonUnstrandedPython,
    ExonSmartUnstrandedPython,
    ExonSmartStrandedPython
)


work_dir = Path("_benchmark_read_counting")
work_dir.mkdir(exist_ok=True)
os.chdir(work_dir)

bam_name = (
    Path("results")
    / "aligned"
    / "STAR_2.6.1d"
    / "Drosophila_melanogaster_94"
    / "ERR2984187"
    / "ERR2984187.bam"
)

if not bam_name.exists():
    # leverage pipeline to get some sample data

    import mbf_align
    import mbf_externals

    ppg.new_pipegraph()

    genome = mbf_genomes.EnsemblGenome("Drosophila_melanogaster", 94)
    aligner = mbf_externals.aligners.STAR()

    # just some random drospohila lane.
    samples = {"ERR2984187": "ERR2984187"}
    raw = {
        name: mbf_align.Sample(
            name,
            mbf_align.strategies.FASTQsFromAccession(err),
            reverse_reads=False,
            pairing="only_first",
        )
        for name, err in samples.items()
    }

    aligned = {lane.name: lane.align(aligner, genome, {}) for lane in raw.values()}
    ppg.run_pipegraph()


def do_time(description, callback):
    p = Path("time_%s.txt" % description)
    rp = Path("result_%s.pickle" % description)
    if not p.exists():
        start = time.time()
        e = None
        try:
            result = callback()
        except Exception as ex:
            print("caught")
            e = ex
            pass
        stop = time.time()
        print(description, "%.2f" % (stop - start), "(new)", sep="\t")
        if not e:
            p.write_text("%.2f" % (stop - start))
            with open(rp, "wb") as op:
                pickle.dump(result, op)
        else:
            raise e
    else:
        print(description, p.read_text(), "(cached)", sep="\t")
    return pickle.load(open(rp, "rb"))


class Dummy:
    name = "lane"
    vid = None

    def __init__(self, genome, bam_name):
        self.genome = genome
        self.bam_name = bam_name

    def get_bam(self):
        return pysam.Samfile(self.bam_name)


genome = mbf_genomes.EnsemblGenome("Drosophila_melanogaster", 94)


def in_python_unstranded():
    import mbf_genomics

    ppg.util.global_pipegraph = None
    dummy = Dummy(genome, bam_name)
    genes = mbf_genomics.genes.Genes(genome)
    ac = GeneUnstrandedPython(dummy)
    genes += ac
    return genes.df.set_index("gene_stable_id")[ac.columns[0]]


def in_rust_unstranded():
    import mbf_genomics

    ppg.util.global_pipegraph = None
    dummy = Dummy(genome, bam_name)
    genes = mbf_genomics.genes.Genes(genome)
    ac = mbf_genomics.genes.anno_tag_counts.GeneUnstrandedRust(dummy)
    genes += ac
    return genes.df.set_index("gene_stable_id")[ac.columns[0]]


def in_python_stranded():
    import mbf_genomics

    ppg.util.global_pipegraph = None
    dummy = Dummy(genome, bam_name)
    genes = mbf_genomics.genes.Genes(genome)
    ac = GeneStrandedPython(dummy)
    ac.count_strategy.disable_sanity_check = True
    genes += ac
    return genes.df.set_index("gene_stable_id")[ac.columns[0]]


def in_rust_stranded():
    import mbf_genomics

    ppg.util.global_pipegraph = None
    dummy = Dummy(genome, bam_name)
    genes = mbf_genomics.genes.Genes(genome)
    ac = mbf_genomics.genes.anno_tag_counts.GeneStrandedRust(dummy)
    ac.count_strategy.disable_sanity_check = True
    genes += ac
    return genes.df.set_index("gene_stable_id")[ac.columns[0]]


def in_python_exon_stranded():
    import mbf_genomics

    ppg.util.global_pipegraph = None
    dummy = Dummy(genome, bam_name)
    genes = mbf_genomics.genes.Genes(genome)
    ac = ExonStrandedPython(dummy)
    ac.count_strategy.disable_sanity_check = True
    genes += ac
    return genes.df.set_index("gene_stable_id")[ac.columns[0]]


def in_rust_exon_stranded():
    import mbf_genomics

    ppg.util.global_pipegraph = None
    dummy = Dummy(genome, bam_name)
    genes = mbf_genomics.genes.Genes(genome)
    ac = mbf_genomics.genes.anno_tag_counts.ExonStrandedRust(dummy)
    ac.count_strategy.disable_sanity_check = True
    genes += ac
    return genes.df.set_index("gene_stable_id")[ac.columns[0]]


def in_python_exon_smart_stranded():
    import mbf_genomics

    ppg.util.global_pipegraph = None
    dummy = Dummy(genome, bam_name)
    genes = mbf_genomics.genes.Genes(genome)
    ac = ExonSmartStrandedPython(dummy)
    ac.count_strategy.disable_sanity_check = True
    genes += ac
    return genes.df.set_index("gene_stable_id")[ac.columns[0]]


def in_rust_exon_smart_stranded():
    import mbf_genomics

    ppg.util.global_pipegraph = None
    dummy = Dummy(genome, bam_name)
    genes = mbf_genomics.genes.Genes(genome)
    ac = mbf_genomics.genes.anno_tag_counts.ExonSmartStrandedRust(dummy)
    ac.count_strategy.disable_sanity_check = True
    genes += ac
    return genes.df.set_index("gene_stable_id")[ac.columns[0]]


def check_results(py, rust):
    if (rust == py).all():
        print("OK")
    ok = 0
    err = 0
    # for stable_id in genome.df_genes[genome.df_genes.chr == '3L'].index:
    # for stable_id in ['FBgn0035181']:
    for stable_id in genome.df_genes.index:
        if py[stable_id] == rust[stable_id]:
            ok += 1
        elif rust[stable_id] > 0:
            # print(stable_id, py[stable_id], rust[stable_id])
            err += 1
    print("ok", ok, "err", err)


for p in (
    "result_rust_exon_smart_stranded.pickle",
    "time_rust_exon_smart_stranded.txt",
):
    p = Path(p)
    if p.exists():
        p.unlink()

py = do_time("python_exon_smart_stranded", in_python_exon_smart_stranded)
rust = do_time("rust_exon_smart_stranded", in_rust_exon_smart_stranded)
check_results(py, rust)


py = do_time("python_exon_stranded", in_python_exon_stranded)
rust = do_time("rust_exon_stranded", in_rust_exon_stranded)
check_results(py, rust)

py = do_time("python_gene_stranded", in_python_stranded)
rust = do_time("rust_gene_stranded", in_rust_stranded)
check_results(py, rust)


py = do_time("python_gene_unstranded", in_python_unstranded)
rust = do_time("rust_gene_unstranded", in_rust_unstranded)
check_results(py, rust)
