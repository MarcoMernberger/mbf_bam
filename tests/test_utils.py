from mbf_bam import reheader_and_rename_chromosomes, job_reheader_and_rename_chromosomes
from pathlib import Path
import pysam
from mbf_sampledata import get_sample_path
import pypipegraph as ppg
import pytest


class TestReheader:
    def test_rename(self, new_pipegraph):
        ppg.util.global_pipegraph.quiet = False
        input = get_sample_path("mbf_align/ex2.bam")
        output = "out.bam"
        job_reheader_and_rename_chromosomes(
            input, output, {"chr1": "shu", "chr2": "sha"}
        )
        ppg.run_pipegraph()
        assert Path("out.bam").exists()
        f = pysam.Samfile("out.bam")
        assert set(f.references) == set(["shu", "sha"])

    def test_rename_raises_on_no_replacement(self, new_pipegraph):
        ppg.util.global_pipegraph.quiet = False
        input = get_sample_path("mbf_align/ex2.bam")
        output = "out.bam"
        j = job_reheader_and_rename_chromosomes(input, output, {})
        with pytest.raises(ppg.RuntimeError):
            ppg.run_pipegraph()
        assert not Path("out.bam").exists()
        assert "No replacement happened" in str(j.exception)

class TestSubtract:

    def test_subtract_subset(self, new_pipegraph):
        from mbf_sampledata import get_sample_path
        from mbf_bam import subtract_bam

        input = get_sample_path("mbf_align/chipseq_chr22.bam")
        minued = get_sample_path("mbf_align/chipseq_chr22_subset_plus_unmapped.bam")
        output = "output.bam"
        print(input, input.exists())
        print(minued, minued.exists())
        subtract_bam(str(output), str(input.absolute()), str(minued.absolute()))
        f = pysam.Samfile(output)
        should = 80495
        total = sum((x.total for x in f.get_index_statistics()))
        assert should == total
