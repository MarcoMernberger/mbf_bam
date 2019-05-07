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
