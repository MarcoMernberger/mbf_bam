from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:  # pragma: no cover
    # package is not installed
    __version__ = 'unknown'
    pass

try:  # we need to ignore the import error (module not build) for poetry to be able to determine the version
    from .mbf_bam import *
    from pathlib import Path
    import pypipegraph as ppg
    import pysam
    import tempfile
except ImportError:
    pass


def reheader_and_rename_chromosomes(in_bam_file, out_bam_file, replacements):
    with pysam.Samfile(in_bam_file) as f:
        h = str(f.header)
        org_header = h
        for a, b in replacements.items():
            h = h.replace(f"SN:{a}", f"SN:{b}")
        if h == org_header:
            raise ValueError("No replacement happened")
        tf = tempfile.NamedTemporaryFile()
        tf.write(h.encode("utf-8"))
        tf.flush()
        out_bam_file.write_text("")  # must be there for save_stdout to work..
        pysam.reheader(
            tf.name,
            str(in_bam_file.absolute()),
            save_stdout=str(out_bam_file.absolute()).encode("utf-8"),
        )
        pysam.index(str(out_bam_file))


def job_reheader_and_rename_chromosomes(input_bam_path, output_bam_path, replacements):
    input_path_bam = Path(input_bam_path)
    output_bam_path = Path(output_bam_path)

    def do_replace(replacements=replacements):
        reheader_and_rename_chromosomes(input_bam_path, output_bam_path, replacements)

    output_bam_path.parent.mkdir(exist_ok=True, parents=True)
    return ppg.MultiFileGeneratingJob(
        [output_bam_path, output_bam_path.with_suffix(".bam.bai")], do_replace
    ).depends_on(
        ppg.FileInvariant(input_bam_path),
        ppg.FunctionInvariant(
            "mbf_bam.reheader_and_rename_chromosomes", reheader_and_rename_chromosomes
        ),
    )
