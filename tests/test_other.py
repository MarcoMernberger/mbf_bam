def test_version_is_correct():
    import configparser
    from pathlib import Path
    import mbf_bam
    import tomlkit

    # this is the truth
    c = configparser.ConfigParser()
    c.read(Path(__file__).parent.parent / "setup.cfg")
    version = c["metadata"]["version"]

    assert version == mbf_bam.__version__

    cargo_toml = tomlkit.loads((Path(__file__).parent.parent / "Cargo.toml").read_text())
    cargo_version = cargo_toml["package"]["version"]
    assert cargo_version == version
