import sys
import subprocess
from pathlib import Path
import os


os.chdir(str(Path(__file__).parent.absolute()))

sys.path.insert(0, "../src")
try:
    import mbf_bam

    print("Building release for", mbf_bam.__version__)
except (ImportError, AttributeError):
    print("please install mbf_bam before creating a release")

if not Path('.git').exists():
    raise ValueError("Must be run from a git, not a mercurial repo")

cmd = [
    "docker",
    "run",
    "--rm",
    "-v",
    str(Path('..').resolve().absolute()) + ":/io",
    "quay.io/pypa/manylinux1_x86_64",
    "bash", "/io/release/build_wheels.sh"
]
import pprint
for x in cmd:
    print(x, end=' \\\n' if x != '-v' else ' ')
print('')
subprocess.check_call(cmd)
print("now upload wit twine upload dist/*.tar.gz")
