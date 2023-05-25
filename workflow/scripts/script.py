import subprocess
import sys

subprocess.call(
    ["cat"] + sys.argv[1:],
)
