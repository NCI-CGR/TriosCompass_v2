# TriosCompass_v2
This is a trios analysis workflow written in Snakemake.


## Installation

The workflow has been tested under biowulf with snakemake version 7.3.7.  Snakemake is installed under conda.

---

## Quick start 

sbatch -J snakemake -t 200:00:00 --export=ALL --mem=12g -p norm  --wrap='./run_it_gpu.sh '
