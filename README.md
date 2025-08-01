# Accurate and Reproducible De Novo Variants Calling with TriosCompass

## Introduction

This manuscript branch and its documentation are dedicated to the manuscript titled "Accurate and Reproducible De Novo Variants Calling with TriosCompass".  Please use the main branch to follow the latest release of TriosCompass.

## Methods

### Run TriosCompass for the GIAB trio
This *manuscript* branch, along with its dedicated documentation, supports the manuscript titled "Accurate and Reproducible De Novo Variants Calling with TriosCompass." For the latest stable release and ongoing development of TriosCompass, please refer to [the *main* branch](https://github.com/NCI-CGR/TriosCompass_v2).

#### TiosCompass installation
+ Detailed installation instructions for TriosCompass are available in the [main branch's README](https://github.com/NCI-CGR/TriosCompass_v2/tree/main?tab=readme-ov-file#i-installation).

#### Workspace folder structure 
The recommended folder structure for your analysis workspace is as follows:
<code>
[Workspace]
├── bam/...
├── input
│   ├── <a href="./data/input/40X_bam.csv">40X_bam.csv</a>
│   ├── <a href="./data/input/40X_pep.yaml">40X_pep.yaml</a>
│   ├── <a href="./data/input/80X_bam.csv">80X_bam.csv</a>
│   ├── <a href="./data/input/80X_pep.yaml">80X_pep.yaml</a>
│   └── ped
├── <a href="./config/GIAB_40X.yaml">GIAB_40X.yaml</a>
├── <a href="./config/GIAB_80X.yaml">GIAB_80X.yaml</a>
├── <a href="./data/launch_40X.sh">launch_40X.sh</a>
├── <a href="./data/launch_80X.sh">launch_80X.sh</a>
├── ref/...
└── TriosCompass_v2/...
</code>

#### Launch TriosCompass
After setting up your workspace, execute TriosCompass using the following commands. These examples assume activation of the *TriosCompassV2* Conda environment.
```bash
# Activate the TriosCompass environment
conda activate TriosCompassV2

# Launch TriosCompass for 40X coverage
# (Refer to launch_40X.sh for full details)
snakemake --profile TriosCompass_v2/workflow/profiles/slurm --configfile GIAB_40X.yaml --conda-frontend mamba 

# Launch TriosCompass for 80X coverage
# (Refer to launch_80X.sh for full details)
snakemake --profile TriosCompass_v2/workflow/profiles/slurm --configfile GIAB_80X.yaml --conda-frontend mamba
```

---

### Benchmarking with DeepTrio
Details on running DeepTrio for benchmarking purposes are provided in the supplementary methods document: [TriosCompass_Supp_Methods_deeptrio_benchmark.md](./TriosCompass_Supp_Methods_deeptrio_benchmark.md).

