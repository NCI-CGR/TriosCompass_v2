<!-- omit in TOC -->
# Snakemake workflow: TriosCompass

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.3.7-brightgreen.svg)](https://snakemake.github.io)


A Snakemake workflow for DNM (de novo mutation) calling.

---
- [Snakemake workflow: TriosCompass](#snakemake-workflow-trioscompass)
  - [Overview](#overview)
    - [I. Introduction](#i-introduction)
    - [II. Dependencies](#ii-dependencies)
  - [User's guides](#users-guides)
    - [Installation](#installation)
    - [Inputs](#inputs)
    - [Outputs](#outputs)
    - [Run TrisCompass](#run-triscompass)


---

## Overview

### I. Introduction
TriosCompass consists of 4 functional components:
+ Call DNMs (de novo mutations) using *DeepVariant* and *GATK HaplotypeCaller*
+ Phase DNMs using *whatshap*
+ Call dnSTR (de novo simple tandem repeats) using HipSTR and MonSTR.
+ Call dnSV (de novo structural variants) using *Manta*, *GraphType2* and *smoove*

---

### II. Dependencies

All required bioinformatics tools are wrapped as conda, container and etc, so there is no need to for users to install any of them.

Nevertheless, there are still some dependencies required to start Snakemake workflow, which have been specified in [environment.yaml](environment.yaml).  Users can create a new *conda* env for TriosCompassV2
```bash

mamba env create -f environment.yaml

conda activate TriosCompassV2
```

Besides, singularity needs to be installed globally. 

---

## User's guides

### Installation 

### Inputs

### Outputs

### Run TrisCompass