
### Table of Contents:
[Introduction](#introduction)
[Pipeline](#Pipeline)
[Instillation](#Installation)
[Usage](#usage)
[Output](#output)

## Introduction
This nextflow pipeline is designed for fast exploration of many similar genome assembly methods for complex eukayotic genomes using mixed illumina and nanopore data.

This is a more generalized version of the assembly pipeline used for <insert publication>. After finding a well performing method, additional optimization of assembly parameters and analysis may need to be performed.

## Pipeline 

###### 1. Read processing and analysis
- Nanopore Reads
    - Guppy basecalling
    - Nanoplot
    - Nanofilt
    - Genome size estimation
    - Fastqc
    - phylogenetic analysis (subsample)
- Illumina WGS and RNAseq
    - Fastp

###### 2. Assembly
- Pomoxis
- Raven
- Shasta
- Shasta (haplotype mode)
- Flye
- Flye (high quality read mode)
- Canu
- Nextdenovo
- Smartdenovo
    
###### 3. Assembly polishing (where appropriate)
- racon (illumina read polishing)
- racon (nanopore read polishing)
- medaka
- nextpolish
- pilon

###### 4. Analysis (for each assembly)
- Quality
    - quast
    - busco
    - read coverage
    - merqury

- Contamination 
    - blobtools 
    - krona 
    - bacterial contig prediction 

###### 5. Annotation (single genome)
- repeat masker
- Star alignment
- Braker (RNAseq)
- Braker (conserved proteins)
- TSEBRA 
- eggnog functional annotation
- interproscan functional annotation

###### 6. Organeller Assembly (optional)
- Targeted read recruitment from evolutionily related organeller genomes
- Unicycler assembly
- get_organelle


## Installation

Most dependencies can be installed with the supplied conda yml file

```
    conda create env -f requirements.yml
```

A small handful of dependencies need to be install manually because they don't play nicly with conda.
1. [braker](https://github.com/Gaius-Augustus/BRAKER#installation) see special instructions, update parameters in config file for executable paths
2. [masurca](https://github.com/alekseyzimin/masurca#2-installation-instructions) needs to be installed manually
3. [Diamond](https://github.com/bbuchfink/diamond) database must be downloaded and built before use
4. [Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview) needs to be downloaded from ONT nanopore

Sometime, conda installs are not perfectly reproducible. It may be better to adopt a containerized approach.

## Usage
 To run the pipeline, first update the config file (example config provided at configs/example.config), then run the nextflow run_all.nx script
 ```
    nextflow run_all.nx -c <config.txt> [ -with-report <report.html> -w <work_directory> ]
 ```

 Depending on available compute resources and genome complexity, compute resouces in the nextflow.config may need to be updated. Currently, the nextflow.config is set to run on a large slurm cluster with ample resources

## Output

Results are published to the results directory provided in the configuration file. An overview of the results structure is provided below


```
Results_dir  
└─── Documentation
│   └─── Copy of all code used to analysis
└─── Reads
│   └─── Nanopore reads
│   └─── Illumina reads
└─── Assembly method 1
│   └─── Assembly
│   └─── Assembly analysis
│       ...      
└─── Assembly method N
```
