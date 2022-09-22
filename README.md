# NanoGAP


The Nanopore Genome Assembly Pipeline (`NanoGAP`) utilises a variety of open-source tools and software to automate a non-hybrid (long reads only) genome assembly and self-correction of error prone nanopore reads.

<!-- ---

## Updates
 -->

---

## Features

- Developed specifically for bacterial genomes.
- Requires raw/uncorrected Oxford Nanopore Technologies (ONT) long reads only (`.fastq`).
- Assembly of unknown isolates (no predicted genome size needed).
- Identify nearest bacterial organism.

---

## Installation

Conda is required for installation. You can download and install conda for Linux [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

NanoGAP uses a number of open source projects:
- [flye](https://github.com/fenderglass/Flye)
- [medaka](https://github.com/nanoporetech/medaka)
- [barrnap](https://github.com/tseemann/barrnap)
- [blast](https://github.com/ncbi/blast_plus_docs)


```sh
git clone https://github.com/escasinas/nanogap.git
cd ngap
source install.sh
```

---

## Usage

**Input**
- A single `.fastq` file.

or

- A directory containing multiple `.fastq` files.

**Output**

- Genome assembly in `.fasta` format.
- BLAST output of the genome's 16S rRNA.
- CSV output containing assembly information.

**Command**

For a single fastq file
```sh
conda activate nanogap
python nanogap.py path/to/file.fastq [options]
```

For a directory containing fastq files
```sh
conda activate nanogap
python nanogap.py path/to/reads_directory [options]
```

---

## Options

`-h | --help` Show help message and exit.

`-t | --threads` Number of threads/CPUs to run Minimap2, Flye, Racon, Medaka and Barrnap (**default: 4**).

`-o | --outdir` Name of output directory (**default: ngap_output**).

`-m | --model` Medaka model. Please see the [Medaka repo](https://github.com/nanoporetech/medaka#models) for more information (**default: r941_min_high_g360**).
