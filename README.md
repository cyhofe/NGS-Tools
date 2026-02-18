# NGS-Tools

Lightweight, reproducible wrapper scripts for common **Illumina and Nanopore sequencing workflows**.

These scripts standardize frequently used short- and long-read processing pipelines while keeping execution transparent and HPC-friendly.

---

## Installation

```bash
git clone https://github.com/cyhofe/NGS-Tools.git
cd NGS-Tools
chmod +x scripts/*.pl
```

Optional (recommended):

```bash
export PATH="$PWD/scripts:$PATH"
```

---

## How these scripts work

Each Perl wrapper script generates a self-contained bash workflow (`*.sh`) that contains all commands that will be executed.

### Typical flow

1. Validate input arguments  
2. Create output directory  
3. Write a reproducible bash script (e.g. `PREFIX.qualcheck.sh`)  
4. Either:
   - run it immediately with `--run`, or  
   - leave it for manual execution / SLURM submission  

---

## Example

```bash
Nanopore-quality-filter.pl -f reads.fq.gz -p sampleA -o sampleA.qc --run
```

This will:

- Create `sampleA.qc/`
- Generate `sampleA.qualcheck.sh` (full reproducible workflow)
- Execute it (because `--run` was set)

This design is ideal for HPC because you can generate workflows for many samples and run them using GNU parallel or your scheduler.

---

## Included workflows

### Illumina

#### Illumina-quality-filter.pl

Paired-end preprocessing:

- Interleave reads
- Quality trimming
- Vector removal
- Adapter trimming
- FastQC reports

---

### Nanopore (Metagenomics)

#### Nanopore-quality-filter.pl

Quality + length filtering and adapter trimming.

#### Nanopore-flye-assembly.pl

Flye metagenome assembly + contig filtering/renaming.

#### Nanopore-racon-polishing.pl

Long-read polishing using racon.

#### Hybrid-pilon-polishing.pl

Hybrid polishing using Illumina reads (BWA + samtools + Pilon).

#### Nanopore-short-read-mapping.pl

Map Illumina reads to a long-read assembly (sorted BAM + stats).

---

### Nanopore (16S Amplicon)

#### Nanopore-quality-filter-16srRNA.pl

16S amplicon preprocessing with size-window filtering.

#### Nanopore-reverse-complement-16srRNA.pl

Strand/orientation correction using minimap2 + seqtk.

#### Nanopore-chimera-checker-16srRNA.pl

Chimera filtering using minimap2 overlaps + yacrd + NanoStat.

---

## Designed for

- Reproducible microbial genomics workflows  
- HPC / cluster environments  
- GNU parallel batch processing  
- SLURM submission (generate scripts, then submit `*.sh`)  
- Transparent logs and easy debugging  

---

## Dependencies

Depending on workflow:

- BBMap suite (`reformat.sh`, `bbduk.sh`, `rename.sh`, `bbmerge.sh`)
- FastQC
- NanoStat (NanoPack)
- minimap2
- seqtk
- chopper
- porechop_abi
- flye
- racon
- bwa
- samtools
- java + `pilon.jar`
- yacrd

Executables are resolved from `$PATH` or via optional `*_dir` parameters in the scripts.

---

## Maintainer

**Cyrill Hofer**  
GitHub: https://github.com/cyhofe
