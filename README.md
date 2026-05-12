# MaMISA — Manage Misassemblies

A comprehensive toolkit for metagenomic assembly quality control and filtering.

## Overview

MaMISA provides five commands that cover the full workflow from raw assembly to taxonomically classified, quality-filtered genomes:

| Command | Purpose |
|---|---|
| `process-large-contigs` | Extract, QC, and classify large contigs before binning |
| `filter-misassemblies` | Split or remove contigs with misassemblies detected by anvi'o |
| `remove-hq-contigs` | Remove HQ genome contigs from an assembly |
| `filter-checkm2` | Organize genomes into quality tiers from CheckM2 results |
| `run-gtdbtk` | Run GTDB-Tk taxonomy classification |

---

## Installation

```bash
git clone https://github.com/yourusername/mamisa.git
cd mamisa
pip install -e .

# Verify
mamisa --version
mamisa --help
```

### External tools (required per command)

| Command | Tool |
|---|---|
| `process-large-contigs` | CheckM2 |
| `filter-misassemblies` | anvi'o |
| `filter-checkm2` | CheckM2 |
| `run-gtdbtk` | GTDB-Tk |

```bash
# anvi'o
conda create -n anvio-8 -c conda-forge -c bioconda anvio=8

# CheckM2
conda create -n checkm2 -c conda-forge -c bioconda checkm2

# GTDB-Tk
conda create -n gtdbtk-2 -c conda-forge -c bioconda gtdbtk=2.3.2
export GTDBTK_DATA_PATH=/path/to/gtdbtk_data
```

---

## Complete Workflow

```
assembly.fa
    │
    ▼
[1] process-large-contigs       → identifies very large contigs, runs CheckM2 on them,
    │                             extracts clean HQ ones, returns updated assembly
    │
    ▼
[2] anvi-script-find-misassemblies  (external — anvi'o)
    │
    ▼
[3] filter-misassemblies        → splits/removes misassembled contigs,
    │                             preserves HQ genomes with recoverable misassemblies
    │
    ▼
[4] binning                     (external — MetaBAT2, MaxBin2, etc.)
    │
    ▼
[5] checkm2 predict             (external — CheckM2)
    │
    ▼
[6] filter-checkm2              → organizes bins into HQ / MQ / LQ tiers
    │
    ▼
[7] run-gtdbtk                  → taxonomic classification of selected genomes
```

---

## Step-by-Step Guide

### Step 1 — Process large contigs

Very long contigs (> 300 kbp by default) often represent single complete genomes and
should be handled separately before binning. This command:

1. Separates large and regular contigs
2. Runs CheckM2 on each large contig individually
3. Classifies each one as: *extract HQ*, *keep for splitting*, or *keep low quality*
4. Returns an updated assembly ready for misassembly filtering

```bash
mamisa process-large-contigs \
    --assembly assembly.fa \
    --misassemblies misasm_dir/ \
    --output-dir 01_large_contigs/ \
    --max-length 300000 \
    --min-completeness 50 \
    --max-contamination 10 \
    --threads 40
```

Output:
```
01_large_contigs/
├── 01_extracted/
│   ├── large_contigs.fa          # contigs > max-length
│   └── assembly_regular.fa       # contigs <= max-length
├── 02_individual/                 # one .fa per large contig (CheckM2 input)
├── 03_checkm2/                    # CheckM2 results
│   └── quality_report.tsv
├── HQ_extracted/                  # clean HQ genomes extracted here
├── filtering_decisions.tsv        # per-contig decision log
└── assembly_for_filtering.fa      # updated assembly → use in Step 3
```

To skip CheckM2 if already run:
```bash
mamisa process-large-contigs \
    --assembly assembly.fa \
    --misassemblies misasm_dir/ \
    --output-dir 01_large_contigs/ \
    --skip-checkm2 \
    --checkm2-results 01_large_contigs/03_checkm2/
```

---

### Step 2 — Detect misassemblies with anvi'o (external)

```bash
anvi-script-find-misassemblies \
    -b mapping.bam \
    -o misassemblies/MisAsm \
    -T 40
```

This produces `MisAsm-clipping.txt` and `MisAsm-zero_cov.txt` in `misassemblies/`.
MaMISA uses only the `-clipping.txt` files.

---

### Step 3 — Filter misassemblies

Processes the assembly from Step 1, using misassembly data from Step 2.

```bash
# Preview first (dry-run)
mamisa filter-misassemblies \
    --assembly 01_large_contigs/assembly_for_filtering.fa \
    --misassemblies misassemblies/ \
    --hq-genomes 01_large_contigs/HQ_extracted/ \
    --dry-run

# Apply
mamisa filter-misassemblies \
    --assembly 01_large_contigs/assembly_for_filtering.fa \
    --misassemblies misassemblies/ \
    --hq-genomes 01_large_contigs/HQ_extracted/ \
    --output 02_clean_assembly/assembly_clean.fa \
    --mode split \
    --min-length 2500 \
    --preserve-hq-with-issues \
    --stats 02_clean_assembly/stats.tsv
```

**HQ preservation logic with `--preserve-hq-with-issues`:**

| Contig type | Action |
|---|---|
| HQ + no misassemblies | Removed (already extracted in Step 1) |
| HQ + misassemblies | Split at clipping positions, fragments kept |
| Non-HQ + misassemblies | Split at clipping positions |
| Non-HQ + clean | Kept intact |

Without `--preserve-hq-with-issues`, all HQ contigs are removed regardless.

---

### Step 4 — Binning (external)

Run your preferred binner on the clean assembly. Example with MetaBAT2:

```bash
# Generate depth file first
jgi_summarize_bam_contig_depths --outputDepth depth.txt mapping.bam

# Bin
metabat2 \
    -i 02_clean_assembly/assembly_clean.fa \
    -a depth.txt \
    -o 03_bins/bin \
    -t 40
```

---

### Step 5 — Quality assessment with CheckM2 (external)

```bash
checkm2 predict \
    --threads 40 \
    --input 03_bins/ \
    --output-directory 04_checkm2/ \
    -x fa
```

---

### Step 6 — Filter genomes by quality

Organizes bins into quality tiers following MiMAG standards.

```bash
mamisa filter-checkm2 \
    --checkm2-root 04_checkm2/ \
    --genomes-dir 03_bins/ \
    --output 05_filtered_genomes/ \
    --tiers HQ,MQ \
    --hq-comp-min 90 \
    --hq-cont-max 5 \
    --mq-comp-min 70 \
    --mq-cont-max 10 \
    --symlink
```

Output:
```
05_filtered_genomes/
├── merged_quality.tsv
└── Selected/
    ├── HQ/
    └── MQ/
```

**MiMAG quality thresholds (defaults):**

| Tier | Completeness | Contamination |
|---|---|---|
| HQ | ≥ 90% | ≤ 5% |
| MQ | ≥ 70% | ≤ 10% |
| LQ | ≥ 50% | ≤ 10% |
| Fail | < 50% or > 10% contamination | |

---

### Step 7 — Remove HQ contigs from assembly (optional)

If you need a version of the assembly without any HQ genome contigs
(e.g. for a second round of binning):

```bash
mamisa remove-hq-contigs \
    --assembly 02_clean_assembly/assembly_clean.fa \
    --hq-dir 05_filtered_genomes/Selected/HQ/ \
    --output 06_assembly_no_hq/assembly_no_hq.fa \
    --min-length 1000 \
    --stats 06_assembly_no_hq/stats.tsv
```

---

### Step 8 — Taxonomic classification with GTDB-Tk

Classify all selected genomes by quality tier:

```bash
mamisa run-gtdbtk \
    --selected-dir 05_filtered_genomes/Selected/ \
    --output 07_taxonomy/ \
    --extension fa \
    --cpus 40 \
    --tiers HQ,MQ
```

Or classify a single directory:

```bash
mamisa run-gtdbtk \
    --genome-dir 05_filtered_genomes/Selected/HQ/ \
    --output 07_taxonomy/HQ/ \
    --extension fa \
    --cpus 40
```

---

## Command Reference

### process-large-contigs

```
mamisa process-large-contigs [OPTIONS]

Required:
  -a, --assembly PATH         Input assembly FASTA file
  -m, --misassemblies PATH    Directory with *-clipping.txt files
  -o, --output-dir PATH       Output directory

Thresholds:
  --max-length INT            Contigs longer than this are large (default: 300000)
  --min-completeness FLOAT    Min completeness to classify as HQ (default: 50)
  --max-contamination FLOAT   Max contamination to classify as HQ (default: 10)

CheckM2:
  --threads INT               Threads for CheckM2 (default: 1)
  --skip-checkm2              Skip CheckM2 (requires --checkm2-results)
  --checkm2-results PATH      Existing CheckM2 output directory
  --checkm2-env STR           Conda env name for CheckM2 (default: checkm2)
```

### filter-misassemblies

```
mamisa filter-misassemblies [OPTIONS]

Required:
  -a, --assembly PATH         Input assembly FASTA file
  -m, --misassemblies PATH    Directory with *-clipping.txt files

Optional:
  -g, --hq-genomes PATH       Directory with HQ genome files
  -o, --output PATH           Output filtered assembly
  -l, --min-length INT        Minimum contig/fragment length (default: 2500)
  --mode {remove,split}       Misassembly handling (default: split)
  --preserve-hq-with-issues   Keep HQ genomes with misassemblies (split them)
  --dry-run                   Preview without writing output
  --stats PATH                Save statistics to TSV
```

### remove-hq-contigs

```
mamisa remove-hq-contigs [OPTIONS]

Required:
  -a, --assembly PATH         Input assembly FASTA file
  --hq-dir PATH               Directory with HQ genome files
    OR
  --hq-list PATH              Text file with HQ contig IDs (one per line)

Optional:
  -o, --output PATH           Output filtered assembly
  -l, --min-length INT        Minimum contig length (default: 0)
  --dry-run                   Preview without writing output
  --stats PATH                Save statistics to TSV
```

### filter-checkm2

```
mamisa filter-checkm2 [OPTIONS]

Required:
  --checkm2-root PATH         Root directory with CheckM2 results
  --genomes-dir PATH          Directory with genome files
  -o, --output PATH           Output directory

Quality thresholds:
  --hq-comp-min FLOAT         HQ min completeness (default: 90)
  --hq-cont-max FLOAT         HQ max contamination (default: 5)
  --mq-comp-min FLOAT         MQ min completeness (default: 70)
  --mq-cont-max FLOAT         MQ max contamination (default: 10)
  --lq-comp-min FLOAT         LQ min completeness (default: 50)
  --lq-cont-max FLOAT         LQ max contamination (default: 10)

Other:
  --tiers LIST                Tiers to select, comma-separated (default: HQ,MQ,LQ)
  --extensions LIST           Genome file extensions (default: fa,fasta,fna,fa.gz,...)
  --symlink | --copy          Link or copy files (default: symlink)
  --name-map PATH             TSV mapping report names to filenames
  --strip-prefix STR          Remove prefix from genome names
  --strip-suffix STR          Remove suffix from genome names
  --name-prefix STR           Add prefix to genome names
  --name-suffix STR           Add suffix to genome names
  --dry-run                   Preview without writing output
```

### run-gtdbtk

```
mamisa run-gtdbtk [OPTIONS]

Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --extension STR             Genome file extension (default: fa)
  --cpus INT                  CPUs for GTDB-Tk (default: 1)
  --mash-db PATH              GTDB-Tk mash database
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --gtdbtk-args STR           Extra arguments passed to gtdbtk
```

---

## Understanding Misassembly Detection

MaMISA uses clipping information from read mapping to detect misassemblies. A
**clipping position** is a genomic location where most or all mapped reads are
soft-clipped — strong evidence that two unrelated sequences were joined during assembly.

A clipping ratio of 1.0 (100% of reads clipped at a position) is the most reliable
signal. Zero-coverage regions are also indicative but less specific; MaMISA currently
uses only clipping data.

When a contig is split, the clipping position itself is the cut point. Fragments shorter
than `--min-length` are discarded.

---

## Citation

If you use MaMISA in your research, please cite:

```
[citation here]
```

## License

MIT License — see LICENSE for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

Issues and feature requests: https://github.com/yourusername/mamisa/issues
