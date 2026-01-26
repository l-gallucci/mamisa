# MaMISA - Manage Misassemblies

A comprehensive toolkit for metagenomic assembly quality control and filtering.

## Features

- **Filter Misassemblies**: Detect and handle misassemblies using anvi'o clipping data
- **HQ Genome Management**: Remove or preserve high-quality genomes with intelligent misassembly detection
- **CheckM2 Integration**: Filter and organize genomes by quality tiers (HQ/MQ/LQ)
- **GTDB-Tk Wrapper**: Streamlined taxonomy classification workflow

## Installation

### From source

```bash
git clone https://github.com/yourusername/mamisa.git
cd mamisa
pip install -e .
```

### Requirements

- Python >= 3.7
- For specific commands:
  - `filter-misassemblies`: anvi'o (for generating misassembly files)
  - `filter-checkm2`: CheckM2
  - `run-gtdbtk`: GTDB-Tk

## Quick Start

### 1. Filter Misassemblies

Remove or split contigs with misassemblies detected by anvi'o:

```bash
# Dry-run to see what would happen
mamisa filter-misassemblies \
    --assembly assembly.fa \
    --misassemblies misasm_dir/ \
    --hq-genomes hq_dir/ \
    --dry-run

# Process with HQ preservation (recommended)
mamisa filter-misassemblies \
    --assembly assembly.fa \
    --misassemblies misasm_dir/ \
    --hq-genomes hq_dir/ \
    --output clean.fa \
    --mode split \
    --min-length 2500 \
    --preserve-hq-with-issues
```

**Key features:**
- Intelligently preserves HQ genomes that have recoverable misassemblies
- Splits contigs at misassembly positions instead of discarding entirely
- Filters fragments by minimum length
- Provides detailed statistics

### 2. Remove HQ Contigs

Remove high-quality genome contigs from assembly:

```bash
mamisa remove-hq-contigs \
    --assembly assembly.fa \
    --hq-dir hq_genomes/ \
    --output filtered.fa \
    --min-length 1000
```

### 3. Filter by CheckM2 Quality

Organize genomes by quality tiers based on CheckM2 results:

```bash
mamisa filter-checkm2 \
    --checkm2-root checkm2_results/ \
    --genomes-dir genomes/ \
    --output filtered/ \
    --tiers HQ,MQ \
    --hq-comp-min 90 \
    --hq-cont-max 5
```

**Output structure:**
```
filtered/
├── merged_quality.tsv
└── Selected/
    ├── HQ/
    ├── MQ/
    └── LQ/
```

### 4. Run GTDB-Tk

Classify genomes taxonomically:

```bash
# On tier-organized genomes
mamisa run-gtdbtk \
    --selected-dir filtered/Selected/ \
    --output gtdbtk_results/ \
    --extension fa \
    --cpus 40

# On single directory
mamisa run-gtdbtk \
    --genome-dir my_genomes/ \
    --output gtdbtk_results/ \
    --extension fa \
    --cpus 40
```

## Workflow Example

Complete workflow from assembly to classified genomes:

```bash
# 1. Run anvi'o misassembly detection (external)
anvi-script-find-misassemblies \
    -b mapping.bam \
    -o misassemblies/prefix

# 2. Filter assembly preserving HQ genomes without issues
mamisa filter-misassemblies \
    --assembly assembly.fa \
    --misassemblies misassemblies/ \
    --hq-genomes initial_hq/ \
    --output clean_assembly.fa \
    --mode split \
    --min-length 2500 \
    --preserve-hq-with-issues

# 3. Bin genomes with your favorite binner
# ... binning step ...

# 4. Run CheckM2 (external)
checkm2 predict --input bins/ --output-directory checkm2_results/

# 5. Filter genomes by quality
mamisa filter-checkm2 \
    --checkm2-root checkm2_results/ \
    --genomes-dir bins/ \
    --output filtered_genomes/ \
    --tiers HQ,MQ

# 6. Classify with GTDB-Tk
mamisa run-gtdbtk \
    --selected-dir filtered_genomes/Selected/ \
    --output taxonomy/ \
    --extension fa \
    --cpus 40
```

## Command Reference

### filter-misassemblies

```bash
mamisa filter-misassemblies [OPTIONS]

Required arguments:
  -a, --assembly PATH          Input assembly FASTA file
  -m, --misassemblies PATH     Directory with *-clipping.txt files

Optional arguments:
  -g, --hq-genomes PATH        Directory with HQ genome files
  -o, --output PATH            Output filtered assembly file
  -l, --min-length INT         Minimum contig/fragment length (default: 2500)
  --mode {remove,split}        Misassembly handling mode (default: split)
  --preserve-hq-with-issues    Keep HQ genomes with misassemblies
  --dry-run                    Show statistics without generating output
  --stats PATH                 Save statistics to TSV file
```

### remove-hq-contigs

```bash
mamisa remove-hq-contigs [OPTIONS]

Required arguments:
  -a, --assembly PATH          Input assembly FASTA file
  --hq-dir PATH | --hq-list PATH   HQ genomes directory or list file

Optional arguments:
  -o, --output PATH            Output filtered assembly file
  -l, --min-length INT         Minimum contig length (default: 0)
  --dry-run                    Show statistics without generating output
  --stats PATH                 Save statistics to TSV file
```

### filter-checkm2

```bash
mamisa filter-checkm2 [OPTIONS]

Required arguments:
  --checkm2-root PATH          Root directory with CheckM2 results
  --genomes-dir PATH           Directory containing genome files
  -o, --output PATH            Output directory

Quality thresholds:
  --hq-comp-min FLOAT          HQ minimum completeness (default: 90)
  --hq-cont-max FLOAT          HQ maximum contamination (default: 5)
  --mq-comp-min FLOAT          MQ minimum completeness (default: 70)
  --mq-cont-max FLOAT          MQ maximum contamination (default: 10)
  --lq-comp-min FLOAT          LQ minimum completeness (default: 50)
  --lq-cont-max FLOAT          LQ maximum contamination (default: 10)

Other options:
  --tiers LIST                 Comma-separated tiers to select (default: HQ,MQ,LQ)
  --extensions LIST            Genome file extensions (default: fa,fasta,fna)
  --symlink | --copy           Link or copy files (default: symlink)
  --dry-run                    Show what would be done
```

### run-gtdbtk

```bash
mamisa run-gtdbtk [OPTIONS]

Required arguments:
  --selected-dir PATH | --genome-dir PATH    Input genomes
  -o, --output PATH            Output directory

Optional arguments:
  --extension STR              Genome file extension (default: fa)
  --cpus INT                   Number of CPUs (default: 1)
  --mash-db PATH               GTDB-Tk mash database
  --tiers LIST                 Tiers to process (default: HQ,MQ,LQ)
  --gtdbtk-args STR            Additional gtdbtk arguments
```

## Understanding Misassembly Detection

MaMISA uses clipping information from read mapping to detect misassemblies:

- **Clipping positions**: Locations where most/all reads are clipped
- **High clipping ratio**: Indicates misjoined sequences
- **Zero coverage**: Regions with no read support (less reliable)

**Key insight**: A clipping ratio of 1.0 (100% of reads clipped) is strong evidence of misassembly.

## HQ Genome Preservation Logic

When using `--preserve-hq-with-issues`:

1. **HQ Clean**: No misassemblies → Removed from assembly (they're already binned)
2. **HQ Problematic**: Has misassemblies → Kept and split (recoverable information)
3. **Non-HQ with misassemblies**: Split at misassembly positions
4. **Non-HQ clean**: Kept intact

This strategy maximizes information retention while maintaining quality.

## Quality Tiers (MiMAG Standard)

| Tier | Completeness | Contamination |
|------|--------------|---------------|
| HQ   | ≥90%         | ≤5%           |
| MQ   | ≥70%         | ≤10%          |
| LQ   | ≥50%         | ≤10%          |
| Fail | <50% or >10% contamination |

## Citation

If you use MaMISA in your research, please cite:

```
[Your citation here]
```

## License

MIT License - see LICENSE file for details

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## Support

- **Issues**: https://github.com/yourusername/mamisa/issues
- **Documentation**: https://github.com/yourusername/mamisa/wiki

## Authors

- Your Name (@yourusername)

## Acknowledgments

- anvi'o for misassembly detection methods
- CheckM2 for genome quality assessment
- GTDB-Tk for taxonomic classification
