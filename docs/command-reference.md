# Command reference

Every subcommand also has `-h/--help`. Flags below are the stable ones.

## process-large-contigs

```
Required:
  -a, --assembly PATH         Input assembly FASTA
  -m, --misassemblies PATH    Directory with *-clipping.txt files
  -o, --output-dir PATH       Output directory

Thresholds:
  --max-length INT            Large contig threshold (default: 300000)
  --min-completeness FLOAT    Min completeness for HQ (default: 50)
  --max-contamination FLOAT   Max contamination for HQ (default: 10)

CheckM2:
  --threads INT               Threads (default: 1)
  --skip-checkm2              Skip CheckM2 (requires --checkm2-results)
  --checkm2-results PATH      Existing CheckM2 output directory
```

## check-read-chimeras

```
Required:
  --bam PATH                  Sorted, indexed BAM file
  -o, --output-dir PATH       Output directory
  one classifier (mutually exclusive):
    --kraken2-output PATH     Kraken2 per-read classification output, or
    --kaiju-output PATH       Kaiju per-read output (protein-level, more
                              sensitive on divergent taxa)

Optional:
  --kraken2-report PATH       Kraken2 report (adds taxon names; with Kraken2)
  --kaiju-names PATH          NCBI names.dmp (adds taxon names; with Kaiju)
  --assembly PATH             Assembly FASTA (for contig lengths if not in BAM)
  --min-mapq INT              Min mapping quality (default: 20)
  --min-reads INT             Min reads per contig to report (default: 10)
  --window INT                Sliding window size in bp (default: 10000)
  --window-step INT           Window step in bp (default: 5000)
  --window-threshold INT      Min contig length for windowed analysis (default: 50000)
  --exclude-unclassified      Exclude taxid=0 reads from diversity calculations
  --dry-run
```

## check-chimeras

Combines GC heterogeneity, windowed GC (circular contigs), CheckM2 contamination,
GTDB-Tk warnings and GUNC into a per-bin High/Medium/Low/Clean risk. With
`--gunc-dir` the summary also reports the aggregate `pass.GUNC` / chimera count.

```
Required:
  --bins-dir PATH             Directory containing bin FASTA files

Optional:
  -o, --output PATH           Output TSV (default: chimera_report.tsv)
  --gtdbtk-dir PATH           GTDB-Tk output directory (adds taxonomy signals)
  --checkm2-report PATH       CheckM2 report (adds contamination signal)
  --gunc-dir PATH             GUNC output directory (adds gene-level clade-consistency signal)
  --gc-window INT             GC window size in bp (default: 5000)
  --gc-step INT               GC step in bp (default: 2500)
  --taxonomy-level STR        Taxonomy level for comparison (default: phylum)
  --extensions LIST           Genome extensions (default: fa,fasta,fna)
  --dry-run
```

## classify-clipping

```
Required:
  --bam PATH                  Sorted, indexed BAM file
  -m, --misassemblies PATH    Directory with *-clipping.txt files
  -o, --output PATH           Output classification TSV

Optional:
  --min-mapq INT              Min mapping quality (default: 20)
  --window INT                Read window around each position in bp (default: 500)
  --min-clip-coverage N       Only classify positions with >=N clipped reads
  --taxonomy-windows TSV      chimera_read_windows.tsv from check-read-chimeras
  --taxonomy-flank BP         Extend shift windows by this many bp (default: 5000)
  --assembly FASTA            Assembly FASTA - required for --self-blast
  --self-blast                Self-BLAST repeat detection (needs --assembly + BLAST+)
  --dry-run
```

The `chimeric_join` label: when the soft-clipped parts of reads at a position map
(via their `SA:Z:` supplementary-alignment tag) consistently to one other contig,
the two contigs are joined there. This read-level join signal (Trigodet et al.
2025) is computed automatically from the BAM; see output columns
`sa_partner_contig` / `sa_partner_fraction`.

## filter-misassemblies

```
Required:
  -a, --assembly PATH         Input assembly FASTA
  -m, --misassemblies PATH    Directory with *-clipping.txt files

Optional:
  -g, --hq-genomes PATH       Directory with HQ genome files
  -o, --output PATH           Output filtered assembly
  -l, --min-length INT        Minimum fragment length (default: 2500)
  --mode {remove,split}       Misassembly handling (default: split)
  --preserve-hq-with-issues   Split HQ contigs with clipping instead of removing

Clipping selection (mutually exclusive):
  --safe                      Use ALL anvi'o-reported positions (default)
  --min-clip-coverage N       Only split where >=N reads are clipped

Chimera awareness:
  --chimera-report PATH       TSV from check-chimeras or check-read-chimeras
  --chimera-risk-threshold    Min risk level to act on (default: Medium)

HQ circular handling:
  --split-hq-circular         Force-split HQ circular contigs with clipping
  --hq-circular-min-length    Length threshold for circular detection (default: 200000)

Other:
  --dry-run
  --stats PATH                Save statistics to TSV
```

## remove-hq-contigs

```
Required:
  -a, --assembly PATH         Input assembly FASTA
  --hq-dir PATH               Directory with HQ genome files
    OR
  --hq-list PATH              Text file with HQ contig IDs (one per line)

Optional:
  -o, --output PATH           Output filtered assembly
  -l, --min-length INT        Minimum contig length (default: 0)
  --dry-run
  --stats PATH                Save statistics to TSV
```

## organize-mags  (formerly filter-checkm2)

Splits genomes into `Selected/HQ`, `Selected/MQ`, `Selected/LQ` by MiMAG tier, and
can rename each output `<sample>__<taxon>__<original>` using the sample parsed from
the filename and the GTDB-Tk taxonomy. `filter-checkm2` still works as a deprecated
alias.

```
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

Renaming (sample + taxonomy):
  --rename                    Rename outputs as <sample>__<taxon>__<original>
  --gtdbtk-dir PATH           GTDB-Tk summaries, for the <taxon> part of the name
  --sample-regex RE           Capture group for the sample id
                              (default: ^([^._]+) = token before first . or _)
  --tax-level LEVEL           Rank used in name: domain..species (default: genus;
                              falls back up ranks, NoTax if no GTDB hit)

MIMAG rRNA/tRNA HQ gate:
  --mimag-rna-dir PATH        run-mimag-rna output (mimag_rna_summary.tsv);
                              demotes HQ genomes lacking 5S/16S/23S rRNA or
                              >=18 tRNA amino acids to MQ

Other:
  --tiers LIST                Tiers to select (default: HQ,MQ,LQ)
  --extensions LIST           Genome extensions (default: fa,fasta,fna,...)
  --symlink | --copy          Link or copy files (default: symlink)
  --name-map / --strip-prefix / --strip-suffix / --name-prefix / --name-suffix
  --dry-run
```

MiMAG note: the defaults split the MIMAG medium band into MQ/LQ. For strict MIMAG
tiers (HQ comp>90/cont<5, MQ comp>=50/cont<10, LQ comp<50) run with
`--mq-comp-min 50 --mq-cont-max 10 --lq-comp-min 0 --lq-cont-max 10`. Full MIMAG HQ
also needs rRNA/tRNA (`--mimag-rna-dir`).

## run-checkm2

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --extension STR             Genome file extension (default: fa)
  --threads INT               Threads for CheckM2 (default: 1)
  --database PATH             CheckM2 diamond database (.dmnd FILE)
  --force                     Overwrite existing CheckM2 output
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --checkm2-args STR          Extra arguments passed to checkm2
```

## run-gtdbtk

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --extension STR             Genome file extension (default: fa)
  --cpus INT                  CPUs for GTDB-Tk (default: 1)
  --place-species             Pass --place_species (GTDB-Tk >=2.7): place in the
                              pplacer tree even when skani classifies
  --skip-existing             Skip classify_wf if a gtdbtk.*.summary.tsv already
                              exists in the output (reuse previous results)
  --mash-db PATH              DEPRECATED/ignored - GTDB-Tk removed Mash in v2.5.0
                              (now uses skani)
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --gtdbtk-args STR           Extra arguments passed to gtdbtk
```

GTDB-Tk version note (latest 2.7.2): since v2.5.0 GTDB-Tk screens with skani and no
longer uses Mash, so `--mash_db` is gone. The old `--skip_ani_screen` was replaced
by `--place_species` in v2.7.0. MaMISA parses both the old (`fastani_*`) and new
(`closest_genome_*`) summary column names, so reports from any recent GTDB-Tk
version work.

## run-mimag-rna

Checks the MIMAG rRNA (5S/16S/23S) + tRNA (>=18 amino acids) high-quality criterion
that CheckM2 cannot evaluate, via barrnap + tRNAscan-SE. Writes
`mimag_rna_summary.tsv`; feed it to `organize-mags --mimag-rna-dir` to gate HQ.

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --extension STR             Genome file extension (default: fa)
  --kingdom {bac,arc,euk}     Kingdom for barrnap/tRNAscan-SE (default: bac;
                              with --split-by-domain, the fallback only)
  --split-by-domain           Pick bac/arc per genome from GTDB-Tk (needs
                              --gtdbtk-dir); handles mixed bacteria+archaea sets
  --gtdbtk-dir PATH           GTDB-Tk output dir for --split-by-domain
  --threads INT               Threads for barrnap (default: 1)
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --keep-intermediate         Keep per-genome barrnap/tRNAscan output files
```

Mixed bacteria + archaea: run once with `--split-by-domain --gtdbtk-dir <taxonomy>`.
GTDB-Tk writes one summary per domain (`gtdbtk.bac120.summary.tsv`,
`gtdbtk.ar53.summary.tsv`), so each genome gets barrnap/tRNAscan with the correct
kingdom automatically; genomes absent from GTDB-Tk fall back to `--kingdom`.

## run-gunc

Runs GUNC gene-level chimerism/contamination detection, complementing
`check-chimeras`. Feed the output back with `check-chimeras --gunc-dir`. The
command runs an env-health check and refuses to launch a GUNC broken by an
incompatible pandas/numpy/python (GUNC is old - see installation.md).

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --file-suffix STR           Genome file suffix (default: .fa)
  --threads INT               Threads (default: 1)
  --db-file PATH              GUNC diamond database (.dmnd); else uses $GUNC_DB
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --gunc-args STR             Extra arguments passed to gunc run
```

## setup-workflow

Configures the Snakemake workflow: points `workflow/envs/*.yaml` at this checkout,
writes paths into `workflow/config.yaml`, validates database kinds, and (when env
names are given) sanity-checks them. `--install` builds every per-rule conda env.

```
Optional:
  --repo PATH                 MaMISA checkout (default: auto-detected)
  --workflow-dir PATH         Workflow dir (default: <repo>/workflow)
  --genomes-dir PATH          Input genomes directory (config)
  --genome-ext STR            Genome extension (config)
  --outdir PATH / --threads INT
  --gtdbtk-data DIR           GTDB-Tk data DIRECTORY
  --gunc-db FILE              GUNC .dmnd FILE
  --checkm2-db FILE           CheckM2 .dmnd FILE
  --sample-regex RE / --tax-level LEVEL
  --mamisa-env / --checkm2-env / --gtdbtk-env / --gunc-env / --mimag-env NAME
                              Use these existing envs (sets use_named_envs)
  --install                   Build per-rule conda envs (snakemake)
  --no-check                  Skip the env sanity check
  --dry-run
```

## fetch-databases

Registers databases you already have (validated, written to config) or downloads
missing ones with each tool's own downloader. Path KIND is enforced: GTDB-Tk =
directory, GUNC and CheckM2 = `.dmnd` file.

```
Optional:
  --config PATH               config.yaml to update (default: workflow/config.yaml)
  --db-dir PATH               Where to download (default: ./databases)
  --tools LIST                Which DBs (default: gtdbtk,gunc,checkm2)
  --download                  Download any DB not supplied as an existing path
  --gtdbtk-data DIR           Existing GTDB-Tk data DIRECTORY
  --gunc-db FILE              Existing GUNC .dmnd FILE
  --gunc-db-source {progenomes_2.1,progenomes_3,gtdb_95,gtdb_214}
                              Which GUNC DB to download (GUNC >=1.1.1; default:
                              progenomes_2.1; use gtdb_214 for environmental MAGs)
  --checkm2-db FILE           Existing CheckM2 .dmnd FILE
  --gtdbtk-env / --gunc-env / --checkm2-env NAME   Env to run each downloader in
  --dry-run
```

## check-envs

Sanity-checks existing conda envs: reports each tool's version, whether MaMISA is
importable there, and (for GUNC) pandas/numpy/python compatibility. Prints the
exact install/fix command for anything that fails.

```
Optional (at least one):
  --mamisa-env NAME           Env with the light CLI tools + MaMISA
  --checkm2-env NAME          Env with CheckM2
  --gtdbtk-env NAME           Env with GTDB-Tk
  --gunc-env NAME             Env with GUNC
  --mimag-env NAME            Env with barrnap + tRNAscan-SE
  --repo PATH                 MaMISA checkout used in fix hints
```
