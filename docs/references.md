# References

MaMISA is an orchestrator: it drives established tools and signals. Please cite
the underlying methods you actually used.

## Misassembly detection (anvi'o)

The misassembly branch of MaMISA (`classify-clipping`, `filter-misassemblies`,
`process-large-contigs`, `check-zero-coverage`) consumes the output of
`anvi-script-find-misassemblies`, which reports soft-clipping and zero-coverage
positions from long reads mapped back to their own assembly. The `chimeric_join`
label in `classify-clipping` uses the read-level supplementary-alignment (`SA:Z:`)
join signal described in the same work.

- Trigodet F. et al. (2025). *Assemblies of long-read metagenomes suffer from
  diverse forms of errors.* bioRxiv. doi:10.1101/2025.04.22.649783 —
  https://doi.org/10.1101/2025.04.22.649783
- Published version: *Troubleshooting common errors in assemblies of long-read
  metagenomes.* Nature Biotechnology (2025). doi:10.1038/s41587-025-02971-8
- Reproducible workflow: https://merenlab.org/data/benchmarking-long-read-assemblers/
- Program: `anvi-script-find-misassemblies` —
  https://anvio.org/help/9/programs/anvi-script-find-misassemblies/
- anvi'o: Eren A.M. et al. (2021). *Community-led, integrated, reproducible
  multi-omics with anvi'o.* Nature Microbiology 6:3-6.
  doi:10.1038/s41564-020-00834-3

## Quality tiers (MIMAG)

HQ/MQ/LQ tiers in `organize-mags` follow the MIMAG standard. Full HQ additionally
requires 5S/16S/23S rRNA and tRNAs for at least 18 amino acids, which
`run-mimag-rna` evaluates (CheckM2 does not).

- Bowers R.M. et al. (2017). *Minimum information about a single amplified genome
  (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea.*
  Nature Biotechnology 35:725-731. doi:10.1038/nbt.3893

## Tools driven by MaMISA

- CheckM2: Chklovski A. et al. (2023). Nature Methods 20:1203-1212.
  doi:10.1038/s41592-023-01940-w
- GTDB-Tk: Chaumeil P.-A. et al. (2022). Bioinformatics 38:5315-5316.
  doi:10.1093/bioinformatics/btac672 — v2.5+ screens with skani (below)
- skani: Shaw J. & Yu Y.W. (2023). Nature Methods 20:1661-1665.
  doi:10.1038/s41592-023-02018-3
- GUNC: Orakov A. et al. (2021). Genome Biology 22:178.
  doi:10.1186/s13059-021-02393-0
- barrnap: Seemann T. — https://github.com/tseemann/barrnap
- tRNAscan-SE 2.0: Chan P.P. et al. (2021). Nucleic Acids Research 49:9077-9096.
  doi:10.1093/nar/gkab688
- Kraken2: Wood D.E. et al. (2019). Genome Biology 20:257.
  doi:10.1186/s13059-019-1891-0
- Kaiju: Menzel P. et al. (2016). *Fast and sensitive taxonomic classification for
  metagenomics with Kaiju.* Nature Communications 7:11257.
  doi:10.1038/ncomms11257
