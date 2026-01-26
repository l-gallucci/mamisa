# Changelog

All notable changes to MaMISA will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-XX-XX

### Added
- Initial release of MaMISA
- `filter-misassemblies` command for intelligent misassembly handling
- `remove-hq-contigs` command for HQ genome removal
- `filter-checkm2` command for quality-based genome filtering
- `run-gtdbtk` command for taxonomy classification
- Comprehensive documentation and examples
- Streaming FASTA processing for memory efficiency
- Dry-run mode for all commands
- Detailed statistics reporting

### Features
- HQ genome preservation with misassembly detection
- Contig splitting at misassembly positions
- CheckM2 integration with MiMAG-compliant tiers
- GTDB-Tk workflow automation
- Multi-format ID extraction and matching
- Progress reporting and logging

### Documentation
- Complete README with workflow examples
- Installation guide
- Command reference
- Citation information

## [Unreleased]

### Planned
- HTML report generation
- Multi-threading support for large assemblies
- Integration with additional QC tools
- Docker container
- Test suite
- Continuous integration
