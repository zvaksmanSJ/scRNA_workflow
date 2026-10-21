# RNA workflow setup

`setup.sh` provisions the conda environments, installs checksum-pinned R
packages, downloads a checksum-verified Reactome pathway archive and Nextflow,
and optionally downloads public Ensembl 104 references and builds the STAR and
Salmon indexes. It writes runtime paths to
`wdl_runtime/runtime_config.json`.

## Requirements

- Linux, conda or mamba, curl, gzip, tar, `cksum`, and `sha256sum`.
- Approximately 100 GB free disk for environments, references, and indexes.
- At least 48 GB RAM for STAR index generation (64 GB recommended).
- Cromwell and Java are separate prerequisites for `main.wdl`; they are not
  installed by this setup.

## Commands

```bash
# Install software only.
bash setup/setup.sh --skip-references

# Install software and build Ensembl 104 STAR/Salmon references.
bash setup/setup.sh --threads 8 --memory-mb 64000

# Install outside the source tree.
bash setup/setup.sh --runtime-dir /path/with/space/wdl_runtime
```

Reactome is downloaded to the ignored `analysis_resources/` asset path expected
by the downstream analysis; source control contains only its README and hash.
The reference builder verifies downloaded Ensembl files against the release
`CHECKSUMS` entries before unpacking, derives annotation helper files from the
same GTF, and refuses to overwrite an incomplete Salmon index. Existing
references and environments are retained; inspect incomplete installations
before removing them.

## Distribution notes

The private BaseJumper Nextflow source is deliberately not installed or
redistributed by these scripts. The WDL runner requires a separately obtained,
authorized checkout. The downloaded package archives and environment
manifests are stored under the selected runtime directory, not in source
control. A clean install and index build have not yet been executed in this
release-preparation pass.
