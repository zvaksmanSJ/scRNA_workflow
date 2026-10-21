# BaseJumper/Nextflow stage implementation

This directory contains the adapters used by the Cromwell WDL stages. They
prepare an isolated copy of the separately authorized BaseJumper source,
preserve its scientific process definitions, and expose stage-specific
Nextflow entry points and output audits. The WDL run tree receives only the
small scripts in `stage_adapters/`; its runtime paths are supplied through
environment variables and the WDL run manifest.

This is not a standalone public BaseJumper distribution. The upstream
workflow source is excluded from source control because redistribution
permission has not been confirmed. Do not run historical launchers against
completed production results. See `../wdl/README.md` for the supported
orchestration path and its current validation limits.
