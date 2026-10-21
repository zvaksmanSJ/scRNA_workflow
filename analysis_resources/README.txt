Analysis resources
==================

setup/setup.sh downloads ReactomePathways.gmt.zip from the official Reactome
release endpoint (https://reactome.org/download/current/ReactomePathways.gmt.zip)
and verifies SHA-256
8c1dbc8578431da5d2d5118262718c60b553a9be3398e93658daa069e4a9afd4 before use.
The archive is excluded from source control. Reactome content is distributed
under CC BY 4.0. The analysis records the source URL, archive hash, and
pathways used in final_results.

Institutional logos are intentionally not included in the public source
package. A user may provide an approved logo separately with RNA_LOGO; without
one, the dashboard uses its text-only design.

The Salmon analysis environment is isolated from the completed ResolveOME RNA
environment and pins Salmon 1.6.0 to match the version recorded in the shared
Ensembl GRCh38 Salmon index.
