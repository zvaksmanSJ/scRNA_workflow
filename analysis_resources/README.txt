Analysis resources
==================

For new runs, setup/setup.sh downloads Reactome release 98 from
https://reactome.org/download/98/ReactomePathways.gmt.zip
as ReactomePathways.v98.gmt.zip
and verifies SHA-256
6056fcd95dc8f09740c8bb64c2765bea1386ce91b9b19cd10bcd45d30451805c before use.
The older unversioned archive and completed results are not overwritten.
Pathway results from different releases may differ.
The archive is excluded from source control. Reactome content is distributed
under CC BY 4.0. The analysis records the source URL, archive hash, and
pathways used in final_results.

Institutional logos are intentionally not included in the public source
package. A user may provide an approved logo separately with RNA_LOGO; without
one, the dashboard uses its text-only design.

The Salmon analysis environment is isolated from the completed ResolveOME RNA
environment and pins Salmon 1.6.0 to match the version recorded in the shared
Ensembl GRCh38 Salmon index.
