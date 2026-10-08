#!/usr/bin/env bash
set -euo pipefail

# Exact IntaRNA commands reported as used for the NE-NRVM project.
# These are preserved as provenance. The repository maintainer may parameterize
# paths/sample inputs, but should preserve the analytical settings:
#   --mode=M
#   --outMode=C
#
# Downstream manuscript networks retain interactions with MFE <= -20 kcal/mol.

# circRNA-miRNA
IntaRNA --mode=M -t circRNA_UP_6h_all_records.fa -q miRNAs_DOWN_6hr.fa --outMode=C --out Interact_circUP_miDOWN_6h_all_modeM.csv
IntaRNA --mode=M -t circRNA_DOWN_6h_all_records.fa -q miRNAs_UP_6hr.fa --outMode=C --out Interact_circDOWN_miUP_6h_all_modeM.csv
IntaRNA --mode=M -t circRNA_UP_24h_all_records.fa -q miRNAs_DOWN_24hr.fa --outMode=C --out Interact_circUP_miDOWN_24h_all_modeM.csv
IntaRNA --mode=M -t circRNA_DOWN_24h_all_records.fa -q miRNAs_UP_24hr.fa --outMode=C --out Interact_circDOWN_miUP_24h_all_modeM.csv

# lncRNA-miRNA
IntaRNA --mode=M -t Candidates_lncRNA_6h_UP.fa -q miRNAs_DOWN_6hr.fa --outMode=C --out Interact_LncUP_miDOWN_6h_modeM.csv
IntaRNA --mode=M -t Candidates_lncRNA_6h_DOWN.fa -q miRNAs_UP_6hr.fa --outMode=C --out Interact_LncDOWN_miUP_6h_modeM.csv
IntaRNA --mode=M -t Candidates_lncRNA_24h_UP.fa -q miRNAs_DOWN_24hr.fa --outMode=C --out Interact_LncUP_miDOWN_24h_modeM.csv
IntaRNA --mode=M -t Candidates_lncRNA_24h_DOWN.fa -q miRNAs_UP_24hr.fa --outMode=C --out Interact_LncDOWN_miUP_24h_modeM.csv
