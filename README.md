# PanCancer_Protein_Networks

*Correspondence: Won-Min Song (won-min.song@mssm.edu or wonmin1984@gmail.com)

### Data download
The github page houses companion R codes to reproduce results and figures for pan-cancer protein network study. The input data can be dowloaded at https://www.synapse.org/#!Synapse:syn30554232, or using the following R codes via synapser package: 

library(synapser)

library(synapserutils)

synLogin("synapse_username","password")
files <- synapserutils::syncFromSynapse("syn30554232")

The .tar.gz file must be extracted under "Data" folder.

### Data folder organization: 

- CrossCancer_Module_Overlap: Fisher's Exact Test (FET) results to to evaluate significance of overlap between modules from different cancer types
- CrossCancer_Module_Presevation: Module preservation statistics when modules from different cancer types are compared. 
- DEG: PanCancer TCGA transcriptome-based differential expression analysis between tumor and matched normal (limma results)
- DEP: Differentially expressed protein results between tumor and matched normal samples (limma results) 
- MEGENA: MEGENA outputs for modules and networks. Modules and network edge list are concatenated across different cancer types, and are marked for their respective cancer types. 
- MSigDB: MSigDB v6.2 signatures 
- MutationalDrivers: List of mutational driver genes from TCGA PanCancer studies. 
- QCed_Proteome: Batch, age and gender adjusted protein expression data that were used for network constructions. 
- PRO_Clinical_Signatures.RDS: Pre-compiled protein signatures in .RDS format for quick loading in R workspace


