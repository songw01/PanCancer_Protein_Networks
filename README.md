# PanCancer_Protein_Networks

*Correspondence: Won-Min Song (won-min.song@mssm.edu or wonmin1984@gmail.com)

### Data download
The github page houses companion R codes to reproduce results and figures for pan-cancer protein network study. The input data can be dowloaded at https://www.synapse.org/#!Synapse:syn30595684, or using the following R codes via synapser package: 

library(synapser)

library(synapserutils)

synLogin("synapse_username","password")

files <- synapserutils::syncFromSynapse("syn30595684")

The .tar.gz file must be extracted under "Data" folder.

### Data folder organization: 

- **CrossCancer_Module_Overlap**: Fisher's Exact Test (FET) results to to evaluate significance of overlap between modules from different cancer types
- **CrossCancer_Module_Presevation**: Module preservation statistics when modules from different cancer types are compared. 
- **DEG**: PanCancer TCGA transcriptome-based differential expression analysis between tumor and matched normal (limma results)
- **DEP**: Differentially expressed protein results between tumor and matched normal samples (limma results) 
- **MEGENA**: MEGENA outputs for modules and networks. Modules and network edge list are concatenated across different cancer types, and are marked for their respective cancer types. Subfolders include *ModuleRanking* (ranked list of modules by enrichment of DEP signatures per cancer), *Individual_MEGENA_Proteome* (MEGENA output per cancer type from proteome data), *hubs* (hub statistics per cancer type). 
- **KD_stratification**: holds enrichment of CRISPRi signature from Archilles data set (FDR < 0.05), or LINCS gene perturbation signatures in network hubs from proteome and transcriptome data (*KD_Stratification_Evaluaton.Transcriptome/Proteome.txt*). Also holds broader list of network key driver list (*summarized_KDA_results.txt*).  
- **LINCS_validation**: holds .RDS files for GSEA analysis results (from fgsea R package) to evaluate enrichments of LINCS gene perturbation signatures in respective target gene neighborhoods. Proteome/Transcriptome subfolders hold the results by using proteome-/transcriptome-network neighborhoods. 
- **MSigDB**: MSigDB v6.2 signatures 
- **MutationalDrivers**: List of mutational driver genes from TCGA PanCancer studies. 
- **QCed_Proteome**: Batch, age and gender adjusted protein expression data that were used for network constructions. 
- **PRO_Clinical_Signatures.RDS**: Pre-compiled protein signatures in .RDS format for quick loading in R workspace



