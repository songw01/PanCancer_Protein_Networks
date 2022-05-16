# PanCancer_Protein_Networks

*Correspondence: Won-Min Song (won-min.song@mssm.edu or wonmin1984@gmail.com)

########## 
The github page houses companion R codes to reproduce results and figures for pan-cancer protein network study. The input data can be dowloaded at https://www.synapse.org/#!Synapse:syn30554232, or using the following R codes via synapser package: 

library(synapser)
library(synapserutils)

synLogin("synapse_username","password")
files <- synapserutils::syncFromSynapse("syn30554232")

The .tar.gz file must be extracted under "Data" folder.


