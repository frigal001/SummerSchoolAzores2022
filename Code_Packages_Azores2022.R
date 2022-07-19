######################################################################
######################################################################
######---------------FD summer school Azores 2022---------------######
######################################################################
######################################################################

install.packages(c("ade4", "BAT", "vegan", "tidyverse", "gawdis", "devtools", 
                   "ggrepel", "hypervolume", "picante", "GGally", "ggpubr", "nlme", 
                   "emmeans", "ggcorrplot", "devtools", "BiocManager", 
                   "ggnewscale", "ggstar", "FD")) 

BiocManager::install("ggtree")
devtools::install_github("xiangpin/ggtreeExtra")

lapply(c("ade4", "BAT", "vegan", "tidyverse", "gawdis", "devtools", 
         "ggrepel", "hypervolume", "picante", "GGally", "ggpubr", "nlme", 
         "emmeans", "ggcorrplot", "devtools", "BiocManager", 
         "ggnewscale", "ggstar", "ggtreeExtra", "ggtree"), require, character.only = TRUE)

######################################################################
######################################################################
######################################################################