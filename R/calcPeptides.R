library(OrgMassSpecR)
library(cleaver)
library(seqinr)
library(tidyverse)
source("R/mapPeptides.R")

prot <- unlist(read.fasta("data/proteins.fasta", as.string = T))
prot[1] <- substring(prot[1], 1763) # for triple KS
#prot[1] <- substring(prot[1], 500) # for debugging

tripleKS <- mapPeptides(prot[1], antigen = substring(prot[1], nchar(prot[1])-20), minLength = 100)
write_csv(tripleKS, "results/triple_KS.csv")

BurA <- mapPeptides(prot[2], antigen = substring(prot[2], nchar(prot[2])-20), minLength = 100)
write_csv(BurA, "results/BurA.csv")

ZmaK <- mapPeptides(prot[3], antigen = substring(prot[3], nchar(prot[3])-20), minLength = 100)
write_csv(ZmaK, "results/ZmaK.csv")

PKS_NRPS_Cysteine <- mapPeptides(prot[4], antigen = substring(prot[4], nchar(prot[4])-20), minLength = 100)
write_csv(PKS_NRPS_Cysteine, "results/PKS_NRPS_Cysteine")