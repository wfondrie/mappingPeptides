library(OrgMassSpecR)
library(cleaver)
library(seqinr)
library(tidyverse)

prot <- toupper(unlist(read.fasta("data/triple_KS_Sequenc.fasta", as.string = T)))

dat <- tibble(pep = unlist(cleave(prot, custom = ".", missedCleavages = 100:nchar(prot))))

dat <- dat %>%
    filter(pep != "") %>%
    group_by(pep) %>%
    mutate(mass = MolecularWeight(ConvertPeptide(pep)))