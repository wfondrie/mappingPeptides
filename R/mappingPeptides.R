library(OrgMassSpecR)
library(cleaver)
library(seqinr)
library(tidyverse)

prot <- toupper(unlist(read.fasta("data/triple_KS_Sequence.fasta", as.string = T)))
prot <- substring(prot, 1763)

dat <- tibble(pep = unlist(cleave(prot, custom = ".", missedCleavages = 100:nchar(prot))))

dat <- dat %>%
    filter(pep != "") %>%
    group_by(pep) %>%
    mutate(mass = MolecularWeight(ConvertPeptide(pep)))
