library(OrgMassSpecR)
library(cleaver)
library(seqinr)
library(tidyverse)

mapPeptides <- function(fasta, antigen = NULL, minLength = 1, sub = NULL){
    prot <- toupper(unlist(read.fasta(fasta), as.string = T))[1]
    
    mc <- minLength - length(antigen)
    
    if(mx < 0) mc <- 0
    
    dat <- tibble(pep = unlist(cleave(prot, custom = ".", missedCleavages = mc:nchar(prot))))
    
    dat <- dat %>%
        filter(pep != "") %>%
        group_by(pep) %>%
        mutate(mass = MolecularWeight(ConvertPeptide(pep)))
}


prot <- toupper(unlist(read.fasta("data/triple_KS_Sequence.fasta", as.string = T)))
prot <- substring(prot, 1763)

prot <- toupper("aklhind")

dat <- tibble(pep = unlist(cleave(prot, custom = "hind", missedCleavages = 0:nchar(prot))))

dat <- dat %>%
    filter(pep != "") %>%
    group_by(pep) %>%
    mutate(mass = MolecularWeight(ConvertPeptide(pep)))
