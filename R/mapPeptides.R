library(OrgMassSpecR)
library(cleaver)
library(seqinr)
library(tidyverse)

mapPeptides <- function(seq, antigen = NULL, minLength = 1){
    seq <- toupper(seq)
    
    mc <- minLength - nchar(antigen) - 1
    
    if(mc < 0) mc <- 0
    
    if(is.null(antigen)) antigen <- "."
    
    dat <- tibble(peptide = unlist(cleave(seq, custom = ".", missedCleavages = mc:nchar(seq))))
    
    dat <- dat %>%
        filter(peptide != "", grepl(antigen, peptide)) %>%
        group_by(peptide) %>%
        mutate(mass = MolecularWeight(ConvertPeptide(peptide)))
    
    return(dat)
}
