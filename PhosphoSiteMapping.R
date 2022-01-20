# 21.09.2021 TM ################################################################


library(tidyverse)
library(seqinr)
# library(biomaRt) # only needed if FASTA is not available


##### load file ################################################################

RES = read_delim("PrecInt_PPep_WideFormat.csv", delim=",")




#######  Get protein phosphosites  ##############################################

# use this code if FASTA is unavailable
# ensembl = useDataset( "hsapiens_gene_ensembl", mart = useMart("ensembl") )
# SequencesAllIsoforms = getSequence(id = RES$uniprot,
#                          type = "uniprotswissprot",
#                          seqType = "peptide",
#                          mart = ensembl)


FASTA = read.fasta( file="mouse_20211121_17090entires.fasta", seqtype="AA", as.string=TRUE) %>% flatten()
FASTA = do.call(rbind, lapply(FASTA, data.frame, stringsAsFactors=FALSE)) %>% rownames_to_column(var="FastaEntryName") %>% as_tibble()
FASTA = mutate(FASTA, FastaUniProt = sapply( str_split(FASTA$FastaEntryName,fixed("|") ), "[", 3 ) )
names(FASTA)[2] = c("FastaProtSeq")

# Add Phospho only sequence, plain sequence and protein sequence to the RES table

ModList <- c("\\(UniMod\\:4\\)")
ReplacementText <- rep(c(""), times  = length(ModList))
names(ReplacementText) <- ModList

library(plyr)

RES["PhosSeq"] = str_replace_all(RES$Modified.Sequence, ReplacementText)
RES["PlainSeq"] = str_replace_all(RES$Modified.Sequence, "[\\(UniMod\\:4\\)  \\(UniMod\\:21\\)]", "")
RES["ProtSeq"] = mapvalues(RES$prot.names, 
                           from=FASTA$FastaUniProt, 
                           to=FASTA$FastaProtSeq)
RES["PhosProtSeq"] = str_replace(RES$ProtSeq, RES$PlainSeq, RES$PhosSeq)

detach("package:plyr")

RES_pSites = RES %>% 
  group_by(PhosProtSeq, Modified.Sequence) %>% 
  do(as.data.frame(str_locate_all(.$PhosProtSeq, "\\(UniMod\\:21\\)"))) %>%
  mutate(start_mod = start - 0:(n()-1)*8 )%>%
  mutate(PhosphoSite = paste(str_sub(PhosProtSeq, start=start-1, end=start-1), as.character(start_mod-1), sep="") ) %>% 
  group_by(Modified.Sequence) %>%
  mutate(PhosphoSite_All = paste0(PhosphoSite, collapse = "|") ) %>%
  filter(!duplicated(Modified.Sequence)) %>%
  ungroup()

# add phosphosite information to RES table

library(plyr)

RES["pSites"] = mapvalues(RES$Modified.Sequence, from=RES_pSites$Modified.Sequence, to=RES_pSites$PhosphoSite_All )

detach("package:plyr")


######## export results (checkpoint) ############################################



write_csv(RES, "RES_wSites.csv")
