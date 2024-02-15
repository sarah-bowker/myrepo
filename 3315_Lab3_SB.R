#1:
library(Biostrings)

read_fa_files <- function(fileToRead) {
  # Read fasta files using the Biostrings reading function
  read_file <- readBStringSet( fileToRead)


  # Find the number of sequences in the file
  num_seq <- length(read_file)

  # Find the length of each sequence in the file
  len_seq <- width(read_file)

  # Store the result of both variables in a list so they could be returned
  seq_result <- list(Number_Seq_Present = num_seq, Length_Each_Seq = len_seq)

  return(seq_result)
}


#2:
library(dplyr)
library(stringr)
rev_comp <- function(seq)
{
  seq <- toupper(seq)
  sqc <- str_split(seq, "")[[1]]
  
  is_RNA <- grepl("U", seq)
 
  
  if(is_RNA)
  {
    sqc_RNA <- if_else(sqc == "U", "A",
                       if_else(sqc =="A","T", 
                               if_else(sqc =="G","C",if_else(sqc =="C","G",sqc))))
    seqrc <- rev(sqc_RNA)
    return (paste(seqrc, collapse = ""))
  }
  
  sqc_DNA <- if_else(sqc == "T", "A",
                  if_else(sqc =="A","T", 
                          if_else(sqc =="G","C",if_else(sqc =="C","G",sqc))))
  
  seqrc <- rev(sqc_DNA)
  return (paste(seqrc, collapse = ""))
  
}

# Check the rev_comp funxtion:
sq_DNA_rc<- rev_comp("AGCTGCAATGGTCACGTAC")
sq_DNA_rc

sq_RNA_rc <- rev_comp("GCUCCAUCUUUAACGCCAU")
sq_RNA_rc

#3:
library(Biostrings)

my_fave_proteins <- character(0)

i = 1
for(i in 1:5)
{
  set.seed(124+i)
  
  AA_seq <- AAString(paste(sample(AA_STANDARD, 500, replace = T),collapse = ""))
  
  my_fave_proteins <- c(my_fave_proteins, as.character(AA_seq))
  
}

#4:
my_fave_proteins_set <- AAStringSet(my_fave_proteins)
my_fave_proteins_set

#5:

i=1
for(i in seq_along(my_fave_proteins_set))
{
  my_fave_proteins_set[[i]][1] <- "M"
}

names(my_fave_proteins_set) <- c("Curious_Protein", "Bored_Protein", "Miserable_Protein", "Amused_Protein", "Loving_Protein")
my_fave_proteins_set

#6:
library(ggplot2)
library(tidyr)

# Calculate AA frequencies for each protein
my_fave_proteins_set_frq <- lapply(my_fave_proteins_set, function(seq) {
  alphabetFrequency(seq, as.prob = TRUE)[alphabetFrequency(seq) > 0] * 100
})

my_fave_proteins_set_frq <- do.call(rbind, my_fave_proteins_set_frq)
my_fave_proteins_set_frq <- as.data.frame(my_fave_proteins_set_frq)

my_fave_proteins_set_frq$AminoAcids <- rownames(my_fave_proteins_set_frq)

my_fave_proteins_set_frq <- tidyr::gather(
  my_fave_proteins_set_frq,
  key = "Protein",
  value = "Frequency",
  -AminoAcids
)

#Box Plot
Plot1 <- ggplot(my_fave_proteins_set_frq, aes(x = Protein, y = Frequency)) +
  geom_boxplot() +
  labs(title = "                    Amino Acid Frequencies in All The Proteins",
       x = "Amino Acids",
       y = "Frequency (%)") +
  theme_minimal()

Plot1

ggsave(plot = Plot1, filename = "Lab_3_plot.pdf", 
       units = "in", width = 11.69, height = 8.27 )


#7:
library(Biostrings)
writeXStringSet(my_fave_proteins_set,"Sarah_fave_prots.fasta")

#8:
library(seqinr)
ripos <- seqinr::choosebank()

bank <- seqinr::choosebank(bank=ripos[str_which(ripos, "16s")], infobank = T)

bank
bacteria <- seqinr::query( listname = "Bacteria", query="TID=2")
bacteria
bacteria2 <- unlist(seqinr::getSequence(bacteria$req[1:200],as.string = T))

#9:
my_fave_bac_set <- DNAStringSet(bacteria2)
names(my_fave_bac_set) <- unlist(lapply(bacteria$req[1:200], '[[',1))
my_fave_bac_set
my_fave_bac_set_rc <- reverseComplement(my_fave_bac_set)
my_fave_bac_set_rc

#10: 
Longest_seq <- max(width(my_fave_bac_set))
Longest_seq

#11:
library(Biostrings)
writeXStringSet(my_fave_bac_set,"Sarah__fave_bac.fa")

