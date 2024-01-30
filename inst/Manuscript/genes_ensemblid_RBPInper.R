#_______________________________________
#FUNCTION TO EXTRACT POSSIBLE GENE NAMES====

library(qdapRegex)

extract_genes <- function(x,pt, remove = character(0))
{
  #pre-defined pattern, what I want to subset with space for easier text manipulation
  for (i in pt) {
    x <- gsub(i, " ", x)
  }
  #splitting sentences into separate words
  words <- unlist(strsplit(x, split = " "))
  #getting words in only capital letters or capital letters and numbers
  gene_words <- words[str_detect(words, "\\b[A-Z]+\\b|\\b[A-Z0-9]+\\b")]

  #removes period at the end of the words
  #I put it separately because if it was together (with the line below),
  #when there is a word in parentheses at the end of the sentence it would just remove the period and not the () itself
  gene_words_a <- gsub("\\.$", "", gene_words)
  #removes () and []
  cleaned <- gsub("^\\(|\\)$|^\\[|\\]$", "", gene_words_a)


  #removing pre-defined words
  filtered_words <- cleaned[!cleaned %in% remove]
  #removes numbers (no numbers)
  nono <- rm_number(filtered_words)
  #removes words with less than two characters
  end <- nono[nchar(nono) > 2]

  return(paste(end, collapse = " "))
}


#________________
#________________

#applying the the function and extracting possible gene names
res_all_genes$Gene_Words <- apply(res_all_genes, 1, function(column) {
  extract_genes(column["Title"], pt=c("/", "-", ":", "·", ",", '"'), remove = c("RNA", "DNA", "SFPQ", "PSF", "3D"))
})

#creating a table with frequency of filtered words
gene_table <- table(unlist(strsplit(res_all_genes$Gene_Words, split = " ")))
unique_gene_table <- data.frame(symbol = names(gene_table), Freq = as.vector(gene_table))

#_genes to ensemblID___====
unique_gene_table$ensemblID <- mapIds(org.Hs.eg.db,
                                      keys = unique_gene_table$symbol,
                                      keytype = "SYMBOL",
                                      column = "ENSEMBL")

lit_genes <- unique_gene_table[!is.na(unique_gene_table$ensemblID),]

write.csv(lit_genes, file="lit_extract_genes.csv", row.names = F)
