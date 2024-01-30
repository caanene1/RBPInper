# Install
# BiocManager::install("rtracklayer")
# BiocManager::install("GenomicRanges")

#' @title Load gff, gtf, gff3 into GRanges object
#'
#' @description Read a gtf, gff, gtf3 into GRanges for processing
#'
#' @param path Path to the gff file to load
#'
#' @return GRanges object of the gff, gtf
#'
#' @keywords RBPInper, bed annotation
#'
#' @examples See the manuscript.
#'
#' @export
#'
gffload <- function(path){
  gff <- rtracklayer::import(path)
  return(gff)
}


#' @title Load and sort a narrow peaks bed file into dataframe
#'
#' @description Read and sort a narrow peaks bed file into a data frame
#'
#' @param path Path to the bed file to load
#'
#' @return Data frame of the bed file
#'
#' @keywords RBPInper, bed load
#'
#' @examples See the manuscript.
#'
#' @export
#'
bedload <- function(path) {
  # Note pValue - is -log10 in narrow peak bed
  bed <- read.table(path, col.names = c("chrom", "start", "end", "name", "score",
                                        "strand", "signal", "pvalue", "qvalue",
                                        "peak"), header = FALSE)
  bed <- dplyr::arrange(bed, chrom, start)
  bed$chrom <- gsub("chr", "", bed$chrom)
  return(bed)
}


#' @title P-value fusion on data frame
#'
#' @description Apply fusion on data frame
#'
#' @param x Data frame of p-values to combine
#'
#' @param alpha Alpha method to use for the Binomial method
#'
#' @param type Method to use for merging, fisher, binomial, bonferroni, harmonic
#'
#' @return Data frame of fused p-values with threshold labels
#'
#' @keywords RBPInper, P-value combine
#'
#' @examples See the manuscript.
#'
#' @export
#'
fs_matrix <- function(x, alpha = 0.05, type="fisher"){
  # Check that data is a matrix
  # Matrix is used here to force ID to row names, reduce errors

  oup <- apply(x, 1, RBPInper::fusion, alpha=alpha, type=type)
  oup <- as.data.frame(t(oup))
  label <- ifelse(oup <= alpha, "Hit", "NA")

  return(oup)
}


#' @title Prepare macs2 bed file and save it to file
#'
#' @description Prepare macs2 peak call file with different column arrangement
#'
#' @param fie Data frame of the bed columns
#'
#' @param nam The name of the saved bed file
#'
#' @param soure Indicate the source of the bed file to allow column changes
#'
#' @return Bed file saved to the current working directory
#'
#' @keywords RBPInper, bed save
#'
#' @examples See the manuscript.
#'
#' @export
#'
premac2bed <- function(fie, nam="fie.bed", soure=c("base", "encode")){
  if(soure=="encode"){
    fie$V6 <- "*"
    # Using the q-values as the evidence
    fie$V8 <- fie$V9

  } else if(soure=="base"){
    # This numbers are used to allow easy fail not match
    fie <- fie[c(1, 2, 3, 10, 6, 4, 8, 7, 9, 5)]
    fie$length <- "*"
  } else{
    stop("ERROR: This function makes sense with soure set")
  }


  write.table(fie,nam, sep="\t",
              col.names = FALSE,
              row.names = FALSE, quote = FALSE)
}

