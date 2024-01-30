#' @title Prepare bed file for RBPInper
#'
#' @description Annotate narrow peak bed files using gtf,gff,gff3 annotation file
#'
#' @param bed Path to bed file to annotate and summaries "can take a vector of paths"
#'
#' @param gtf Path to gtf, gff, gff3 annotation file "same genome version as bed"
#'
#' @param fTarget The gtf feature to use for summarising the peaks, default "gene".
#' gtf used must have the type column
#'
#' @param islog10 Is bed file p-value -log10 transformed, default TRUE
#'
#' @param mMerge Method to merge peak P-values to gene p-values
#'
#' @return Data frame of summaries p-value
#'
#' @keywords RBPInper, bed annotate
#'
#' @examples See the manuscript.
#'
#' @export
#'
prebed <- function(bed, gtf, fTarget="gene", islog10=TRUE, mMerge="Bonferroni"){
  library(GenomicRanges)

 # Load gtf once for speed
 gtf <- RBPInper::gffload(gtf)

 # Output hold
 output <- list()

 # Process each bed in turn
 for(bb in bed){
   print(paste("Processing....", bb, sep = ""))
   # Set the name
   bb_name <- basename(bb)
   bb_name <- substr(bb_name, 1, nchar(bb_name) - 4)

   # Load the ith bed into GRnages directly
   bedgr <- makeGRangesFromDataFrame(RBPInper::bedload(bb),
                                     keep.extra.columns=TRUE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     seqnames.field="chrom",
                                     start.field="start",
                                     end.field="end",
                                     strand.field="strand",
                                     starts.in.df.are.0based=TRUE)

  # Search overlap and get the index
  bedgtf <- findOverlaps(bedgr, gtf)

  # Bind the subject info with the query location
  gtfindex <- cbind(as.data.frame(gtf[subjectHits(bedgtf)]),
                    as.data.frame(bedgtf))

  # Filter by the required feature
  gtfindex <- gtfindex[gtfindex$type == fTarget, ]

  # Remove duplicates in the query index (same gene if any)
  gtfindex <- gtfindex[!duplicated(gtfindex$queryHits), ]

  # Get columns to merge
  gtfindex <- gtfindex[c("gene_id", "queryHits")]

  # Convert the bed file GRanges to table with meta data
  # None overlaps are unimportant
  bedgr <- as.data.frame(bedgr)
  bedgr$index <- as.integer(rownames(bedgr))
  bedgr <- merge(bedgr, gtfindex, by.x = "index",
                  by.y = "queryHits", all.x=F)

  # Merge Peak P-values
  if(islog10){
    bedgr$pvalue <- 10^(-bedgr$pvalue)
  }

  # Get the output, it is merged on Ensemble ID for compatibility
  out <- aggregate(bedgr["pvalue"], by = list(bedgr[["gene_id"]]),
                   FUN = RBPInper::fusion, type=mMerge)
  names(out) <- c("gene_id", bb_name)

  # Final append
  output[[bb_name]] <- out

 }

 # Merge into one data frame
 output <- Reduce(function(x, y) merge(x, y,
                                       by = "gene_id",
                                       all = TRUE),
                  output)

 # Set NA values to 1
 output[is.na(output)] <- 1

 # Next get the gene symbol column or anything else
 mda <- mcols(gtf)
 mda <- as.data.frame(mda)
 mda <- mda[mda$type == fTarget, ][c("gene_id", "gene_name")]
 output <- merge(mda, output, by="gene_id")

 return(output)

}



