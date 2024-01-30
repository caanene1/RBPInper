# Class
base.class = "data.frame"

# Set the class
setClass("RBPInper", slots = list(evidence = base.class,
                                  meta = base.class,
                                  cells = "character",
                                  gene_id = "character",
                                  id = "character",
                                  information = "list",
                                 L1.result = base.class,
                                 L2.result = base.class))


# Generics
setGeneric("load.inper", function(x, evi, info, gene.col="gene_id", info.id="ID",
                            cell.col="cell") standardGeneric("load.inper"))

setGeneric("integrate.inper", function(x, L1="bonfe", L2="fishe",
                                 p=0.05) standardGeneric("integrate.inper"))

# Methods
setMethod("load.inper", "RBPInper", function(x, evi, info, gene.col, info.id,
                                       cell.col) {
  # Sanity check
  if(!is(object = x, class2 = "RBPInper")) {
    stop('Object must be of class "RBPInper"\n')
  }

  # Load the evidence
  pvals <- names(evi)[sapply(evi, function(x) is.numeric(x) || is.integer(x))]

  x@meta <- evi
  rownames(x@meta) <- paste(rownames(evi), evi[[gene.col]], sep = "_")
  x@evidence <- x@meta[pvals]

  # Load the info
  x@cells <- unique(info[[cell.col]])
  x@information <- lapply(x@cells, function(cc){
    subinfo <- info[info[[cell.col]] %in% cc, ]
  })
  names(x@information) <- x@cells
  x@id <- info.id
  x@gene_id <- gene.col

  x
})


setMethod("integrate.inper", "RBPInper", function(x, L1, L2, p) {
  # Sanity check
  if(!is(object = x, class2 = "RBPInper")) {
    stop('Object must be of class "RBPInper"\n')
  }

  inteed <- lapply(x@cells, function(ci){
  idinteed <- x@evidence[x@information[[ci]][[x@id]]]
  idinteed[[ci]] <- apply(idinteed, 1, RBPInper::fusion, alpha=p, type=L1)
  idinteed <- idinteed[ci]
  })


  x@L1.result <- do.call(cbind, inteed)
  res <- x@L1.result
  x@L2.result <- x@L1.result
  x@L2.result$global <- apply(x@L2.result, 1, RBPInper::fusion,
                              alpha=p, type=L2)

  x@L2.result <- cbind(x@meta[x@gene_id], x@L2.result["global"])
  x@L2.result$adjP <- p.adjust(x@L2.result$global, method = "BH")
  x@L2.result$call <- ifelse(x@L2.result$adjP <= p, "Hit", "")

  x
})


#' @title Integrate RBP interaction profiles
#'
#' @description Meta-analysis of RNA binding protein interaction profiles
#'
#' @param evi Data frame of p-value evidence, genes in row and data in columns
#'
#' @param info Data frame of information file, at least data ID and cell type include
#'
#' @param gene.col Column name in "evi" p-value matrix with gene IDs
#'
#' @param info.id Column name in "info" information matrix with matching sample IDs
#'
#' @param cell.col Column name in "info" information matrix with cell type name
#'
#' @param L1 Method for p-value merging at cell type level, default "bonf"
#'
#' @param L2 Method for p-value merging at global level, default "fisher"
#'
#' @param p P-value threshold for calls
#'
#' @return RBPInper object with cell type level results and global results
#'
#' @keywords RBPInper, p-value merging, RNA binding protein
#'
#' @examples
#'
#' @export
#'
rbpinper.run <- function(evi, info, gene.col="gene_id",
                         info.id="ID", cell.col="cell",
                         L1="bonfe", L2="fishe", p=0.05){

  res <- new("RBPInper")
  res <- load.inper(res, evi=evi, info=info, gene.col=gene.col,
              info.id=info.id, cell.col)

  res <- integrate.inper(res, L1=L1, L2=L2, p=p)

  print("Integration complete")
 return(res)
}
