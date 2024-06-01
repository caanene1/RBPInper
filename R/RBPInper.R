# Class
base.class = "data.frame"

# Set the class
setClass("RBPInper", slots = list(evidence = base.class,
                                  meta = base.class,
                                  cell.col = "character",
                                  cells = "character",
                                  gene_id = "character",
                                  id = "character",
                                  pen.group="character",
                                  method.column=base.class,
                                  penalised = "list",
                                  information = "list",
                                 L1.result = base.class,
                                 L2.result = base.class))


# Generics
setGeneric("load.inper", function(x, evi, info, gene.col="gene_id", info.id="ID",
                            cell.col="cell", pen="no",
                            method.col ="info.method") standardGeneric("load.inper"))

setGeneric("integrate.inper", function(x, L1="bonfe", L2="fishe",
                                 p=0.05, penalise=F) standardGeneric("integrate.inper"))

# Methods
setMethod("load.inper", "RBPInper", function(x, evi, info, gene.col, info.id,
                                       cell.col, pen, method.col) {
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
  x@cell.col <- cell.col

  # Load the method to penalize
  if(pen != "no"){
    x@method.column <- info[info[[method.col]] %in% pen, ]
    x@pen.group <- unique(x@method.column[[x@cell.col]])
    x@penalised <- lapply(x@pen.group, function(cp){
      subinfo <- x@method.column[x@method.column[[x@cell.col]] %in% cp, ]
    })
    names(x@penalised) <- x@pen.group
  }

  x
})


setMethod("integrate.inper", "RBPInper", function(x, L1, L2, p, penalise) {
  # Sanity check
  if(!is(object = x, class2 = "RBPInper")) {
    stop('Object must be of class "RBPInper"\n')
  }

  if(penalise){

  x@information <- lapply(x@information, function(tab){
    tabb <- tab[!tab[[x@id]] %in% x@method.column[[x@id]], ]
    if(nrow(tabb) >= 1){
      return(tabb)
    }else{
    }
  })

  x@information <- x@information[sapply(x@information, function(x) !is.null(x))]

  x@cells <- unique(unlist(sapply(x@information, function(tab2){
    tab2[[x@cell.col]] })))
  }

  ##
  inteed <- lapply(x@cells, function(ci){
  idinteed <- x@evidence[x@information[[ci]][[x@id]]]
  idinteed[[ci]] <- apply(idinteed, 1, RBPInper::fusion, alpha=p, type=L1)
  idinteed <- idinteed[ci]
  })

  x@L1.result <- do.call(cbind, inteed)

 # Do penalty
  if(penalise){
    print("Runing the penalised version")
    inteedP <- lapply(x@pen.group, function(ci){
      idinteedP <- x@evidence[x@penalised[[ci]][[x@id]]]

      idinteedP[[ci]] <- apply(idinteedP, 1, RBPInper::fusion, alpha=p, type=L1)
      idinteedP <- idinteedP[ci]
    })
   L1.pen <- do.call(cbind, inteedP)
   L1.pen$pen <- apply(L1.pen, 1, RBPInper::fusion,
                       alpha=p, type=L2)
   x@L1.result$pen <- L1.pen$pen
  }



  # res <- x@L1.result
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
#' @param penalise Indicate if a method should be penalized. Sub-process for DE RNA-Seq.
#'
#' @param pen If penalise TRUE, then indicate the method to be penalized found in method.col below.
#'            This must match information in info file.
#'
#' @param method.col If penalise TRUE and pen given, then indicate the column in the info file that match pen.
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
                         L1="bonfe", L2="fishe", p=0.05, penalise=F,
                         pen="no", method.col ="info.method"){

  res <- new("RBPInper")
  res <- load.inper(res, evi=evi, info=info, gene.col=gene.col,
              info.id=info.id, cell.col=cell.col, pen=pen, method.col=method.col)

  res <- integrate.inper(res, L1=L1, L2=L2, p=p, penalise=penalise)

  print("Integration complete")
 return(res)
}
