### Run the final analysis
predata <- read.csv("~/Documents/GitHub/RBPInper/inst/data/pre_integration.csv")
infodata <- read.csv("~/Documents/GitHub/RBPInper/inst/data/info.csv")

infoRNA <- infodata[infodata$method %in% c("eCLIP", "RIP-Seq", "RNA-Seq"), ]
infoDNA <- infodata[infodata$method %in% c("ChIP", "RNA-Seq"), ]

## First get the normal calculations
resultRNA <- RBPInper::rbpinper.run(evi=predata, info = infoRNA)
resultRNAP <- RBPInper::rbpinper.run(evi=predata, info = infoRNA,
                                     penalise=T, pen="RNA-Seq",
                                     method.col="method")
resultRNA <- resultRNA@L2.result[c(1, 4)]
names(resultRNAP)
resultRNAP <- resultRNAP@L2.result[c(1, 4)]
names(resultRNAP) <- c("gene_id.pen", "call.pen")
both <- cbind(resultRNA, resultRNAP)
both$side <- ifelse(both$call == "Hit" & both$call.pen == "Hit", "both",
                    ifelse(both$call != "Hit" & both$call.pen == "Hit",
                           "pen",
                           ifelse(both$call == "Hit" & both$call.pen != "Hit", "unpen", "")))
## Load the meta file
meta <- read.csv("~/Desktop/RBPInper/raw files/meta_results_indviduals.csv")
meta <- meta[c(1, 2, 13:15)]
both <- merge(both, meta, by="gene_id")
both2 <- both[both$side %in% c("unpen","pen"),]
write.csv(both, file = "effect_penalisation.csv", row.names = F)


## Get and count the interactors
mergdf <- resultRNA@L2.result[c(1, 4)]
nae <- infoRNA$ID
mergls <- cbind(predata[1], predata[nae])

for(i in nae){
  mergls[i] <- ifelse(mergls[[i]] <= 0.05, "Hit", "")
}
##
mergls$Union <- rowSums(mergls == "Hit")
mergls$Union <- ifelse(mergls$Union >= 1, "Hit", "")
mergls <- mergls[c("gene_id", "Union")]
mergdf <- merge(mergdf, mergls, by="gene_id")


### Do simulations
sim <- names(predata[3:17])

null <- function(nam=infoRNA$ID,
                 id=names(predata[1:2]),
                 data = predata,
                 percent = c(10, 20, 30, 40, 50)){

# Set seed
set.seed(123)

# Create the numbers for the percentage
percent_to_n <- lapply(percent, function(p){
    return(ceiling((p/100) * length(nam)))
  })

# Create the combinations
combinations <- lapply(percent_to_n, function(n){
  combn(nam, n, simplify = FALSE)})


 shuffle_run <- lapply(combinations, function(ll){

   elmres <- list()

   counter <- 0

     for (lll in ll) {

       shffuled <- data[c(id, nam)]

       # Shuffle the targets
       for (ii in lll) {
         shffuled[[ii]] <- sample(data[[ii]])
       }

       # Run RBPInper
       res <- RBPInper::rbpinper.run(evi=shffuled, info = infoRNA)
       res <- res@L2.result
       res <- res[c(names(res[1]), "call")]
       names(res)[2] <- paste("RBPInper_", counter)


       # Run Union
       for(i in nam){
         shffuled[i] <- ifelse(shffuled[[i]] <= 0.05, "Hit", "")
       }

       shffuled[paste("Union_", counter)] <- rowSums(shffuled == "Hit")
       shffuled[paste("Union_", counter)] <- ifelse(shffuled$Union >= 1, "Hit", "")
       shffuled <- shffuled[c(names(shffuled[1]), paste("Union_", counter))]

       out <- merge(res, shffuled, by="gene_id")

       ## Merge
       # Gene id is hard coded here because of the data
       elmres[[paste("run_", counter)]] <- out
       counter = counter+1

     }

   return(elmres)
 })


 save(shuffle_run, mergdf, file="shuffle_run.RData")
 #
 load("shuffle_run.RData")

 bind <- lapply(shuffle_run, function(llt){
 out_ls <- lapply(llt, function(gg){
   base <- merge(mergdf, gg, by = "gene_id")


   rbpinper <- rowSums(base[c(2, 4)] == "Hit")
   union <- rowSums(base[c(3, 5)] == "Hit")
   ###
   TP_r <- length(rbpinper[rbpinper == 2])
   TP_u <- length(union[union == 2])

   ###
   TN_r <- length(rbpinper[rbpinper == 0])
   TN_u <- length(rbpinper[rbpinper == 0])

   ##
   FP_r <- sum(ifelse(base[[c(4)]] =="Hit" & base[[c(2)]] !="Hit",
                      1, 0))
   FP_u <- sum(ifelse(base[[c(5)]] =="Hit" & base[[c(3)]] !="Hit",
                      1, 0))

   ##
   FN_r <- sum(ifelse(base[[c(2)]] =="Hit" & base[[c(4)]] !="Hit",
                      1, 0))
   FN_u <- sum(ifelse(base[[c(3)]] =="Hit" & base[[c(5)]] !="Hit",
                      1, 0))

   ##
   sens_r <- TP_r / (TP_r + FN_r)
   sens_u <- TP_u / (TP_u + FN_u)

   ##
   spec_r <- TN_r/ (TN_r + FP_r)
   spec_u <- TN_u/ (TN_u + FP_u)


   out <- data.frame(Sensitivity = c(sens_r, sens_u),
                     Specificity = c(spec_r, spec_u),
                     Group = c("RBPInper", "Union"))

   return(out)
 } )

 out_ls <- do.call(rbind, out_ls)

 boxplot(Sensitivity~Group, data=out_ls, main="Sensitivity",
         xlab="", ylab="")

 print(aggregate(out_ls$Sensitivity, list(out_ls$Group), mean))


 boxplot(Specificity~Group, data=out_ls, main="Specificity",
         xlab="", ylab="")
 print("Specificity")
 print(aggregate(out_ls$Specificity, list(out_ls$Group), mean))

 return(out_ls)

 })

 bindindex <- lapply(1:length(bind), function(x){
   get <- bind[[x]]
   get$n_null <- x
   get$per_null <- x*10
   return(get)
 })

 bindindex <- do.call(rbind, bindindex)
 write.csv(bindindex, file = "performance_eval_noise.csv", row.names = F)

 library(ggpubr)


 # Use position = position_dodge()

 # Create line plots of means
 pp <- ggline( bindindex, x = "n_null", y = "Specificity",
        add = c("mean_se"),
        color = "Group", palette = c("darkred", "#00AFBB")) +
         geom_hline(yintercept = 0.8, linetype='dotted', col = "grey30")

 ggpar(pp,legend = "right")
 ### Split it up better

 pp2 <- ggline( bindindex, x = "n_null", y = "Sensitivity",
         add = c("mean_se"),
         color = "Group", palette = c("darkred", "#00AFBB")) +
   geom_hline(yintercept = 0.8, linetype='dotted', col = "grey30")
 ggpar(pp2,legend = "right")

 }



