# This function allows multiple differential expression methods to be implemented. 

# DESeq2 ------------------------------------------------------------------

DEA_DESeq2 <- function(counts, designTable, model,fdr = 0.05, lfc  = 1,
                       shrinkage = F,geneNames=NULL){

cat("Starting differential expression analysis with DESeq2 \n")
dds <-DESeqDataSetFromMatrix(countData = countData_dgeList, colData = designTable, design= model)
dds <- DESeq(dds)
results <- DESeq2::results(dds, alpha =  0.05)
if(shrinkage){
cat("Shrinkage of LFCs \n")
coeffContrast <- tail(resultsNames(dds),n =1 )
results <- lfcShrink(dds, coef=coeffContrast,res=results) 
}
if(!is.null(geneNames)){
  cat("Adding gene names\n")
  results <- merge(as.data.frame(results),unique(geneNames[,c(1,5)]), by.x=0,by.y=1,all.x=T)
  names(results)[1] <- "Ensembl.ID"
}
results <- as.data.frame(results) %>% mutate(sig = case_when(abs(results$log2FoldChange) > lfc & results$padj < fdr ~ "sig", 
                                             .default =" not sig")) %>% rename(., all_of(c(FDR = "padj",LFC = "log2FoldChange")))
cat(paste("Number of DEGs: ", table(results$sig)[2], "\n"))
cat("Done!\n")
return(results)
}


# Limma -------------------------------------------------------------------

DEA_limma <- function(dge, designTable, model, fdr = 0.05, lfc = 1, 
                      voom=F, voomWeights=F, randomEffect=F, norm = "TMM", 
                      plots=F, quantile=F, geneNames=NULL){
  cat("Starting differential expression analysis with Limma \n")
  
  design <- model.matrix(model, data =designTable)
  
  dge <- calcNormFactors(dge, method = norm)
  logCPM <- cpm(dge,log =T)
  
  if(voom){
    cat("Begin voom \n")
    if(randomEffect){
      cat("Random effect included \n") # random effect must be listed first 
      modelR <-  update(model, paste(" ~.-",all.vars(model)[1],sep=""))
      design <- model.matrix(modelR, data=designTable)
      v <- voom(dge, design , plot= T)
      corfit <- duplicateCorrelation(v,design,block=designTable$patient)
      v <- voom(dge,  design, plot =T,block=batch$patient,correlation=corfit$consensus.correlation)
      corfit <- duplicateCorrelation(v,design,block=designTable$patient)
      fit <- lmFit(v,design,block=designTable$patient,correlation=corfit$consensus.correlation)
    }else{
      voom <- voom(dge, design, plot =T)
     fit <- lmFit(voom, design)
    }}
    if(voomWeights){
      cat("Begin voom with quality weights \n")
      if(randomEffect){
        cat("Random effect included \n") # random effect must be listed first 
        modelR <-  update(model, paste(" ~.-",all.vars(model)[1],sep=""))
        design <- model.matrix(modelR, data=designTable)
      vWeights <- voomWithQualityWeights(dge, design, plot = T)
      corfitWeights <- duplicateCorrelation(vWeights,design,block=designTable$patient)
      vWeights <- voomWithQualityWeights(dge, design, plot = T, correlation= corfitWeights$consensus.correlation)
      corfitWeights <- duplicateCorrelation(vWeights,design,block=designTable$patient)
      fit <- lmFit(vWeights, designR, plot=T, correlation=corfitWeights$consensus.correlation)
      }else{
      vWeights <- voomWithQualityWeights(dge, design, plot =T)  
      fit <- lmFit(vWeights, design)
      }}
  else{ 
    fit <- lmFit(logCPM, design)
  }
  fit <- eBayes(fit, trend=T)

   if(!is.null(geneNames)){
    cat("Adding gene names\n")
    results <- merge(as.data.frame(results),unique(geneNames[,c(1,5)]), by.x=0,by.y=1,all.x=T)
    names(results)[1] <- "Ensembl.ID"
  } 
  results <- topTable(fit,n=Inf, coef=ncol(designR)) %>% 
    mutate(sig = case_when(abs(logFC) > lfc & adj.P.Val < fdr ~ "sig", .default =" not sig")) %>%
    rename(., all_of(c(pval = "P.Value", FDR = "adj.P.Val",LFC = "logFC")))
  cat(paste("Number of DEGs: ", table(results$sig)[2], "\n"))
  cat("Done!")
  return(results)
  
}



# edgeR -------------------------------------------------------------------

DEA_edgeR <- function(dge, designTable, model, fdr = 0.05, lfc = 1, 
                      norm = "TMM", geneNames=NULL){
  cat("Starting differential expression analysis with EdgeR \n")
  
  design <- model.matrix(model, data =designTable)
  
  dge <- calcNormFactors(dge, method = norm)
  # estimate dispersion
  dge <- estimateDisp(dge, design, robust=T)
  dge$common.dispersion
  fit <- glmFit(dge,design)
  lrt <- glmLRT(fit)
  results <- topTags(lrt, n =Inf) %>% as.data.frame() %>%
    mutate(sig = case_when(abs(logFC) > 1 & FDR < fdr~ "sig", .default =" not sig"))  %>% select(-1) %>%
    rename(., pval =  PValue )
  if(!is.null(geneNames)){
    cat("Adding gene names\n")
    results <- merge(as.data.frame(results),unique(geneNames[,c(1,5)]), by.x=0,by.y=1,all.x=T)
    names(results)[1] <- "Ensembl.ID"
  }
  cat(paste("Number of DEGs: ", table(results$sig)[2], "\n"))
  cat("Done!")
  return(results)
}


# GLMM --------------------------------------------------------------------
DEA_glmm <- function(dge, designTable, model , dfe = 0.05, lfc =1, 
                     glmm = T, lmm= F, cores = 4,
                     norm = "TMM", geneNames= NULL){
  
  # Required at least 4 cores 
  cl <- makeCluster(cores)
  # create meta table, only have variables that are in the model 
  selectVars <- c("sample",all.vars(model))
  meta <- designTable %>%  select(selectVars) %>% as.data.frame()
  rownames(meta) <- meta$sample
  meta <- meta[,-1]
  
  cpm <- cpm(dge)
  
  if(lmm){ # Gaussian mixed model (lmmSeq)
    lcpm <- log2(cpm + 1)
    results <- lmmSeq( ~ tissue + (1| patient) ,lcpm,metadata = meta)
  }
  if(glmm){ # General linear mixed model (glmmSeq)
    disp <- apply(cpm, 1, function(x) { (var(x)-mean(x))/(mean(x))^2}) ## margin 1 to apply over rows
    sizeFactors <- dge$samples$norm.factors # Estimate size factors
    results <- glmmSeq(model , countdata= cpm,metadata = meta,
                             dispersion = disp, sizeFactors = sizeFactors)
  }
  stopCluster(cl)
  results <- data.frame(FDR  = p.adjust(results@stats$pvals, method="BH"), LFC = results@stats$coef[,2]) %>% 
    mutate(sig = case_when(abs(LFC) > lfc & FDR < fdr ~ "sig", .default =" not sig")) 
  
  if(!is.null(geneNames)){
    cat("Adding gene names\n")
    results <- merge(as.data.frame(results),unique(geneNames[,c(1,5)]), by.x=0,by.y=1,all.x=T)
    names(results)[1] <- "Ensembl.ID"
  }
  cat(paste("\nNumber of DEGs: ", table(results$sig)[2], "\n"))
  cat("Done!")
  return(results)
  
}







