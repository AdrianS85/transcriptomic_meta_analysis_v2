get_full_topTable <- function(series_, platform_, group_names_, col_names_, p_cutoff_, AnnotGPL_)
{
  # load series and platform data from GEO
  gset <- getGEO(series_, GSEMatrix =TRUE, AnnotGPL = AnnotGPL_)
  if (length(gset) > 1) idx <- grep(platform_, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  
  
  
  # make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  
  sml <- c()
  for (i in 1:nchar(group_names_)) { sml[i] <- substr(group_names_,i,i) }
  
  # eliminate samples marked as "X"
  sel <- which(sml != "X")
  sml <- sml[sel]
  gset <- gset[ ,sel]
  
  # log2 transform
  ex <- exprs(gset)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
  
  # set up the data and proceed with analysis
  sml <- paste("G", sml, sep="")    # set group names
  fl <- as.factor(sml)
  gset$description <- fl
  design <- model.matrix(~ description + 0, gset)
  colnames(design) <- levels(fl)
  fit <- lmFit(gset, design)
  cont.matrix <- makeContrasts(G1-G0, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=length(fit2$genes[[1]]))
  
  tT <- subset(tT, tT$P.Value < p_cutoff_, select= col_names_)
  
  return(tT)
}


get_control_and_drug_pairs <- function(subset_to_use_df_)
{
  control_ <- subset(subset_to_use_df_, subset = subset_to_use_df_$`Perturbation type` == 'vehicle')
  
  control_return_ <- stringr::str_c(control_[['Accession']], collapse = ', ')
  
  drugs_ <- subset(subset_to_use_df_, subset = subset_to_use_df_$`Perturbation type` == 'test')
  
  if (length(subset_to_use_df_[[1]]) != (length(control_[[1]]) + length(drugs_[[1]]))) {
    warning('Input is of different lenght than control + drugs')
    return(NULL)
  }
  
  return_ <- list()
  row_nb <- 1
  for (drug in unique(drugs_[['Perturbagen']])) {
    drug_instances <- subset(drugs_, subset = drugs_[['Perturbagen']] == drug)
    
    drug_instances_return_ <- stringr::str_c(drug_instances[['Accession']], collapse = ', ')
    
    return_[[row_nb]] <- paste0(drug_instances_return_, ' vs ', control_return_)
    
    names(return_)[[row_nb]] <- drug
    
    row_nb <- row_nb + 1
  }

  return(return_)
  
}