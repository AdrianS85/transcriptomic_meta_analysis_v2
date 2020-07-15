library(Biobase)
library(GEOquery)
library(limma)
source('functions_GEO_download.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
# rm(cont.matrix, design, ex, fit, fit2, gset, tT, fl, gsms, i, idx, lofC, qx, sml)
# rm(list = ls(pattern = 'temp.*|test.*'))


general_opts <- list('dir' = 'geo/geo_r',
                     'Pub.' = 71,
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
                     'p_cutoff' = 0.05,
                     'series' = "GSE59895")




### 2 ###
analysis_specific_opts <- list(
  'Exp.' = 2,
  'sample_comp' = paste0('paroxetine_1d'),
  'platform' = 'GPL5425',
  'group_design' = paste0("0XXXXXXXXX111XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, analysis_specific_opts$platform, sep = '_')

analysis_71_2 <- get_full_topTable(series_ = general_opts$series, platform_ = analysis_specific_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_71_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###


### 3 ###
analysis_specific_opts <- list(
  'Exp.' = 3,
  'sample_comp' = paste0('paroxetine_3d'),
  'platform' = 'GPL5425',
  'group_design' = paste0("X000XXXXXXXXX11XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, analysis_specific_opts$platform, sep = '_')

analysis_71_3 <- get_full_topTable(series_ = general_opts$series, platform_ = analysis_specific_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_71_3, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 3 ###


### 4 ###
analysis_specific_opts <- list(
  'Exp.' = 4,
  'sample_comp' = paste0('paroxetine_5d'),
  'platform' = 'GPL5425',
  'group_design' = paste0("XXXX000XXXXXXXX111XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, analysis_specific_opts$platform, sep = '_')

analysis_71_4 <- get_full_topTable(series_ = general_opts$series, platform_ = analysis_specific_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_71_4, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 4 ###


### 6 ###
analysis_specific_opts <- list(
  'Exp.' = 6,
  'sample_comp' = paste0('phenelzine_1d'),
  'platform' = 'GPL5425',
  'group_design' = paste0("0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX111XXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, analysis_specific_opts$platform, sep = '_')

analysis_71_6 <- get_full_topTable(series_ = general_opts$series, platform_ = analysis_specific_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_71_6, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 6 ###


### 7 ###
analysis_specific_opts <- list(
  'Exp.' = 7,
  'sample_comp' = paste0('phenelzine_3d'),
  'platform' = 'GPL5425',
  'group_design' = paste0("X000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX111XXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, analysis_specific_opts$platform, sep = '_')

analysis_71_7 <- get_full_topTable(series_ = general_opts$series, platform_ = analysis_specific_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_71_7, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 7 ###


### 8 ###
analysis_specific_opts <- list(
  'Exp.' = 8,
  'sample_comp' = paste0('phenelzine_5d'),
  'platform' = 'GPL5425',
  'group_design' = paste0("XXXX000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX111XXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, analysis_specific_opts$platform, sep = '_')

analysis_71_8 <- get_full_topTable(series_ = general_opts$series, platform_ = analysis_specific_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_71_8, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 8 ###




### 14 ###
analysis_specific_opts <- list(
  'Exp.' = 14,
  'sample_comp' = paste0('tramadol_1d'),
  'platform' = 'GPL5425',
  'group_design' = paste0("0XXXXXXXXXXXXXXXXXXXX11XXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, analysis_specific_opts$platform, sep = '_')

analysis_71_14 <- get_full_topTable(series_ = general_opts$series, platform_ = analysis_specific_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_71_14, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 14 ###


### 15 ###
analysis_specific_opts <- list(
  'Exp.' = 15,
  'sample_comp' = paste0('tramadol_3d'),
  'platform' = 'GPL5425',
  'group_design' = paste0("X000XXXXXXXXXXXXXXXXXXX111XXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, analysis_specific_opts$platform, sep = '_')

analysis_71_15 <- get_full_topTable(series_ = general_opts$series, platform_ = analysis_specific_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_71_15, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 15 ###



### 16 ###
analysis_specific_opts <- list(
  'Exp.' = 16,
  'sample_comp' = paste0('tramadol_5d'),
  'platform' = 'GPL5425',
  'group_design' = paste0("XXXX000XXXXXXXXXXXXXXXXXXX111XXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, analysis_specific_opts$platform, sep = '_')

analysis_71_16 <- get_full_topTable(series_ = general_opts$series, platform_ = analysis_specific_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_71_16, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 16 ###













############
### TEST ###
############

gset <- getGEO("GSE59895", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL5425", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX111XXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

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
# 
tT <- subset(tT, tT$P.Value < 0.05)
