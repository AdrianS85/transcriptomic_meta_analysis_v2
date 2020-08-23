library(Biobase)
library(GEOquery)
library(limma)
source('functions_GEO_download.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
# rm(cont.matrix, design, ex, fit, fit2, gset, tT, fl, gsms, i, idx, lofC, qx, sml)
# rm(list = ls(pattern = 'temp.*|test.*'))



general_opts <- list('platform' = "",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"))


### 1 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 1,
  'sample_comp' = paste0('mianserin_1h'),
  'group_design' = paste0("XXXXXX0XXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX0XXXXXXXXX",
                          "XXXXXXXXXXXXXX0XXXXXXXXX1XXXXXXXXXXXXX0XXXXXXXXX1X",
                          "XXXXXXXXXXXXXXXX0XXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXX0XXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXX1XXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_1 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_1, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 1 ###


### 2 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 2,
  'sample_comp' = paste0('imipramine_1h'),
  'group_design' = paste0("XX1XXX0XXXXXXXXXXXXX0XXXXXXX1XXXXXXXXXXX0XXXXXXX1X",
                          "XXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX0XXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXX0XXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_2 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_2, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###



### 3 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 3,
  'sample_comp' = paste0('fluoxetine_1h'),
  'group_design' = paste0("XXXXXX0XXX1XXXXXXXXX0XXX1XXXXXXXXXXXXXXX0XXX1XXXXX",
                          "XXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX0XXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXX0XXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_3 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_3, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 3 ###



### 4 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 4,
  'sample_comp' = paste0('bupropion_1h'),
  'group_design' = paste0("XXXXXX0XXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX0XXXXXXXXX",
                          "XXXXXXXXXX1XXX0XXXXXXXXXXXXXXXXXXX1XXX0XXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX0XXXXXXX0XXXXXXXXXXXXX1XXXXXXXXXXX",
                          "XXXXXXXXXXXXXXX0XXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_4 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_4, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 4 ###

### 5 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 5,
  'sample_comp' = paste0('tianeptine_1h'),
  'group_design' = paste0("XXXXXX0XXXXX1XXXXXXX0XXXXXXXXXXX1XXXXXXX0XXXXXXXXX",
                          "XX1XXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX0XXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXX0XXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_5 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_5, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 5 ###

### 6 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 6,
  'sample_comp' = paste0('tranylcypromine_1h'),
  'group_design' = paste0("XXXXXX0XXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX0XXXXXXXXX",
                          "XXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX0XXXXXXX0XXXXXXXXXXXXXXXXX1XXXXX1X",
                          "XXXXXXXXXXXXXXX0XXXXXXXXXXX0XXXXXXXXXXXXXXXXXX1XXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_6 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_6, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 6 ###

### 7 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 7,
  'sample_comp' = paste0('mianserin_2h'),
  'group_design' = paste0("XXXXXXXXXXX0XXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX0XXXX",
                          "XXXXXXXXXXXXXXXXXXXX0XXXXXXXXX1XXXXXXXXXXXXX0XXXXX",
                          "XXXXXXX0X0XXXXXXXXXXXXXXXXXXXXXXXX1XXXXXXXXXXXXXXX",
                          "XXXXXXXXXXX1XXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXX0XXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_7 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_7, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 7 ###

### 8 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 8,
  'sample_comp' = paste0('imipramine_2h'),
  'group_design' = paste0("XXXXXXXXXXX0X1XXXXXXXXXXX0XXXXXXX1XXXXXXXXXXX0XXXX",
                          "XXX1XXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX0XXXXX",
                          "XXXXXXX0X0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXX0XXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_8 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_8, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 8 ###

### 9 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 9,
  'sample_comp' = paste0('fluoxetine_2h'),
  'group_design' = paste0("XXX1XXXXXXX0XXXXXXXXXXXXX0XXX1XXXXXXXXXXXXXXX0XXX1",
                          "XXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX0XXXXX",
                          "XXXXXXX0X0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXX0XXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_9 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_9, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 9 ###


### 10 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 10,
  'sample_comp' = paste0('bupropion_2h'),
  'group_design' = paste0("XXXXXXXXXXX0XXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX0XXXX",
                          "XXXXXXXXXXXXXXXX1XXX0XXXXXXXXXXXXXXXXXXX1XXX0XXXXX",
                          "XXXXXXX0X0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "1XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXX0XXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_10 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_10, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 10 ###


### 11 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 11,
  'sample_comp' = paste0('tianeptine_2h'),
  'group_design' = paste0("XXXXXXXXXXX0XXXXX1XXXXXXX0XXXXXXXXXXX1XXXXXXX0XXXX",
                          "XXXXXXX1XXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX0XXXXX",
                          "XXXXXXX0X0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXX0XXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_11 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_11, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 11 ###


### 12 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 12,
  'sample_comp' = paste0('tranylcypromine_2h'),
  'group_design' = paste0("XXXXXXXXXXX0XXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX0XXXX",
                          "XXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX0XXXXX",
                          "XXXXXXX0X0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXX1XXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXX1XXXX0XXXX",
                          "XXXXXXXX1XXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_12 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_12, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 12 ###


### 13 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 13,
  'sample_comp' = paste0('mianserin_4h'),
  'group_design' = paste0("XXXX0XXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX",
                          "0XXXXXXXXXXX1XXXXXXXXXXXXX0XXXXXXXXX1XXXXXXXXXXXXX",
                          "0XXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXX1XXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXX0XXXXX0XX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_13 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_13, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 13 ###


### 14 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 14,
  'sample_comp' = paste0('imipramine_4h'),
  'group_design' = paste0("XXXX0XXXXXXXXXXXXX1XXXXXXXXXXX0XXXXXXX1XXXXXXXXXXX",
                          "0XXXXXXX1XXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX",
                          "0XXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXX0XXXXX0XX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_14 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_14, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 14 ###


### 15 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 15,
  'sample_comp' = paste0('fluoxetine_4h'),
  'group_design' = paste0("XXXX0XXXXXXXXX1XXXXXXXXXXXXXXX0XXX1XXXXXXXXXXXXXXX",
                          "0XXX1XXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX",
                          "0XXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXX0XXXXX0XX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_15 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_15, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 15 ###


### 16 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 16,
  'sample_comp' = paste0('bupropion_4h'),
  'group_design' = paste0("XXXX0XXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX",
                          "0XXXXXXXXXXXXXXXXXXXXX1XXX0XXXXXXXXXXXXXXXXXXX1XXX",
                          "0XXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXX1XXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXX0XXXXX0XX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_16 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_16, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 16 ###


### 17 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 17,
  'sample_comp' = paste0('tianeptine_4h'),
  'group_design' = paste0("XXXX0XXX1XXXXXXXXXXXXX1XXXXXXX0XXXXXXXXXXX1XXXXXXX",
                          "0XXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXX",
                          "0XXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXX0XXXXX0XX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_17 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_17, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 17 ###


### 18 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 18,
  'sample_comp' = paste0('tranylcypromine_4h'),
  'group_design' = paste0("XXXX0XXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX",
                          "0XXXXXXXXXXXXXXXXXXXX1XXXX0XXXXXXXXXXXXXXXXXXXXXXX",
                          "0XXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXX1XXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                          "XX1XXXX0XXXXX0XX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_18 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_18, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 18 ###


### 19 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 19,
  'sample_comp' = paste0('mianserin_8h'),
  'group_design' = paste0("XXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXX",
                          "XXXXX0XXXXXXXXXXXX1XXXXXXXXXXXXX0XXXXXXXXX1XXXXXXX",
                          "XXX0XXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXX0XXXXXX0XXXXXX",
                          "XX1XXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_19 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_19, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 19 ###


### 20 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 20,
  'sample_comp' = paste0('imipramine_8h'),
  'group_design' = paste0("XXXXXXXXX1XXXXX0XXXXXXX1XXXXXXXXXXX0XXXXXXX1XXXXXX",
                          "XXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXX",
                          "XXX0XXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXX0XXXXXX0XXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_20 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_20, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### x ###


### 21 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 21,
  'sample_comp' = paste0('fluoxetine_8h'),
  'group_design' = paste0("XXXXXXXXXXXXXXX0XXX1XXXXXXXXXXXXXXX0XXX1XXXXXXXXXX",
                          "XXXXX0XXX1XXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXX",
                          "XXX0XXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXX0XXXXXX0XXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_21 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_21, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 21 ###


### 22 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 22,
  'sample_comp' = paste0('bupropion_8h'),
  'group_design' = paste0("XXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXX",
                          "XXXXX0XXXXXXXXXXXXXXXXXXXXXX1XXX0XXXXXXXXXXXXXXXXX",
                          "XXX0XXXXXXXXXXXXXXXX0XXXXX1XXXXX1XXX0XXXXXX0XXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_22 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_22, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 22 ###


### 23 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 23,
  'sample_comp' = paste0('tianeptine_8h'),
  'group_design' = paste0("X1XXXXXXXXXXXXX0XXXXXXXXXXX1XXXXXXX0XXXXXXXXXXX1XX",
                          "XXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXX",
                          "XXX0XXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXX0XXXXXX0XXXXXX",
                          "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_23 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_23, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 23 ###


### 24 ###
analysis_specific_opts <- list(
  'Pub.' = 33,
  'Exp.' = 24,
  'sample_comp' = paste0('tranylcypromine_8h'),
  'group_design' = paste0("XXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXX0XXXXXXXXXXXXXX",
                          "XXXXX0XXXXXXXXX1XXXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXX",
                          "XXX0XXXXXXXXXXXXXXXX0XXXXXXXXX1XXXXX0XXXXXX0XXXXXX",
                          "XXXXXXXXXXXXXXXX1XXXXXXXXXXXX0XXXXXXXXXXXXXXXXXXXX",
                          "XXXXXXXXXXXXXXXX")
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_33_24 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_33_24, file = paste0(opts$dir_r_downloaded_data, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 24 ###










############
### TEST ###
############

# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=length(fit2$genes[[1]]))
# # 
# tT <- subset(tT, tT$P.Value < 0.05)