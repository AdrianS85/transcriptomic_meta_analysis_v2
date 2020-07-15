library(Biobase)
library(GEOquery)
library(limma)
source('functions_GEO_download.R')





### PUBLICATION 36 / GSE47541 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL17223",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Search_key","ProbeId","Gid","Transcript","GB_ACC","Symbol","Definition"),
                     'p_cutoff' = 0.05)

### 1 ###

analysis_specific_opts <- list(
  'Pub.' = 36,
  'Exp.' = 1,
  'sample_comp' = paste0('cortex'),
  'series' = "GSE47541",
  'group_design' = paste0('XXXXXXX11110000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_36_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_36_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')

### 1 ###



### 2 ###

analysis_specific_opts <- list(
  'Pub.' = 36,
  'Exp.' = 2,
  'sample_comp' = paste0('hippocampus'),
  'series' = "GSE47541",
  'group_design' = paste0('XXXXXXXXXXXXXXXXXXXXXXX11110000XXXXXXXXXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_36_2 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_36_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')

### 2 ###



### 3 ###

analysis_specific_opts <- list(
  'Pub.' = 36,
  'Exp.' = 3,
  'sample_comp' = paste0('dorsal_raphe'),
  'series' = "GSE47541",
  'group_design' = paste0('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX11110000')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_36_3 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_36_3, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')

### 3 ###
### PUBLICATION 36 / GSE47541 ###





### PUBLICATION 66 / GSE6476 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL1261",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
                     'p_cutoff' = 0.05)

analysis_specific_opts <- list(
  'Pub.' = 66,
  'Exp.' = 1,
  'sample_comp' = paste0('none'),
  'series' = "GSE6476",
  'group_design' = paste0('0011')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_66_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_66_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### PUBLICATION 66 / GSE6476 ###





### PUBLICATION 68 / GSE26836 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL6885",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","GI","GenBank.Accession"),
                     'p_cutoff' = 0.05)

### 1 ###
analysis_specific_opts <- list(
  'Pub.' = 68,
  'Exp.' = 1,
  'sample_comp' = paste0('cortex'),
  'series' = "GSE26836",
  'group_design' = paste0('XXXXXX111000')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_68_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_68_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')

### 1 ###


### 2 ###
analysis_specific_opts <- list(
  'Pub.' = 68,
  'Exp.' = 2,
  'sample_comp' = paste0('hippocampus'),
  'series' = "GSE26836",
  'group_design' = paste0('000111XXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_68_2 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_68_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###
### PUBLICATION 68 / GSE26836 ###





### PUBLICATION 67 / GSE26364 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL6246",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
                     'p_cutoff' = 0.05)

analysis_specific_opts <- list(
  'Pub.' = 67,
  'Exp.' = 1,
  'sample_comp' = paste0('none'),
  'series' = "GSE26364",
  'group_design' = paste0('0011')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_67_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_67_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')

### PUBLICATION 67 / GSE26364 ###





### PUBLICATION 70 / GSE35761-GSE35763-GSE35765 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL1261",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
                     'p_cutoff' = 0.05)

### 1 ###
analysis_specific_opts <- list(
  'Pub.' = 70,
  'Exp.' = 1,
  'sample_comp' = paste0('S100a10_bacTRAP'),
  'series' = "GSE35761",
  'group_design' = paste0('101010')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_70_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_70_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 1 ###

### 2 ###
analysis_specific_opts <- list(
  'Pub.' = 70,
  'Exp.' = 2,
  'sample_comp' = paste0('Glt25d2_bacTRAP'),
  'series' = "GSE35763",
  'group_design' = paste0('101010')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_70_2 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_70_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###

### 3 ###
analysis_specific_opts <- list(
  'Pub.' = 70,
  'Exp.' = 3,
  'sample_comp' = paste0('p11_KO_S100a10_bacTRAP'),
  'series' = "GSE35765",
  'group_design' = paste0('110010')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_70_3 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_70_3, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 3 ###

### PUBLICATION 70 / GSE35761-GSE35763-GSE35765 ###





### PUBLICATION 77 / GSE63469 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL1261",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
                     'p_cutoff' = 0.05)

### 1 ###
analysis_specific_opts <- list(
  'Pub.' = 77,
  'Exp.' = 1,
  'sample_comp' = paste0('30_mg'),
  'series' = "GSE63469",
  'group_design' = paste0('011X0X')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_77_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_77_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 1 ###

### 2 ###
analysis_specific_opts <- list(
  'Pub.' = 77,
  'Exp.' = 2,
  'sample_comp' = paste0('100_mg'),
  'series' = "GSE63469",
  'group_design' = paste0('0XX101')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_77_2 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_77_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###

### PUBLICATION 77 / GSE63469 ###





### PUBLICATION 78 / GSE118668-GSE118669 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL1261",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI"),
                     'p_cutoff' = 0.05)

### 1 ###
analysis_specific_opts <- list(
  'Pub.' = 78,
  'Exp.' = 1,
  'sample_comp' = paste0('cortex'),
  'series' = "GSE118668",
  'group_design' = paste0('1111111100000000')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_78_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_78_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 1 ###

### 2 ###
analysis_specific_opts <- list(
  'Pub.' = 78,
  'Exp.' = 2,
  'sample_comp' = paste0('hippocampus'),
  'series' = "GSE118669",
  'group_design' = paste0('1111111100000000')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_78_2 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_78_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###

### PUBLICATION 78 / GSE118668-GSE118669 ###





### PUBLICATION 28 / GSE27532 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL6885",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
                     'p_cutoff' = 0.05)

### 1 ###
analysis_specific_opts <- list(
  'Pub.' = 28,
  'Exp.' = 1,
  'sample_comp' = paste0('HA'),
  'series' = "GSE27532",
  'group_design' = paste0('XXXXXXXX11110000')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_28_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_28_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 1 ###

### 2 ###
analysis_specific_opts <- list(
  'Pub.' = 28,
  'Exp.' = 2,
  'sample_comp' = paste0('LA'),
  'series' = "GSE27532",
  'group_design' = paste0('11110000XXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_28_2 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_28_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###

### PUBLICATION 28 / GSE27532 ###





### PUBLICATION 72 / GSE73798-GSE73799 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL6887",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
                     'p_cutoff' = 0.05)

### 1 ###
analysis_specific_opts <- list(
  'Pub.' = 72,
  'Exp.' = 1,
  'sample_comp' = paste0('hippocampus_1h'),
  'series' = "GSE73798",
  'group_design' = paste0('1XXXX0XXXXXX1XXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXX0X1XXXXXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_72_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_72_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 1 ###

### 2 ###
analysis_specific_opts <- list(
  'Pub.' = 72,
  'Exp.' = 2,
  'sample_comp' = paste0('hippocampus_2h'),
  'series' = "GSE73798",
  'group_design' = paste0('XXXXXXXXXXXXXXXXX01XXXX0XXXXXXXXXXXX1XXXXXXXXXXXXXXXX01XXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_72_2 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_72_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###

### 3 ###
analysis_specific_opts <- list(
  'Pub.' = 72,
  'Exp.' = 3,
  'sample_comp' = paste0('hippocampus_4h'),
  'series' = "GSE73798",
  'group_design' = paste0('XXXXXX1XXXXXXXXXXXXXXXXXXXXXXX1XXXX0XXXXX01XXXXXXXXXXXXXXXX0')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_72_3 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_72_3, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 3 ###

### 4 ###
analysis_specific_opts <- list(
  'Pub.' = 72,
  'Exp.' = 4,
  'sample_comp' = paste0('hippocampus_8h'),
  'series' = "GSE73798",
  'group_design' = paste0('X1XXXXXXXXX0XXXXXXXXXXXX1XXXX0XXXXXXXXXXXXXXXXX0X1XXXXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_72_4 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_72_4, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 4 ###

### 5 ###
analysis_specific_opts <- list(
  'Pub.' = 72,
  'Exp.' = 5,
  'sample_comp' = paste0('striatum_1h'),
  'series' = "GSE73799",
  'group_design' = paste0('XXXXXX1XXXX0XXXX0X1XXXXXXXXXXXXXXXXXXXXX0X1XXXXXXXXXXXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_72_5 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_72_5, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 5 ###

### 6 ###
analysis_specific_opts <- list(
  'Pub.' = 72,
  'Exp.' = 6,
  'sample_comp' = paste0('striatum_2h'),
  'series' = "GSE73799",
  'group_design' = paste0('XXXXXXXXXXXXXXXXXXXXXXX01XXXX01XXXXXXXXXXXXXXXX01XXXXXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_72_6 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_72_6, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 6 ###

### 7 ###
analysis_specific_opts <- list(
  'Pub.' = 72,
  'Exp.' = 7,
  'sample_comp' = paste0('striatum_4h'),
  'series' = "GSE73799",
  'group_design' = paste0('1XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX01XXXXXXXXXXXXXXXX01XXXX0')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_72_7 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_72_7, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 7 ###

### 8 ###
analysis_specific_opts <- list(
  'Pub.' = 72,
  'Exp.' = 8,
  'sample_comp' = paste0('striatum_8h'),
  'series' = "GSE73799",
  'group_design' = paste0('XXXXX0X1XXXX1XXXX0XXXXXXXXXXXXXXXXXXXXXXX0X1XXXXXXXXXXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_72_8 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = T)

write.table(x = analysis_72_8, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 8 ###

### PUBLICATION 72 / GSE73798-GSE73799 ###





### PUBLICATION 60 / GSE97916-GSE82016 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL10427",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","accessions","ProbeUID","ProbeName","GeneName","SystematicName","Description","GB_LIST"),
                     'p_cutoff' = 0.05)

### 1 ###
analysis_specific_opts <- list(
  'Pub.' = 60,
  'Exp.' = 1,
  'sample_comp' = paste0('hippocampus'),
  'series' = "GSE97916",
  'group_design' = paste0('0001111XXXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_60_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_60_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 1 ###


### 2 ###
analysis_specific_opts <- list(
  'Pub.' = 60,
  'Exp.' = 2,
  'sample_comp' = paste0('hypothalamus'),
  'series' = "GSE82016",
  'group_design' = paste0('0001111XXXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_60_2 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_60_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###


### PUBLICATION 60 / GSE97916-GSE82016 ###





### PUBLICATION 45 / GSE76110 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL10427",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","accessions","ProbeUID","ProbeName","GeneName","SystematicName","Description","GB_LIST"),
                     'p_cutoff' = 0.05)

analysis_specific_opts <- list(
  'Pub.' = 45,
  'Exp.' = 1,
  'sample_comp' = paste0('none'),
  'series' = "GSE76110",
  'group_design' = paste0('00011111XXXXXXXX')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_45_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_45_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### PUBLICATION 45 / GSE76110 ###





### PUBLICATION 34 / GSE42940 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL8160",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC"),
                     'p_cutoff' = 0.05)

analysis_specific_opts <- list(
  'Pub.' = 34,
  'Exp.' = 1,
  'sample_comp' = paste0('none'),
  'series' = "GSE42940",
  'group_design' = paste0('00001111')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_34_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_34_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### PUBLICATION 34 / GSE42940 ###





### PUBLICATION 56 / GSE93041 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL16570",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","probeset_id","gene_assignment","mrna_assignment","swissprot","unigene"),
                     'p_cutoff' = 0.05)

analysis_specific_opts <- list(
  'Pub.' = 56,
  'Exp.' = 1,
  'sample_comp' = paste0('none'),
  'series' = "GSE93041",
  'group_design' = paste0('000XXX111')
)
analysis_specific_opts$name <- paste(analysis_specific_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, analysis_specific_opts$series, general_opts$platform, sep = '_')

analysis_56_1 <- get_full_topTable(series_ = analysis_specific_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_56_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### PUBLICATION 56 / GSE93041 ###





### PUBLICATION 50 / GSE63005 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL16570",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","probeset_id","gene_assignment","mrna_assignment","swissprot","unigene"),
                     'Pub.' = 50,
                     'series' = "GSE63005",
                     'p_cutoff' = 0.05)


### 1 ###
analysis_specific_opts <- list(
  'Exp.' = 1,
  'sample_comp' = paste0('amitriptyline'),
  'group_design' = paste0('11XX00XXXXXX')
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_50_1 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_50_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 1 ###



### 2 ###
analysis_specific_opts <- list(
  'Exp.' = 2,
  'sample_comp' = paste0('citalopram'),
  'group_design' = paste0('XX1100XXXXXX')
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_50_2 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_50_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###


### 3 ###
analysis_specific_opts <- list(
  'Exp.' = 3,
  'sample_comp' = paste0('duloxetine'),
  'group_design' = paste0('XXXX0011XXXX')
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_50_3 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_50_3, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 3 ###


### 4 ###
analysis_specific_opts <- list(
  'Exp.' = 4,
  'sample_comp' = paste0('imipramine'),
  'group_design' = paste0('XXXX00XX11XX')
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_50_4 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_50_4, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 4 ###


### 5 ###
analysis_specific_opts <- list(
  'Exp.' = 5,
  'sample_comp' = paste0('mirtazapine'),
  'group_design' = paste0('XXXX00XXXX11')
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_50_5 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_50_5, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 5 ###
### PUBLICATION 50 / GSE63005 ###






### PUBLICATION 52 / GSE84183 ###
general_opts <- list('dir' = 'geo/geo_r',
                     'platform' = "GPL13912",
                     'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","NAME","GB_ACC","GENE_ID","GENE_SYMBOL","GENE_NAME","UNIGENE_ID","ENSEMBL_ID","ACCESSION_STRING","DESCRIPTION"),
                     'series' = "GSE84183",
                     'Pub.' = 52,
                     'p_cutoff' = 0.05)

### 1 ###
analysis_specific_opts <- list(
  'Exp.' = 1,
  'sample_comp' = paste0('hippocampus'),
  'group_design' = paste0('00000000XXXXXXXXXXXXXXXX11111111XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_52_1 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_52_1, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 1 ###



### 2 ###
analysis_specific_opts <- list(
  'Exp.' = 2,
  'sample_comp' = paste0('cortex'),
  'group_design' = paste0('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX00000000XXXXXXXXXXXXXXXX11111111')
)
analysis_specific_opts$name <- paste(general_opts$Pub., analysis_specific_opts$Exp., analysis_specific_opts$sample_comp, general_opts$series, general_opts$platform, sep = '_')

analysis_52_2 <- get_full_topTable(series_ = general_opts$series, platform_ = general_opts$platform, group_names_ = analysis_specific_opts$group_design, col_names_ = general_opts$col_names, p_cutoff_ = general_opts$p_cutoff, AnnotGPL_ = F)

write.table(x = analysis_52_2, file = paste0(general_opts$dir, '/', analysis_specific_opts$name, '.tsv'), sep = '\t', row.names = F, dec = ',')
### 2 ###
### PUBLICATION 52 / GSE84183 ###




##############################
### ADDITIONAL EXPERIMENTS ###
##############################

##############################
### ADDITIONAL EXPERIMENTS ###
##############################







###############
### TESTING ###
###############





# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=length(fit2$genes[[1]]))
# 
# tT <- subset(tT, tT$P.Value < 0.05)
# 
# rm(cont.matrix, design, ex, fit, fit2, gset, tT, fl, gsms, i, idx, lofC, qx, sml)


