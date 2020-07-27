# BiocManager::install("DESeq2")
# rm(list = ls(pattern = 'temp.*|test.*'))
# source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
# source('functions_pre_preparation.R')




#################
### GSE110256 ###
#################
#counts
seq_GSE110256 <- list(
  'dir' = paste0(opts$dir_seq_input, '/GSE110256')
)

# These should be done only once
# dir.create(seq_GSE110256$dir)
# GEOquery::getGEOSuppFiles(GEO = 'GSE110256', makeDirectory = F, baseDir = seq_GSE110256$dir)

purrr::walk(.x = list.files(path = seq_GSE110256$dir, pattern = '.*gz'), .f = function(x) {GEOquery::gunzip(filename = paste0(seq_GSE110256$dir, '/', x))})

seq_GSE110256$design <- as.data.frame(readr::read_tsv(file = paste0(seq_GSE110256$dir, '/GSE110256_sample_list.txt'), col_types = 'cccccc'))

#Deseq2 uses rownames for sample names I think
rownames(seq_GSE110256$design) <- seq_GSE110256$design$sample_name_readable

#Ensembl IDs, but not symbols are unique thus data ready for analysis
seq_GSE110256$input <- as.data.frame(readr::read_tsv(file = paste0(seq_GSE110256$dir, '/GSE110256_combined.count.txt')))

#This should be count matrix for use in deseq2
seq_GSE110256$input_only_samples_of_interest <- subset(seq_GSE110256$input, select = c('ensembl_id', 'symbol', seq_GSE110256$design$sample_name))

#Here we recode column names for count matrix so that its the same as in design table
seq_GSE110256$design$colnames_in_count_mat <- colnames(seq_GSE110256$input_only_samples_of_interest)[3:length(colnames(seq_GSE110256$input_only_samples_of_interest))]

colnames(seq_GSE110256$input_only_samples_of_interest)[3:length(colnames(seq_GSE110256$input_only_samples_of_interest))] <- as.character(
  recode_values_based_on_key(
    to_recode_chrvec = colnames(seq_GSE110256$input_only_samples_of_interest)[3:length(colnames(seq_GSE110256$input_only_samples_of_interest))],
    replace_this_chrvec = seq_GSE110256$design$sample_name,
    with_this_chrvec = seq_GSE110256$design$sample_name_readable))

#Count matrix should be actual matrix, with rownames as feature identifiers
seq_GSE110256$count_matrix_of_intr <- seq_GSE110256$input_only_samples_of_interest 
rownames(seq_GSE110256$count_matrix_of_intr) <- paste0(seq_GSE110256$count_matrix_of_intr$ensembl_id, '___', seq_GSE110256$count_matrix_of_intr$symbol)
seq_GSE110256$count_matrix_of_intr <- as.matrix(seq_GSE110256$count_matrix_of_intr[-c(1,2)])



#For Deseq2, order of samples in design and matrix needs to be the same. So is the case
all(colnames(seq_GSE110256$count_matrix_of_intr) == seq_GSE110256$design$sample_name_readable)

#Get analysis sets (experimental group pairs)
seq_GSE110256$analyses_to_do <- get_list_of_group_pairs_for_deseq2_analysis(
  cnt_group_name = 'control',
  exp_groups_name_vec = seq_GSE110256$design$treatment[10:23], 
  design_df = seq_GSE110256$design, 
  col_in_design_df_with_group_names = 'treatment', 
  count_matrix = seq_GSE110256$count_matrix_of_intr,
  sample_name_col_in_design_df = 'sample_name_readable')

seq_GSE110256$analyses_done <- purrr::map(
  .x = seq_GSE110256$analyses_to_do, 
  .f = function(x){
    calculate_deseq2(count_matrix_ = x$count_matrix, design_ = x$design, col_with_comparisons_in_design = 'treatment', reference_name_in_design = 'control')
  })

seq_GSE110256$analyses_done_reformated <- purrr::map2(
  .x = seq_GSE110256$analyses_done,
  .y = names(seq_GSE110256$analyses_done),
  .f = function(x,y){
    
    x <- subset(x = x, subset = x$pvalue < 0.05)
    
    properly_formated_data <- set_df_for_append(check_cols_ = opts$check_cols_for_input, length = length(x[[1]]))
    
    properly_formated_data[[opts$check_cols_for_input[['Exp.']] ]] <- seq_GSE110256$design[['Exp.']][[1]]
    properly_formated_data[[opts$check_cols_for_input[['Comp.']] ]] <- as.character(recode_values_based_on_key(to_recode_chrvec = y, replace_this_chrvec = seq_GSE110256$design$treatment[10:23], with_this_chrvec = seq_GSE110256$design[['Comp.']][10:23]))
    properly_formated_data[[opts$check_cols_for_input$`P,Value`]] <- x$pvalue
    properly_formated_data[[opts$check_cols_for_input$`adj,P,Val`]] <- x$padj
    properly_formated_data[[opts$check_cols_for_input$logFC]] <- x$log2FoldChange
    temp_ENSG_ <- stringr::str_remove(string = rownames(x), pattern = '___.*')
    properly_formated_data[[opts$check_cols_for_input$ENSG_]] <- stringr::str_remove(string = temp_ENSG_, pattern = '\\..*')
    properly_formated_data[[opts$check_cols_for_input$Symbol]] <- stringr::str_remove(string = rownames(x), pattern = '.*___')

    return(properly_formated_data)
    
  })

seq_GSE110256$output <- rlist::list.rbind(.data = seq_GSE110256$analyses_done_reformated)

write.table(x = seq_GSE110256$output, file = paste0(opts$dir_r_downloaded_data, '/GSE110256_prepared.tsv'), sep = '\t', row.names = F, dec = ',')
#################
### GSE110256 ###
#################



























#################
### GSE117174 ###
#################
#fpkm normalization - rnaseqGene/cufflinks
seq_GSE117174 <- list(
  'dir' = paste0(opts$dir_seq_input, '/GSE117174')
)

# dir.create(seq_GSE117174$dir)

seq_GSE117174$design <- readr::read_tsv(file = paste0(seq_GSE117174$dir, '/GSE117174_sample_list.txt'))

seq_GSE117174$samples_to_get <- paste(seq_GSE117174$design$sample, sep = ', ')

purrr::walk(.x = seq_GSE117174$samples_to_get, .f = function(x) {GEOquery::getGEOSuppFiles(GEO = x, baseDir = seq_GSE117174$dir, makeDirectory = F)})

purrr::walk(.x = list.files(path = seq_GSE117174$dir, pattern = '.*gz'), .f = function(x) {GEOquery::gunzip(filename = paste0(seq_GSE117174$dir, '/', x))})

seq_GSE117174$samples_to_get <- list.files(path = seq_GSE117174$dir, pattern = 'GSM.*txt')

seq_GSE117174$input_list <- purrr::map(
  .x = seq_GSE117174$samples_to_get, 
  .f = function(x) {readr::read_tsv(paste0(seq_GSE117174$dir, '/', x))})

names(seq_GSE117174$input_list) <- stringr::str_remove(string = seq_GSE117174$samples_to_get, pattern = '_.*')

# extract only these columns from the dataframes, then merge them
# c("tracking_id", "gene_id", "gene_short_name", "locus", "coverage", "FPKM")

# length(unique(seq_GSE110256$input_only_samples_of_interest$ensembl_id))
# length(seq_GSE110256$input_only_samples_of_interest$ensembl_id)
#################
### GSE117174 ###
#################


unique(data_annotation$input$data_NOT_for_probe_annotation$Pub.)



#################
### GSE129143 ###
#################
# DESeq normalization
# counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene - https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

seq_GSE129143 <- list(
  'dir' = paste0(opts$dir_seq_input, '/GSE129143')
)

dir.create(seq_GSE129143$dir)

GEOquery::getGEOSuppFiles(GEO = 'GSE129143', makeDirectory = F, baseDir = seq_GSE129143$dir)

purrr::walk(.x = list.files(path = seq_GSE129143$dir, pattern = '.*gz'), .f = function(x) {GEOquery::gunzip(filename = paste0(seq_GSE129143$dir, '/', x))})

seq_GSE129143$design <- readr::read_tsv(file = paste0(seq_GSE129143$dir, '/GSE129143_sample_list.txt'))

seq_GSE129143$input <- readr::read_tsv(file = paste0(seq_GSE129143$dir, '/GSE129143_RNA_Seq_Normalized_Abundance.txt'))

seq_GSE129143$input_only_samples_of_interest <- subset(seq_GSE129143$input, select = c('GENE', 'NAME', seq_GSE129143$design$sample_name))

# length(unique(seq_GSE110256$input_only_samples_of_interest$ensembl_id))
# length(seq_GSE110256$input_only_samples_of_interest$ensembl_id)
#################
### GSE129143 ###
#################








# https://bioconductor.org/packages/3.11/bioc/html/NOISeq.html
# https://bioconductor.org/packages/3.11/bioc/html/PROPER.html
# https://bioconductor.org/packages/3.11/bioc/html/R4RNA.html
# https://bioconductor.org/packages/3.11/bioc/html/RBM.html
# https://bioconductor.org/packages/3.11/bioc/html/rnaSeqMap.html
# https://bioconductor.org/packages/3.11/bioc/html/RNASeqPower.html
# https://bioconductor.org/packages/3.11/bioc/html/RNASeqR.html
# https://bioconductor.org/packages/3.11/bioc/html/SCBN.html
# https://bioconductor.org/packages/3.11/bioc/html/SCnorm.html
# https://bioconductor.org/packages/3.11/bioc/html/XBSeq.html
# https://bioconductor.org/packages/3.11/bioc/html/ABSSeq.html
# https://bioconductor.org/packages/3.11/bioc/html/APAlyzer.html
# https://bioconductor.org/packages/3.11/bioc/html/BADER.html
# https://bioconductor.org/packages/3.11/bioc/html/BitSeq.html
# https://bioconductor.org/packages/3.11/bioc/html/consensusDE.html
# https://bioconductor.org/packages/3.11/bioc/html/DaMiRseq.html
# https://bioconductor.org/packages/3.11/bioc/html/dearseq.html
# https://bioconductor.org/packages/3.11/bioc/html/DEGseq.html
# https://bioconductor.org/packages/3.11/bioc/html/easyRNASeq.html
# https://bioconductor.org/packages/3.11/bioc/html/EDASeq.html
# https://bioconductor.org/packages/3.11/bioc/html/edge.html
# https://bioconductor.org/packages/3.11/bioc/html/metaSeq.html
# https://bioconductor.org/packages/3.11/bioc/html/metaseqR2.html
# https://bioconductor.org/packages/3.11/bioc/html/MLSeq.html





















# GSM3701233	Nortriptyline replicate 1
# GSM3701234	Nortriptyline replicate 2
# GSM3701235	Nortriptyline replicate 3
# vs
# GSM3701203	Q111SST replicate 1
# GSM3701204	Q111SST replicate 2
# GSM3701205	Q111SST replicate 3








# 75	1	venlafaxine 	PFC
# GSM3272951	PFC_SALIIp1		
# GSM3272952	PFC_SALIIp2		
# GSM3272953	PFC_SALIIp3		
# GSM3272954	PFC_SALIIp4		
# GSM3272955	PFC_SALIp1		
# GSM3272956	PFC_SALIp2		
# GSM3272957	PFC_SALIp3		
# GSM3272958	PFC_SALIp4		
# vs			
# GSM3272959	PFC_WENp1		
# GSM3272960	PFC_WENp2		
# GSM3272961	PFC_WENp3		
# GSM3272962	PFC_WENp4		



# 75	2	venlafaxine 	NAc(septi)
# GSM3272919	NAc_SALIIp1		
# GSM3272920	NAc_SALIIp2		
# GSM3272921	NAc_SALIIp3		
# GSM3272922	NAc_SALIIp4		
# GSM3272923	NAc_SALIp1		
# GSM3272924	NAc_SALIp2		
# GSM3272925	NAc_SALIp3		
# GSM3272926	NAc_SALIp4		
# vs			
# GSM3272927	NAc_WENp1		
# GSM3272928	NAc_WENp2		
# GSM3272929	NAc_WENp3		
# GSM3272930	NAc_WENp4		



# 75	3	mianserin	PFC
# GSM3272951	PFC_SALIIp1		
# GSM3272952	PFC_SALIIp2		
# GSM3272953	PFC_SALIIp3		
# GSM3272954	PFC_SALIIp4		
# GSM3272955	PFC_SALIp1		
# GSM3272956	PFC_SALIp2		
# GSM3272957	PFC_SALIp3		
# GSM3272958	PFC_SALIp4		
# vs			
# GSM3272943	PFC_MIAp1		
# GSM3272944	PFC_MIAp2		
# GSM3272945	PFC_MIAp3		
# GSM3272946	PFC_MIAp4		



# 75	4	mianserin	NAc(septi)
# GSM3272919	NAc_SALIIp1		
# GSM3272920	NAc_SALIIp2		
# GSM3272921	NAc_SALIIp3		
# GSM3272922	NAc_SALIIp4		
# GSM3272923	NAc_SALIp1		
# GSM3272924	NAc_SALIp2		
# GSM3272925	NAc_SALIp3		
# GSM3272926	NAc_SALIp4		
# vs			
# GSM3272911	NAc_MIAp1		
# GSM3272912	NAc_MIAp2		
# GSM3272913	NAc_MIAp3		
# GSM3272914	NAc_MIAp4		



# 75	5	ketamine	PFC
# GSM3272951	PFC_SALIIp1		
# GSM3272952	PFC_SALIIp2		
# GSM3272953	PFC_SALIIp3		
# GSM3272954	PFC_SALIIp4		
# GSM3272955	PFC_SALIp1		
# GSM3272956	PFC_SALIp2		
# GSM3272957	PFC_SALIp3		
# GSM3272958	PFC_SALIp4		
# vs			
# GSM3272935	PFC_KETp1		
# GSM3272936	PFC_KETp2		
# GSM3272937	PFC_KETp3		
# GSM3272938	PFC_KETp4		



# 75	6	ketamine	NAc(septi)
# GSM3272919	NAc_SALIIp1		
# GSM3272920	NAc_SALIIp2		
# GSM3272921	NAc_SALIIp3		
# GSM3272922	NAc_SALIIp4		
# GSM3272923	NAc_SALIp1		
# GSM3272924	NAc_SALIp2		
# GSM3272925	NAc_SALIp3		
# GSM3272926	NAc_SALIp4		
# vs			
# GSM3272903	NAc_KETp1		
# GSM3272904	NAc_KETp2		
# GSM3272905	NAc_KETp3		
# GSM3272906	NAc_KETp4




# #Create DESeqDataSet object
# seq_GSE110256$deseq_test <- DESeq2::DESeqDataSetFromMatrix(
#   countData = seq_GSE110256$analyses_to_do$bupropion$count_matrix, 
#   colData =  seq_GSE110256$analyses_to_do$bupropion$design, 
#   design = ~treatment)
# 
# #Filter low count rows
# # seq_GSE110256$deseq_test_rows_to_keep <- rowSums(counts(seq_GSE110256$deseq_test)) >= 10
# # seq_GSE110256$deseq_test <- seq_GSE110256$deseq_test[seq_GSE110256$deseq_test$rows_to_keep,]
# 
# #Set reference to proper factor
# seq_GSE110256$deseq_test$treatment <- relevel(seq_GSE110256$deseq_test$treatment, ref = "control")
# 
# #Calculate results and get result table
# seq_GSE110256$deseq_test <- DESeq2::DESeq(seq_GSE110256$deseq_test)
# seq_GSE110256$deseq_test_results <- DESeq2::results(seq_GSE110256$deseq_test)
# seq_GSE110256$deseq_test_results_df <- as.data.frame(seq_GSE110256$deseq_test_results)