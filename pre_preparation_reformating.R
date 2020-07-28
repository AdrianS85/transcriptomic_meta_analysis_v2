# rm(list = ls(pattern = 'temp.*|test.*'))
source('opts.R')
library(loggit)

set_logfile(paste0(opts$pre_prep_reform$folder, "/loggit.log"))





######################
### LOAD PRE-INPUT ###
######################
pre_prep_reform <- list()

pre_prep_reform <- purrr::map(
  .x = opts$pre_prep_reform[names(opts$pre_prep_reform) %in% opts$pre_prep_reform[['platform_to_get_exp_and_comp_nb_for']]], 
  .f = function(platform){
    list('pre-input' = read_files_and_add_pubExp(
      opts_pre_prep_ = platform, 
      directory = opts$dir_r_downloaded_data, 
      input_col_types = opts$pre_prep_reform$input_col_types, 
      decimal = opts$decimal, 
      check_cols = opts$check_cols_for_input))
  })
######################
### LOAD PRE-INPUT ###
######################











#####################
### PREPARE INPUT ###
#####################
opts$pre_prep_reform[['pass_nb']] <- 0
for (platform in opts$pre_prep_reform[['platform_to_get_exp_and_comp_nb_for']]) {
  message(sprintf('processing platform %s...', platform))
  
  message('binding pre-input into single df...')
  pre_prep_reform[[platform]][['input']] <- rlist::list.rbind(pre_prep_reform[[platform]][['pre-input']])
  
  
  
  message('performing qa...')
  pre_prep_reform[[platform]][['qa_input']] <- verify_df(
    df_ = pre_prep_reform[[platform]][['input']], 
    only_qa = T)
  

  
  message('setting output df to be appended with reformated identifers...')

  pre_prep_reform[[platform]][['output']] <- set_df_for_append(
    check_cols_ = opts$check_cols_for_input, 
    length = length(pre_prep_reform[[platform]][['input']][[1]])) ### !!! Here all columns are returned as character columns!
  
  
  
  message('copying standard columns returned by geo2r into output df...')
  
  if (opts$pre_prep_reform[['pass_nb']] == 0) { 
    message(sprintf('NOTE!!! columns copied are: %s', copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = opts[['check_cols_for_input']], return_copied_col_names = T))) 
  }
  
  pre_prep_reform[[platform]][['output']] <- copy_out_of_box_ready_columns_into_final_output_wrapper(
    check_cols_ = opts$check_cols_for_input, 
    final_output_ = pre_prep_reform[[platform]][['output']], 
    input_ = pre_prep_reform[[platform]][['input']])
  
  
  
  opts$pre_prep_reform[['pass_nb']] <- opts$pre_prep_reform[['pass_nb']] + 1
  message(sprintf('platform %s processed.', platform))
}
#####################
### PREPARE INPUT ###
#####################






############################################
### Platform-specific column reformating ###
############################################
### GPL1261 - 66, 70, 77, 78 ###
pre_prep_reform$GPL1261$output$Probe <- pre_prep_reform$GPL1261$input$ID
pre_prep_reform$GPL1261$output$Description <- pre_prep_reform$GPL1261$input$`Gene.title`
pre_prep_reform$GPL1261$output$Gene_ID <- stringr::str_replace_all(string = pre_prep_reform$GPL1261$input$`Gene.ID`, pattern = '///', replacement = ', ')
pre_prep_reform$GPL1261$output$Unigene <- pre_prep_reform$GPL1261$input$`UniGene.ID`
pre_prep_reform$GPL1261$output$Nucleotide <- pre_prep_reform$GPL1261$input$GI
### GPL1261 - 66, 70, 77, 78 ###

### GPL5425 - 71 ###
pre_prep_reform$GPL5425$output$Description <- pre_prep_reform$GPL5425$input$`Gene.title`
pre_prep_reform$GPL5425$output$Gene_ID <- stringr::str_replace_all(string = pre_prep_reform$GPL5425$input$`Gene.ID`, pattern = '///', replacement = ', ')
pre_prep_reform$GPL5425$output$Unigene <- pre_prep_reform$GPL5425$input$`UniGene.ID`
pre_prep_reform$GPL5425$output$Nucleotide <- pre_prep_reform$GPL5425$input$GI
### GPL5425 - 71 ###

### GPL6246 - 67 ###
pre_prep_reform$GPL6246$output$Probe <- pre_prep_reform$GPL6246$input$ID
pre_prep_reform$GPL6246$output$Description <- pre_prep_reform$GPL6246$input$`Gene.title`
pre_prep_reform$GPL6246$output$Gene_ID <- stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$`Gene.ID`, pattern = '///', replacement = ', ')
pre_prep_reform$GPL6246$output$Unigene <- stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$`UniGene.ID`, pattern = '///', replacement = ', ')
pre_prep_reform$GPL6246$output$Nucleotide <- stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$GI, pattern = '///', replacement = ', ')
### GPL6246 - 67 ###

### GPL6885 - 28, 68 ###
pre_prep_reform$GPL6885$output$Probe <- pre_prep_reform$GPL6885$input$ID
pre_prep_reform$GPL6885$output$Description <- pre_prep_reform$GPL6885$input$`Gene.title`
pre_prep_reform$GPL6885$output$Gene_ID <- pre_prep_reform$GPL6885$input$`Gene.ID`
pre_prep_reform$GPL6885$output$Unigene <- pre_prep_reform$GPL6885$input$`UniGene.ID`
pre_prep_reform$GPL6885$output$Nucleotide <- pre_prep_reform$GPL6885$input$GI
### GPL6885 - 28, 68 ###

### GPL6887 - 33, 72 ###
pre_prep_reform$GPL6887$output$Probe <- pre_prep_reform$GPL6887$input$ID
pre_prep_reform$GPL6887$output$Description <- pre_prep_reform$GPL6887$input$`Gene.title`
pre_prep_reform$GPL6887$output$Gene_ID <- stringr::str_replace_all(string = pre_prep_reform$GPL6887$input$`Gene.ID`, pattern = '///', replacement = ', ')
pre_prep_reform$GPL6887$output$Unigene <- pre_prep_reform$GPL6887$input$`UniGene.ID`
pre_prep_reform$GPL6887$output$Nucleotide <- pre_prep_reform$GPL6887$input$GI
### GPL6887 - 33, 72 ###

### GPL8160 - 34 ###
pre_prep_reform$GPL8160$output$Probe <- pre_prep_reform$GPL8160$input$ID
### GPL8160 - 34 ###

### GPL10427 - 60, 45 ###
pre_prep_reform$GPL10427$output$Probe <- pre_prep_reform$GPL10427$input$ProbeName
pre_prep_reform$GPL10427$output$Description <- pre_prep_reform$GPL10427$input$Description
### GPL10427 - 60, 45 ###

### GPL13912 - 52 ###
pre_prep_reform$GPL13912$output$Description <- pre_prep_reform$GPL13912$input$GENE_NAME
pre_prep_reform$GPL13912$output$Probe <- pre_prep_reform$GPL13912$input$NAME
pre_prep_reform$GPL13912$output$Gene_ID <- pre_prep_reform$GPL13912$input$GENE_ID
pre_prep_reform$GPL13912$output$Unigene <- pre_prep_reform$GPL13912$input$UNIGENE_ID
### GPL13912 - 52 ###

### GPL16570 - 56, 50 ###
pre_prep_reform$GPL16570$output$Probe <- pre_prep_reform$GPL16570$input$probeset_id
### GPL16570 - 56, 50 ###

### GPL17223 - 36 ###
pre_prep_reform$GPL17223$output$Probe <- pre_prep_reform$GPL17223$input$ProbeId
pre_prep_reform$GPL17223$output$Nucleotide <- pre_prep_reform$GPL17223$input$Gid
pre_prep_reform$GPL17223$output$Description <- pre_prep_reform$GPL17223$input$Definition
### GPL17223 - 36 ###

### GPL25480 - 76 ###
pre_prep_reform$GPL25480$output$Probe <- pre_prep_reform$GPL25480$input$SPOT_ID
### GPL25480 - 76 ###
############################################
### Platform-specific column reformating ###
############################################







  





################################################################################
### Prepare specific, composite bind_for_extraction column for each platform ###
################################################################################
### GPL1261 - 66, 70, 77, 78 ###
pre_prep_reform$GPL1261$input$bind_for_extraction <- paste(
  stringr::str_replace_all(string = pre_prep_reform$GPL1261$input$`Gene.symbol`, pattern = '///', replacement = opts$pre_prep_reform$bind_sep), 
  pre_prep_reform$GPL1261$input$`UniGene.symbol`, 
  pre_prep_reform$GPL1261$input$`GenBank.Accession`, 
  sep = opts$pre_prep_reform$bind_sep) 

### GPL5425 - 71 ###
pre_prep_reform$GPL5425$input$bind_for_extraction <- paste(
  stringr::str_remove(string = pre_prep_reform$GPL5425$input$ID, pattern = '_P.*'), 
  stringr::str_replace_all(string = pre_prep_reform$GPL5425$input$`Gene.symbol`, pattern = '///', replacement = opts$pre_prep_reform$bind_sep), 
  pre_prep_reform$GPL5425$input$`GenBank.Accession`, 
  pre_prep_reform$GPL5425$input$`UniGene.symbol`, 
  sep = opts$pre_prep_reform$bind_sep)

### GPL6246 - 67 ###
pre_prep_reform$GPL6246$input$bind_for_extraction <- paste(
  stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$`Gene.symbol`, pattern = '///', replacement = opts$pre_prep_reform$bind_sep), 
  stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$`UniGene.symbol`, pattern = '///', replacement = opts$pre_prep_reform$bind_sep), 
  stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$`GenBank.Accession`, pattern = '///', replacement = opts$pre_prep_reform$bind_sep), 
  sep = opts$pre_prep_reform$bind_sep)

### GPL6885 - 28, 68 ###
pre_prep_reform$GPL6885$input$bind_for_extraction <- paste(
  pre_prep_reform$GPL6885$input$`Gene.symbol`, 
  pre_prep_reform$GPL6885$input$`UniGene.symbol`, 
  pre_prep_reform$GPL6885$input$`GenBank.Accession`, 
  sep = opts$pre_prep_reform$bind_sep)

### GPL6887 - 33, 72 ###
pre_prep_reform$GPL6887$input$bind_for_extraction <- paste(
  stringr::str_replace_all(string = pre_prep_reform$GPL6887$input$`Gene.symbol`, pattern = '///', replacement = opts$pre_prep_reform$bind_sep), 
  pre_prep_reform$GPL6887$input$`UniGene.symbol`, 
  pre_prep_reform$GPL6887$input$`GenBank.Accession`, 
  sep = opts$pre_prep_reform$bind_sep)

### GPL8160 - 34 ###
pre_prep_reform$GPL8160$input$bind_for_extraction <- pre_prep_reform$GPL8160$input$GB_ACC # 

### GPL10427 - 60, 45
pre_prep_reform$GPL10427$input$bind_for_extraction <- paste(
  prepare_uncorrupted_GPL10427_accessions_wrapper(GPL10427_accessions = pre_prep_reform$GPL10427$input$accessions), 
  pre_prep_reform$GPL10427$input$GeneName, 
  pre_prep_reform$GPL10427$input$SystematicName, 
  stringr::str_replace_all(string = pre_prep_reform$GPL10427$input$GB_LIST, pattern = ', ', replacement = opts$pre_prep_reform$bind_sep),
  sep = opts$pre_prep_reform$bind_sep)

### GPL13912 - 52 ### ### !!! REMOVING 'NAP012158-001' NAMES, A_55_P2114333 ### !!!
pre_prep_reform$GPL13912$input$bind_for_extraction <- paste(
  GPL13912_cols = pre_prep_reform$GPL13912$input$GB_ACC,
  pre_prep_reform$GPL13912$input$GENE_SYMBOL,
  pre_prep_reform$GPL13912$input$ENSEMBL_ID,
  prepare_uncorrupted_GPL13912_ACCESSION_STRING_wrapper(pre_prep_reform$GPL13912$input$ACCESSION_STRING),
  sep = opts$pre_prep_reform$bind_sep)

### GPL16570 - 56, 50 ###
pre_prep_reform$GPL16570$input$bind_for_extraction <- paste(
  prepare_noncorrupted_GPL16570_cols(GPL16570_cols = pre_prep_reform$GPL16570$input$gene_assignment),
  sep = opts$pre_prep_reform$bind_sep)

### GPL17223 - 36 
pre_prep_reform$GPL17223$input$bind_for_extraction <- paste(
  pre_prep_reform$GPL17223$input$Symbol, 
  pre_prep_reform$GPL17223$input$Transcript, 
  pre_prep_reform$GPL17223$input$GB_ACC, 
  sep = opts$pre_prep_reform$bind_sep)

### GPL25480 - 76 ###
pre_prep_reform$GPL25480$input$bind_for_extraction <- NULL
################################################################################
### Prepare specific, composite bind_for_extraction column for each platform ###
################################################################################









####################################################
### PREPARE FOR EXTRACTION FROM COMPOSITE COLUMN ###
####################################################
for (platform in opts$pre_prep_reform[['platform_to_get_exp_and_comp_nb_for']])
{
  message(sprintf('processing platform %s...', platform))
  
  if (!is.null(pre_prep_reform[[platform]][['input']][['bind_for_extraction']])) {

    col_with_ids_for_extraction <- length(pre_prep_reform[[platform]][['input']])
    
    message(sprintf('column bind_for_extraction set to position %i', col_with_ids_for_extraction))
    
    
    
    message('separating identifiers from single column...')

    pre_prep_reform[[platform]][['sep_all_identifiers']] <- as.data.frame(
      stringr::str_split_fixed(
        string = pre_prep_reform[[platform]][['input']][[col_with_ids_for_extraction]],
        pattern = opts$pre_prep_reform$bind_sep,
        n = Inf))
  }

  message(sprintf('platform %s processed', platform))
}
####################################################
### PREPARE FOR EXTRACTION FROM COMPOSITE COLUMN ###
####################################################




#########################################
### EXTRACT IDS FROM COMPOSITE COLUMN ###
#########################################
for (platform in opts$pre_prep_reform[['platform_to_get_exp_and_comp_nb_for']]) 
{
  if (!is.null(pre_prep_reform[[platform]][['input']][['bind_for_extraction']])) {
    
    pre_prep_reform[[platform]][['output']] <- extract_identifers(
      sep_all_identifiers_ = pre_prep_reform[[platform]][['sep_all_identifiers']], 
      search_by = opts$pre_prep_reform[['search_bys']][['search_by']], 
      output_ = pre_prep_reform[[platform]][['output']], 
      opts_pre_check_cols = opts[['check_cols_for_input']],
      range_of_cols_with_IDs_start = opts[['pre_prep_reform']][['range_of_cols_with_IDs_in_output_start']], 
      range_of_cols_with_IDs_end = opts[['pre_prep_reform']][['range_of_cols_with_IDs_in_output_end']]) ### !!! does unigene search_by really identifies Mms?
    
    
    
    message('performing qa...')
    pre_prep_reform[[platform]][['output_qa']] <- verify_df(df_ = pre_prep_reform[[platform]][['output']], only_qa = T)
    
    
    
    message(sprintf('writing data to %s', opts$dir_r_downloaded_data))
    write.table(
      x = pre_prep_reform[[platform]][['output']], 
      file = paste0(opts$pre_prep_reform$folder, '/', platform,  '_prepared.tsv'), 
      sep = '\t', 
      row.names = F, 
      dec = ',')
  }
}
#########################################
### EXTRACT IDS FROM COMPOSITE COLUMN ###
#########################################








################################################
### MANUAL - 3 5 7 8 12 14 15 18 26 58 38 30 ### 
#### 19 27 16 41 46 11 20 80 81 83 82 17 42 ####
################################################
pre_prep_reform[['manual']][['input']] <- openxlsx::read.xlsx(xlsxFile = 'Antidepressants_metadata_140520.xlsx', sheet = 'data')[1:25]


pre_prep_reform[['manual']][['input']] <- verify_df(
  df_ = pre_prep_reform[['manual']][['input']],
  repeated_spaces_ = T,
  trailing_spaces_ = T, 
  character_NAs_ = T,
  change_to_lower_ = F, 
  to_ascii_ = F)[['df']] # This helps with errors during further manual df processing


colnames(pre_prep_reform[['manual']][['input']]) == as.character(opts$check_cols_for_input)


colnames(pre_prep_reform[['manual']][['input']]) <- as.character(opts$check_cols_for_input)


pre_prep_reform[['manual']][['input']][opts$pre_prep_reform$range_of_cols_with_IDs_in_output_start:opts$pre_prep_reform$range_of_cols_with_IDs_in_output_end] <- purrr::map_dfc(
  .x = pre_prep_reform[['manual']][['input']][opts$pre_prep_reform$range_of_cols_with_IDs_in_output_start:opts$pre_prep_reform$range_of_cols_with_IDs_in_output_end], 
  .f = as.character)

pre_prep_reform[['manual']][['input']][['Exp.']] <- as.character(pre_prep_reform[['manual']][['input']][['Exp.']])

pre_prep_reform[['manual']][['input']][['Comp']] <- as.character(pre_prep_reform[['manual']][['input']][['Comp']])

pre_prep_reform[['manual']][['input']][['ID']] <- as.character(pre_prep_reform[['manual']][['input']][['ID']])

# pre_prep_reform[['manual']][['input']]$Symbol[pre_prep_reform[['manual']][['input']]$Symbol == 'Ggh ‡'] <- 'Ggh'

pre_prep_reform[['manual']][['input']][['Symbol']] <- stringr::str_replace_all(
  string = prepare_noncorrupted_symbol_col(symbol_col_ = pre_prep_reform[['manual']][['input']][['Symbol']]), 
  pattern = ', ', 
  replacement = opts$pre_prep_reform$bind_sep)

pre_prep_reform[['manual']][['input_QA_checklist']] <- verify_df(
  df_ = pre_prep_reform[['manual']][['input']], 
  only_qa = T)

pre_prep_reform[['manual']][['output']] <- pre_prep_reform[['manual']][['input']]


pre_prep_reform[['manual']][['sep_all_identifiers']] <- as.data.frame(
  stringr::str_split_fixed(
    string = pre_prep_reform[['manual']][['input']][['Symbol']], ### !!! change this into colname
    pattern = opts$pre_prep_reform$bind_sep,
    n = Inf))

pre_prep_reform[['manual']][['output']][['Symbol']] <- NA
pre_prep_reform[['manual']][['output']][['Symbol']] <- as.character(pre_prep_reform[['manual']][['output']][['Symbol']])

pre_prep_reform[['manual']][['output']] <- extract_identifers(
  sep_all_identifiers_ = pre_prep_reform[['manual']][['sep_all_identifiers']], 
  search_by = opts$pre_prep_reform[['search_bys']][['search_by']], 
  output_ = pre_prep_reform[['manual']][['output']], 
  opts_pre_check_cols = opts[['check_cols_for_input']], 
  range_of_cols_with_IDs_start = opts[['pre_prep_reform']][['range_of_cols_with_IDs_in_output_start']], 
  range_of_cols_with_IDs_end = opts[['pre_prep_reform']][['range_of_cols_with_IDs_in_output_end']])


pre_prep_reform[['manual']][['output_QA_checklist']] <- verify_df(df_ = pre_prep_reform[['manual']][['output']], only_qa = T)


write.table(
  x = pre_prep_reform[['manual']][['output']], 
  file = paste0(opts$pre_prep_reform$folder, '/manual_prepared.tsv'),
  sep = '\t', 
  row.names = F, 
  dec = ',')
################################################
### MANUAL - 3 5 7 8 12 14 15 18 26 58 38 30 ### 
#### 19 27 16 41 46 11 20 80 81 83 82 17 42 ####
################################################










############
### TEST ###
############
pre_prep_reform[['check_i/o']] <- purrr::map2(
  .x = pre_prep_reform, 
  .y = names(pre_prep_reform), 
  .f = function(platform, pl_names){
    tryCatch({
      compare_sample_rows_in_i_vs_o(
        in_df = platform[['input']], 
        out_df = platform[['output']], 
        analysis_name = pl_names, 
        write_to_folder = opts$pre_prep_reform$folder,
        drop_o_cols = c("adj,P,Val", "P,Value", "t", "B", "Description", "Protein", "Unknown"))
    }, error = function(cond) { return(NA) }
    )})
############
### TEST ###
############





#########################
### CONSTRUCT DATASET ###
#########################
pre_prep_reform$combined$data <- rlist::list.rbind(
  .data = pre_prep_reform[names(pre_prep_reform) %in% c(opts$pre_prep_reform[['which_platform_to_get_exp_and_comp_nb_for_regex']], 'manual')])


### !!! check can You bind data with non-compatible colnames
test2 <- test
colnames(test2)[[1]] <- 'ass'

rlist::list.rbind(list(test, test))
rlist::list.rbind(list(test, test2))
### !!! check can You bind data with non-compatible colnames


pre_prep_reform$combined$data <- verify_df(
  df_ = pre_prep_reform$combined$data, 
  sort_by_col = F, 
  repeated_spaces_ = T, 
  trailing_spaces_ = T, 
  character_NAs_ = T, 
  change_to_lower_ = F, 
  to_ascii_ = F)

pre_prep_reform$combined$qa <- verify_df(df_ = pre_prep_reform$combined$data, only_qa = T)

#Remove columns without single value
pre_prep_reform$combined$data[['NP_']] <- NULL
pre_prep_reform$combined$data[['XP_']] <- NULL
pre_prep_reform$combined$data[['Protein']] <- NULL
pre_prep_reform$combined$data[['Unknown']] <- NULL


# Checked - ID, Probe, NR_, Nucleotide, Unigene, ENSG_, ENST_, Gene_ID, Accession
pre_prep_reform$combined$data[['NM_']] <- stringr::str_remove_all(string = pre_prep_reform$combined$data[['NM_']], pattern = '\\.[0-9]{1,}')
pre_prep_reform$combined$data[['XM_']] <- stringr::str_remove_all(string = pre_prep_reform$combined$data[['XM_']], pattern = '\\.[0-9]{1,}')
pre_prep_reform$combined$data[['XR_']] <- stringr::str_remove_all(string = pre_prep_reform$combined$data[['XR_']], pattern = '\\.[0-9]{1,}')

pre_prep_reform$combined$data[['Symbol']] <- input_data_symbol_uncorruption_wrapper(manual_input_symbol = pre_prep_reform$combined$data[['Symbol']])

pre_prep_reform$combined$data <- pre_prep_reform$combined$data[order(pre_prep_reform$combined$data_reformated[['Pub.']], pre_prep_reform$combined$data_[['Exp.']]),]

pre_prep_reform$combined$qa_2 <- verify_df(df_ = pre_prep_reform$combined$data, only_qa = T)

write.table(
  x = pre_prep_reform$combined$data, 
  file = 'input_data_for_analysis.tsv', 
  sep = '\t', 
  dec = ',', 
  row.names = F)

#########################
### CONSTRUCT DATASET ###
#########################

save(pre_prep_reform, file = paste0(opts$pre_prep_reform$folder, '/pre_prep_reform'))
save(opts$pre_prep_reform, file = paste0(opts$pre_prep_reform$folder, '/opts_pre_prep_reform'))


########################
### GPL6885 - 28, 68 ###
########################
# colnames(pre_prep_reform$GPL6885$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6885$input$ID) # pre_prep_reform$GPL6885$temp_output$Probe - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6885$input$`Gene.symbol`) # + - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6885$input$`Gene.title`) # pre_prep_reform$GPL6885$temp_output$Description - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6885$input$`Gene.ID`) # pre_prep_reform$GPL6885$temp_output$Gene_ID (entrez) - ok
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6885$input$`UniGene.title`) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6885$input$`UniGene.symbol`) # + - ok
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6885$input$`UniGene.ID`) # pre_prep_reform$GPL6885$temp_output$Unigene - ok
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6885$input$GI) # pre_prep_reform$GPL6885$temp_output$Nucleotide - ok
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6885$input$`GenBank.Accession`) # + - ok
########################
### GPL6885 - 28, 68 ###
########################

########################
### GPL6887 - 33, 72 ###
########################
# colnames(pre_prep_reform$GPL6887$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6887$input$ID) # pre_prep_reform$GPL6887$temp_output$Probe - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6887$input$`Gene.symbol`) # ### !!! additional divider - '///' + - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6887$input$`Gene.title`) # pre_prep_reform$GPL6887$temp_output$Description
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6887$input$`Gene.ID`) #  additional divider - '///'  pre_prep_reform$GPL6887$temp_output$Gene_ID - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6887$input$`UniGene.title`) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6887$input$`UniGene.symbol`) # + - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6887$input$`UniGene.ID`) # pre_prep_reform$GPL6887$temp_output$Unigene - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6887$input$`Nucleotide.Title`) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6887$input$GI) # pre_prep_reform$GPL6887$temp_output$Nucleotide - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6887$input$`GenBank.Accession`) # + - OK
########################
### GPL6887 - 33, 72 ###
########################

#####################
### GPL17223 - 36 ###
#####################
# colnames(pre_prep_reform$GPL17223$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL17223$input$ID) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL17223$input$Search_key) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL17223$input$ProbeId) # pre_prep_reform$GPL17223$temp_output$Probe - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL17223$input$Gid) # pre_prep_reform$GPL17223$temp_output$Nucleotide - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL17223$input$Transcript) # + - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL17223$input$GB_ACC) # + - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL17223$input$Symbol) # + - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL17223$input$Definition) # pre_prep_reform$GPL17223$temp_output$Description - OK
#####################
### GPL17223 - 36 ###
#####################

################################
### GPL1261 - 66, 70, 77, 78 ###
################################
# colnames(pre_prep_reform$GPL1261$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL1261$input$ID) # pre_prep_reform$GPL1261$temp_output$Probe # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL1261$input$`Gene.symbol`) # ### !!! additional divider - '///' + # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL1261$input$`Gene.title`) # pre_prep_reform$GPL1261$temp_output$Description # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL1261$input$`Gene.ID`) #  additional divider - '///'  pre_prep_reform$GPL1261$temp_output$Gene_ID # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL1261$input$`UniGene.symbol`) # +
# get_all_symbols_in_chrvec(pre_prep_reform$GPL1261$input$`UniGene.ID`) # pre_prep_reform$GPL1261$temp_output$Unigene# - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL1261$input$`Nucleotide.Title`) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL1261$input$GI) # pre_prep_reform$GPL1261$temp_output$Nucleotide # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL1261$input$`GenBank.Accession`) # + # - OK
################################
### GPL1261 - 66, 70, 77, 78 ###
################################

####################
### GPL6246 - 67 ###
####################
# colnames(pre_prep_reform$GPL6246$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6246$input$ID) # pre_prep_reform$GPL6246$temp_output$Probe # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6246$input$`Gene.symbol`) # ### !!! additional divider - '///' + # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6246$input$`Gene.title`) # pre_prep_reform$GPL6246$temp_output$Description # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6246$input$`Gene.ID`) # ### !!! additional divider - '///' pre_prep_reform$GPL6246$output$Gene_ID # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6246$input$`UniGene.symbol`) # ### !!! additional divider - '///' + # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6246$input$`UniGene.ID`) # ### !!! additional divider - '///' pre_prep_reform$GPL6246$temp_output$Unigene # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6246$input$`Nucleotide.Title`) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6246$input$GI) # ### !!! additional divider - '///' pre_prep_reform$GPL6246$temp_output$Nucleotide # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL6246$input$`GenBank.Accession`) # ### !!! additional divider - '///' +
####################
### GPL6246 - 67 ###
####################

####################
### GPL5425 - 71 ###
####################
# colnames(pre_prep_reform$GPL5425$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL5425$input$ID) # remove _P.* and + # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL5425$input$`Gene.symbol`) # ### !!! additional divider - '///' + # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL5425$input$`Gene.title`) # pre_prep_reform$GPL5425$temp_output$Description # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL5425$input$`Gene.ID`) # ### !!! additional divider - '///' pre_prep_reform$GPL5425$temp_output$Gene_ID
# get_all_symbols_in_chrvec(pre_prep_reform$GPL5425$input$`UniGene.title`) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL5425$input$`UniGene.symbol`) # +
# get_all_symbols_in_chrvec(pre_prep_reform$GPL5425$input$`UniGene.ID`) # pre_prep_reform$GPL5425$temp_output$Unigene # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL5425$input$`Nucleotide.Title`) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL5425$input$GI) # pre_prep_reform$GPL5425$temp_output$Nucleotide # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL5425$input$`GenBank.Accession`) # + # - OK
####################
### GPL5425 - 71 ###
####################

#####################
### GPL25480 - 76 ###
#####################
# colnames(pre_prep_reform$GPL25480$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL25480$input$ID) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL25480$input$SPOT_ID) # - OK
#####################
### GPL25480 - 76 ###
#####################

#########################
### GPL10427 - 60, 45 ###
#########################
# colnames(pre_prep_reform$GPL10427$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL10427$input$ID) # OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL10427$input$accessions) # for extraction, | as devider, remove also 'ug|ref'
# get_all_symbols_in_chrvec(pre_prep_reform$GPL10427$input$ProbeUID) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL10427$input$ProbeName) # pre_prep_reform$GPL10427$temp_output$Probe
# get_all_symbols_in_chrvec(pre_prep_reform$GPL10427$input$GeneName) # for extraction
# get_all_symbols_in_chrvec(pre_prep_reform$GPL10427$input$SystematicName) # for extraction
# get_all_symbols_in_chrvec(pre_prep_reform$GPL10427$input$Description) # pre_prep_reform$GPL10427$temp_output$Description
# get_all_symbols_in_chrvec(pre_prep_reform$GPL10427$input$GB_LIST) # for extraction,  to 
#########################
### GPL10427 - 60, 45 ###
#########################

####################
### GPL8160 - 34 ###
####################
# colnames(pre_prep_reform$GPL8160$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL8160$input$ID) # pre_prep_reform$GPL8160$temp_output$Probe
# get_all_symbols_in_chrvec(pre_prep_reform$GPL8160$input$GB_ACC) # for extraction
####################
### GPL8160 - 34 ###
####################

#########################
### GPL16570 - 56, 50 ###
#########################
# unique(pre_prep_reform$GPL16570$input$Pub.)
# colnames(pre_prep_reform$GPL16570$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL16570$input$ID) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL16570$input$probeset_id) # pre_prep_reform$GPL16570$temp_output$Probe
# get_all_symbols_in_chrvec(pre_prep_reform$GPL16570$input$gene_assignment) # 
# get_all_symbols_in_chrvec(pre_prep_reform$GPL16570$input$mrna_assignment) # IGNORE(Too many wierd symbols, uncortrolable)
# get_all_symbols_in_chrvec(pre_prep_reform$GPL16570$input$swissprot) # IGNORE(Too many wierd symbols, uncortrolable)
# get_all_symbols_in_chrvec(pre_prep_reform$GPL16570$input$unigene) # IGNORE(Too many wierd symbols, uncortrolable)
#########################
### GPL16570 - 56, 50 ###
#########################

#####################
### GPL13912 - 52 ###
#####################
# colnames(pre_prep_reform$GPL13912$input) # - OK
# get_all_symbols_in_chrvec(pre_prep_reform$GPL13912$input$ID) # IGNORE
# get_all_symbols_in_chrvec(pre_prep_reform$GPL13912$input$NAME) # pre_prep_reform$GPL13912$temp_output$Probe 
# get_all_symbols_in_chrvec(pre_prep_reform$GPL13912$input$GB_ACC) # extraction
# get_all_symbols_in_chrvec(pre_prep_reform$GPL13912$input$GENE_ID) # pre_prep_reform$GPL13912$temp_output$Gene_ID
# get_all_symbols_in_chrvec(pre_prep_reform$GPL13912$input$GENE_SYMBOL) # both pre_prep_reform$GPL13912$temp_output$Symbol and extraction for LOCs
# get_all_symbols_in_chrvec(pre_prep_reform$GPL13912$input$GENE_NAME) # pre_prep_reform$GPL13912$temp_output$Description
# get_all_symbols_in_chrvec(pre_prep_reform$GPL13912$input$UNIGENE_ID) # pre_prep_reform$GPL13912$temp_output$Unigene
# get_all_symbols_in_chrvec(pre_prep_reform$GPL13912$input$ENSEMBL_ID) # pre_prep_reform$GPL13912$temp_output$ENST_
# get_all_symbols_in_chrvec(pre_prep_reform$GPL13912$input$ACCESSION_STRING) # extraction, but needs symbol removal
# get_all_symbols_in_chrvec(pre_prep_reform$GPL13912$input$DESCRIPTION) # IGNORE
#####################
### GPL13912 - 52 ###
#####################

# get_all_symbols_in_chrvec(manual$input$ID) # OK
# get_all_symbols_in_chrvec(manual$input$Probe) # OK
# get_all_symbols_in_chrvec(manual$input$Symbol) #
# get_all_symbols_in_chrvec(manual$input$Description) # OK
# get_all_symbols_in_chrvec(manual$input$NM_) # OK
# get_all_symbols_in_chrvec(manual$input$Accession) # OK
# get_all_symbols_in_chrvec(manual$input$Gene_ID) # OK
# get_all_symbols_in_chrvec(manual$input$ENST_) # OK
# get_all_symbols_in_chrvec(manual$input$ENSG_) # OK
# get_all_symbols_in_chrvec(manual$input$XM_) # remove
# get_all_symbols_in_chrvec(manual$input$XR_) # OK





#####################################################
### Platform-specific steps before input creation ###
#####################################################
### GPL6885 - 28, 68 ###
# for (nb in seq_along(pre_prep_reform[['GPL6885']][['pre-input']])) { 
#   pre_prep_reform[['GPL6885']][['pre-input']][[nb]][['Nucleotide.Title']] <- NULL }
### GPL6885 - 28, 68 ###


### GPL1261 - 66, 70, 77, 78 ###
# for (nb in seq_along(pre_prep_reform[['GPL1261']][['pre-input']])) { 
#   pre_prep_reform[['GPL1261']][['pre-input']][[nb]][['UniGene.title']] <- NULL }

# pre_prep_reform[['GPL1261']][['pre-input']][[6]] <- tibble::add_column(.data = pre_prep_reform[['GPL1261']][['pre-input']][[6]], 'GenBank.Accession' = NA, .after = 'GI')
### GPL1261 - 66, 70, 77, 78 ###
#####################################################
### Platform-specific steps before input creation ###
#####################################################