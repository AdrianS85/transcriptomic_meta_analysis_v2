# rm(list = ls(pattern = 'temp.*|test.*'))
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
source('functions_pre_preparation.R')
source('opts.R')
library(loggit)

opts_pre_prep <- list()
pre_prep_reform <- list()
dir.create(opts$pre_prep_reform$folder)
set_logfile(paste0(opts$pre_prep_reform$folder, "/loggit.log"))


###################
### SET OPTIONS ###
###################
opts_pre_prep <- list(
  'GPL6885' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^28.*tsv|^68.*tsv')),
  'GPL6887' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^33.*tsv|^72.*tsv')),
  'GPL17223' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^36.*tsv')),
  'GPL1261' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^66.*tsv|^70.*tsv|^77.*tsv|^78.*tsv')),
  'GPL6246' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^67.*tsv')),
  'GPL5425' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^71.*tsv')),
  'GPL25480' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^76.*tsv')),
  'GPL10427' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^45.*tsv|^60.*tsv')),
  'GPL8160' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^34.*tsv')),
  'GPL16570' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^56.*tsv|^50.*tsv')),
  'GPL13912' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '^52.*tsv')),
  'which_platform_to_get_pap_and_exp_nb_for_regex' = '^GPL.*',
  'exp_names_extraction_regex' = '^[[0-9]]{1,3}',
  'comp_names_extraction_regex' = '^[[0-9]]{1,3}_[[0-9]]{1,3}',
  'exp_comp_divider_symbol' = '_',
  'manual' = list('input' = openxlsx::read.xlsx(xlsxFile = 'Antidepressants_metadata_140520.xlsx', sheet = 'data')[1:25])
  )

opts_pre_prep[['platform_to_get_pap_and_exp_nb_for']] <- names(opts_pre_prep[stringr::str_detect(string = names(opts_pre_prep), pattern = opts_pre_prep[['which_platform_to_get_pap_and_exp_nb_for_regex']])])



for (platform in opts_pre_prep[['platform_to_get_pap_and_exp_nb_for']]) 
{
  opts_pre_prep[[platform]][['exp_names']] <- stringr::str_extract(
    string = opts_pre_prep[[platform]][['file_names']], 
    pattern = opts_pre_prep[['exp_names_extraction_regex']])
  
  opts_pre_prep[[platform]][['comp_names']] <- stringr::str_remove(
    string = stringr::str_extract(
      string = opts_pre_prep[[platform]][['file_names']],
      pattern = opts_pre_prep[['comp_names_extraction_regex']]),
    pattern = paste0('.*', opts_pre_prep[['exp_comp_divider_symbol']]))
}
###################
### SET OPTIONS ###
###################





######################
### LOAD PRE-INPUT ###
######################
pre_prep_reform <- purrr::map(
  .x = opts_pre_prep[names(opts_pre_prep) %in% opts_pre_prep[['platform_to_get_pap_and_exp_nb_for']]], 
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





#####################################################
### Platform-specific steps before input creation ###
#####################################################
### GPL6885 - 28, 68 ###
for (nb in seq_along(pre_prep_reform[['GPL6885']][['pre-input']])) { 
  pre_prep_reform[['GPL6885']][['pre-input']][[nb]][['Nucleotide.Title']] <- NULL }
### GPL6885 - 28, 68 ###


### GPL1261 - 66, 70, 77, 78 ###
# for (nb in seq_along(pre_prep_reform[['GPL1261']][['pre-input']])) { 
#   pre_prep_reform[['GPL1261']][['pre-input']][[nb]][['UniGene.title']] <- NULL }

pre_prep_reform[['GPL1261']][['pre-input']][[6]] <- tibble::add_column(.data = pre_prep_reform[['GPL1261']][['pre-input']][[6]], 'GenBank.Accession' = NA, .after = 'GI')
### GPL1261 - 66, 70, 77, 78 ###
#####################################################
### Platform-specific steps before input creation ###
#####################################################





#####################
### PREPARE INPUT ###
#####################
opts_pre_prep[['pass_nb']] <- 0
for (platform in opts_pre_prep[['platform_to_get_pap_and_exp_nb_for']]) {
  message(sprintf('processing platform %s...', platform))
  
  message('binding pre-input into single df...')
  pre_prep_reform[[platform]][['input']] <- rlist::list.rbind(pre_prep_reform[[platform]][['pre-input']])
  
  
  
  message('performing qa...')
  pre_prep_reform[[platform]][['qa_input']] <- verify_df(
    df_ = pre_prep_reform[[platform]][['input']], 
    only_qa = T)
  

  
  message('setting output df to be appended with reformated identifers...')
  
  if (opts_pre_prep[['pass_nb']] == 0) { message('NOTE!!! Need to set -temp_output-, cause set_search_bys function builds its own -output- df from input, which in case of data downloaded from GEO via R is not compatible with further functions. Lol. Poor software engineering again.') }
  
  pre_prep_reform[[platform]][['temp_output']] <- set_df_for_append(
    check_cols_ = opts$check_cols_for_input, 
    length = length(pre_prep_reform[[platform]][['input']][[1]]))
  
  
  
  message('copying standard columns returned by geo2r into output df...')
  
  if (opts_pre_prep[['pass_nb']] == 0) { 
    message(sprintf('NOTE!!! columns copied are: %s', copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = opts[['check_cols_for_input']], return_copied_col_names = T))) 
  }
  
  pre_prep_reform[[platform]][['temp_output']] <- copy_out_of_box_ready_columns_into_final_output_wrapper(
    check_cols_ = opts$check_cols_for_input, 
    final_output_ = pre_prep_reform[[platform]][['temp_output']], 
    input_ = pre_prep_reform[[platform]][['input']])
  
  
  
  opts_pre_prep[['pass_nb']] <- opts_pre_prep[['pass_nb']] + 1
  message(sprintf('platform %s processed.', platform))
}
#####################
### PREPARE INPUT ###
#####################





############################################
### Platform-specific column reformating ###
############################################
### GPL6885 - 28, 68 ###
pre_prep_reform$GPL6885$temp_output$Probe <- pre_prep_reform$GPL6885$input$ID # - OK
pre_prep_reform$GPL6885$temp_output$Description <- pre_prep_reform$GPL6885$input$`Gene.title` # - OK
pre_prep_reform$GPL6885$temp_output$Gene_ID <- pre_prep_reform$GPL6885$input$`Gene.ID` # - OK
pre_prep_reform$GPL6885$temp_output$Unigene <- pre_prep_reform$GPL6885$input$`UniGene.ID` # - OK
pre_prep_reform$GPL6885$temp_output$Nucleotide <- pre_prep_reform$GPL6885$input$GI # - OK

pre_prep_reform$GPL6885$input$bind_for_extraction <- paste(pre_prep_reform$GPL6885$input$`Gene.symbol`, pre_prep_reform$GPL6885$input$`UniGene.symbol`, pre_prep_reform$GPL6885$input$`GenBank.Accession`, sep = '___') # - OK
### GPL6885 - 28, 68 ###

### GPL6887 - 33, 72 ###
pre_prep_reform$GPL6887$temp_output$Probe <- pre_prep_reform$GPL6887$input$ID # - OK
pre_prep_reform$GPL6887$temp_output$Description <- pre_prep_reform$GPL6887$input$`Gene.title` # - OK
pre_prep_reform$GPL6887$temp_output$Gene_ID <- stringr::str_replace_all(string = pre_prep_reform$GPL6887$input$`Gene.ID`, pattern = '///', replacement = ', ') # - OK
pre_prep_reform$GPL6887$temp_output$Unigene <- pre_prep_reform$GPL6887$input$`UniGene.ID` # - OK
pre_prep_reform$GPL6887$temp_output$Nucleotide <- pre_prep_reform$GPL6887$input$GI # - OK

pre_prep_reform$GPL6887$input$bind_for_extraction <- paste(stringr::str_replace_all(string = pre_prep_reform$GPL6887$input$`Gene.symbol`, pattern = '///', replacement = '___'), pre_prep_reform$GPL6887$input$`UniGene.symbol`, pre_prep_reform$GPL6887$input$`GenBank.Accession`, sep = '___') # - OK
### GPL6887 - 33, 72 ###

### GPL17223 - 36 ###
pre_prep_reform$GPL17223$temp_output$Probe <- pre_prep_reform$GPL17223$input$ProbeId # - OK
pre_prep_reform$GPL17223$temp_output$Nucleotide <- pre_prep_reform$GPL17223$input$Gid # - OK
pre_prep_reform$GPL17223$temp_output$Description <- pre_prep_reform$GPL17223$input$Definition # - OK

pre_prep_reform$GPL17223$input$bind_for_extraction <- paste(pre_prep_reform$GPL17223$input$Symbol, pre_prep_reform$GPL17223$input$Transcript, pre_prep_reform$GPL17223$input$GB_ACC, sep = '___') # - OK
### GPL17223 - 36 ###

### GPL1261 - 66, 70, 77, 78 ###
pre_prep_reform$GPL1261$temp_output$Probe <- pre_prep_reform$GPL1261$input$ID # - OK
pre_prep_reform$GPL1261$temp_output$Description <- pre_prep_reform$GPL1261$input$`Gene.title` # - OK
pre_prep_reform$GPL1261$temp_output$Gene_ID <- stringr::str_replace_all(string = pre_prep_reform$GPL1261$input$`Gene.ID`, pattern = '///', replacement = ', ') # - OK
pre_prep_reform$GPL1261$temp_output$Unigene <- pre_prep_reform$GPL1261$input$`UniGene.ID` # - OK
pre_prep_reform$GPL1261$temp_output$Nucleotide <- pre_prep_reform$GPL1261$input$GI # - OK

pre_prep_reform$GPL1261$input$bind_for_extraction <- paste(stringr::str_replace_all(string = pre_prep_reform$GPL1261$input$`Gene.symbol`, pattern = '///', replacement = '___'), pre_prep_reform$GPL1261$input$`UniGene.symbol`, pre_prep_reform$GPL1261$input$`GenBank.Accession`, sep = '___') # - OK
### GPL1261 - 66, 70, 77, 78 ###

### GPL6246 - 67 ###
pre_prep_reform$GPL6246$temp_output$Probe <- pre_prep_reform$GPL6246$input$ID # - OK
pre_prep_reform$GPL6246$temp_output$Description <- pre_prep_reform$GPL6246$input$`Gene.title` # - OK
pre_prep_reform$GPL6246$temp_output$Gene_ID <- stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$`Gene.ID`, pattern = '///', replacement = ', ') # - OK
pre_prep_reform$GPL6246$temp_output$Unigene <- stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$`UniGene.ID`, pattern = '///', replacement = ', ') # - OK
pre_prep_reform$GPL6246$temp_output$Nucleotide <- stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$GI, pattern = '///', replacement = ', ') # - OK

pre_prep_reform$GPL6246$input$bind_for_extraction <- paste(stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$`Gene.symbol`, pattern = '///', replacement = '___'), stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$`UniGene.symbol`, pattern = '///', replacement = '___'), stringr::str_replace_all(string = pre_prep_reform$GPL6246$input$`GenBank.Accession`, pattern = '///', replacement = '___'), sep = '___') # - OK
### GPL6246 - 67 ###

### GPL5425 - 71 ###
pre_prep_reform$GPL5425$temp_output$Description <- pre_prep_reform$GPL5425$input$`Gene.title` # - OK
pre_prep_reform$GPL5425$temp_output$Gene_ID <- stringr::str_replace_all(string = pre_prep_reform$GPL5425$input$`Gene.ID`, pattern = '///', replacement = ', ') # - OK
pre_prep_reform$GPL5425$temp_output$Unigene <- pre_prep_reform$GPL5425$input$`UniGene.ID` # - OK
pre_prep_reform$GPL5425$temp_output$Nucleotide <- pre_prep_reform$GPL5425$input$GI # - OK

pre_prep_reform$GPL5425$input$bind_for_extraction <- paste(stringr::str_remove(string = pre_prep_reform$GPL5425$input$ID, pattern = '_P.*'), stringr::str_replace_all(string = pre_prep_reform$GPL5425$input$`Gene.symbol`, pattern = '///', replacement = '___'), pre_prep_reform$GPL5425$input$`GenBank.Accession`, pre_prep_reform$GPL5425$input$`UniGene.symbol`, sep = '___')
### GPL5425 - 71 ###

### GPL25480 - 76 ###
pre_prep_reform$GPL25480$temp_output$Probe <- pre_prep_reform$GPL25480$input$SPOT_ID # - OK
### GPL25480 - 76 ###

### GPL10427 - 60, 45 ###
pre_prep_reform$GPL10427$temp_output$Probe <- pre_prep_reform$GPL10427$input$ProbeName # 
pre_prep_reform$GPL10427$temp_output$Description <- pre_prep_reform$GPL10427$input$Description #

pre_prep_reform$GPL10427$input$bind_for_extraction <- paste(
  prepare_uncorrupted_GPL10427_accessions_wrapper(GPL10427_accessions = pre_prep_reform$GPL10427$input$accessions), 
  pre_prep_reform$GPL10427$input$GeneName, 
  pre_prep_reform$GPL10427$input$SystematicName, 
  stringr::str_replace_all(string = pre_prep_reform$GPL10427$input$GB_LIST, pattern = ', ', replacement = '___'),
  sep = '___') # - OK
### GPL10427 - 60, 45 ###

### GPL8160 - 34 ###
pre_prep_reform$GPL8160$temp_output$Probe <- pre_prep_reform$GPL8160$input$ID # 

pre_prep_reform$GPL8160$input$bind_for_extraction <- pre_prep_reform$GPL8160$input$GB_ACC # 
### GPL8160 - 34 ###

### GPL16570 - 56, 50 ###
pre_prep_reform$GPL16570$temp_output$Probe <- pre_prep_reform$GPL16570$input$probeset_id #

pre_prep_reform$GPL16570$input$bind_for_extraction <- paste(
  prepare_noncorrupted_GPL16570_cols(GPL16570_cols = pre_prep_reform$GPL16570$input$gene_assignment),
  sep = '___')
### GPL16570 - 56, 50 ###

### GPL13912 - 52 ###
pre_prep_reform$GPL13912$temp_output$Description <- pre_prep_reform$GPL13912$input$GENE_NAME # 
pre_prep_reform$GPL13912$temp_output$Probe <- pre_prep_reform$GPL13912$input$NAME # 
pre_prep_reform$GPL13912$temp_output$Gene_ID <- pre_prep_reform$GPL13912$input$GENE_ID
pre_prep_reform$GPL13912$temp_output$Symbol <- pre_prep_reform$GPL13912$input$GENE_SYMBOL
pre_prep_reform$GPL13912$temp_output$Unigene <- pre_prep_reform$GPL13912$input$UNIGENE_ID
pre_prep_reform$GPL13912$temp_output$ENST_ <- pre_prep_reform$GPL13912$input$ENSEMBL_ID

pre_prep_reform$GPL13912$input$bind_for_extraction <- paste(
  GPL13912_cols = pre_prep_reform$GPL13912$input$GB_ACC,
  pre_prep_reform$GPL13912$input$GENE_SYMBOL,
  prepare_uncorrupted_GPL13912_ACCESSION_STRING_wrapper(pre_prep_reform$GPL13912$input$ACCESSION_STRING),
  sep = '___')
### GPL13912 - 52 ###
############################################
### Platform-specific column reformating ###
############################################





####################################################
### PREPARE FOR EXTRACTION FROM COMPOSITE COLUMN ###
####################################################
for (platform in opts_pre_prep[['platform_to_get_pap_and_exp_nb_for']]) 
{
  message(sprintf('processing platform %s...', platform))
  
  if (!is.null(pre_prep_reform[[platform]][['input']][['bind_for_extraction']])) {

    col_with_ids_for_extraction <- length(pre_prep_reform[[platform]][['input']])
    
    message(sprintf('column bind_for_extraction set to position %i', col_with_ids_for_extraction))
    
    
    
    message('separating identifiers from single column...')
    pre_prep_reform[[platform]][['sep_all_identifiers']] <- get_sep_all_identifiers(
      input_df_name = pre_prep_reform[[platform]][['input']], 
      bad_cols_to_process_seq = seq(col_with_ids_for_extraction,col_with_ids_for_extraction), 
      GPL_ID = 'R')
    
    
    
    message('merging...')
    pre_prep_reform[[platform]] <- rlist::list.merge(
      pre_prep_reform[[platform]], 
      set_search_bys(
        input_df_ = pre_prep_reform[[platform]][['input']], 
        search_by_p1_ = "stringr::str_detect(string = ", 
        search_by_n_ = 'identifer', 
        check_cols_ = opts$check_cols_for_input))
  }
  
  message('copying temporary output to output...')
  pre_prep_reform[[platform]][['output']] <- pre_prep_reform[[platform]][['temp_output']]
  
  message(sprintf('platform %s processed', platform))
}
####################################################
### PREPARE FOR EXTRACTION FROM COMPOSITE COLUMN ###
####################################################





#########################################
### EXTRACT IDS FROM COMPOSITE COLUMN ###
#########################################
for (platform in opts_pre_prep[['platform_to_get_pap_and_exp_nb_for']]) 
{
  if (!is.null(pre_prep_reform[[platform]][['input']][['bind_for_extraction']])) {
    
    pre_prep_reform[[platform]]output <- extract_identifers(
      sep_all_identifiers_ = pre_prep_reform[[platform]][['sep_all_identifiers']], 
      search_by_ = pre_prep_reform[[platform]][['search_by']], 
      output_ = pre_prep_reform[[platform]][['output']], 
      opts_pre_check_cols = opts[['check_cols_for_input']], 
      range_of_cols_with_IDs_start = opts[['pre_prep_reform']][['range_of_cols_with_IDs_in_output_start']], 
      range_of_cols_with_IDs_end = opts[['pre_prep_reform']][['range_of_cols_with_IDs_in_output_end']])
    
    
    
    message('performing qa...')
    pre_prep_reform[[platform]][['output_qa']] <- verify_df(df_ = pre_prep_reform[[platform]][['output']], only_qa = T)
    
    
    
    message(sprintf('writing data to %s', opts$dir_r_downloaded_data))
    write.table(
      x = pre_prep_reform[[platform]][['output']], 
      file = paste0(opts$dir_r_downloaded_data, '/', platform,  '_prepared.tsv'), 
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
colnames(manual$input) <- as.character(opts$check_cols_for_input)

manual$input[opts$pre_prep_reform$range_of_cols_with_IDs_in_output_start:opts$pre_prep_reform$range_of_cols_with_IDs_in_output_end] <- purrr::map_dfc(.x = manual$input[opts$pre_prep_reform$range_of_cols_with_IDs_in_output_start:opts$pre_prep_reform$range_of_cols_with_IDs_in_output_end], .f = as.character)

manual$input$Symbol[manual$input$Symbol == 'Ggh ‡'] <- 'Ggh'

manual$input_QA_checklist <- get_input_QA_checklist(input_table_ = manual$input, Pub._col_ = 'Pub.')

manual$input$NM_ <- prepare_noncorrupted_NXMR_col(nxmr_col_ = manual$input$NM_)
manual$input$XM_ <- prepare_noncorrupted_NXMR_col(nxmr_col_ = manual$input$XM_)
manual$input$Accession <- stringr::str_remove(string = manual$input$Accession, pattern = ' $')
manual$input$Symbol <- stringr::str_replace_all(string = prepare_noncorrupted_symbol_col(symbol_col_ = manual$input$Symbol), pattern = ', ', replacement = '___')

manual$temp_output <- manual$input

manual$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = manual$input, bad_cols_to_process_seq = seq(10,10), GPL_ID = 'R')

manual <- rlist::list.merge(manual, set_search_bys(input_df_ = manual$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = opts$check_cols_for_input))

manual$output <- manual$temp_output

manual$output <- extract_identifers(sep_all_identifiers_ = manual$sep_all_identifiers, search_by_ = manual$search_by, output_ = manual$output, opts_pre_check_cols = opts$check_cols_for_input, range_of_cols_with_IDs_start = opts$pre_prep_reform$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = opts$pre_prep_reform$range_of_cols_with_IDs_in_output_end)

manual$output_cols_check <- verify_df(df_ = manual$output, only_qa = T)

write.table(
  x = manual$output, 
  file = paste0(opts$dir_r_downloaded_data, '/manual_prepared.tsv'), 
  sep = '\t', 
  row.names = F, 
  dec = ',')


################################################
### MANUAL - 3 5 7 8 12 14 15 18 26 58 38 30 ### 
#### 19 27 16 41 46 11 20 80 81 83 82 17 42 ####
################################################




save(pre_prep_reform, file = paste0(opts$pre_prep_reform$folder, '/pre_prep_reform'))
save(opts_pre_prep, file = paste0(opts$pre_prep_reform$folder, '/opts_pre_prep'))


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