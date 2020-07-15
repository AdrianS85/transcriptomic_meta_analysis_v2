# rm(list = ls(pattern = 'temp.*|test.*'))
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
source('functions_pre_preparation.R')
gen_opts_pre <- list(
  'dir' = 'geo/geo_r',
  'check_cols' = set_col_check(), 
  'decimal' = ',',
  'input_col_types' = 'cnnnnncccccccccccccccccccccccccccccccc',
  # 'numeric_col_types' = c(2:6),
  'range_of_cols_with_IDs_in_output_start' = 9,
  'range_of_cols_with_IDs_in_output_end' = 25
)





########################
### GPL6885 - 28, 68 ###
########################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^28.*tsv|^68.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL6885 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

#some files have this, some dont
for (nb in seq_along(GPL6885$`pre-input`)) { 
  GPL6885$`pre-input`[[nb]]$Nucleotide.Title <- NULL }

GPL6885$input <- rlist::list.rbind(GPL6885$`pre-input`)
# GPL6885$input <- GPL6885$input[1:10,]
#Need to set 'temp_output', cause set_search_bys function builds its own 'output' df from input, which in case of data downloaded from GEO via R is not compatible with further functions. Lol. Poor software engineering again.
GPL6885$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL6885$input[[1]]))

GPL6885$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL6885$temp_output, input_ = GPL6885$input)

GPL6885$temp_output$Probe <- GPL6885$input$ID # - OK
GPL6885$temp_output$Description <- GPL6885$input$`Gene.title` # - OK
GPL6885$temp_output$Gene_ID <- GPL6885$input$`Gene.ID` # - OK
GPL6885$temp_output$Unigene <- GPL6885$input$`UniGene.ID` # - OK
GPL6885$temp_output$Nucleotide <- GPL6885$input$GI # - OK

GPL6885$input$bind_for_extraction <- paste(GPL6885$input$`Gene.symbol`, GPL6885$input$`UniGene.symbol`, GPL6885$input$`GenBank.Accession`, sep = '___') # - OK

GPL6885$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = GPL6885$input, bad_cols_to_process_seq = seq(17,17), GPL_ID = 'R')

GPL6885 <- rlist::list.merge(GPL6885, set_search_bys(input_df_ = GPL6885$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

GPL6885$output <- GPL6885$temp_output

GPL6885$output <- extract_identifers(sep_all_identifiers_ = GPL6885$sep_all_identifiers, search_by_ = GPL6885$search_by, output_ = GPL6885$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)

GPL6885$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL6885$output)

write.table(x = GPL6885$output, file = paste0(gen_opts_pre$dir, '/GPL6885_prepared.tsv'), sep = '\t', row.names = F, dec = ',')

# colnames(GPL6885$input) # - OK
# get_all_symbols_in_chrvec(GPL6885$input$ID) # GPL6885$temp_output$Probe - OK
# get_all_symbols_in_chrvec(GPL6885$input$`Gene.symbol`) # + - OK
# get_all_symbols_in_chrvec(GPL6885$input$`Gene.title`) # GPL6885$temp_output$Description - OK
# get_all_symbols_in_chrvec(GPL6885$input$`Gene.ID`) # GPL6885$temp_output$Gene_ID (entrez) - ok
# get_all_symbols_in_chrvec(GPL6885$input$`UniGene.title`) # IGNORE
# get_all_symbols_in_chrvec(GPL6885$input$`UniGene.symbol`) # + - ok
# get_all_symbols_in_chrvec(GPL6885$input$`UniGene.ID`) # GPL6885$temp_output$Unigene - ok
# get_all_symbols_in_chrvec(GPL6885$input$GI) # GPL6885$temp_output$Nucleotide - ok
# get_all_symbols_in_chrvec(GPL6885$input$`GenBank.Accession`) # + - ok
########################
### GPL6885 - 28, 68 ###
########################





########################
### GPL6887 - 33, 72 ###
########################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^33.*tsv|^72.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL6887 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

GPL6887$input <- rlist::list.rbind(GPL6887$`pre-input`)

GPL6887$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL6887$input[[1]]))

GPL6887$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL6887$temp_output, input_ = GPL6887$input)

GPL6887$temp_output$Probe <- GPL6887$input$ID # - OK
GPL6887$temp_output$Description <- GPL6887$input$`Gene.title` # - OK
GPL6887$temp_output$Gene_ID <- stringr::str_replace_all(string = GPL6887$input$`Gene.ID`, pattern = '///', replacement = ', ') # - OK
GPL6887$temp_output$Unigene <- GPL6887$input$`UniGene.ID` # - OK
GPL6887$temp_output$Nucleotide <- GPL6887$input$GI # - OK

GPL6887$input$bind_for_extraction <- paste(stringr::str_replace_all(string = GPL6887$input$`Gene.symbol`, pattern = '///', replacement = '___'), GPL6887$input$`UniGene.symbol`, GPL6887$input$`GenBank.Accession`, sep = '___') # - OK

GPL6887$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = GPL6887$input, bad_cols_to_process_seq = seq(18,18), GPL_ID = 'R')

GPL6887 <- rlist::list.merge(GPL6887, set_search_bys(input_df_ = GPL6887$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

GPL6887$output <- GPL6887$temp_output

GPL6887$output <- extract_identifers(sep_all_identifiers_ = GPL6887$sep_all_identifiers, search_by_ = GPL6887$search_by, output_ = GPL6887$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)

GPL6887$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL6887$output)

write.table(x = GPL6887$output, file = paste0(gen_opts_pre$dir, '/GPL6887_prepared.tsv'), sep = '\t', row.names = F, dec = ',')


# colnames(GPL6887$input) # - OK
# get_all_symbols_in_chrvec(GPL6887$input$ID) # GPL6887$temp_output$Probe - OK
# get_all_symbols_in_chrvec(GPL6887$input$`Gene.symbol`) # ### !!! additional divider - '///' + - OK
# get_all_symbols_in_chrvec(GPL6887$input$`Gene.title`) # GPL6887$temp_output$Description
# get_all_symbols_in_chrvec(GPL6887$input$`Gene.ID`) #  additional divider - '///'  GPL6887$temp_output$Gene_ID - OK
# get_all_symbols_in_chrvec(GPL6887$input$`UniGene.title`) # IGNORE
# get_all_symbols_in_chrvec(GPL6887$input$`UniGene.symbol`) # + - OK
# get_all_symbols_in_chrvec(GPL6887$input$`UniGene.ID`) # GPL6887$temp_output$Unigene - OK
# get_all_symbols_in_chrvec(GPL6887$input$`Nucleotide.Title`) # IGNORE
# get_all_symbols_in_chrvec(GPL6887$input$GI) # GPL6887$temp_output$Nucleotide - OK
# get_all_symbols_in_chrvec(GPL6887$input$`GenBank.Accession`) # + - OK
########################
### GPL6887 - 33, 72 ###
########################





#####################
### GPL17223 - 36 ###
#####################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^36.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL17223 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

GPL17223$input <- rlist::list.rbind(GPL17223$`pre-input`)

GPL17223$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL17223$input[[1]]))

GPL17223$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL17223$temp_output, input_ = GPL17223$input)

GPL17223$temp_output$Probe <- GPL17223$input$ProbeId # - OK
GPL17223$temp_output$Nucleotide <- GPL17223$input$Gid # - OK
GPL17223$temp_output$Description <- GPL17223$input$Definition # - OK

GPL17223$input$bind_for_extraction <- paste(GPL17223$input$Symbol, GPL17223$input$Transcript, GPL17223$input$GB_ACC, sep = '___') # - OK

GPL17223$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = GPL17223$input, bad_cols_to_process_seq = seq(16,16), GPL_ID = 'R')

GPL17223 <- rlist::list.merge(GPL17223, set_search_bys(input_df_ = GPL17223$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

GPL17223$output <- GPL17223$temp_output

GPL17223$output <- extract_identifers(sep_all_identifiers_ = GPL17223$sep_all_identifiers, search_by_ = GPL17223$search_by, output_ = GPL17223$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)

GPL17223$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL17223$output)

write.table(x = GPL17223$output, file = paste0(gen_opts_pre$dir, '/GPL17223_prepared.tsv'), sep = '\t', row.names = F, dec = ',')


# colnames(GPL17223$input) # - OK
# get_all_symbols_in_chrvec(GPL17223$input$ID) # IGNORE
# get_all_symbols_in_chrvec(GPL17223$input$Search_key) # IGNORE
# get_all_symbols_in_chrvec(GPL17223$input$ProbeId) # GPL17223$temp_output$Probe - OK
# get_all_symbols_in_chrvec(GPL17223$input$Gid) # GPL17223$temp_output$Nucleotide - OK
# get_all_symbols_in_chrvec(GPL17223$input$Transcript) # + - OK
# get_all_symbols_in_chrvec(GPL17223$input$GB_ACC) # + - OK
# get_all_symbols_in_chrvec(GPL17223$input$Symbol) # + - OK
# get_all_symbols_in_chrvec(GPL17223$input$Definition) # GPL17223$temp_output$Description - OK
#####################
### GPL17223 - 36 ###
#####################




################################
### GPL1261 - 66, 70, 77, 78 ###
################################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^66.*tsv|^70.*tsv|^77.*tsv|^78.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL1261 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

for (nb in seq_along(GPL1261$`pre-input`)) { 
  GPL1261$`pre-input`[[nb]]$UniGene.title <- NULL }

GPL1261$`pre-input`[[6]] <- tibble::add_column(.data = GPL1261$`pre-input`[[6]], GenBank.Accession = NA, .after = 'GI')

GPL1261$input <- rlist::list.rbind(GPL1261$`pre-input`)

GPL1261$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL1261$input[[1]]))

GPL1261$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL1261$temp_output, input_ = GPL1261$input)

GPL1261$temp_output$Probe <- GPL1261$input$ID # - OK
GPL1261$temp_output$Description <- GPL1261$input$`Gene.title` # - OK
GPL1261$temp_output$Gene_ID <- stringr::str_replace_all(string = GPL1261$input$`Gene.ID`, pattern = '///', replacement = ', ') # - OK
GPL1261$temp_output$Unigene <- GPL1261$input$`UniGene.ID` # - OK
GPL1261$temp_output$Nucleotide <- GPL1261$input$GI # - OK

GPL1261$input$bind_for_extraction <- paste(stringr::str_replace_all(string = GPL1261$input$`Gene.symbol`, pattern = '///', replacement = '___'), GPL1261$input$`UniGene.symbol`, GPL1261$input$`GenBank.Accession`, sep = '___') # - OK

GPL1261$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = GPL1261$input, bad_cols_to_process_seq = seq(17,17), GPL_ID = 'R')

GPL1261 <- rlist::list.merge(GPL1261, set_search_bys(input_df_ = GPL1261$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

GPL1261$output <- GPL1261$temp_output

GPL1261$output <- extract_identifers(sep_all_identifiers_ = GPL1261$sep_all_identifiers, search_by_ = GPL1261$search_by, output_ = GPL1261$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)

GPL1261$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL1261$output)

write.table(x = GPL1261$output, file = paste0(gen_opts_pre$dir, '/GPL1261_prepared.tsv'), sep = '\t', row.names = F, dec = ',')


# colnames(GPL1261$input) # - OK
# get_all_symbols_in_chrvec(GPL1261$input$ID) # GPL1261$temp_output$Probe # - OK
# get_all_symbols_in_chrvec(GPL1261$input$`Gene.symbol`) # ### !!! additional divider - '///' + # - OK
# get_all_symbols_in_chrvec(GPL1261$input$`Gene.title`) # GPL1261$temp_output$Description # - OK
# get_all_symbols_in_chrvec(GPL1261$input$`Gene.ID`) #  additional divider - '///'  GPL1261$temp_output$Gene_ID # - OK
# get_all_symbols_in_chrvec(GPL1261$input$`UniGene.symbol`) # +
# get_all_symbols_in_chrvec(GPL1261$input$`UniGene.ID`) # GPL1261$temp_output$Unigene# - OK
# get_all_symbols_in_chrvec(GPL1261$input$`Nucleotide.Title`) # IGNORE
# get_all_symbols_in_chrvec(GPL1261$input$GI) # GPL1261$temp_output$Nucleotide # - OK
# get_all_symbols_in_chrvec(GPL1261$input$`GenBank.Accession`) # + # - OK
################################
### GPL1261 - 66, 70, 77, 78 ###
################################





####################
### GPL6246 - 67 ###
####################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^67.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL6246 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

GPL6246$input <- rlist::list.rbind(GPL6246$`pre-input`)

GPL6246$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL6246$input[[1]]))

GPL6246$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL6246$temp_output, input_ = GPL6246$input)

GPL6246$temp_output$Probe <- GPL6246$input$ID # - OK
GPL6246$temp_output$Description <- GPL6246$input$`Gene.title` # - OK
GPL6246$temp_output$Gene_ID <- stringr::str_replace_all(string = GPL6246$input$`Gene.ID`, pattern = '///', replacement = ', ') # - OK
GPL6246$temp_output$Unigene <- stringr::str_replace_all(string = GPL6246$input$`UniGene.ID`, pattern = '///', replacement = ', ') # - OK
GPL6246$temp_output$Nucleotide <- stringr::str_replace_all(string = GPL6246$input$GI, pattern = '///', replacement = ', ') # - OK

GPL6246$input$bind_for_extraction <- paste(stringr::str_replace_all(string = GPL6246$input$`Gene.symbol`, pattern = '///', replacement = '___'), stringr::str_replace_all(string = GPL6246$input$`UniGene.symbol`, pattern = '///', replacement = '___'), stringr::str_replace_all(string = GPL6246$input$`GenBank.Accession`, pattern = '///', replacement = '___'), sep = '___') # - OK

GPL6246$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = GPL6246$input, bad_cols_to_process_seq = seq(17,17), GPL_ID = 'R')

GPL6246 <- rlist::list.merge(GPL6246, set_search_bys(input_df_ = GPL6246$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

GPL6246$output <- GPL6246$temp_output

GPL6246$output <- extract_identifers(sep_all_identifiers_ = GPL6246$sep_all_identifiers, search_by_ = GPL6246$search_by, output_ = GPL6246$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)
 
GPL6246$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL6246$output)

write.table(x = GPL6246$output, file = paste0(gen_opts_pre$dir, '/GPL6246_prepared.tsv'), sep = '\t', row.names = F, dec = ',')


# colnames(GPL6246$input) # - OK
# get_all_symbols_in_chrvec(GPL6246$input$ID) # GPL6246$temp_output$Probe # - OK
# get_all_symbols_in_chrvec(GPL6246$input$`Gene.symbol`) # ### !!! additional divider - '///' + # - OK
# get_all_symbols_in_chrvec(GPL6246$input$`Gene.title`) # GPL6246$temp_output$Description # - OK
# get_all_symbols_in_chrvec(GPL6246$input$`Gene.ID`) # ### !!! additional divider - '///' GPL6246$output$Gene_ID # - OK
# get_all_symbols_in_chrvec(GPL6246$input$`UniGene.symbol`) # ### !!! additional divider - '///' + # - OK
# get_all_symbols_in_chrvec(GPL6246$input$`UniGene.ID`) # ### !!! additional divider - '///' GPL6246$temp_output$Unigene # - OK
# get_all_symbols_in_chrvec(GPL6246$input$`Nucleotide.Title`) # IGNORE
# get_all_symbols_in_chrvec(GPL6246$input$GI) # ### !!! additional divider - '///' GPL6246$temp_output$Nucleotide # - OK
# get_all_symbols_in_chrvec(GPL6246$input$`GenBank.Accession`) # ### !!! additional divider - '///' +
####################
### GPL6246 - 67 ###
####################





####################
### GPL5425 - 71 ###
####################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^71.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL5425 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

GPL5425$input <- rlist::list.rbind(GPL5425$`pre-input`)

GPL5425$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL5425$input[[1]]))

GPL5425$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL5425$temp_output, input_ = GPL5425$input)

GPL5425$temp_output$Description <- GPL5425$input$`Gene.title` # - OK
GPL5425$temp_output$Gene_ID <- stringr::str_replace_all(string = GPL5425$input$`Gene.ID`, pattern = '///', replacement = ', ') # - OK
GPL5425$temp_output$Unigene <- GPL5425$input$`UniGene.ID` # - OK
GPL5425$temp_output$Nucleotide <- GPL5425$input$GI # - OK

GPL5425$input$bind_for_extraction <- paste(stringr::str_remove(string = GPL5425$input$ID, pattern = '_P.*'), stringr::str_replace_all(string = GPL5425$input$`Gene.symbol`, pattern = '///', replacement = '___'), GPL5425$input$`GenBank.Accession`, GPL5425$input$`UniGene.symbol`, sep = '___')

GPL5425$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = GPL5425$input, bad_cols_to_process_seq = seq(18,18), GPL_ID = 'R')

GPL5425 <- rlist::list.merge(GPL5425, set_search_bys(input_df_ = GPL5425$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

GPL5425$output <- GPL5425$temp_output

GPL5425$output <- extract_identifers(sep_all_identifiers_ = GPL5425$sep_all_identifiers, search_by_ = GPL5425$search_by, output_ = GPL5425$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)

GPL5425$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL5425$output)

write.table(x = GPL5425$output, file = paste0(gen_opts_pre$dir, '/GPL5425_prepared.tsv'), sep = '\t', row.names = F, dec = ',')


# colnames(GPL5425$input) # - OK
# get_all_symbols_in_chrvec(GPL5425$input$ID) # remove _P.* and + # - OK
# get_all_symbols_in_chrvec(GPL5425$input$`Gene.symbol`) # ### !!! additional divider - '///' + # - OK
# get_all_symbols_in_chrvec(GPL5425$input$`Gene.title`) # GPL5425$temp_output$Description # - OK
# get_all_symbols_in_chrvec(GPL5425$input$`Gene.ID`) # ### !!! additional divider - '///' GPL5425$temp_output$Gene_ID
# get_all_symbols_in_chrvec(GPL5425$input$`UniGene.title`) # IGNORE
# get_all_symbols_in_chrvec(GPL5425$input$`UniGene.symbol`) # +
# get_all_symbols_in_chrvec(GPL5425$input$`UniGene.ID`) # GPL5425$temp_output$Unigene # - OK
# get_all_symbols_in_chrvec(GPL5425$input$`Nucleotide.Title`) # IGNORE
# get_all_symbols_in_chrvec(GPL5425$input$GI) # GPL5425$temp_output$Nucleotide # - OK
# get_all_symbols_in_chrvec(GPL5425$input$`GenBank.Accession`) # + # - OK
####################
### GPL5425 - 71 ###
####################








#####################
### GPL25480 - 76 ###
#####################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^76.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL25480 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

GPL25480$input <- rlist::list.rbind(GPL25480$`pre-input`)

GPL25480$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL25480$input[[1]]))

GPL25480$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL25480$temp_output, input_ = GPL25480$input)

GPL25480$temp_output$Probe <- GPL25480$input$SPOT_ID # - OK

GPL25480$output <- GPL25480$temp_output

GPL25480$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL25480$output)

write.table(x = GPL25480$output, file = paste0(gen_opts_pre$dir, '/GPL25480_prepared.tsv'), sep = '\t', row.names = F, dec = ',')

# colnames(GPL25480$input) # - OK
# get_all_symbols_in_chrvec(GPL25480$input$ID) # IGNORE
# get_all_symbols_in_chrvec(GPL25480$input$SPOT_ID) # - OK
#####################
### GPL25480 - 76 ###
#####################





#########################
### GPL10427 - 60, 45 ###
#########################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^45.*tsv|^60.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL10427 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

GPL10427$input <- rlist::list.rbind(GPL10427$`pre-input`)

GPL10427$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL10427$input[[1]]))

GPL10427$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL10427$temp_output, input_ = GPL10427$input)

# colnames(GPL10427$input) # - OK
# get_all_symbols_in_chrvec(GPL10427$input$ID) # OK
# get_all_symbols_in_chrvec(GPL10427$input$accessions) # for extraction, | as devider, remove also 'ug|ref'
# get_all_symbols_in_chrvec(GPL10427$input$ProbeUID) # IGNORE
# get_all_symbols_in_chrvec(GPL10427$input$ProbeName) # GPL10427$temp_output$Probe
# get_all_symbols_in_chrvec(GPL10427$input$GeneName) # for extraction
# get_all_symbols_in_chrvec(GPL10427$input$SystematicName) # for extraction
# get_all_symbols_in_chrvec(GPL10427$input$Description) # GPL10427$temp_output$Description
# get_all_symbols_in_chrvec(GPL10427$input$GB_LIST) # for extraction,  to 

GPL10427$temp_output$Probe <- GPL10427$input$ProbeName # 
GPL10427$temp_output$Description <- GPL10427$input$Description # 

GPL10427$input$bind_for_extraction <- paste(
  prepare_uncorrupted_GPL10427_accessions_wrapper(GPL10427_accessions = GPL10427$input$accessions), 
  GPL10427$input$GeneName, 
  GPL10427$input$SystematicName, 
  stringr::str_replace_all(string = GPL10427$input$GB_LIST, pattern = ', ', replacement = '___'),
  sep = '___') # - OK

GPL10427$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = GPL10427$input, bad_cols_to_process_seq = seq(16,16), GPL_ID = 'R')

GPL10427 <- rlist::list.merge(GPL10427, set_search_bys(input_df_ = GPL10427$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

GPL10427$output <- GPL10427$temp_output

GPL10427$output <- extract_identifers(sep_all_identifiers_ = GPL10427$sep_all_identifiers, search_by_ = GPL10427$search_by, output_ = GPL10427$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)

GPL10427$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL10427$output)

write.table(x = GPL10427$output, file = paste0(gen_opts_pre$dir, '/GPL10427_prepared.tsv'), sep = '\t', row.names = F, dec = ',')
#########################
### GPL10427 - 60, 45 ###
#########################




####################
### GPL8160 - 34 ###
####################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^34.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL8160 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

GPL8160$input <- rlist::list.rbind(GPL8160$`pre-input`)

GPL8160$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL8160$input[[1]]))

GPL8160$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL8160$temp_output, input_ = GPL8160$input)

# colnames(GPL8160$input) # - OK
# get_all_symbols_in_chrvec(GPL8160$input$ID) # GPL8160$temp_output$Probe
# get_all_symbols_in_chrvec(GPL8160$input$GB_ACC) # for extraction

GPL8160$temp_output$Probe <- GPL8160$input$ID # 

GPL8160$input$bind_for_extraction <- GPL8160$input$GB_ACC # 

GPL8160$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = GPL8160$input, bad_cols_to_process_seq = seq(10,10), GPL_ID = 'R')

GPL8160 <- rlist::list.merge(GPL8160, set_search_bys(input_df_ = GPL8160$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

GPL8160$output <- GPL8160$temp_output

GPL8160$output <- extract_identifers(sep_all_identifiers_ = GPL8160$sep_all_identifiers, search_by_ = GPL8160$search_by, output_ = GPL8160$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)

GPL8160$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL8160$output)

write.table(x = GPL8160$output, file = paste0(gen_opts_pre$dir, '/GPL8160_prepared.tsv'), sep = '\t', row.names = F, dec = ',')
####################
### GPL8160 - 34 ###
####################




#########################
### GPL16570 - 56, 50 ###
#########################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^56.*tsv|^50.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL16570 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

GPL16570$input <- rlist::list.rbind(GPL16570$`pre-input`)

GPL16570$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL16570$input[[1]]))

GPL16570$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL16570$temp_output, input_ = GPL16570$input)

# unique(GPL16570$input$Pub.)
# colnames(GPL16570$input) # - OK
# get_all_symbols_in_chrvec(GPL16570$input$ID) # IGNORE
# get_all_symbols_in_chrvec(GPL16570$input$probeset_id) # GPL16570$temp_output$Probe
# get_all_symbols_in_chrvec(GPL16570$input$gene_assignment) # 
# get_all_symbols_in_chrvec(GPL16570$input$mrna_assignment) # IGNORE(Too many wierd symbols, uncortrolable)
# get_all_symbols_in_chrvec(GPL16570$input$swissprot) # IGNORE(Too many wierd symbols, uncortrolable)
# get_all_symbols_in_chrvec(GPL16570$input$unigene) # IGNORE(Too many wierd symbols, uncortrolable)

GPL16570$temp_output$Probe <- GPL16570$input$probeset_id # 

GPL16570$input$bind_for_extraction <- paste(
  prepare_noncorrupted_GPL16570_cols(GPL16570_cols = GPL16570$input$gene_assignment),
  sep = '___')

GPL16570$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = GPL16570$input, bad_cols_to_process_seq = seq(14,14), GPL_ID = 'R')

GPL16570 <- rlist::list.merge(GPL16570, set_search_bys(input_df_ = GPL16570$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

GPL16570$output <- GPL16570$temp_output

GPL16570$output <- extract_identifers(sep_all_identifiers_ = GPL16570$sep_all_identifiers, search_by_ = GPL16570$search_by, output_ = GPL16570$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)

GPL16570$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL16570$output)



write.table(x = GPL16570$output, file = paste0(gen_opts_pre$dir, '/GPL16570_prepared.tsv'), sep = '\t', row.names = F, dec = ',')
#########################
### GPL16570 - 56, 50 ###
#########################




#####################
### GPL13912 - 52 ###
#####################
an_opts_pre <- list('file_names'  = list.files(path = gen_opts_pre$dir, pattern = '^52.*tsv'))
an_opts_pre$pub_names <- stringr::str_extract(string = an_opts_pre$file_names, pattern = '^[[0-9]]{1,3}')
an_opts_pre$exp_names <- stringr::str_remove(
  string = stringr::str_extract(string = an_opts_pre$file_names, 
                                pattern = '^[[0-9]]{1,3}_[[0-9]]{1,3}'),
  pattern = '.*_')

GPL13912 <- list('pre-input' = read_files_and_add_pubExp(an_opts_pre_ = an_opts_pre, gen_opts_pre_ = gen_opts_pre))

GPL13912$input <- rlist::list.rbind(GPL13912$`pre-input`)

GPL13912$temp_output <- set_df_for_append(check_cols_ = gen_opts_pre$check_cols, length = length(GPL13912$input[[1]]))

GPL13912$temp_output <- copy_out_of_box_ready_columns_into_final_output_wrapper(check_cols_ = gen_opts_pre$check_cols, final_output_ = GPL13912$temp_output, input_ = GPL13912$input)


# colnames(GPL13912$input) # - OK
# get_all_symbols_in_chrvec(GPL13912$input$ID) # IGNORE
# get_all_symbols_in_chrvec(GPL13912$input$NAME) # GPL13912$temp_output$Probe 
# get_all_symbols_in_chrvec(GPL13912$input$GB_ACC) # extraction
# get_all_symbols_in_chrvec(GPL13912$input$GENE_ID) # GPL13912$temp_output$Gene_ID
# get_all_symbols_in_chrvec(GPL13912$input$GENE_SYMBOL) # both GPL13912$temp_output$Symbol and extraction for LOCs
# get_all_symbols_in_chrvec(GPL13912$input$GENE_NAME) # GPL13912$temp_output$Description
# get_all_symbols_in_chrvec(GPL13912$input$UNIGENE_ID) # GPL13912$temp_output$Unigene
# get_all_symbols_in_chrvec(GPL13912$input$ENSEMBL_ID) # GPL13912$temp_output$ENST_
# get_all_symbols_in_chrvec(GPL13912$input$ACCESSION_STRING) # extraction, but needs symbol removal
# get_all_symbols_in_chrvec(GPL13912$input$DESCRIPTION) # IGNORE

GPL13912$temp_output$Description <- GPL13912$input$GENE_NAME # 
GPL13912$temp_output$Probe <- GPL13912$input$NAME # 
GPL13912$temp_output$Gene_ID <- GPL13912$input$GENE_ID
GPL13912$temp_output$Symbol <- GPL13912$input$GENE_SYMBOL
GPL13912$temp_output$Unigene <- GPL13912$input$UNIGENE_ID
GPL13912$temp_output$ENST_ <- GPL13912$input$ENSEMBL_ID

GPL13912$input$bind_for_extraction <- paste(
  GPL13912_cols = GPL13912$input$GB_ACC,
  GPL13912$input$GENE_SYMBOL,
  prepare_uncorrupted_GPL13912_ACCESSION_STRING_wrapper(GPL13912$input$ACCESSION_STRING),
  sep = '___')

GPL13912$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = GPL13912$input, bad_cols_to_process_seq = seq(18,18), GPL_ID = 'R')

GPL13912 <- rlist::list.merge(GPL13912, set_search_bys(input_df_ = GPL13912$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

GPL13912$output <- GPL13912$temp_output

GPL13912$output <- extract_identifers(sep_all_identifiers_ = GPL13912$sep_all_identifiers, search_by_ = GPL13912$search_by, output_ = GPL13912$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)

GPL13912$output_cols_check <- get_all_symbols_in_df_per_column(df_ = GPL13912$output)

write.table(x = GPL13912$output, file = paste0(gen_opts_pre$dir, '/GPL13912_prepared.tsv'), sep = '\t', row.names = F, dec = ',')
#####################
### GPL13912 - 52 ###
#####################





#######################################################################################
### MANUAL - 3 5 7 8 12 14 15 18 26 58 38 30 19 27 16 41 46 11 20 80 81 83 82 17 42 ###
#######################################################################################
manual <- list('input' = openxlsx::read.xlsx(xlsxFile = 'Antidepressants_metadata_140520.xlsx', sheet = 'data')[1:25])

colnames(manual$input) <- as.character(gen_opts_pre$check_cols)

manual$input[gen_opts_pre$range_of_cols_with_IDs_in_output_start:gen_opts_pre$range_of_cols_with_IDs_in_output_end] <- purrr::map_dfc(.x = manual$input[gen_opts_pre$range_of_cols_with_IDs_in_output_start:gen_opts_pre$range_of_cols_with_IDs_in_output_end], .f = as.character)

manual$input$Symbol[manual$input$Symbol == 'Ggh ‡'] <- 'Ggh'

manual$input_QA_checklist <- get_input_QA_checklist(input_table_ = manual$input, Pub._col_ = 'Pub.')

manual$input$NM_ <- prepare_noncorrupted_NXMR_col(nxmr_col_ = manual$input$NM_)
manual$input$XM_ <- prepare_noncorrupted_NXMR_col(nxmr_col_ = manual$input$XM_)
manual$input$Accession <- stringr::str_remove(string = manual$input$Accession, pattern = ' $')
manual$input$Symbol <- stringr::str_replace_all(string = prepare_noncorrupted_symbol_col(symbol_col_ = manual$input$Symbol), pattern = ', ', replacement = '___')

manual$temp_output <- manual$input

manual$sep_all_identifiers <- get_sep_all_identifiers(input_df_name = manual$input, bad_cols_to_process_seq = seq(10,10), GPL_ID = 'R')

manual <- rlist::list.merge(manual, set_search_bys(input_df_ = manual$input, search_by_p1_ = "stringr::str_detect(string = ", search_by_n_ = 'identifer', check_cols_ = gen_opts_pre$check_cols))

manual$output <- manual$temp_output

manual$output <- extract_identifers(sep_all_identifiers_ = manual$sep_all_identifiers, search_by_ = manual$search_by, output_ = manual$output, opts_pre_check_cols = gen_opts_pre$check_cols, range_of_cols_with_IDs_start = gen_opts_pre$range_of_cols_with_IDs_in_output_start, range_of_cols_with_IDs_end = gen_opts_pre$range_of_cols_with_IDs_in_output_end)

manual$output_cols_check <- get_all_symbols_in_df_per_column(df_ = manual$output)

write.table(x = manual$output, file = paste0(gen_opts_pre$dir, '/manual_prepared.tsv'), sep = '\t', row.names = F, dec = ',')

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
#######################################################################################
### MANUAL - 3 5 7 8 12 14 15 18 26 58 38 30 19 27 16 41 46 11 20 80 81 83 82 17 42 ###
#######################################################################################


pre_prep_r_down_all_lists <- list('gen_opts_pre' = gen_opts_pre, "GPL1261" = GPL1261, "GPL17223" = GPL17223, "GPL25480" = GPL25480, "GPL5425" = GPL5425,  "GPL6246" = GPL6246,  "GPL6885" = GPL6885,  "GPL6887" = GPL6887, 'GPL10427' = GPL10427, 'GPL8160' = GPL8160, 'GPL16570' = GPL16570, 'GPL13912' = GPL13912, 'manual' = manual)
save(pre_prep_r_down_all_lists, file = 'pre_prep_r_down_all_lists')
rm(pre_prep_r_down_all_lists)

# test_name <- 'NR_'
# test_2 <- tibble::tibble(GPL16570$output[[test_name]])
# get_all_symbols_in_chrvec(test_2[[1]]) #
# test_subset <- stringr::str_detect(string = test_2[[1]], pattern = "\\(")
# test <- subset(x = test_2, subset = test_subset)