source('opts.R')

set_logfile(paste0(opts$ann$folder, "/loggit.log"))

## !! I need to fork and secure devtools::install_github('PavlidisLab/gemmaAPI.R')
# rm(list = ls(pattern = 'temp.*|test.*'))
# experiment 11 was queried against wrong species. It has names for Mus musculus, while the rat was studied in the experiment



##############################
### SET ADDITIONAL OPTIONS ###
##############################
### Read in the table with description of experiments. They need to have the same identifiers as data in PRE_DATA, that is: Paper(int), Experiment(chr). Also needs Species(chr) column, if probe_id identification function or ncbi query tool is to be used
data_annotation <- list(
  'descriptions' = NA,
  'pre_input' = NA,
  'input' = NA,
  'probes_direct' = NA,
  'gemma' = NA,
  'identifers' = NA,
  'symbols' = NA,
  'post_annot' = NA
)

# load(paste0(opts$ann$folder, '/data_annotation'))

load(paste0(opts$pre_prep_reform$folder, '/pre_prep_reform'))
load(paste0(opts$dir_descriptions, '/descriptions_for_analysis'))

data_annotation$pre_input$data <- pre_prep_reform$combined[['data']]
data_annotation$descriptions$desc <- descriptions[['output']]

rm(pre_prep_reform, descriptions)

#data_annotation$pre_input[['sample_data']] <- dplyr::sample_n(tbl = data_annotation$pre_input[['data']], size = 2500)




### Check if every exp has unique platform.
data_annotation$descriptions[['unique_platforms']] <- subset(
  data_annotation$descriptions[['desc']], 
  select = c(opts$ann[['exp_id_col']], opts$ann[['des_platform_col']])) %>%
  dplyr::group_by( deparse(substitute(opts$ann[['exp_id_col']])) ) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    'nb_of_platforms_per_exps' = purrr::map(.x = data, .f = function(x){length(unique(x))}))
  


assert_desc_and_input_have_the_same_exp_col_wrapper(
  data_ = data_annotation$pre_input[['data']], 
  desc_ = data_annotation$descriptions[['desc']], 
  exp_col = opts$ann[['exp_id_col']])
##############################
### SET ADDITIONAL OPTIONS ###
##############################







##############################################
### IDENTIFY PLATFORM FOR PROBE ANNOTATION ###
##############################################
data_annotation$pre_input[['input_exp_with_probes']] <- unique(
  subset(
    x = data_annotation$pre_input[['data']], 
    subset = !is.na(data_annotation$pre_input[['data']][[opts$check_cols_for_input[['Probe']]]])
    )[[opts$check_cols_for_input[['Exp.']]]]
  )


data_annotation$pre_input[['exp_for_probe_analysis']][['input']] <- subset(
  x = data_annotation$pre_input[['data']], 
  subset = data_annotation$pre_input[['data']][[opts$check_cols_for_input[['Exp.']]]] %in% data_annotation$pre_input[['input_exp_with_probes']]
  )


data_annotation$pre_input[['exp_NOT_for_probe_analysis']][['input']] <- subset(
  x = data_annotation$pre_input[['data']], 
  subset = data_annotation$pre_input[['data']][[opts$check_cols_for_input[['Exp.']]]] %nin% data_annotation$pre_input[['input_exp_with_probes']]
  )


length(data_annotation$pre_input[['data']][[1]]) == 
  length(data_annotation$pre_input[['exp_for_probe_analysis']][['input']][[1]]) + 
  length(data_annotation$pre_input[['exp_NOT_for_probe_analysis']][['input']][[1]])


data_annotation$pre_input$exp_for_probe_analysis[['highest_hit_analysis']] <- master_annotator(
  descriptions_df = data_annotation[['descriptions']][['desc']], 
  exp_id_col = opts$ann[['exp_id_col']], 
  des_species_col = opts$ann[['des_species_col']], 
  input_df = data_annotation$pre_input[['exp_for_probe_analysis']][['input']], 
  input_id_col = 'Probe',
  PERFORM_B_get_the_highest_hit_returning_id_type = T,
  B_highest_hit_int_Probe_IDs_to_test = 200)




### REMOVE EXPERIMENTS WITH IMPROPERLY DETECTED PLATFORMS ###
data_annotation$pre_input$exp_for_probe_analysis$improperly_detected_exp <- 22


data_annotation$pre_input$exp_for_probe_analysis$platform_ids_for_probe_annotation <- subset(
  x = data_annotation$pre_input$exp_for_probe_analysis[['highest_hit_analysis']][['best_ID']],
  subset = data_annotation$pre_input$exp_for_probe_analysis[['highest_hit_analysis']][['best_ID']][['Exp_ID']] %nin%
    data_annotation$pre_input$exp_for_probe_analysis[['improperly_detected_exp']])

data_annotation$pre_input$exp_for_probe_analysis$platform_ids_for_probe_annotation$platform_to_use <- as.character(data_annotation$pre_input$exp_for_probe_analysis$platform_ids_for_probe_annotation[['platform_to_use']])


data_annotation$pre_input$exp_for_probe_analysis$input_proper <- subset(
  x = data_annotation$pre_input$exp_for_probe_analysis[['input']], 
  subset = data_annotation$pre_input$exp_for_probe_analysis[['input']][[ opts$check_cols_for_input[['Exp.']] ]] %in% 
    data_annotation$pre_input$exp_for_probe_analysis[['platform_ids_for_probe_annotation']][['Exp_ID']])


data_annotation$pre_input$exp_NOT_for_probe_analysis$input_proper <- rbind(
  data_annotation$pre_input$exp_NOT_for_probe_analysis$input,
  subset(
    x = data_annotation$pre_input$exp_for_probe_analysis[['input']], 
    subset = data_annotation$pre_input$exp_for_probe_analysis[['input']][[ opts$check_cols_for_input[['Exp.']] ]] %nin% 
      data_annotation$pre_input$exp_for_probe_analysis[['platform_ids_for_probe_annotation']][['Exp_ID']])
)

data_annotation$pre_input$was_properization_ok <- length(data_annotation$pre_input$exp_NOT_for_probe_analysis$input[[1]]) + length(data_annotation$pre_input$exp_for_probe_analysis$input[[1]]) ==
  length(data_annotation$pre_input$exp_NOT_for_probe_analysis$input_proper[[1]]) + length(data_annotation$pre_input$exp_for_probe_analysis$input_proper[[1]])





data_annotation$input <- list(
  'data_for_probe_annotation' = data_annotation$pre_input[['exp_for_probe_analysis']][['input_proper']],
  'platform_ids_for_probe_annotation' = data_annotation$pre_input$exp_for_probe_analysis[['platform_ids_for_probe_annotation']][['platform_to_use']],
  'data_NOT_for_probe_annotation' = data_annotation$pre_input[['exp_NOT_for_probe_analysis']][['input_proper']]
)

data_annotation$input$qa$data_for_probe_annotation <- verify_df(df_ = data_annotation$input[['data_for_probe_annotation']], only_qa = T)
data_annotation$input$qa$data_NOT_for_probe_annotation <- verify_df(df_ = data_annotation$input[['data_NOT_for_probe_annotation']], only_qa = T)


message('Saving platform identification as pre_input_analysis...')
pre_input <- data_annotation$pre_input
save(pre_input, file = paste0(opts$ann[['folder']], '/pre_input_analysis'))
message('... saved. Please remember to check quality of data - use "qa" element of saved list.')

##############################################
### IDENTIFY PLATFORM FOR PROBE ANNOTATION ###
##############################################









#######################
### Annotate probes ###
#######################
message('*** STARTING PROBE-CENTERED ANNOTATION... ***')

message('Starting annotation...')
data_annotation$probes_direct$annotations <- master_annotator(
  descriptions_df = data_annotation[['descriptions']][['desc']],
  exp_id_col = opts$ann[['exp_id_col']],
  des_species_col = opts$ann[['des_species_col']],
  input_df = data_annotation$input[['data_for_probe_annotation']],
  input_id_col = 'Probe',
  PERFORM_C_annotation = T,
  C_legacy_D_str_identifier_type__ = data_annotation$input[['platform_ids_for_probe_annotation']],
  A_D_C_string_separator__ = opts$str_sep)
message('... annotated.')



message('Subsetting annotation into finalized and leftover data...')
data_annotation$probes_direct$finalized <- subset(x = data_annotation$probes_direct$annotations, subset = !is.na(data_annotation$probes_direct[['annotations']][['external_gene_name']]))

data_annotation$probes_direct$leftover <- subset(x = data_annotation$probes_direct$annotations, subset = is.na(data_annotation$probes_direct[['annotations']][['external_gene_name']]))



message('Performing QA...')
length(data_annotation$input[['data_for_probe_annotation']][[1]]) == length(data_annotation$probes_direct[['annotations']][[1]]) 

was_splitting_ok <- check_was_the_spliting_of_df_by_filtering_ok(
  df_original = data_annotation$probes_direct[['annotations']], 
  list_df_splited = list(
    data_annotation$probes_direct[['finalized']], 
    data_annotation$probes_direct[['leftover']]))

if (was_splitting_ok) {
  message(sprintf(
    opts$ann[['splitting_message_ok']], 
    deparse(substitute(data_annotation$probes_direct[['annotations']])),
    deparse(substitute(data_annotation$probes_direct[['finalized']])),
    deparse(substitute(data_annotation$probes_direct[['leftover']]))))
} else stop(opts$ann[['splitting_message_bad']])

message('... subsetting done.')

data_annotation$probes_direct$qa$finalized <- verify_df(df_ = data_annotation$probes_direct[['finalized']], only_qa = T)
data_annotation$probes_direct$qa$leftover <- verify_df(df_ = data_annotation$probes_direct[['leftover']], only_qa = T)
message('... performed.')



message('Saving probe-centered analysis as probe_centered_analysis...')
probes_direct <- data_annotation$probes_direct
save(probes_direct, file = paste0(opts$ann[['folder']], '/probe_centered_analysis'))
message('... saved. Please remember to check quality of data - use "qa" element of saved list.')

rm(probes_direct, was_splitting_ok)

message('*** PROBE-CENTERED ANNOTATION FINISHED! ***')
#######################
### Annotate probes ###
#######################





##################################
### Annotate probes with GEMMA ###
##################################
message('*** STARTING GEMMA PRE-ANNOTATION OF EXPERIMENTS/PLATFORMS NOT FOUND IN ENSEMBL... ***')

data_annotation$gemma$input_from_probes_leftover <- dplyr::select(data_annotation$probes_direct[['leftover']], -'external_gene_name')


message('Setting experiments to check against gemma platform base (these experiments that had 0 hits in probe-centered analysis)... ')
data_annotation$gemma$candidate_exp_for_gemma <- 
  unique(
    data_annotation$gemma[['input_from_probes_leftover']][[opts$ann[['exp_id_col']]]]
    )[unique(data_annotation$gemma[['input_from_probes_leftover']][[opts$ann[['exp_id_col']]]]) %nin% unique(data_annotation$probes_direct[['finalized']][[opts$ann[['exp_id_col']]]])]

assertthat::assert_that(opts$ann[['exp_id_col']] == opts$ann[['exp_id_col']], msg = 'exp_id_col and exp_id_col need to be the same in current version of tool. see below and above commands for reason.')

data_annotation$gemma$candidate_exp_for_gemma <- 
  unique(subset(
    x = data_annotation$descriptions[['desc']], 
    subset = data_annotation$descriptions[['desc']][[opts$ann[['exp_id_col']]]] %in% data_annotation$gemma[['candidate_exp_for_gemma']], 
    select = c(opts$ann[['exp_id_col']], opts$ann[['des_platform_col']])))

message(sprintf('Set platforms %s for experiments: %s.', paste(data_annotation$gemma[['candidate_exp_for_gemma']][[opts$ann[['des_platform_col']]]], collapse = ', '), paste(data_annotation$gemma[['candidate_exp_for_gemma']][[opts$ann[['exp_id_col']]]], collapse = ', ')))



message('Downloading gemma annotation files...')
data_annotation$gemma$gemma_platforms <- download_platforms_from_gemma(unique(data_annotation$gemma[['candidate_exp_for_gemma']][[opts$ann[['des_platform_col']]]]))
message('... downloaded.')



data_annotation$gemma$exp_found_in_gemma <- subset(
  x = data_annotation$gemma[['candidate_exp_for_gemma']], 
  subset = toupper(data_annotation$gemma[['candidate_exp_for_gemma']][[opts$ann[['des_platform_col']]]]) %in% names(data_annotation$gemma[['gemma_platforms']][!is.na(data_annotation$gemma[['gemma_platforms']])]))

data_annotation$gemma$exp_NOT_found_in_gemma <-  subset(
  x = data_annotation$gemma[['candidate_exp_for_gemma']], 
  subset = toupper(data_annotation$gemma[['candidate_exp_for_gemma']][[opts$ann[['des_platform_col']]]]) %in% names(data_annotation$gemma[['gemma_platforms']][is.na(data_annotation$gemma[['gemma_platforms']])]))


length(data_annotation$gemma[['candidate_exp_for_gemma']][[1]]) == length(data_annotation$gemma$exp_found_in_gemma[[1]]) + length(data_annotation$gemma[['exp_NOT_found_in_gemma']][[1]])



message(sprintf('Platforms to be pre-annotated in gemma are: %s for expications: %s.', paste(data_annotation$gemma[['exp_found_in_gemma']][[opts$ann[['des_platform_col']]]], collapse = ', '), paste(data_annotation$gemma[['exp_found_in_gemma']][[opts$ann[['exp_id_col']]]], collapse = ', ')))



data_annotation$gemma$leftover$leftovers_not_run_through_gemma_cause_these_exp_were_analyzed_in_probes_step <- subset(
  x = data_annotation$gemma[['input_from_probes_leftover']], 
  subset = data_annotation$gemma[['input_from_probes_leftover']][[opts$ann[['exp_id_col']]]] %nin% data_annotation$gemma[['candidate_exp_for_gemma']][[opts$ann[['exp_id_col']]]])

data_annotation$gemma$leftover$leftovers_not_run_through_gemma_cause_platform_was_not_found <- subset(
  x = data_annotation$gemma[['input_from_probes_leftover']], 
  subset = data_annotation$gemma[['input_from_probes_leftover']][[opts$ann[['exp_id_col']]]] %in% data_annotation$gemma[['exp_NOT_found_in_gemma']][[opts$ann[['exp_id_col']]]])



message('Pre-annotating data with various identfiers found in gemma platform description files (to be annotated in next, identifer-centered step)... ')
data_annotation$gemma$annotation_raw <- get_gemma_annotations_for_data(
  descriptions_df = data_annotation[['descriptions']][['desc']], 
  des_platform_col = opts$ann[['des_platform_col']], 
  exp_id_col = opts$ann[['exp_id_col']], 
  named_list_gemma_platforms = data_annotation$gemma[['gemma_platforms']], 
  gemma_probe_col = 'ProbeName', 
  input_df = data_annotation$gemma[['input_from_probes_leftover']],
  input_probe_col = 'Probe')




data_annotation$gemma$annotation$GPL8160 <- GPL8160_post_gemma_wrapper(
  annotation_raw = data_annotation$gemma[['annotation_raw']][['GPL8160']], 
  divider_of_values_in_serie_str_ = opts$str_sep, 
  old_divider_to_be_replaced = '\\|', 
  uniqualize_ = T, 
  Gene_ID_col = 'Gene_ID', 
  gene_description_col = 'Description', 
  gene_symbol_col = 'Symbol', 
  cols_to_nullify = c('GeneSymbols', 'NCBIids', 'GeneNames', 'GOTerms', 'GemmaIDs'))

message('... following platforms were annotated: ', paste(names(data_annotation$gemma[['annotation']]), collapse = ', '))



data_annotation$gemma$finalized_full_gemma_input_pre_annotated <- rbind(
  data_annotation$gemma$leftover[['leftovers_not_run_through_gemma_cause_these_exp_were_analyzed_in_probes_step']], 
  data_annotation$gemma$leftover[['leftovers_not_run_through_gemma_cause_platform_was_not_found']], 
  data_annotation$gemma$annotation[['GPL8160']])



message('Performing QA...')
length(data_annotation$gemma[['finalized_full_gemma_input_pre_annotated']][[1]]) == length(data_annotation$gemma[['input_from_probes_leftover']][[1]]) 

i_and_o_has_the_same_nb_of_rows <- length(data_annotation$gemma[['finalized_full_gemma_input_pre_annotated']][[1]]) == length(data_annotation$gemma[['input_from_probes_leftover']][[1]])

if (i_and_o_has_the_same_nb_of_rows) {
  message(sprintf(
    'Ok. Input (%s) and output (%s) of this step have the same length.', 
    deparse(substitute(data_annotation$gemma[['input_from_probes_leftover']][[1]])), 
    deparse(substitute(data_annotation$gemma[['finalized_full_gemma_input_pre_annotated']][[1]]))))
} else stop('Error! Input and output of this step have different lengths.')

data_annotation$gemma$qa$finalized_full_gemma_input_pre_annotated <- verify_df(df_ = data_annotation$gemma[['finalized_full_gemma_input_pre_annotated']], only_qa = T)
message('... performed.')



message('Saving gemma-centered pre-analysis as gemma_centered_pre_analysis...')
gemma <- data_annotation$gemma
save(gemma, file = paste0(opts$ann$folder, '/gemma_centered_pre_analysis'))
message('... saved. Please remember to check quality of data - use "qa" element of saved list.')

rm(i_and_o_has_the_same_nb_of_rows, gemma)

message('*** ... GEMMA PRE-ANNOTATION FINISHED ***')
##################################
### Annotate probes with GEMMA ###
##################################



message('*** MERGING DATA FOR NON-PROBE INDENTIFER-CENTERED ANNOTATION... ***')
message('Merging data originally without probe identifers and unannotated entries from probe-centered analysis (both passed through gemma pre-annotation and not)... ')

if (all(colnames(data_annotation$gemma[['finalized_full_gemma_input_pre_annotated']]) == colnames(data_annotation$input[['data_NOT_for_probe_annotation']]))) {
  message('Datasets have the same column names')
} else stop('Error! Datasets have different column names')

if (length(colnames(data_annotation$gemma[['finalized_full_gemma_input_pre_annotated']])) == length(colnames(data_annotation$input[['data_NOT_for_probe_annotation']]))) {
  message('Datasets have the same number of columns')
} else stop('Error! Datasets have different number of columns')

data_annotation$identifers$input <- rbind(data_annotation$gemma[['finalized_full_gemma_input_pre_annotated']], data_annotation$input[['data_NOT_for_probe_annotation']])

message('*** ... MERGED! ***')





#####################################
### Pseudo-memoization non-probes ###
#####################################
message('*** STARTING NON-PROBE INDENTIFER-CENTERED ANNOTATION... ***')

message('Creating database for annotation...')
data_annotation$identifers$pseudo_memoized_annotation_db <- master_annotator(
  descriptions_df = data_annotation[['descriptions']][['desc']],
  exp_id_col = opts$ann[['exp_id_col']],
  des_species_col = opts$ann[['des_species_col']],
  input_df = data_annotation$identifers[['input']],
  input_id_col = 'Probe',
  A_C_all_ENSG_col = 'ENSG_', A_C_all_ENST_col = 'ENST_', A_C_all_Gene_ID_col = 'Gene_ID', A_C_all_NM_col = 'NM_', A_C_all_Accession_col = 'Accession', A_C_all_NR_col = 'NR_', A_C_all_XM_col = 'XM_', A_C_all_XR_col = 'XR_',
  PERFORM_A_should_i_prepare_dbs_for_pseudomemoization = T,
  A_return_qa_of_pseudomemoization = F,
  A_D_C_string_separator__ = opts$str_sep)
message('... created.')
# actually, mgi_id are not Mm., that is why we are not using Unigene column for annotation




message('Starting annotation...')
data_annotation$identifers$annotations$list <- master_annotator(
  descriptions_df = data_annotation$descriptions[['desc']],
  exp_id_col = opts$ann[['exp_id_col']],
  des_species_col = opts$ann[['des_species_col']],
  input_df = data_annotation$identifers[['input']],
  input_id_col = 'Probe',
  PERFORM_C_annotation = T,
  A_C_all_ENSG_col = 'ENSG_', A_C_all_ENST_col = 'ENST_', A_C_all_Gene_ID_col = 'Gene_ID', A_C_all_NM_col = 'NM_', A_C_all_Accession_col = 'Accession', A_C_all_NR_col = 'NR_', A_C_all_XM_col = 'XM_', A_C_all_XR_col = 'XR_',
  A_D_C_string_separator__ = opts$str_sep,
  C_pseudo_memoized_db = data_annotation$identifers[['pseudo_memoized_annotation_db']])
message('... annotated.')



data_annotation$identifers[['annotations']][['df']] <- rlist::list.rbind(data_annotation$identifers[['annotations']][['list']])

message('Subsetting annotation into finalized and leftover data...')

data_annotation$identifers$finalized <- subset(x = data_annotation$identifers[['annotations']][['df']], subset = !is.na(data_annotation$identifers[['annotations']][['df']][['external_gene_name']]))

data_annotation$identifers$leftover <- subset(x = data_annotation$identifers[['annotations']][['df']], subset = is.na(data_annotation$identifers[['annotations']][['df']][['external_gene_name']]))



message('Performing QA...')
length(data_annotation$identifers[['annotations']][['df']][[1]]) == length(data_annotation$identifers$input[[1]]) 

was_splitting_ok <- check_was_the_spliting_of_df_by_filtering_ok(
  df_original = data_annotation$identifers$annotations[['df']], 
  list_df_splited = list(data_annotation$identifers[['finalized']], data_annotation$identifers[['leftover']]))

if (was_splitting_ok) {
  message(sprintf(
    opts$ann[['splitting_message_ok']], 
    deparse(substitute(data_annotation$identifers[['annotations']][['df']])),
    deparse(substitute(data_annotation$identifers[['finalized']])),
    deparse(substitute(data_annotation$identifers[['leftover']]))))
} else stop(opts$ann[['splitting_message_bad']])

message('... subsetting done.')

data_annotation$identifers$qa$finalized <- verify_df(df_ = data_annotation$identifers[['finalized']], only_qa = T)
data_annotation$identifers$qa$leftover <- verify_df(df_ = data_annotation$identifers[['leftover']], only_qa = T)
message('... performed.')



message('Saving non-probe identifier-centered analysis as identifier_centered_analysis...')
identifers <- data_annotation$identifers
save(identifers, file = paste0(opts$ann[['folder']], '/identifier_centered_analysis'))
message('... saved. Please remember to check quality of data - use "qa" element of saved list.')

rm(identifers, was_splitting_ok)

message('*** ... NON-PROBE INDENTIFER-CENTERED ANNOTATION FINISHED ***')
#####################################
### Pseudo-memoization non-probes ###
#####################################



#####################################
### Annotate gene names with ncbi ###
#####################################
message('*** STARTING GENE SYMBOL-CENTERED ANNOTATION... ***')

data_annotation$symbols$input_from_non_probes_leftovers <- data_annotation$identifers[['leftover']]



message('Translating all gene names in the data to gene ids using ncbi... ')
data_annotation$symbols$ncbi_annotation_of_symbols_to_gene_ids <- master_annotator(
  descriptions_df = data_annotation[['descriptions']][['desc']],
  exp_id_col = opts$ann[['exp_id_col']],
  des_species_col = opts$ann[['des_species_col']],
  input_df = data_annotation$symbols[['input_from_non_probes_leftovers']],
  input_id_col = opts$check_cols_for_input[['Symbol']],
  C_legacy_D_str_identifier_type__ = 'Gene name',
  PERFORM_D_ncbi_annotation = T, 
  A_D_C_string_separator__ = opts$str_sep)
message('... translated.')



message('Adding ncbi-derived gene ids to original data... ')
data_annotation$symbols$input_plus_new_gene_id_from_ncbi_column <- master_annotator(
  descriptions_df = data_annotation[['descriptions']][['desc']],
  exp_id_col = opts$ann[['exp_id_col']],
  des_species_col = opts$ann[['des_species_col']],
  input_df = data_annotation$symbols[['input_from_non_probes_leftovers']],
  input_id_col = opts$check_cols_for_input[['Symbol']],
  PERFORM_E_add_new_gene_id_col_originating_from_ncbi = T, 
  E_PERFORM_D_output = data_annotation$symbols[['ncbi_annotation_of_symbols_to_gene_ids']],
  A_D_C_string_separator__ = opts$str_sep)
message('.. added.')



message('Creating database for annotation...')
data_annotation$symbols$pseudomemoize_external_gene_names_to_gene_id <- master_annotator(
  descriptions_df = data_annotation$symbols[['ncbi_annotation_of_symbols_to_gene_ids']],
  exp_id_col = 'dummy_exp',
  des_species_col = 'organism',
  input_df = data_annotation$symbols[['ncbi_annotation_of_symbols_to_gene_ids']],
  input_id_col = opts$check_cols_for_input[['Gene_ID']],
  A_C_all_Gene_ID_col = opts$check_cols_for_input[['Gene_ID']],
  PERFORM_A_should_i_prepare_dbs_for_pseudomemoization = T,
  A_return_qa_of_pseudomemoization = F,
  A_D_C_string_separator__ = opts$str_sep)
message('... created.')



message('Starting annotation...')
data_annotation$symbols$annotations$list <- master_annotator(
  descriptions_df = data_annotation[['descriptions']][['desc']],
  exp_id_col = opts$ann[['exp_id_col']],
  des_species_col = opts$ann[['des_species_col']],
  input_df = data_annotation$symbols[['input_plus_new_gene_id_from_ncbi_column']],
  input_id_col = opts$check_cols_for_input[['Probe']],
  PERFORM_C_annotation = T,
  C_legacy_D_str_identifier_type__ = opts$check_cols_for_input[['Gene_ID']],
  A_C_all_Gene_ID_col = 'Gene_ID_from_ncbi',
  C_pseudo_memoized_db = data_annotation$symbols[['pseudomemoize_external_gene_names_to_gene_id']],
  A_D_C_string_separator__ = opts$str_sep)
message('... annotated.')



data_annotation$symbols$annotations$df <- rlist::list.rbind(data_annotation$symbols[['annotations']][['list']])

message('Subsetting annotation into finalized and leftover data...')
data_annotation$symbols$finalized <- subset(x = data_annotation$symbols[['annotations']][['df']], subset = !is.na(data_annotation$symbols[['annotations']][['df']][['external_gene_name']]))

data_annotation$symbols$leftover <- subset(x = data_annotation$symbols[['annotations']][['df']], subset = is.na(data_annotation$symbols[['annotations']][['df']][['external_gene_name']]))



message('Performing QA...')
length(data_annotation$symbols[['annotations']][['df']][[1]]) == length(data_annotation$symbols[['input_from_non_probes_leftovers']][[1]])

was_splitting_ok <- check_was_the_spliting_of_df_by_filtering_ok(
  df_original = data_annotation$symbols[['annotations']][['df']], 
  list_df_splited = list(
    data_annotation$symbols[['finalized']], 
    data_annotation$symbols[['leftover']]))

if (was_splitting_ok) {
  message(sprintf(
    opts$ann[['splitting_message_ok']], 
    deparse(substitute(data_annotation$symbols[['annotations']][['df']])),
    deparse(substitute(data_annotation$symbols[['finalized']])),
    deparse(substitute(data_annotation$symbols[['leftover']]))))
} else stop(opts$ann[['splitting_message_bad']])

message('... subsetting done.')

data_annotation$symbols$qa$finalized <- verify_df(df_ = data_annotation$symbols[['finalized']], only_qa = T)
data_annotation$symbols$qa$leftover <- verify_df(df_ = data_annotation$symbols[['leftover']], only_qa = T)
message('... performed.')



message('Saving ncbi-centered analysis as ncbi_centered_analysis...')
symbols <- data_annotation[['symbols']]
save(symbols, file = paste0(opts$ann[['folder']], '/ncbi_centered_analysis'))
message('... saved. Please remember to check quality of data - use "qa" element of saved list.')

rm(symbols, was_splitting_ok)

message('*** ... GENE SYMBOL-CENTERED ANNOTATION FINISHED! ***')
#####################################
### Annotate gene names with ncbi ###
#####################################










#######################
### Post-annotation ###
#######################
message('*** MERGING ALL FINALIZED DATA... ***')

message('Normalizing columns...')
data_annotation$probes_direct[['finalized']][['identifed_identifer']] <- data_annotation[['probes_direct']][['finalized']][['Probe']]
data_annotation$probes_direct[['finalized']][['unidentifed_identifers']] <- NA
data_annotation$probes_direct[['finalized']][['Gene_ID_from_ncbi']] <- NA

data_annotation$identifers[['finalized']][['Gene_ID_from_ncbi']] <- NA


data_annotation$full_finalized$data <- rbind(
  data_annotation$probes_direct[['finalized']], 
  data_annotation$identifers[['finalized']], 
  data_annotation$symbols[['finalized']])

data_annotation$true_leftovers$data <- data_annotation[['symbols']][['leftover']]

if ((length(data_annotation$input[['data_for_probe_annotation']][[1]]) + length(data_annotation$input[['data_NOT_for_probe_annotation']][[1]])) == (length(data_annotation[['full_finalized']][['data']][[1]]) + length(data_annotation[['true_leftovers']][['data']][[1]]))) {
  message('Annotation input and out datasets have equal number of entries.')
} else stop('Error! Annotation input and out datasets have different number of entries!')

data_annotation$full_finalized$qa <- verify_df(df_ = data_annotation[['full_finalized']][['data']], only_qa = T)
data_annotation$true_leftovers$qa <- verify_df(df_ = data_annotation[['true_leftovers']][['data']], only_qa = T)



# Take first gene symbol from true_leftovers and put it into external_gene_name
data_annotation$post_annot$usable_leftovers_list <- get_usable_symbols_subset_of_leftovers(
  true_leftovers_ = data_annotation$true_leftovers[['data']], 
  input_symbol_col = opts$check_cols_for_input[['Symbol']],
  input_geneid_col =  opts$check_cols_for_input[['Gene_ID']],
  other_input_cols = c(
    opts$check_cols_for_input[['ENSG_']],
    opts$check_cols_for_input[['ENST_']],
    opts$check_cols_for_input[['NM_']],
    opts$check_cols_for_input[['Accession']],
    opts$check_cols_for_input[['XM_']],
    opts$check_cols_for_input[['XR_']],
    opts$check_cols_for_input[['NR_']],
    opts$check_cols_for_input[['Unigene']]
    ),
  output_ext_gene_name_col = opts[['ext_gene_name_col']]
  )



data_annotation$post_annot$final_dataset <- rbind(
  data_annotation$full_finalized[['data']], 
  data_annotation$post_annot$usable_leftovers_list[['usable']])



data_annotation$post_annot$final_dataset <- data_annotation$post_annot$final_dataset[order(data_annotation$post_annot$final_dataset[[ opts$check_cols_for_input[['Exp.']] ]]),]



data_annotation$post_annot$final_dataset_qa <- verify_df(df_ = data_annotation$post_annot[['final_dataset']], only_qa = T)



data_annotation$post_annot$are_are_input_and_final_dataset_the_same_lenght <- length(data_annotation$input$data_for_probe_annotation[[1]]) +
  length(data_annotation$input$data_NOT_for_probe_annotation[[1]]) ==
  length(data_annotation$post_annot$final_dataset[[1]]) +
  length(data_annotation$post_annot$usable_leftovers_list$unusable[[1]])






message('Saving whole annotation data as data_annotation')
save(data_annotation, file = paste0(opts$ann[['folder']], '/data_annotation'))
message('... saved.')


message('*** ... FINALIZED DATA MERGED! ***')
#######################
### Post-annotation ###
#######################














# test <- dplyr::sample_n(tbl = data_annotation$post_annot$final_dataset, size = 100)
# test <- subset(test, select = c(1,2,3,8,9,10,12:20,22:25))

test_1 <- 12

test <- subset(
  data_annotation$post_annot$final_dataset, 
  subset = data_annotation$post_annot$final_dataset$Exp. == test_1, 
  select = c(1,2,3,8,9,10,12:20,22:25))
test_left <- subset(
  data_annotation$post_annot$usable_leftovers_list$unusable, 
  subset = data_annotation$post_annot$usable_leftovers_list$unusable$Exp. == test_1, 
  select = c(1,2,3,8,9,10,12:20,22:25))


test_input <- rbind(data_annotation$input$data_for_probe_annotation, data_annotation$input$data_NOT_for_probe_annotation)
test_input <- subset(test_input, select = c(1,2,3,8,9,10,12:20))




