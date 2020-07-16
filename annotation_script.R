source('functions_for_annotation.R')

library(loggit)
default::default(message) <- list(echo = F)

### !!! need to add step which checks wheter all papers have unique platforms
### !!! acutally, we need to use experiment id, instead of paper id, which has the same platform... species is implied then. perhaps lets call it dataset?
### !!! I need to fork and secure devtools::install_github('PavlidisLab/gemmaAPI.R')


### Prepare working list ###
# load('data_prepared_for_annotation')

data_annotation <- list(
  'descriptions' = data_prepared_for_annotation$descriptions,
  'input' = NA,
  'probes_direct' = NA,
  'gemma' = NA,
  'identifers' = NA,
  'symbols' = NA
)

data_annotation$input <- list(
  'data_for_probe_annotation' = data_prepared_for_annotation$data_for_probe_annotation,
  'platform_ids_for_probe_annotation' = data_prepared_for_annotation$platform_ids_for_probe_annotation,
  'data_NOT_for_probe_annotation' = data_prepared_for_annotation$data_NOT_for_probe_annotation
)

data_annotation$input$qa$data_for_probe_annotation <- verify_df(df_ = data_annotation$input$data_for_probe_annotation, only_qa = T)
data_annotation$input$qa$platform_ids_for_probe_annotation <- verify_df(df_ = data_annotation$input$platform_ids_for_probe_annotation, only_qa = T)
data_annotation$input$qa$data_NOT_for_probe_annotation <- verify_df(df_ = data_annotation$input$data_NOT_for_probe_annotation, only_qa = T)

### Prepare working list ###

opts_da <- list('des_paper_id_col' = 'Pub.',
                'des_species_col' = 'Species',
                'input_paper_col' = 'Pub.',
                'des_platform_col' = 'Assay',
                'annotation_folder' = 'annotation_qc',
                'splitting_message_ok' = 'Splitting of %s into %s ad %s was ok.',
                'splitting_message_bad' = 'Annotated data do not have the same number of rows as finalized and leftover data. Splitting failed!'
                )

dir.create(opts_da$annotation_folder)

set_logfile(paste0(opts_da$annotation_folder, "/loggit.log"))

message('Saving analysis options...')
save(opts_da, file = 'opts_da_used_for_analysis')
message('... options saved.')

#######################
### Annotate probes ###
#######################
message('*** STARTING PROBE-CENTERED ANNOTATION... ***')

message('Starting annotation...')
data_annotation$probes_direct$annotations <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = opts_da$des_paper_id_col,
  des_species_col = opts_da$des_species_col,
  input_df = data_annotation$input$data_for_probe_annotation,
  input_paper_col = opts_da$input_paper_col,
  input_id_col = 'Probe',
  PERFORM_C_annotation = T,
  C_legacy_D_str_identifier_type__ = data_annotation$input$platform_ids_for_probe_annotation)
message('... annotated.')



message('Subsetting annotation into finalized and leftover data...')
data_annotation$probes_direct$finalized <- subset(x = data_annotation$probes_direct$annotations, subset = !is.na(data_annotation$probes_direct$annotations$external_gene_name))

data_annotation$probes_direct$leftover <- subset(x = data_annotation$probes_direct$annotations, subset = is.na(data_annotation$probes_direct$annotations$external_gene_name))


### !!! What about annotated vs input???
was_splitting_ok <- check_was_the_spliting_of_df_by_filtering_ok(
  df_original = data_annotation$probes_direct$annotations, 
  list_df_splited = list(
    data_annotation$probes_direct$finalized, 
    data_annotation$probes_direct$leftover))

if (was_splitting_ok) {
  message(sprintf(
    opts_da$splitting_message_ok, 
    deparse(substitute(data_annotation$probes_direct$annotations)),
            deparse(substitute(data_annotation$probes_direct$finalized)),
                    deparse(substitute(data_annotation$probes_direct$leftover))))
} else stop(opts_da$splitting_message_bad)

message('... subsetting done.')



data_annotation$probes_direct$qa$finalized <- verify_df(df_ = data_annotation$probes_direct$finalized, only_qa = T)
data_annotation$probes_direct$qa$leftover <- verify_df(df_ = data_annotation$probes_direct$leftover, only_qa = T)



message('Saving probe-centered analysis as probe_centered_analysis...')
probes_direct <- data_annotation$probes_direct
save(probes_direct, file = paste0(opts_da$annotation_folder, '/probe_centered_analysis'))
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

data_annotation$gemma$input_from_probes_leftover <- dplyr::select(data_annotation$probes_direct$leftover, -'external_gene_name')

### !!! Are platforms with NA returned in getting_highest_hit also here?
message('Setting publications to check against gemma platform base (these publications that had 0 hits in probe-centered analysis)... ')
data_annotation$gemma$candidate_publ_for_gemma <- 
  unique(data_annotation$probes_direct$leftover[[opts_da$input_paper_col]])[unique(data_annotation$probes_direct$leftover[[opts_da$input_paper_col]]) %nin% unique(data_annotation$probes_direct$finalized[[opts_da$input_paper_col]])]

data_annotation$gemma$candidate_publ_for_gemma <- 
  unique(subset(
    x = data_annotation$descriptions, 
    subset = data_annotation$descriptions[[opts_da$des_paper_id_col]] %in% data_annotation$gemma$candidate_publ_for_gemma, 
    select = c(opts_da$des_paper_id_col, opts_da$des_platform_col))) 

message(sprintf('Set platforms %s for publications: %s.', paste(data_annotation$gemma$candidate_publ_for_gemma$Assay, collapse = ', '), paste(data_annotation$gemma$candidate_publ_for_gemma$Pub., collapse = ', ')))



message('Downloading gemma annotation files...')
data_annotation$gemma$gemma_platforms <- download_platforms_from_gemma(unique(data_annotation$gemma$candidate_publ_for_gemma[[opts_da$des_platform_col]])) ### !!! changed from data_annotation$gemma$candidate_publ_for_gemma$Assay
message('... downloaded.')



data_annotation$gemma$pubs_found_in_gemma <- subset(
  x = data_annotation$gemma$candidate_publ_for_gemma, 
  subset = toupper(data_annotation$gemma$candidate_publ_for_gemma[[opts_da$des_platform_col]]) %in% names(data_annotation$gemma$gemma_platforms[!is.na(data_annotation$gemma$gemma_platforms)])) ### !!! changed from data_annotation$gemma$candidate_publ_for_gemma$Assay

data_annotation$gemma$pubs_NOT_found_in_gemma <-  subset(
  x = data_annotation$gemma$candidate_publ_for_gemma, 
  subset = toupper(data_annotation$gemma$candidate_publ_for_gemma[[opts_da$des_platform_col]]) %in% names(data_annotation$gemma$gemma_platforms[is.na(data_annotation$gemma$gemma_platforms)])) ### !!! changed from data_annotation$gemma$candidate_publ_for_gemma$Assay

message(sprintf('Platforms to be pre-annotated in gemma are: %s for publications: %s.', paste(data_annotation$gemma$pubs_found_in_gemma$Assay, collapse = ', '), paste(data_annotation$gemma$pubs_found_in_gemma$Pub., collapse = ', ')))



data_annotation$gemma$leftover$leftovers_not_run_through_gemma_cause_these_pubs_were_analyzed_in_probes_step <- subset(
  x = data_annotation$gemma$input_from_probes_leftover, 
  subset = data_annotation$gemma$input_from_probes_leftover$Pub. %nin% data_annotation$gemma$candidate_publ_for_gemma$Pub.)

data_annotation$gemma$leftover$leftovers_not_run_through_gemma_cause_platform_was_not_found <- subset(x = data_annotation$gemma$input_from_probes_leftover, subset = data_annotation$gemma$input_from_probes_leftover$Pub. %in% data_annotation$gemma$pubs_NOT_found_in_gemma$Pub.)


message('Pre-annotating data with various identfiers found in gemma platform description files (to be annotated in next, identifer-centered step)... ')
data_annotation$gemma$annotation_raw <- get_gemma_annotations_for_data(
  descriptions_df = data_annotation$descriptions, 
  des_platform_col = opts_da$des_platform_col, 
  des_paper_id_col = opts_da$des_paper_id_col, 
  named_list_gemma_platforms = data_annotation$gemma$gemma_platforms, 
  gemma_probe_col = 'ProbeName', 
  input_df = data_annotation$gemma$input_from_probes_leftover, 
  input_paper_id_col = opts_da$input_paper_col, 
  input_probe_col = 'Probe')


data_annotation$gemma$annotation$GPL8160 <- GPL8160_post_gemma_wrapper(annotation_raw = data_annotation$gemma$annotation_raw$GPL8160, divider_of_values_in_serie_str_ = ', ', old_divider_to_be_replaced = '\\|', uniqualize_ = T, Gene_ID_col = 'Gene_ID', gene_description_col = 'Description', gene_symbol_col = 'Symbol', cols_to_nullify = c('GeneSymbols', 'NCBIids', 'GeneNames', 'GOTerms', 'GemmaIDs'))
message('... following platforms were annotated: %s.', paste(names(data_annotation$gemma$annotation), collapse = ', '))


### !!! Now we need to subset finalized and leftovers, and we can move on
data_annotation$gemma$finalized_full_gemma_input_pre_annotated <- rbind(
  data_annotation$gemma$leftover$leftovers_not_run_through_gemma_cause_these_pubs_were_analyzed_in_probes_step, 
  data_annotation$gemma$leftover$leftovers_not_run_through_gemma_cause_platform_was_not_found, 
  data_annotation$gemma$annotation$GPL8160)

i_and_o_has_the_same_nb_of_rows <- length(data_annotation$gemma$finalized_full_gemma_input_pre_annotated[[1]]) == length(data_annotation$gemma$input_from_probes_leftover[[1]])
if (i_and_o_has_the_same_nb_of_rows) {
  message(sprintf(
    'Ok. Input (%s) and output (%s) of this step have the same length.', 
    deparse(substitute(data_annotation$gemma$input_from_probes_leftover[[1]])), 
    deparse(substitute(data_annotation$gemma$finalized_full_gemma_input_pre_annotated[[1]]))))
} else stop('Error! Input and output of this step have different lengths.')



data_annotation$gemma$qa$finalized_full_gemma_input_pre_annotated <- verify_df(df_ = data_annotation$gemma$finalized_full_gemma_input_pre_annotated, only_qa = T)



message('Saving gemma-centered pre-analysis as gemma_centered_pre_analysis...')
gemma <- data_annotation$gemma
save(gemma, file = paste0(opts_da$annotation_folder, '/gemma_centered_pre_analysis'))
message('... saved. Please remember to check quality of data - use "qa" element of saved list.')

rm(i_and_o_has_the_same_nb_of_rows, gemma)

message('*** ... GEMMA PRE-ANNOTATION FINISHED ***')
##################################
### Annotate probes with GEMMA ###
##################################



message('*** MERGING DATA FOR NON-PROBE INDENTIFER-CENTERED ANNOTATION... ***')
message('Merging data originally without probe identifers and unannotated entries from probe-centered analysis (both passed through gemma pre-annotation and not)... ')

if (colnames(data_annotation$gemma$finalized_full_gemma_input_pre_annotated) == colnames(data_annotation$input$data_NOT_for_probe_annotation)) {
  message('Datasets have the same column names')
} else stop('Error! Datasets have different column names')
if (length(colnames(data_annotation$gemma$finalized_full_gemma_input_pre_annotated)) == length(colnames(data_annotation$input$data_NOT_for_probe_annotation))) {
  message('Datasets have the same number of columns')
} else stop('Error! Datasets have different number of columns')

data_annotation$identifers$input <- rbind(data_annotation$gemma$finalized_full_gemma_input_pre_annotated, data_annotation$input$data_NOT_for_probe_annotation)

message('*** ... MERGED! ***')



### !!! Add not searching for identifier if it is NA - unidentifed_identifers - 'AK045385, TC1717724, TC1611597, NA, NA, NA'
#####################################
### Pseudo-memoization non-probes ###
#####################################
message('*** STARTING NON-PROBE INDENTIFER-CENTERED ANNOTATION... ***')

message('Creating database for annotation...')
data_annotation$identifers$pseudo_memoized_annotation_db <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = opts_da$des_paper_id_col,
  des_species_col = opts_da$des_species_col,
  input_df = data_annotation$identifers$input,
  input_paper_col = opts_da$input_paper_col,
  input_id_col = 'Probe',
  A_C_all_ENSG_col = 'ENSG_', A_C_all_ENST_col = 'ENST_', A_C_all_Gene_ID_col = 'Gene_ID', A_C_all_NM_col = 'NM_', A_C_all_Accession_col = 'Accession', A_C_all_Unigene_col = 'Unigene', A_C_all_NR_col = 'NR_', A_C_all_XM_col = 'XM_', A_C_all_XR_col = 'XR_',
  PERFORM_A_prepare_dbs_for_pseudomemoization = T,
  A_return_qa_of_pseudomemoization = F)
message('... created.')



message('Starting annotation...')
data_annotation$identifers$annotations$list <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = opts_da$des_paper_id_col,
  des_species_col = opts_da$des_species_col,
  input_df = data_annotation$identifers$input,
  input_paper_col = opts_da$input_paper_col,
  input_id_col = 'Probe',
  PERFORM_C_annotation = T,
  A_C_all_ENSG_col = 'ENSG_', A_C_all_ENST_col = 'ENST_', A_C_all_Gene_ID_col = 'Gene_ID', A_C_all_NM_col = 'NM_', A_C_all_Accession_col = 'Accession', A_C_all_Unigene_col = 'Unigene', A_C_all_NR_col = 'NR_', A_C_all_XM_col = 'XM_', A_C_all_XR_col = 'XR_', 
  C_pseudo_memoized_db = data_annotation$identifers$pseudo_memoized_annotation_db)
message('... annotated.')



data_annotation$identifers$annotations$df <- rlist::list.rbind(data_annotation$identifers$annotations$list)



message('Subsetting annotation into finalized and leftover data...')

data_annotation$identifers$finalized <- subset(x = data_annotation$identifers$annotations$df, subset = !is.na(data_annotation$identifers$annotations$df$external_gene_name))

data_annotation$identifers$leftover <- subset(x = data_annotation$identifers$annotations$df, subset = is.na(data_annotation$identifers$annotations$df$external_gene_name))



### !!! What about annotated vs input???
was_splitting_ok <- check_was_the_spliting_of_df_by_filtering_ok(
  df_original = data_annotation$identifers$annotations$df, 
  list_df_splited = list(data_annotation$identifers$finalized, data_annotation$identifers$leftover))

if (was_splitting_ok) {
  message(sprintf(
    opts_da$splitting_message_ok, 
    deparse(substitute(data_annotation$identifers$annotations)),
    deparse(substitute(data_annotation$identifers$finalized)),
    deparse(substitute(data_annotation$identifers$leftover))))
} else stop(opts_da$splitting_message_bad)

message('... subsetting done.')



data_annotation$identifers$qa$finalized <- verify_df(df_ = data_annotation$identifers$finalized, only_qa = T)
data_annotation$identifers$qa$leftover <- verify_df(df_ = data_annotation$identifers$leftover, only_qa = T)



message('Saving non-probe identifier-centered analysis as identifier_centered_analysis...')
identifers <- data_annotation$identifers
save(probes_direct, file = paste0(opts_da$annotation_folder, '/identifier_centered_analysis'))
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

data_annotation$symbols$input_from_non_probes_leftovers <- data_annotation$identifers$leftover



message('Translating all gene names in the data to gene ids using ncbi... ')
data_annotation$symbols$ncbi_annotation_of_symbols_to_gene_ids <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = opts_da$des_paper_id_col,
  des_species_col = opts_da$des_species_col,
  input_df = data_annotation$symbols$input_from_non_probes_leftovers,
  input_paper_col = opts_da$input_paper_col,
  input_id_col = 'Symbol',
  C_legacy_D_str_identifier_type__ = 'Gene name',
  PERFORM_D_ncbi_annotation = T, 
  A_D_string_separator__ = ', ')
### !!! add saving midpoints, super important
message('... translated.')




message('Adding ncbi-derived gene ids to original data... ')
data_annotation$symbols$input_plus_new_gene_id_from_ncbi_column <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = opts_da$des_paper_id_col,
  des_species_col = opts_da$des_species_col,
  input_df = data_annotation$symbols$input_from_non_probes_leftovers,
  input_paper_col = opts_da$input_paper_col,
  input_id_col = 'Symbol',
  PERFORM_E_add_new_gene_id_col_originating_from_ncbi = T, 
  E_PERFORM_D_output = data_annotation$symbols$ncbi_annotation_of_symbols_to_gene_ids)
message('.. added.')




message('Creating database for annotation...')
data_annotation$symbols$pseudomemoize_external_gene_names_to_gene_id <- master_annotator(
  descriptions_df = data_annotation$symbols$ncbi_annotation_of_symbols_to_gene_ids,
  des_paper_id_col = 'dummy_paper',
  des_species_col = 'organism',
  input_df = data_annotation$symbols$ncbi_annotation_of_symbols_to_gene_ids,
  input_paper_col = 'dummy_paper',
  input_id_col = 'Gene_ID',
  Gene_ID_col = 'Gene_ID',
  PERFORM_A_prepare_dbs_for_pseudomemoization = T,
  return_qa_of_pseudomemoization = F)
message('... created.')



message('Starting annotation...')
data_annotation$symbols$annotations$list <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = opts_da$des_paper_id_col,
  des_species_col = opts_da$des_species_col,
  input_df = data_annotation$symbols$input_plus_new_gene_id_from_ncbi_column,
  input_paper_col = opts_da$input_paper_col,
  input_id_col = 'Probe',
  PERFORM_C_annotation = T,
  C_legacy_D_str_identifier_type__ = 'Gene_ID',
  Gene_ID_col = 'Gene_ID_from_ncbi',
  C_pseudo_memoized_db = data_annotation$symbols$pseudomemoize_external_gene_names_to_gene_id)
message('... annotated.')


data_annotation$symbols$annotations$df <- rlist::list.rbind(data_annotation$symbols$annotations$list)




message('Subsetting annotation into finalized and leftover data...')
data_annotation$symbols$finalized <- subset(x = data_annotation$symbols$annotations$df, subset = !is.na(data_annotation$symbols$annotations$df$external_gene_name))

data_annotation$symbols$leftover <- subset(x = data_annotation$symbols$annotations$df, subset = is.na(data_annotation$symbols$annotations$df$external_gene_name))

### !!! What about annotated vs input???
was_splitting_ok <- check_was_the_spliting_of_df_by_filtering_ok(
  df_original = data_annotation$symbols$annotations, 
  list_df_splited = list(
    data_annotation$symbols$finalized, 
    data_annotation$symbols$leftover))

if (was_splitting_ok) {
  message(sprintf(
    opts_da$splitting_message_ok, 
    deparse(substitute(data_annotation$symbols$annotations)),
    deparse(substitute(data_annotation$symbols$finalized)),
    deparse(substitute(data_annotation$symbols$leftover))))
} else stop(opts_da$splitting_message_bad)

message('... subsetting done.')


data_annotation$symbols$qa$finalized <- verify_df(df_ = data_annotation$symbols$finalized, only_qa = T)
data_annotation$symbols$qa$leftover <- verify_df(df_ = data_annotation$symbols$leftover, only_qa = T)



message('Saving ncbi-centered analysis as ncbi_centered_analysis...')
symbols <- data_annotation$symbols
save(symbols, file = paste0(opts_da$annotation_folder, '/ncbi_centered_analysis'))
message('... saved. Please remember to check quality of data - use "qa" element of saved list.')

rm(symbols, was_splitting_ok)

message('*** ... GENE SYMBOL-CENTERED ANNOTATION FINISHED! ***')
###########################
### Annotate non-probes ###
###########################









message('*** MERGING ALL FINALIZED DATA... ***')

message('Normalizing columns...')
data_annotation$probes_direct$finalized$identifed_identifer <- data_annotation$probes_direct$finalized$Probe
data_annotation$probes_direct$finalized$unidentifed_identifers <- NA
data_annotation$probes_direct$finalized$Gene_ID_from_ncbi <- NA

data_annotation$identifers$finalized$Gene_ID_from_ncbi <- NA

if (colnames(data_annotation$probes_direct$finalized) == colnames(data_annotation$identifers$finalized) & colnames(data_annotation$identifers$finalized) == colnames(data_annotation$symbols$finalized)) {
  message('Datasets have the same column names')
} else stop('Error! Datasets have different column names')
if (length(colnames(data_annotation$probes_direct$finalized)) == length(colnames(data_annotation$identifers$finalized)) & length(colnames(data_annotation$identifers$finalized)) == length(colnames(data_annotation$symbols$finalized))) {
  message('Datasets have the same number of columns')
} else stop('Error! Datasets have different number of columns')
message('... columns normalized.')

data_annotation$full_finalized$data <- rbind(data_annotation$probes_direct$finalized, data_annotation$identifers$finalized, data_annotation$symbols$finalized)

data_annotation$true_leftovers$data <- data_annotation$symbols$leftover

if ((length(data_annotation$input$data_for_probe_annotation[[1]]) + length(data_annotation$input$data_NOT_for_probe_annotation[[1]])) == (length(data_annotation$full_finalized$data[[1]]) + length(data_annotation$true_leftovers$data[[1]]))) {
  message('Annotation input and out datasets have equal number of entries.')
} else stop('Error! Annotation input and out datasets have different number of entries!')

data_annotation$full_finalized$qa <- verify_df(df_ = data_annotation$full_finalized$data, only_qa = T)
data_annotation$true_leftovers$qa <- verify_df(df_ = data_annotation$true_leftovers$data, only_qa = T)

message('Saving final finalized data as full_finalized_data...')
full_finalized <- data_annotation$full_finalized
save(full_finalized, file = paste0(opts_da$annotation_folder, '/full_finalized_data'))
message('... saved. Please remember to check quality of data - use "qa" element of saved list.')

message('Saving final leftover data as true_leftovers_data...')
true_leftovers <- data_annotation$true_leftovers
save(true_leftovers, file = paste0(opts_da$annotation_folder, '/true_leftovers_data'))
message('... saved. Please remember to check quality of data - use "qa" element of saved list.')


message('*** ... FINALIZED DATA MERGED! ***')

get_all_symbols_in_chrvec()

args(master_annotator)
deparse(sys.call())
Sys.getenv()
message()
paste(ls(), collapse = ', ')
