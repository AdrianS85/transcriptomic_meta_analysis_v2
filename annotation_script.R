source('test.R')
source('functions_for_annotation.R')
load('data_prepared_for_annotation')


### !!! need to add step which checks wheter all papers have unique platforms
### !!! acutally, we need to use experiment id, instead of paper id, which has the same platform... species is implied then. perhaps lets call it dataset?



### Prepare working list ###
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
### Prepare working list ###



#######################
### Annotate probes ###
#######################

# test <- subset(x = data_prepared_for_annotation$data_for_probe_annotation, subset = data_prepared_for_annotation$data_for_probe_annotation$Pub. %in% c(34)) 



Sys.time()
data_annotation$probes_direct$annotations <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = 'Pub.',
  des_species_col = 'Species',
  input_df = data_annotation$input$data_for_probe_annotation,
  input_paper_col = 'Pub.',
  input_id_col = 'Probe',
  str_experiment_name = 'probes_direct',
  str_identifier_type__ = data_annotation$input$platform_ids_for_probe_annotation)
Sys.time()

# length(unique(data_annotation$annotations$probes$Pub.))
# 
# length(unique(data_annotation$data_for_probe_annotation$Pub.))
# 
# unique(data_annotation$data_for_probe_annotation$Pub.)[order(unique(data_annotation$data_for_probe_annotation$Pub.))]
# test <- dplyr::select(data_annotation$finalized$probes, colnames(data_annotation$platform_ids_for_probe_annotation), dplyr::everything())

# Prepare finalized and leftover files
data_annotation$probes_direct$finalized <- subset(x = data_annotation$probes_direct$annotations, subset = !is.na(data_annotation$probes_direct$annotations$external_gene_name))
# save(data_annotation$probes_direct$finalized, file = 'finalized_probes')
# load(file = 'finalized_probes')

data_annotation$probes_direct$leftover <- subset(x = data_annotation$probes_direct$annotations, subset = is.na(data_annotation$probes_direct$annotations$external_gene_name))
# save(data_annotation$probes_direct$leftover, file = 'leftover_probes')
# load(file = 'leftover_probes')

check_was_the_spliting_of_df_by_filtering_ok(df_original = data_annotation$probes_direct$annotations, list_df_splited = list(data_annotation$probes_direct$finalized, data_annotation$probes_direct$leftover))
#######################
### Annotate probes ###
#######################



##################################
### Annotate probes with GEMMA ###
##################################
data_annotation$gemma$input_from_probes_leftover <- dplyr::select(data_annotation$probes_direct$leftover, -'external_gene_name')

data_annotation$gemma$candidate_publ_for_gemma <- unique(data_annotation$probes_direct$leftover$Pub.)[unique(data_annotation$probes_direct$leftover$Pub.) %nin% unique(data_annotation$probes_direct$finalized$Pub.)]

data_annotation$gemma$candidate_publ_for_gemma <- unique(subset(x = data_annotation$descriptions, subset = data_annotation$descriptions$Pub. %in% data_annotation$gemma$candidate_publ_for_gemma, select = c('Pub.', 'Assay')))

### !!! I need to fork and secure devtools::install_github('PavlidisLab/gemmaAPI.R')
data_annotation$gemma$gemma_platforms <- download_platforms_from_gemma(unique(data_annotation$gemma$candidate_publ_for_gemma$Assay))

data_annotation$gemma$pubs_found_in_gemma <- subset(
  x = data_annotation$gemma$candidate_publ_for_gemma, 
  subset = toupper(data_annotation$gemma$candidate_publ_for_gemma$Assay) %in% names(data_annotation$gemma$gemma_platforms[!is.na(data_annotation$gemma$gemma_platforms)]))

data_annotation$gemma$pubs_NOT_found_in_gemma <-  subset(
  x = data_annotation$gemma$candidate_publ_for_gemma, 
  subset = toupper(data_annotation$gemma$candidate_publ_for_gemma$Assay) %in% names(data_annotation$gemma$gemma_platforms[is.na(data_annotation$gemma$gemma_platforms)]))


data_annotation$gemma$leftover$leftovers_not_run_through_gemma_cause_these_pubs_were_analyzed_in_probes_step <- subset(
  x = data_annotation$gemma$input_from_probes_leftover, 
  subset = data_annotation$gemma$input_from_probes_leftover$Pub. %nin% data_annotation$gemma$candidate_publ_for_gemma$Pub.)

data_annotation$gemma$leftover$leftovers_not_run_through_gemma_cause_platform_was_not_found <- subset(x = data_annotation$gemma$input_from_probes_leftover, subset = data_annotation$gemma$input_from_probes_leftover$Pub. %in% data_annotation$gemma$pubs_NOT_found_in_gemma$Pub.)

data_annotation$gemma$annotation_raw <- get_gemma_annotations_for_data(
  descriptions_df = data_annotation$descriptions, 
  des_platform_col = 'Assay', 
  des_paper_id_col = 'Pub.', 
  named_list_gemma_platforms = data_annotation$gemma$gemma_platforms, 
  gemma_probe_col = 'ProbeName', 
  input_df = data_annotation$gemma$input_from_probes_leftover, 
  input_paper_id_col = 'Pub.', 
  input_probe_col = 'Probe')


data_annotation$gemma$annotation$GPL8160 <- GPL8160_post_gemma_wrapper(annotation_raw = data_annotation$gemma$annotation_raw$GPL8160, divider_of_values_in_serie_str_ = ', ', old_divider_to_be_replaced = '\\|', uniqualize_ = T, Gene_ID_col = 'Gene_ID', gene_description_col = 'Description', gene_symbol_col = 'Symbol', cols_to_nullify = c('GeneSymbols', 'NCBIids', 'GeneNames', 'GOTerms', 'GemmaIDs'))

### !!! Now we need to subset finalized and leftovers, and we can move on
data_annotation$gemma$finalized_full_gemma_input_pre_annotated <- rbind(
  data_annotation$gemma$leftover$leftovers_not_run_through_gemma_cause_these_pubs_were_analyzed_in_probes_step, 
  data_annotation$gemma$leftover$leftovers_not_run_through_gemma_cause_platform_was_not_found, 
  data_annotation$gemma$annotation$GPL8160)

length(data_annotation$gemma$finalized_full_gemma_input_pre_annotated[[1]]) == length(data_annotation$gemma$input_from_probes_leftover[[1]])
##################################
### Annotate probes with GEMMA ###
##################################



### Unanotated probes from annotate probes need to be merged with non-probes ###
colnames(data_annotation$gemma$finalized_full_gemma_input_pre_annotated) == colnames(data_annotation$input$data_NOT_for_probe_annotation)

length(colnames(data_annotation$gemma$finalized_full_gemma_input_pre_annotated)) == length(colnames(data_annotation$input$data_NOT_for_probe_annotation))

data_annotation$identifers$input <- rbind(data_annotation$gemma$finalized_full_gemma_input_pre_annotated, data_annotation$input$data_NOT_for_probe_annotation)
### Unanotated probes from annotate probes need to be merged with non-probes ###



### !!! Add not searching for identifier if it is NA - unidentifed_identifers - 'AK045385, TC1717724, TC1611597, NA, NA, NA'
#####################################
### Pseudo-memoization non-probes ###
#####################################
data_annotation$identifers$pseudo_memoized_annotation_db <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = 'Pub.',
  des_species_col = 'Species',
  input_df = data_annotation$identifers$input,
  input_paper_col = 'Pub.',
  input_id_col = 'Probe',
  str_experiment_name = 'pseudo_memoization',
  str_identifier_type__ = 'pseudo_memoization',
  ENSG_col = 'ENSG_', ENST_col = 'ENST_', Gene_ID_col = 'Gene_ID', NM_col = 'NM_', Accession_col = 'Accession', Unigene_col = 'Unigene', NR_col = 'NR_', XM_col = 'XM_', XR_col = 'XR_',
  should_i_prepare_dbs_for_pseudomemoization = T,
  return_qa_of_pseudomemoization = F)



data_annotation$identifers$annotations$list <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = 'Pub.',
  des_species_col = 'Species',
  input_df = data_annotation$identifers$input,
  input_paper_col = 'Pub.',
  input_id_col = 'Probe',
  str_experiment_name = 'annotations',
  str_identifier_type__ = 'all',
  ENSG_col = 'ENSG_', ENST_col = 'ENST_', Gene_ID_col = 'Gene_ID', NM_col = 'NM_', Accession_col = 'Accession', Unigene_col = 'Unigene', NR_col = 'NR_', XM_col = 'XM_', XR_col = 'XR_', 
  pseudo_memoized_db = data_annotation$identifers$pseudo_memoized_annotation_db)




data_annotation$identifers$annotations$df <- rlist::list.rbind(data_annotation$identifers$annotations$list)

# Prepare finalized and leftover files
data_annotation$identifers$finalized <- subset(x = data_annotation$identifers$annotations$df, subset = !is.na(data_annotation$identifers$annotations$df$external_gene_name))
# save(data_annotation$identifers$finalized, file = 'finalized_non_probes')
# load(file = 'finalized_non_probes')

data_annotation$identifers$leftover <- subset(x = data_annotation$identifers$annotations$df, subset = is.na(data_annotation$identifers$annotations$df$external_gene_name))
# save(data_annotation$identifers$leftover, file = 'leftover_non_probes')
# load(file = 'leftover_non_probes')

check_was_the_spliting_of_df_by_filtering_ok(df_original = data_annotation$identifers$annotations$df, list_df_splited = list(data_annotation$identifers$finalized, data_annotation$identifers$leftover))

#####################################
### Pseudo-memoization non-probes ###
#####################################




#####################################
### Annotate gene names with ncbi ###
#####################################
data_annotation$symbols$input_from_non_probes_leftovers <- data_annotation$identifers$leftover




# test_ <- data_annotation$identifers$input[1:20,]

# Get Gene_IDs for each gene Symbol
data_annotation$symbols$ncbi_annotation_of_symbols_to_gene_ids <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = 'Pub.',
  des_species_col = 'Species',
  input_df = data_annotation$symbols$input_from_non_probes_leftovers,
  input_paper_col = 'Pub.',
  input_id_col = 'Symbol',
  str_experiment_name = 'annotations',
  str_identifier_type__ = 'Gene name',
  perform_ncbi_annotation = T, 
  string_separator__ = ', ')
### !!! add saving midpoints, super important



### Concominant steps ###
data_annotation$symbols$input_plus_new_gene_id_from_ncbi_column <- add_new_gene_id_col_originating_from_ncbi_annotation(
  descriptions_df = data_annotation$descriptions, 
  desc_paper_id_col = 'Pub.', 
  desc_organism_col = 'Species',
  input_df = data_annotation$symbols$input_from_non_probes_leftovers, 
  input_paper_id_col = 'Pub.',
  input_input_id_col = 'Symbol', 
  input_organism_col = 'Species', 
  perform_ncbi_annotation_output = data_annotation$symbols$ncbi_annotation_of_symbols_to_gene_ids)


# Get pseudomemoized list containing external_gene_name_for each Gene_ID
data_annotation$symbols$pseudomemoize_external_gene_names_to_gene_id <- master_annotator(
  descriptions_df = data_annotation$symbols$ncbi_annotation_of_symbols_to_gene_ids,
  des_paper_id_col = 'dummy_paper',
  des_species_col = 'organism',
  input_df = data_annotation$symbols$ncbi_annotation_of_symbols_to_gene_ids,
  input_paper_col = 'dummy_paper',
  input_id_col = 'Gene_ID',
  str_experiment_name = 'pseudo_memoization',
  str_identifier_type__ = 'pseudo_memoization',
  Gene_ID_col = 'Gene_ID',
  should_i_prepare_dbs_for_pseudomemoization = T,
  return_qa_of_pseudomemoization = F)
### Concominant steps ###


### !!! Symbols analysis overrides unidentified_idetifiers columns! It happenes on following step:

# Get pseudoannotate external_gene_name to Gene_IDs in actual data
data_annotation$symbols$annotations$list <- master_annotator(
  descriptions_df = data_annotation$descriptions,
  des_paper_id_col = 'Pub.',
  des_species_col = 'Species',
  input_df = data_annotation$symbols$input_plus_new_gene_id_from_ncbi_column,
  input_paper_col = 'Pub.',
  input_id_col = 'Probe',
  str_experiment_name = 'annotations',
  str_identifier_type__ = 'Gene_ID',
  Gene_ID_col = 'Gene_ID_from_ncbi',
  pseudo_memoized_db = data_annotation$symbols$pseudomemoize_external_gene_names_to_gene_id)


data_annotation$symbols$annotations$df <- rlist::list.rbind(data_annotation$symbols$annotations$list)

# Prepare finalized and leftover files
data_annotation$symbols$finalized <- subset(x = data_annotation$symbols$annotations$df, subset = !is.na(data_annotation$symbols$annotations$df$external_gene_name))
# save(data_annotation$symbols$finalized, file = 'finalized_non_probes')
# load(file = 'finalized_non_probes')

data_annotation$symbols$leftover <- subset(x = data_annotation$symbols$annotations$df, subset = is.na(data_annotation$symbols$annotations$df$external_gene_name))
# save(data_annotation$symbols$leftover, file = 'leftover_non_probes')
# load(file = 'leftover_non_probes')

check_was_the_spliting_of_df_by_filtering_ok(df_original = data_annotation$symbols$annotations$df, list_df_splited = list(data_annotation$symbols$finalized, data_annotation$symbols$leftover))












# test_x <- data_annotation$identifers$pseudo_memoized_annotation_db
# 
# 
# 
# class(test$dummy_paper)
# 
# test_3 <- as.data.frame(t(test))
# 
# test2 <- purrr::map(.x = test, .f = function(x) {x$ids})
# 
# 
# data_annotation$descriptions$Species
# 
# data_annotation$symbols$input_from_non_probes_leftovers$Symbol
# 
# 
# test <- change_vector_of_mixed_normal_and_c_geneNames_into_unique_geneName(char_vec = data_annotation$symbols$input_from_non_probes_leftovers$Symbol, regex_pattern_to_split_with = ', ')
# 
# 
# 
# test <- split_string_by_pattern_and_extract_vec_of_unique_values(chr_vec = data_annotation$symbols$input_from_non_probes_leftovers$Symbol, pattern_to_split_individual_strings_with = ', ')
# 
# 
# 
# testx <- merge(x = data_annotation$symbols$input_from_non_probes_leftovers, y = data_annotation$descriptions, by = 'Pub.', all.x = T)
# 
# 
# test2 <- get_vector_of_single_unique_gene_ids_and_species(
#   input_df = testx,
#   identifer_col = 'Symbol', 
#   species_col = 'Species', 
#   string_separator = ', ')

###########################
### Annotate non-probes ###
###########################
















# Merge leftover_Probe_ID z data_prepared_for_annotation$data_NOT_for_probe_annotation:
# 1) kolejność kolumn sié chyba zmienia?
# Make a list, with elements
#ENSG_
#ENST_
#Gene_ID
#NM_
#Accession
#Unigene	XM_	XR_	NR_	Nucleotide
#Symbol
# test <- return_all_usable_id_types(ENSG_ = 'ENSG_', ENST_ = 'ENST_', Gene_ID = 'Gene_ID', NM_ = 'NM_', Accession = 'Accession', Unigene = 'Unigene', NR_ = 'NR_', XM_ = 'XM_', XR_ = 'XR_')
# 
# for (n in test){
#   print(n$filter_name)
# }
# 
# test[[1]]
# test_mem <- data_prepared_for_annotation$data_NOT_for_probe_annotation$NM_
# 
# test_mem2 <- purrr::map(.x = test_mem, .f = function(x){
#   stringr::str_split(string = x, pattern = ', ', simplify = T)
# })
# 
# test_mem3 <- unique(as.character(rlist::list.cbind(test_mem2)))
# class(test_mem3)
# 
# test_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
# 
# test_mem4 <- biomaRt::getBM(
#   attributes = c('refseq_mrna', "external_gene_name"),
#   filters = 'refseq_mrna',
#   values = test_mem3,
#   uniqueRows = F,
#   mart = test_mart
# )
# 
# 
# data_prepared_for_annotation$annotations$data_list_non_probes <- split(data_prepared_for_annotation$data_NOT_for_probe_annotation, f = data_prepared_for_annotation$data_NOT_for_probe_annotation[['Pub.']])
# 
# data_prepared_for_annotation$annotations$output_data_list_non_probes <- list()
# 
# for(pub_n in seq_along(data_prepared_for_annotation$annotations$data_list_non_probes))
# {
#   Sys.time()
#   data_prepared_for_annotation$annotations$output_data_list_non_probes[[pub_n]] <- master_annotator(
#     descriptions_df = data_prepared_for_annotation$descriptions, 
#     des_paper_id_col = 'Pub.', 
#     des_species_col = 'Species', 
#     input_df = data_prepared_for_annotation$data_list_non_probes[[pub_n]], 
#     input_paper_col = 'Pub.', 
#     input_id_col = 'Probe', 
#     str_experiment_name = 'non-probes', 
#     str_identifier_type__ = 'all',
#     ENSG_ = 'ENSG_', ENST_ = 'ENST_', Gene_ID = 'Gene_ID', NM_ = 'NM_', Accession = 'Accession', Unigene = 'Unigene', NR_ = 'NR_', XM_ = 'XM_', XR_ = 'XR_')
#   Sys.time()
# }



# test <- set_mart_to_be_used(species_ = 'mice', int_loop = n, mouse_name = normalized_species_names$mouse, rat_name = normalized_species_names$rat, human_name = normalized_species_names$human, sheep_name = normalized_species_names$sheep, saimiri_name = normalized_species_names$saimiri)

# Sys.time()
# test_2_annotations_Probe_ID <- master_annotator(
#   descriptions_df = data_annotation$descriptions,
#   des_paper_id_col = 'Pub.',
#   des_species_col = 'Species',
#   input_df = test,
#   input_paper_col = 'Pub.',
#   input_id_col = 'Probe',
#   str_experiment_name = 'test_2',
#   str_identifier_type__ = 'all',
#   ENSG_ = 'ENSG_', ENST_ = 'ENST_', Gene_ID = 'Gene_ID', NM_ = 'NM_', Accession = 'Accession', Unigene = 'Unigene', NR_ = 'NR_', XM_ = 'XM_', XR_ = 'XR_')
# Sys.time()
# # refseq_mrna
# ### !!! i think we need to remove properly those ids whish are na
# # save(annotations_Probe_ID, file = 'annotations_Probe_ID')