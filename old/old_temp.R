

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













# test_desc <- subset(x = id_platf$pubs_for_probe_analysis$desc, subset = id_platf$pubs_for_probe_analysis$desc$Pub. == c(17, 18) )
# 
# test_input <- subset(x = id_platf$pubs_for_probe_analysis$input, subset = id_platf$pubs_for_probe_analysis$input$Pub. == c(17, 18))
# 
# test <- master_annotator(
#   descriptions_df = id_platf$desc, 
#   des_paper_id_col = 'Pub.', 
#   des_species_col = 'Species', 
#   input_df = test_input, 
#   input_paper_col = 'Pub.', 
#   input_id_col = 'Probe', 
#   str_experiment_name = 'test',
#   PERFORM_B_get_the_highest_hit_returning_id_type = T,
#   highest_hit_int_Probe_IDs_to_test = 50)

# test <- get_the_highest_hit_returning_id_type(
#   descriptions_df = test_desc,
#   des_paper_id_col = 'Pub.',
#   des_species_col = 'Species',
#   input_df = test_input,
#   input_paper_col = 'Pub.',
#   input_probe_col = 'Probe',
#   int_Probe_IDs_to_test = 200,
#   str_experiment_name = 'test')
#Check if the platforms were selected properly
#Get the actual number of annotated probes (in some platforms single probe can get multiple output rows)
#Add percentage for highest hitting platform - just add highest_annotated_identifier_percentages write to resulting list - OK
#Add species for which we were annotating stuff - OK
#Add exp numbers to output lists - OK
#Number of input genes - OK
#Check if col type validators are working


# testMart <- biomaRt::useMart(
#   "ENSEMBL_MART_ENSEMBL",
#   dataset = "mmusculus_gene_ensembl")
# 
# testMart_filt <- biomaRt::listFilters(mart = testMart) #MouseWG-6
# 
# test <- biomaRt::getBM(
#   attributes = c('', "external_gene_name"),
#   filters = '',
#   values = 'ILMN_1219935',
#   uniqueRows = T,
#   mart = testMart)

# id_platf$highest_hit_analysis$id_to_be_used_for_annotation <- as.character(id_platf$highest_hit_analysis$platform_to_use)