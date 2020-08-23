source('opts.R')
# rm(list = ls(pattern = 'temp.*|test.*'))
set_logfile(paste0(opts$data[['folder']], "/loggit.log"))

### !!! change naming mapped lists to from within function, not outside it.

### !!! beware when using select option of subset - you may need to check if: if (bazar::is.empty(spread_group) | all(colnames(spread_group) == '0'))

### !!! check also operation in main script, not only functions

test <- data

###########################
### GET PRIMARY DATASET ###
###########################
load(paste0(opts$ann[['folder']], '/data_annotation'))
load(paste0(opts$dir_descriptions, '/descriptions_for_analysis'))

data <- list('final_dataset' = data_annotation$post_annot[['final_dataset']],
             'descriptions' = descriptions[['output']]
             )

rm(data_annotation, descriptions)



data$work_dataset <- data[['final_dataset']]

data$work_dataset[[ opts$data[['exp_and_comp_col']] ]] <- paste0(
  data$work_dataset[[ opts$check_cols_for_input[['Exp.']] ]], 
  '_', 
  data$work_dataset[[ opts$check_cols_for_input[['Comp.']] ]])

data$descriptions[[ opts$data[['exp_and_comp_col']] ]] <- paste0(
  data$descriptions[['Exp.']],
  '_',
  data$descriptions[['Comp.']])




data$work_dataset[[ opts[['ext_gene_name_col']] ]] <- as.character(lapply(
  X = data$work_dataset[[ opts[['ext_gene_name_col']] ]],
  FUN = function(x) {
    return(select_best_geneName_wrapper_for_single_string(x))
  } ) )



data$work_dataset <- subset(
  x = data$work_dataset, 
  subset = abs(data$work_dataset[[ opts$check_cols_for_input[['logFC']] ]]) > opts$ann[['cutoff_for_logfc']] )



data$work_dataset <- data$work_dataset %>%
  dplyr::group_by(
    eval(parse(text = opts$data[['exp_and_comp_col']])),
    eval(parse(text = opts[['ext_gene_name_col']]))) %>%
  dplyr::summarize(
    logFC_median = median( eval(parse(text = opts$check_cols_for_input[['logFC']]))) ) # Medianing values for the same gene within Comparison

colnames(data$work_dataset)[[1]] <- opts$data[['exp_and_comp_col']]
colnames(data$work_dataset)[[2]] <- opts[['ext_gene_name_col']]
colnames(data$work_dataset)[[3]] <- opts$data[['logfc_median_col']]


data$work_dataset <- data$work_dataset[order(
  data$work_dataset[[ opts[['ext_gene_name_col']] ]]),]



data$integrated_dataset <- merge(
  x = data$work_dataset, 
  y = data$descriptions, 
  by = opts$data[['exp_and_comp_col']], 
  all.x = T)



data$whole_dataset <- get_spread_data_and_entries_with_n_or_more_values_per_exp(
  dataset_ = data$work_dataset, 
  spread_key = opts$data[['exp_and_comp_col']], 
  spread_value = opts$data[['logfc_median_col']],
  n_ = 3,
  dataset_name = 'whole_data', 
  save = F)
###########################
### GET PRIMARY DATASET ###
###########################





#########################
### PREPARE SUBGROUPS ###
#########################
data$subgroups$categorical_col_to_analyze <- c('Drug', 'Tissue', 'Sex', 'Age_cat', 'Duration_cat', 'Species')



data$subgroups$data <- purrr::map(
  .x = data$subgroups$categorical_col_to_analyze, 
  .f = function(column_name){
    column <- data$descriptions[[column_name]]
    
    value_types <- unique(column)
    value_types <- value_types[!is.na(value_types)]
    
    experiments_to_subset <- purrr::map(
      .x = value_types,
      .f = function(value_type)
      {
        logic_desc_for_value_type <- data$descriptions[[column_name]] == value_type
        
        desc_for_value_type <- subset(
          x = data$descriptions,
          subset = logic_desc_for_value_type,
          select = c(opts$data[['exp_and_comp_col']], column_name))
        
        logic_work_dataset_ <- data[['work_dataset']][[ opts$data[['exp_and_comp_col']] ]] %in% desc_for_value_type[[ opts$data[['exp_and_comp_col']] ]]
        
        work_dataset_ <- subset(
          x = data[['work_dataset']], 
          subset = logic_work_dataset_)
        
        if (!bazar::is.empty(work_dataset_)) {
          whole_dataset <- get_spread_data_and_entries_with_n_or_more_values_per_exp(
            dataset_ = work_dataset_, 
            spread_key = opts$data[['exp_and_comp_col']], 
            spread_value = opts$data[['logfc_median_col']], 
            n_ = 3, 
            dataset_name = paste0(column_name, '_', value_type)
            )
          
          return(whole_dataset)
        } else return(NA)
      })
    
    names(experiments_to_subset) <- value_types
    
    return(experiments_to_subset)
  })

names(data$subgroups$data) <- data$subgroups$categorical_col_to_analyze
#########################
### PREPARE SUBGROUPS ###
#########################





####################################
### PREPARE ALL RELEVANT DATASET ###
####################################
data$dataset_types <- c("work_dataset", "spread_work_dataset", "work_dataset_3_exps_up", "spread_work_dataset_3_exps_up")



data$all_sets_of_datasets <- purrr::map(
  .x = data$dataset_types, 
  .f = function(set_of_datasets)
    {
    get_dataset_sets_wrapper(
      subgroups_ = data[['subgroups']][['data']], 
      dataset_set_name = set_of_datasets,
      whole_dataset = data[['whole_dataset']])
  })

names(data$all_sets_of_datasets) <- data$dataset_types
####################################
### PREPARE ALL RELEVANT DATASET ###
####################################





###################################
### PREPARE DATA FOR CLUSTERING ###
###################################
purrr::walk2(
  .x = data$all_sets_of_datasets[['spread_work_dataset_3_exps_up']],
  .y = names(data$all_sets_of_datasets[['spread_work_dataset_3_exps_up']]),
  .f = function(sp_dataset, names_of_sp_dataset){
    write_table_for_clustering(
      dataset_ = sp_dataset,
      dataset_name = names_of_sp_dataset)
})
###################################
### PREPARE DATA FOR CLUSTERING ###
###################################





###########################
### PERCENTAGE ANALYSIS ###
###########################
data$percentage_analysis <- purrr::map2(
  .x = data$all_sets_of_datasets[['spread_work_dataset']],
  .y = names(data$all_sets_of_datasets[['spread_work_dataset']]),
  .f = function(dataset, name) {
    
    done <- F
    counter <- 0
    
    while (done == F) {
      message(paste0(name, ': ', counter))
      
      tryCatch( {
        return_ <- get_number_and_percentage_of_directionality_of_exp_and_comp_first_column_names(
          spread_df_ = dataset,
          name_of_df_ = name)
        
        done <- T
        
        counter = counter + 1
        
        Sys.sleep(5)
      },
      error=function(e) {
        Sys.sleep(5)
        warning('something went wrong, trying again in 5 s')
      })
    }
    
    return(return_)
  })
###########################
### PERCENTAGE ANALYSIS ###
###########################





###################################################
### NUMBER OF GENES IN ALL EXPERIMENTS - FIGURE ###
###################################################
data$nb_of_genes_detected_in_given_nb_of_exps$data <- get_number_of_genes_detected_in_given_number_of_experiments(
  dataset_ = data$percentage_analysis[['whole_dataset']],
  number_of_times_detected_col = 'number_of_exp',
  cols_to_remove = c(opts[['ext_gene_name_col']], 'no_of_comps', 'perc_of_upregulated'))



data$nb_of_genes_detected_in_given_nb_of_exps$nb_of_exp_to_collapse_higher_data_into = 10



data$nb_of_genes_detected_in_given_nb_of_exps$for_plot <- collapse_highest_numbers_of_experiments_for_plotting_nb_of_genes_detected_in_given_nb_of_exps(
  data_ = data$nb_of_genes_detected_in_given_nb_of_exps[['data']], 
  nb_of_exp_to_collapse_higher_data_into = data$nb_of_genes_detected_in_given_nb_of_exps[['nb_of_exp_to_collapse_higher_data_into']], 
  number_of_times_detected_col_ = 'number_of_exp') ### !!! Check if this works, if there are missing experiment number values



library(ggplot2)
data$nb_of_genes_detected_in_given_nb_of_exps[['for_plot']] %>%
  ggplot(aes(x = number_of_exp, y = nb_of_genes_detected_in_this_nb_of_papers)) +
  
  geom_col() +
  
  geom_text(data = data$nb_of_genes_detected_in_given_nb_of_exps[['for_plot']],
            aes(
              label = paste0(round(100 * data$nb_of_genes_detected_in_given_nb_of_exps[['for_plot']][['percent_of_genes_detected_in_this_nb_of_papers']], digits = 2), ' %'),
              y = 1000,
              angle = 90),
            size = 4)+
  
  scale_x_continuous(
    breaks = seq(1, data$nb_of_genes_detected_in_given_nb_of_exps[['nb_of_exp_to_collapse_higher_data_into']], 1), 
    labels = c(
      as.character(data$nb_of_genes_detected_in_given_nb_of_exps$for_plot$number_of_exp[1:(which(
        data$nb_of_genes_detected_in_given_nb_of_exps$for_plot$number_of_exp == data$nb_of_genes_detected_in_given_nb_of_exps[['nb_of_exp_to_collapse_higher_data_into']]) - 1)]), 
      paste0(data$nb_of_genes_detected_in_given_nb_of_exps[['nb_of_exp_to_collapse_higher_data_into']], ' <')))+
  
  labs(title = 'number of papers in which gene was detected')+
  
  xlab('number of papers')+
  
  ylab('number of genes')
###################################################
### NUMBER OF GENES IN ALL EXPERIMENTS - FIGURE ###
###################################################





##########################################################
### PREPARE GENESETS TO ANALYZE USING PRIMARY DATASETS ###
##########################################################
data$genesets_to_analyze$groups <- list()

data$genesets_to_analyze$groups[['test']] <- list(
  'genes' = c('akr1a1', 'fasn', 'gap43'),
  'uniformalize_gene_names' = F)

data$genesets_to_analyze$groups[['test_2']] <- list(
  'genes' = 'gapdh',
  'uniformalize_gene_names' = F)

# data$genesets_to_analyze$groups[['gc_genes_all']] <- list(
#   'dataset' = readr::read_tsv(file = 'datasets_and_metadata/gc_input/gc_genes_all_input.txt'),
#   'gene_name_col' = 'gene',
#   'uniformalize_gene_names' = T,
#   'uniformalize_species' = 'mouse')

data$genesets_to_analyze$groups[['gc_genes_88_best']] <- list(
  'dataset' = readr::read_tsv(file = 'datasets_and_metadata/gc_input/gc.tsv')[1:25,],
  'gene_name_col' = 'Gene',
  'uniformalize_gene_names' = T,
  'uniformalize_species' = 'mouse')

data$genesets_to_analyze$groups[['gc_genes_extended']] <- list(
  'dataset' = readr::read_tsv(file = 'datasets_and_metadata/gc_input/extended_gc_input.txt')[1:25,],
  'gene_name_col' = 'lower_final_gene_name',
  'uniformalize_gene_names' = T,
  'uniformalize_species' = 'mouse')

# For human genes - compare with whole dataset and pfc
# data$genesets_to_analyze$groups[['human']] <- list(
#   'genes' = data$all_sets_of_datasets$spread_work_dataset$Species_human[['external_gene_name']],
#   'uniformalize_gene_names' = F)



data$genesets_to_analyze$groups <- extract_genesets_from_datasets(groups_to_analyze = data$genesets_to_analyze$groups)



data$genesets_to_analyze$groups <- purrr::map2(
  .x = data$genesets_to_analyze$groups, 
  .y = names(data$genesets_to_analyze$groups), 
  .f = function(group, group_name){
    
    group[['symbols_input']] <- get_all_symbols_in_chrvec(group[['genes']])
    
    if (group[['uniformalize_gene_names']] == T) {
      
      annotated_df <- ncbi_annotation_for_gene_vectors(
        gene_name_vec = group[['genes']], 
        species = group[['uniformalize_species']])
      
      group[['genes']] <- annotated_df[[ opts[['ext_gene_name_col']] ]]
      
      group[['dataset_annotated']] <- merge(
        group[['dataset']], 
        annotated_df, 
        by.x = group[['gene_name_col']], 
        by.y = 'input_id_col', 
        all.x = T)
      
      group[['symbols_output']] <- get_all_symbols_in_chrvec(group[['genes']])
    }

    return(group)
  })

### !!! compare dataset annotated with genes on larger dataset

##########################################################
### PREPARE GENESETS TO ANALYZE USING PRIMARY DATASETS ###
##########################################################





#########################################
### GET PERCENTAGE NUMBERS FOR GROUPS ###
#########################################
data$genesets_to_analyze$groups_and_numbers <- purrr::map(
  .x = data$genesets_to_analyze$groups, 
  .f = function(group){
    
    merge_group_and_numbers_wrapper(
      group_ = group, 
      dataset_to_merge_group_with = data$percentage_analysis$whole_dataset) 
  })

names(data$genesets_to_analyze$groups_and_numbers) <- names(data$genesets_to_analyze$groups)
#########################################
### GET PERCENTAGE NUMBERS FOR GROUPS ###
#########################################





###################################
### GET FIGURES FOR GENE GROUPS ###
###################################
data$genesets_to_analyze$mapped_groups <- map_data_from_list_of_gene_groups(
  gene_groups = data$genesets_to_analyze[['groups']], 
  gather_dataset = data$whole_dataset[['work_dataset']],
  columns_to_analyze_charvec = data$subgroups[['categorical_col_to_analyze']], 
  element_name_with_gene_names_within_gene_group_element = 'genes',
  descriptions = data$descriptions)



print_figures_for_mapped_groups_wrapper(map_data_from_list_of_gene_groups_output = data$genesets_to_analyze$mapped_groups) 
###################################
### GET FIGURES FOR GENE GROUPS ###
###################################





###################################
### GET CLUSTER FOR GENE GROUPS ###
###################################
data$genesets_to_analyze$clustering_spread_matrixi_for_groups  <- purrr::map(
  .x = data$genesets_to_analyze$groups,
  .f = function(group){

    spread_group <- prepare_spread_data_for_getting_preety_clusters(
      genesets_to_analyze_group = group, 
      logfc_cutoff = 0.1, 
      ratio_of_0_values_cutoff = 0.7,
      spread_datasets = data$all_sets_of_datasets[['spread_work_dataset']],
      gene_name_col_in_geneset_to_analyze_group = 'genes')
    
    return(spread_group)
  }) # All of the genes from this cluster and given experimental factor are there, even if they are 0 in every comparison

names(data$genesets_to_analyze[['clustering_spread_matrixi_for_groups']]) <- names(data$genesets_to_analyze[['groups']])



purrr::walk2(
  .x = data$genesets_to_analyze[['clustering_spread_matrixi_for_groups']],
  .y = names(data$genesets_to_analyze[['clustering_spread_matrixi_for_groups']]),
  .f = function(matrix_list, matrix_list_name){
    
    print_tables_and_figures_for_preety_clusters_wrapper(
      spread_matrixi_list = matrix_list, 
      cluster_name = matrix_list_name,
      comparisons_to_print = 'whole_dataset')
  })
###################################
### GET CLUSTER FOR GENE GROUPS ###
###################################





##################################################
### GET GENERAL SPREAD MATRIXI FOR GENE GROUPS ###
##################################################
data$genesets_to_analyze$spread_matrixi_for_groups  <- purrr::map(
  .x = data$genesets_to_analyze$groups,
  .f = function(group){
    
    spread_group <- prepare_spread_data_for_getting_preety_clusters(
      genesets_to_analyze_group = group, 
      logfc_cutoff = 0, 
      ratio_of_0_values_cutoff = (ncol(data$all_sets_of_datasets[['spread_work_dataset']]$whole_dataset) - 0.1)/ncol(data$all_sets_of_datasets[['spread_work_dataset']]$whole_dataset),
      spread_datasets = data$all_sets_of_datasets[['spread_work_dataset']],
      gene_name_col_in_geneset_to_analyze_group = 'genes')
    
    return(spread_group)
  }) # All of the genes from this cluster and given experimental factor are there, even if they are 0 in every comparison

names(data$genesets_to_analyze[['spread_matrixi_for_groups']]) <- names(data$genesets_to_analyze[['groups']])
##################################################
### GET GENERAL SPREAD MATRIXI FOR GENE GROUPS ###
##################################################


### !!! starts here


##############################################################
### ARE GENES IN THE GROUP EXPRESSED IN THE SAME DIRCTIONS ###
##############################################################
data$genesets_to_analyze$expr_stability <- purrr::map2(
  .x = data$genesets_to_analyze$spread_matrixi_for_groups, 
  .y = names(data$genesets_to_analyze$spread_matrixi_for_groups), 
  .f = function(gene_group, gene_group_name){
    
    are_genes_in_group_expressed_in_the_same_direction(
      spread_data_list_from_genegroup = gene_group, 
      name_ = gene_group_name)
    
  })

names(data$genesets_to_analyze$expr_stability) <- names(data$genesets_to_analyze$spread_matrixi_for_groups)
##############################################################
### ARE GENES IN THE GROUP EXPRESSED IN THE SAME DIRCTIONS ###
##############################################################








#############################################################
### ADD INFO ON ENRICHMENT OF CLUSTER IN GIVEN EXPERIMENT ###
#############################################################
data$genesets_to_analyze$individual_genes_in_group_from_all_comps_gathered <- get_individual_genes_in_group_from_all_comps_gathered(
  groups_ = data$genesets_to_analyze[['groups']], 
  full_work_dataset = data$whole_dataset[['work_dataset']])



data$genesets_to_analyze$is_value_from_column_enriched$named_list_of_cols_and_values_to_analyze <- list(
  'Drug' = list(
    'col_name' = 'Drug',
    'values_of_interest' = 'fluoxetine'),
  'Tissue' = list(
    'col_name' = 'Tissue',
    'values_of_interest' = c('hippocampus', 'amygdala'))
  )



data$genesets_to_analyze$is_value_from_column_enriched$primary_data_on_column_value_enrichment <- get_enrichment_for_values_in_columns_for_groups_wrapper(
  named_list_of_cols_and_values_to_analyze_ = data$genesets_to_analyze$is_value_from_column_enriched[['named_list_of_cols_and_values_to_analyze']],
  individual_genes_in_group_from_all_comps_gathered_ = data$genesets_to_analyze[['individual_genes_in_group_from_all_comps_gathered']],
  full_dataset_ = data$whole_dataset[['work_dataset']])



data$genesets_to_analyze$is_value_from_column_enriched$how_many_times_was_the_gene_detected_in_paper <- get_number_of_times_the_gene_detected_in_paper(
  individual_genes_in_group_from_all_comps_gathered_ = data$genesets_to_analyze[['individual_genes_in_group_from_all_comps_gathered']], 
  pattern_to_extract_exp_from_exp_and_comp = '_.*')
  

  
data$genesets_to_analyze$is_value_from_column_enriched$number_of_exps_for_paper_per_gene <- purrr::map(
  .x = data$genesets_to_analyze$is_value_from_column_enriched[['how_many_times_was_the_gene_detected_in_paper']], 
  .f = function(how_many_times_was_the_gene_detected_in_paper_)
{
  get_number_of_exps_for_paper_per_gene_wrapper(
    how_many_times_were_genes_detected_in_paper = how_many_times_was_the_gene_detected_in_paper_, 
    experiment_to_check_number_of_comparisons_in = c('20'))
})



data$genesets_to_analyze$is_value_from_column_enriched$exps_in_which_gene_was_present_for_given_value_of_interest <- purrr::map(
  .x = data$genesets_to_analyze$individual_genes_in_group_from_all_comps_gathered, 
  .f = function(gene_group)
    {
    get_exps_in_which_gene_was_present_for_given_value_of_interest_wrapper(
      gene_group_ = gene_group,
      list_for_subsetting_name_is_colname_value_is_value_of_interest = data$genesets_to_analyze$is_value_from_column_enriched[['named_list_of_cols_and_values_to_analyze']] )
  })
#############################################################
### ADD INFO ON ENRICHMENT OF CLUSTER IN GIVEN EXPERIMENT ###
#############################################################





######################
### OTHER ANALYSES ###
######################
data$other_analyses$what_percentage_of_our_data_are_protein_coding_genes <- get_percentage_of_our_data_are_protein_coding_genes(
  species = 'mouse', 
  spread_gene_names_col = data$whole_dataset[['spread_work_dataset']][['external_gene_name']])





data$other_analyses$compare_human_and_pfc_subset <- merge(
  x = data$all_sets_of_datasets$work_dataset[['Species_human']], 
  y = data$all_sets_of_datasets$work_dataset[['Tissue_frontal cortex']], 
  by = opts[['ext_gene_name_col']],
  all.x = T)

names(data$other_analyses[['compare_human_and_pfc_subset']]) <- c(opts[['ext_gene_name_col']], 'EC_ID_human', 'logFC_median_human', 'EC_ID_pfc', 'logFC_median_pfc')





data$other_analyses$experiments_in_actual_analysis$exps <- unique(stringr::str_remove(
  string = data$whole_dataset[['work_dataset']][['EC_ID']], 
  pattern = '_.*') )

data$other_analyses$experiments_in_actual_analysis$nb_of_exps <- length(data$other_analyses$experiments_in_actual_analysis[['exps']])

######################
### OTHER ANALYSES ###
######################

save(data, file = 'temp_data')






#####################
### ENRICH GROUPS ###
#####################

# INVOKE PYTHON SCRIPT FROM R

#####################
### ENRICH GROUPS ###
#####################







