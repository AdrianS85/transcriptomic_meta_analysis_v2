# test3 <- unique(data$work_dataset_3_exps_up$external_gene_name)
# 
# 
# unique(data$descriptions[['Duration_cat']])
# 
data$subgroups$categorical_col_to_analyze <- c('Drug', 'Tissue', 'Sex', 'Age_cat', 'Duration_cat', 'Species')



test <- function(){
  purrr::map(
  .x = data$subgroups$categorical_col_to_analyze, 
  .f = function(column_name){
    column <- data$descriptions[[column_name]]
    
    value_types <- unique(column)
    value_types <- value_types[!is.na(value_types)]
    
    experiments_to_subset <- purrr::map(
      .x = value_types,
      .f = function(value_type)
      {
        desc_for_value_type <- subset(
          x = data$descriptions,
          subset = data$descriptions[[column_name]] == value_type,
          select = c('EC_ID', column_name))
        
        work_dataset_ <- subset(x = data[['work_dataset']], subset = data[['work_dataset']][['EC_ID']] %in% desc_for_value_type[['EC_ID']])
        
        if (length(work_dataset_[[1]]) != 0) {
          whole_dataset <- get_spread_data_and_entries_with_n_or_more_values_per_exp(
            dataset_ = work_dataset_, 
            spread_key = 'EC_ID', 
            spread_value = 'logFC_median', 
            gene_name_col_ = opts[['ext_gene_name_col']], 
            n_ = 3)
          
          return(whole_dataset)
        } else return(NA)
      })
    
    names(experiments_to_subset) <- value_types
    
    return(experiments_to_subset)
  })

# names(test) <- data$subgroups$categorical_col_to_analyze
}


test_x <- test()




# test <- function(){
#   
#   purrr::map(
#     .x = data$subgroups$categorical_col_to_analyze, 
#     .f = function(column_name){
#       column <- data$descriptions[[column_name]]
#       
#       value_types <- unique(column)
#       value_types <- value_types[!is.na(value_types)]
#       
#       experiments_to_subset <- purrr::map(
#         .x = value_types,
#         .f = function(value_type)
#         {
#           desc_for_value_type <- subset(
#             x = data$descriptions,
#             subset = data$descriptions[[column_name]] == value_type,
#             select = c('EC_ID', column_name))
#           
#           work_dataset_ <- subset(x = data[['work_dataset']], subset = data[['work_dataset']][['EC_ID']] %in% desc_for_value_type[['EC_ID']])
#           
#           spread_work_dataset_ <- tidyr::spread(
#             data = work_dataset_, 
#             key = EC_ID, 
#             value = logFC_median)
#           
#           ### !!!
#           bool_entries_with_3_or_more_values_ <- get_subset_vector_for_entries_with_n_or_more_values_per_exp(
#             work_dataset_ = work_dataset_,
#             exp_col = 'EC_ID',
#             gene_name_col = opts[['ext_gene_name_col']],
#             n = 3)
#           
#           spread_work_dataset_3_exps_up_ <- subset(
#             x = spread_work_dataset_,
#             subset = bool_entries_with_3_or_more_values_)
#           
#           work_dataset_3_exps_up_ <- subset(
#             x = work_dataset_,
#             subset = work_dataset_[[ opts[['ext_gene_name_col']] ]] %in% spread_work_dataset_3_exps_up_[[ opts[['ext_gene_name_col']]  ]] )
#           
#           are_any_cols_empty_in_3_exps_up_data_ <- are_any_cols_empty_in_3_exps_up_data_wrapper(
#             spread_work_dataset_3_exps_up_ = spread_work_dataset_3_exps_up_)
#           
#           
#           
#           return(list('work_dataset' = work_dataset_, 'spread_work_dataset_' = spread_work_dataset_, 'work_dataset_3_exps_up_' = work_dataset_3_exps_up_, 'spread_work_dataset_3_exps_up_' = spread_work_dataset_3_exps_up_, 'are_any_cols_empty_in_3_exps_up_data_' = are_any_cols_empty_in_3_exps_up_data_))
#           ### !!!
#           
#           
#           # sub_datasets <- rlist::list.ungroup(sub_datasets)
#           # 
#           # sub_spread_datasets <- rlist::list.ungroup(sub_spread_datasets)
#           # 
#           # return(list('desc_for_value_type' = desc_for_value_type, 'sub_datasets' = sub_datasets, 'sub_spread_datasets' = sub_spread_datasets))
#         })
#       
#       names(experiments_to_subset) <- value_types
#       
#       return(experiments_to_subset)
#     })
#   
#   names(test) <- data$subgroups$categorical_col_to_analyze
# }