# rm(list = ls(pattern = 'temp.*|test.*'))
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')


set_col_check <- function()
{
  col_check_ <- list('Pub.' = 'Pub.',	'Exp.' = 'Exp.',	'ID' = 'ID',	'adj,P,Val' = 'adj,P,Val',	'P,Value' = 'P,Value',	't' = 't',	'B' = 'B',	'logFC' = 'logFC', 'Probe' = 'Probe',	'Symbol' = 'Symbol',	'Description' = 'Description',	'NM_' = 'NM_',	'Accession' = 'Accession',	'Gene_ID' = 'Gene_ID',	'ENST_' = 'ENST_',	'ENSG_' = 'ENSG_', 'Unigene' = 'Unigene',	'XM_' = 'XM_',	'XR_' = 'XR_',	'NR_' = 'NR_',	'NP_' = 'NP_',	'XP_' = 'XP_',	'Nucleotide' = 'Nucleotide',	'Protein' = 'Protein',	'Unknown' = 'Unknown')
  
  return(col_check_)
}



load_GPL <- function(file_name, decimal, coltypes = NULL)
{
  file <- list()
  file$input <- readr::read_tsv(file = file_name, col_types = coltypes, locale = readr::locale(decimal_mark = decimal))
  return(file)
}



read_files_and_add_pubExp <- function (an_opts_pre_, gen_opts_pre_)
{
  file_list_ <- purrr::pmap(.l = list(an_opts_pre_$file_names, an_opts_pre_$pub_names, an_opts_pre_$exp_names),
                            .f = function(x, y, z) { 
                              temp <- readr::read_tsv(file = paste0(gen_opts_pre_$dir, '/', x), col_types = gen_opts_pre_$input_col_types, locale = readr::locale(decimal_mark = gen_opts_pre_$decimal))
                              temp[[gen_opts_pre_$check_cols$Pub.]] <- y 
                              temp[[gen_opts_pre_$check_cols$Exp.]] <- z
                              return(temp)})
  
  return(file_list_)
}



load_all_exps_from_paper <- function(dir_r_, pattern_, decimal_)
{
  files <- list.files(path = dir_r_, pattern = pattern_)
  
  table_list <- purrr::map(
    .x = files, 
    .f = function(x) {
      table <- readr::read_tsv(file = paste0(dir_r_, '/', x), locale = readr::locale(decimal_mark = decimal_))
      
      pub <- stringr::str_remove(string = x, pattern = '_.*')
      
      exp_ <- stringr::str_extract(string = x, pattern = '^[0-9]{1,3}_[0-9]{1,3}')
      exp_ <- stringr::str_remove(string = exp_, pattern = '.*_')
      
      tibble::add_column(.data = table, 'Pub.' = pub, 'Exp.' = exp_, .before = 1)
    })
  
  tables <- rlist::list.rbind(table_list)
  
  return(tables)
}



get_sep_all_identifiers <- function(input_df_name, bad_cols_to_process_seq, GPL_ID)
{
  temp_list <- list()
  
  if (GPL_ID == 'GPL16570') {
    temp_list <- lapply(
      X = bad_cols_to_process_seq,
      FUN = function(x) {
        temp <- stringr::str_replace_all(string = input_df_name[[x]], pattern = ' /// ', replacement = ' // ')
        
        as.data.frame(stringr::str_split_fixed(string = temp, pattern = ' // ', n = Inf))
      }
    )
  } else if (GPL_ID == 'GPL13912') {
    temp_list <- lapply(
      X = bad_cols_to_process_seq,
      FUN = function(x) {
        as.data.frame(stringr::str_split_fixed(string = input_df_name[[x]], pattern = '\\|', n = Inf))
      } 
    )
  } else if (GPL_ID == 'R') {
    temp_list <- lapply(
      X = bad_cols_to_process_seq,
      FUN = function(x) {
        as.data.frame(stringr::str_split_fixed(string = input_df_name[[x]], pattern = '___', n = Inf))
      } 
    )
  }
  
  sep_all_identifiers <- rlist::list.cbind(temp_list)
  
  return(sep_all_identifiers)
}




extract_identifers <- function(sep_all_identifiers_, search_by_, output_, opts_pre_check_cols, range_of_cols_with_IDs_start, range_of_cols_with_IDs_end)
{
  oc <- opts_pre_check_cols
  
  for (row_nb in seq(length(sep_all_identifiers_[[1]]))) {
    
    row <- unique(t(sep_all_identifiers_[row_nb,]))
    
    for(identifer in row)
    {
      if (eval(parse(text = search_by_$NM_))) {
        output_[row_nb,oc[['NM_']]] <- paste(output_[row_nb,oc[['NM_']]], stringr::str_remove(string = identifer, pattern = ' .*'), sep = ', ')
      } else if (eval(parse(text = search_by_$ENST_))){
        output_[row_nb,oc[['ENST_']]] <- paste(output_[row_nb,oc[['ENST_']]], stringr::str_remove(string = identifer, pattern = ' .*'), sep = ', ')
      } else if (eval(parse(text = search_by_$ENSG_))){
        output_[row_nb,oc[['ENSG_']]] <- paste(output_[row_nb,oc[['ENSG_']]], stringr::str_remove(string = identifer, pattern = ' .*'), sep = ', ')
      } else if (eval(parse(text = search_by_$Gene_ID))){
        output_[row_nb,oc[['Gene_ID']]] <- paste(output_[row_nb,oc[['Gene_ID']]], stringr::str_remove(string = identifer, pattern = ' .*'), sep = ', ')
      } else if (eval(parse(text = search_by_$Accession))){
        output_[row_nb,oc[['Accession']]] <- paste(output_[row_nb,oc[['Accession']]], stringr::str_remove(string = identifer, pattern = ' .*'), sep = ', ')
      } else if (eval(parse(text = search_by_$Unigene))){
        output_[row_nb,oc[['Unigene']]] <- paste(output_[row_nb,oc[['Unigene']]], identifer, sep = ', ')
      } else if (eval(parse(text = search_by_$XM_))){
        output_[row_nb,oc[['XM_']]] <- paste(output_[row_nb,oc[['XM_']]], stringr::str_remove(string = identifer, pattern = ' .*'), sep = ', ')
      } else if (eval(parse(text = search_by_$XR_))){
        output_[row_nb,oc[['XR_']]] <- paste(output_[row_nb,oc[['XR_']]], stringr::str_remove(string = identifer, pattern = ' .*'), sep = ', ')
      } else if (eval(parse(text = search_by_$NR_))){
        output_[row_nb,oc[['NR_']]] <- paste(output_[row_nb,oc[['NR_']]], stringr::str_remove(string = identifer, pattern = ' .*'), sep = ', ')
      } else if (eval(parse(text = search_by_$NP_))){
        output_[row_nb,oc[['NP_']]] <- paste(output_[row_nb,oc[['NP_']]], stringr::str_remove(string = identifer, pattern = ' .*'), sep = ', ')
      } else if (eval(parse(text = search_by_$XP_))){
        output_[row_nb,oc[['XP_']]] <- paste(output_[row_nb,oc[['XP_']]], stringr::str_remove(string = identifer, pattern = ' .*'), sep = ', ')
      }  else if (eval(parse(text = search_by_$LOC))){
        identifer_loc_removed <- stringr::str_remove(string = stringr::str_remove(string = identifer, pattern = ' .*'), pattern = 'LOC')
        output_[row_nb,oc[['Gene_ID']]] <- paste(output_[row_nb,oc[['Gene_ID']]], identifer_loc_removed, sep = ', ')
      } else if (eval(parse(text = search_by_$Symbol))){
        output_[row_nb,oc[['Symbol']]] <- paste(output_[row_nb,oc[['Symbol']]], identifer, sep = ', ')
      }
    }
  }
  
  output_[range_of_cols_with_IDs_start:range_of_cols_with_IDs_end] <- purrr::map_dfc(output_[range_of_cols_with_IDs_start:range_of_cols_with_IDs_end], .f = function(x){
    temp <- stringr::str_remove(string = x, pattern = 'NA, ')
    temp_1 <- stringr::str_remove(string = temp, pattern = ', $')
    temp_2 <- stringr::str_remove(string = temp_1, pattern = ', NA')
    temp_3 <- stringr::str_remove(string = temp_2, pattern = 'NA')
    temp_4 <- stringr::str_remove(string = temp_3, pattern = '^, ')
  })
  output_[output_ == ''] <- NA
  
  return(output_)
}



set_search_bys <- function(input_df_, search_by_p1_, search_by_n_, check_cols_)
{
  results_list <- list(
    'search_by_p1' = search_by_p1_, 
    'search_by_n' = search_by_n_)

  results_list$search_by[check_cols_$NM_] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^NM_')"))
  results_list$search_by[check_cols_$ENST_] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^ENS.*T')"))
  results_list$search_by[check_cols_$ENSG_] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^ENS.*G')"))
  results_list$search_by[check_cols_$Gene_ID] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^[0-9]*$')"))
  results_list$search_by[check_cols_$Accession] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^[A-Z]{1,2}[0-9]*$')"))
  results_list$search_by[check_cols_$Unigene] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^Mm(?![a-zA-Z]|[0-9]|[ ])')"))
  results_list$search_by[check_cols_$XM_] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^XM_')"))
  results_list$search_by[check_cols_$XR_] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^XR_')"))
  results_list$search_by[check_cols_$NR_] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^NR_')"))
  results_list$search_by[check_cols_$NP_] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^NP_')"))
  results_list$search_by[check_cols_$XP_] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^XP_')"))
  results_list$search_by['LOC'] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^LOC')"))
  results_list$search_by[check_cols_$Symbol] = list(paste0(results_list$search_by_p1, results_list$search_by_n, ", pattern = '^[a-zA-Z]{1,6}(?![a-zA-Z]|[ ]|[:])|.*Rik$')"))

  results_list$output <- input_df_
  
  return(results_list)
}





get_proper_NM_and_Accession_NR_XM_cols <- function(loaded_GPL = GPL13912)
{
  source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
  
  stopifnot(exprs = {
    exists(x = 'get_all_symbols_in_chrvec')
    is.na(any(stringr::str_detect(string = loaded_GPL$input$NM_, pattern = ',| |\\.|/|;'))) || !any(stringr::str_detect(string = loaded_GPL$input$NM_, pattern = ',| |\\.|/|;'))
    is.na(any(stringr::str_detect(string = loaded_GPL$input$Accession, pattern = ',| |\\.|/|;'))) || !any(stringr::str_detect(string = loaded_GPL$input$Accession, pattern = ',| |\\.|/|;'))
  })
  
  NM_logical <- stringr::str_detect(string = loaded_GPL$input$NM_, pattern = '^NM_', negate = T)
  Accession_logical <- stringr::str_detect(string = loaded_GPL$input$Accession, pattern = '^[A-Z]{1,2}[0-9]*$', negate = T)
  XM_logical <- stringr::str_detect(string = loaded_GPL$input$XM_, pattern = '^XM_', negate = T)
  NR_logical <- stringr::str_detect(string = loaded_GPL$input$NR_, pattern = '^NR_', negate = T)
  
  loaded_GPL$input$NM_[NM_logical] <- NA
  loaded_GPL$input$Accession[Accession_logical] <- NA
  loaded_GPL$input$XM_[XM_logical] <- NA
  loaded_GPL$input$NR_[NR_logical] <- NA
  
  return(loaded_GPL)
}



remove_non_symbols_from_GPL13912_wrapper <- function(GPL13912_output)
{
  temp1 <- stringr::str_remove_all(string = GPL13912_output[['Symbol']], pattern = '(ref, )|(ens, )|(gb, )|(tc, )|(riken, )|(nap, )')
  temp2 <- stringr::str_remove_all(string = temp1, pattern = '(ref)|(ens)|(gb)|(tc)|(riken)|(nap)')
  temp3 <- stringr::str_remove(string = temp2, pattern = '(, )$')
  temp4 <- stringr::str_remove(string = temp3, pattern = '(chr)[0-9]:.*')
  temp4[temp4 == ''] <- NA

  GPL13912_output[['Symbol']] <- temp4

  return(GPL13912_output)
}




remove_repeated_IDs_from_final_column <- remove_repeated_values_from_string_series_from_final_column 



set_df_for_append <- function(check_cols_, length)
{
  column_ <- rep(NA, length)
  
  df_ <- purrr::map_dfc(
    .x = check_cols_, 
    .f = function(x) { as.character(rep(NA, length)) })
  
  colnames(df_) <- as.character(check_cols_)
  
  return(df_)
}



copy_out_of_box_ready_columns_into_final_output_wrapper <- function(check_cols_, final_output_, input_)
{
  final_output_[[check_cols_$Pub.]] <- input_[[check_cols_$Pub.]]
  final_output_[[check_cols_$Exp.]] <- input_[[check_cols_$Exp.]]
  final_output_[[check_cols_$`adj,P,Val`]] <- input_$adj.P.Val
  final_output_[[check_cols_$`P,Val`]] <- input_$P.Value
  final_output_[[check_cols_$t]] <- input_$t
  final_output_[[check_cols_$B]] <- input_$B
  final_output_[[check_cols_$logFC]] <- input_$logFC
  final_output_[[check_cols_$ID]] <- input_$ID
  
  return(final_output_)
}



prepare_noncorrupted_symbol_col <- function(symbol_col_)
{
  symbol_col_ <- stringr::str_replace(string = symbol_col_, pattern = ", beta-carotene 15,15'-monooxygenase", replacement = '')
  symbol_col_ <- stringr::str_replace_all(string = symbol_col_, pattern = '‐', replacement = '-')
  symbol_col_[symbol_col_ == '---'] <- NA
  symbol_col_[symbol_col_ == '-'] <- NA
  symbol_col_[symbol_col_ == '4-Sep'] <- NA
  symbol_col_[symbol_col_ == 'Cyp3a23/3a1'] <- 'Cyp3a23, Cyp3a1'
  symbol_col_[symbol_col_ == 'Fam164a/Zc2hc1a'] <- 'Fam164a, Zc2hc1a'
  symbol_col_[symbol_col_ == 'X|X'] <- NA
  symbol_col_[symbol_col_ == 'Y|Y'] <- NA
  symbol_col_[symbol_col_ == 'Ggh ‡'] <- 'Ggh'
  symbol_col_[stringr::str_detect(string = symbol_col_, pattern = 'IGHG1_J00453\\$V00793_IG_HEAVY_CONSTANT_GAMMA_1_792')] <- 'IGHG1'
  symbol_col_ <- stringr::str_replace_all(string = symbol_col_, pattern = ' /// ', replacement = ', ')
  symbol_col_ <- stringr::str_replace_all(string = symbol_col_, pattern = '///', replacement = ', ')
  symbol_col_ <- stringr::str_replace(string = symbol_col_, pattern = 'Pcdha@$', replacement = 'Pcdha2')
  symbol_col_ <- stringr::str_remove_all(string = symbol_col_, pattern = '^NA, ')
  symbol_col_ <- stringr::str_remove(string = symbol_col_, pattern = ', NA$')
  symbol_col_[stringr::str_detect(string = symbol_col_, pattern = '^chr.{1,2}:[0-9]')] <- NA
  symbol_col_ <- stringr::str_remove(string = symbol_col_, pattern = ' $')
  symbol_col_ <- stringr::str_remove(string = symbol_col_, pattern = '\\*$')
  symbol_col_ <- stringr::str_remove(string = symbol_col_, pattern = ', X\\|X$|, Y\\|Y')
  
  
  return(symbol_col_)
}



prepare_noncorrupted_NXMR_col <- function(nxmr_col_)
{
  nxmr_col_ <- stringr::str_remove(string = nxmr_col_, pattern = ' $')
  
  nxmr_col_ <- stringr::str_remove_all(string = nxmr_col_, pattern = '\\.[0-9]{1,}')
  
  return(nxmr_col_)
}



get_input_QA_checklist <- function(input_table_, Pub._col_)
{
  input_QA_checklist_ <- list (
    'are_all_pub_present' = unique(input_table_[[Pub._col_]]),
    'are_all_cols_noncorrupted' = get_all_symbols_in_df_per_column(input_table_),
    'col_types' = purrr::map(.x = input_table_, class),
    'numeric_range' = purrr::map(.x = input_table_, .f = function(x) {ifelse(test = class(x) == 'numeric', yes = paste0(min(x[!is.na(x)]), ' - ', max(x[!is.na(x)])), no = NA)}),
    'number_of_NAs' = purrr::map(.x = input_table_, .f = function(x){length(x[is.na(x)])}),
    'percent_of_NAs' = purrr::map(.x = input_table_, .f = function(x){length(x[is.na(x)])/length(x)}))
  
  return(input_QA_checklist_)
}



prepare_uncorrupted_GPL10427_accessions_wrapper <- function(GPL10427_accessions)
{
  temp1 <- stringr::str_replace_all(string = GPL10427_accessions, pattern = '\\|', replacement = '___')
  temp2 <- stringr::str_remove_all(string = temp1, pattern = 'ug|ref')
  temp3 <- stringr::str_remove(string = temp2, pattern = '^___')
  temp4 <- stringr::str_replace_all(string = temp3, pattern = '_{4,}', replacement = '___')
  
  return(temp4)
}




prepare_noncorrupted_GPL16570_cols <- function(GPL16570_cols)
{
  GPL16570_cols[GPL16570_cols == '---'] <- NA
  GPL16570_cols <- stringr::str_replace_all(string = GPL16570_cols, pattern = ' /// ', replacement = '___')
  GPL16570_cols <- stringr::str_replace_all(string = GPL16570_cols, pattern = ' // ', replacement = '___')
  
  return(GPL16570_cols)
}




prepare_uncorrupted_GPL10427_accessions_wrapper <- function(GPL10427_accessions)
{
  temp1 <- stringr::str_replace_all(string = GPL10427_accessions, pattern = '\\|', replacement = '___')
  temp2 <- stringr::str_remove_all(string = temp1, pattern = 'ug|ref')
  temp3 <- stringr::str_remove(string = temp2, pattern = '^___')
  temp4 <- stringr::str_replace_all(string = temp3, pattern = '_{4,}', replacement = '___')
  
  return(temp4)
}

prepare_uncorrupted_GPL13912_ACCESSION_STRING_wrapper <- function(GPL13912_ACCESSION_STRING)
{
  temp1 <- stringr::str_replace_all(string = GPL13912_ACCESSION_STRING, pattern = '\\|', replacement = '___')
  temp2 <- stringr::str_remove_all(string = temp1, pattern = 'ref|ens|gb|tc|riken|nap')
  temp3 <- stringr::str_remove(string = temp2, pattern = '^___')
  temp4 <- stringr::str_replace_all(string = temp3, pattern = '_{4,}', replacement = '___')
  temp5 <- stringr::str_remove(string = temp4, pattern = '(chr)([0-9]{1,2}|x|y|X|Y):.*___')
  temp6 <- stringr::str_remove(string = temp5, pattern = '(chr)([0-9]{1,2}|x|y|X|Y):.*')
  temp7 <- stringr::str_remove(string = temp6, pattern = '(, )$')
  temp7[temp7 == ''] <- NA
  
  return(temp7)
}




get_list_of_group_pairs_for_deseq2_analysis <- function(cnt_group_name, exp_groups_name_vec, design_df, col_in_design_df_with_group_names, count_matrix, sample_name_col_in_design_df)
{
  warning('Beware!!! Name of each sample (column name) has to start with group name: e.g. control (group name) and control_1 (sample name)')
  
  analyses_to_do <- purrr::map(
    .x = exp_groups_name_vec, 
    .f = function(x) {
      temp <- list()
      
      temp$design <-  subset(x = design_df, subset = design_df[[col_in_design_df_with_group_names]] %in% c(cnt_group_name, x))
      
      temp$count_matrix <- subset(x = count_matrix, select = stringr::str_detect(string = colnames(count_matrix), pattern = paste0('^', cnt_group_name, '.*|^', x, '.*')))

      ifelse(
        test = all(temp$design[[sample_name_col_in_design_df]] == colnames(temp$count_matrix)),
        yes = print('...'),
        no = return('colnames in count_matrix differ from sample names in design_df'))

      return(temp)
    })
  
  names(analyses_to_do) <- exp_groups_name_vec
  
  return(analyses_to_do)
}




calculate_deseq2 <- function(count_matrix_, design_, col_with_comparisons_in_design, reference_name_in_design, filter_low_counts_int = NA, return_only_df = T)
{
  result_list <- list()
  #Create DESeqDataSet object
  result_list$deseq_data <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix_, 
    colData =  design_, 
    design = as.formula(paste0('~ ', col_with_comparisons_in_design))
  )
  
  #Filter low count rows
  if (!is.na(filter_low_counts_int)) {
    rows_to_keep <- rowSums(counts(result_list$deseq_data)) >= filter_low_counts_int
    result_list$deseq_data <- result_list$deseq_data[rows_to_keep,]
  }
  
  #Set reference to proper factor
  result_list$deseq_data[[col_with_comparisons_in_design]] <- relevel(result_list$deseq_data[[col_with_comparisons_in_design]], ref = reference_name_in_design)
  
  #Calculate results and get result table
  result_list$deseq_data_post_deseq <- DESeq2::DESeq(result_list$deseq_data)
  result_list$deseq_data_post_deseq_results <- DESeq2::results(result_list$deseq_data_post_deseq)
  result_list$deseq_data_post_deseq_results_df <- as.data.frame(result_list$deseq_data_post_deseq_results)
  
  if (return_only_df == T) {
    return(result_list$deseq_data_post_deseq_results_df)
  } else {
    return(result_list)
  }
  
  return('I am Error')
}



input_data_symbol_uncorruption_wrapper <- function(manual_input_symbol)
{
  return_ <- manual_input_symbol
  
  # vs '
  return_[return_ == "Bcmo1, beta-carotene 15,15'-monooxygenase"] <- 'Bcmo1'
  
  # vs @
  return_[return_ == 'Gm37013, Pcdhga12, Pcdhga11, Pcdhga10, Pcdhga9, Pcdhga8, Pcdhga7, Pcdhga6, Pcdhga5, Pcdhga4, Pcdhga3, Pcdhga2, Pcdhga1, Pcdhgc5, Pcdhgc4, Pcdhgc3, Pcdhgb8, Pcdhgb7, Pcdhgb6, Pcdhgb5, Pcdhgb4, Pcdhgb2, Pcdhgb1, Pcdha@']  <- 'Gm37013, Pcdhga12, Pcdhga11, Pcdhga10, Pcdhga9, Pcdhga8, Pcdhga7, Pcdhga6, Pcdhga5, Pcdhga4, Pcdhga3, Pcdhga2, Pcdhga1, Pcdhgc5, Pcdhgc4, Pcdhgc3, Pcdhgb8, Pcdhgb7, Pcdhgb6, Pcdhgb5, Pcdhgb4, Pcdhgb2, Pcdhgb1, Pcdha2'
  
  # vs +
  return_ <- stringr::str_remove_all(string = return_, pattern = '(Na\\+/K\\+)|(Na\\+/H\\+)|(ATPase, Ca\\+\\+)|(ATPase, H\\+)|(Ca2\\+-)|(ATPase, Cu\\+\\+)|(UDP-Gal:)|(UDP-N-).*')
  
  # vs |
  return_ <- stringr::str_remove_all(string = return_, pattern = '(X\\|X)|(Y\\|Y).*(, ){0,1}')
  return_ <- stringr::str_remove(string = return_, pattern = '(, XC3\\|X)|(, XE3\\|X)|(, XA1.1\\|X)')
  
  
  return_ <- stringr::str_replace_all(string = return_, pattern = 'Cyp3a23/3a1', replacement = 'Cyp3a, Cyp3a23')
  
  return_ <- stringr::str_replace_all(string = return_, pattern = '(, ){2,}', replacement = ', ')
  
  return_ <- stringr::str_remove(string = return_, pattern = ', $')
  
  return(return_)
}



#Works for mouse/rat (days), rat (weight), and embrionic (- days)
age_category_wrapper <- function(age_col, age_days_col, species_col, return_test)
{
  age <- tibble::tibble('age_col' = age_col, 'age_days_col' = age_days_col, 'species_col' = species_col, 'age_cat' = NA)
  
  for (row in seq_along(age_days_col)) 
  {
    if (
      !is.na(age_days_col[[row]]) &&
      age_days_col[[row]] < 0)
    {
      age$age_cat[[row]] <- 'fetus'
    } else if (
      is.na(age_days_col[[row]]) &&
      !is.na(species_col[[row]]) &&
      !is.na(age_col[[row]]) &&
      species_col[[row]] == 'rat' &&
      stringr::str_detect(string = age_col[[row]], pattern = '[0-9]{1,} g.*')) 
    {
      grams <- stringr::str_remove(string = age_col[[row]], pattern = ' g.*')
      
      age$age_cat[[row]] <- dplyr::case_when(
        dplyr::between(grams, 0, 49) ~ 'juvenile',
        dplyr::between(grams, 50, 320) ~ 'adolescent',
        dplyr::between(grams, 321, 480) ~ 'adult'
      )
    } else if (
      !is.na(age_days_col[[row]]) &&
      !is.na(species_col[[row]]) &&
      !is.na(age_col[[row]]) &&
      species_col[[row]] == 'rat')
    {
      days <- age_days_col[[row]]
      
      age$age_cat[[row]] <- dplyr::case_when(
        dplyr::between(days, 0, 20) ~ 'juvenile',
        dplyr::between(days, 21, 70) ~ 'adolescent',
        dplyr::between(days, 71, 600) ~ 'adult',
        days > 600 ~ 'old'
      )
    } else if (
      !is.na(age_days_col[[row]]) &&
      !is.na(species_col[[row]]) &&
      !is.na(age_col[[row]]) &&
      species_col[[row]] == 'mouse')
    {
      days <- age_days_col[[row]]
      
      age$age_cat[[row]] <- dplyr::case_when(
        dplyr::between(days, 0, 22) ~ 'juvenile',
        dplyr::between(days, 23, 60) ~ 'adolescent',
        dplyr::between(days, 61, 720) ~ 'adult',
        days > 721 ~ 'old'
      )
    } 
    
  }
  
  if (return_test == T) {
    return(age)
  } else return(age$age_cat)
  
}
