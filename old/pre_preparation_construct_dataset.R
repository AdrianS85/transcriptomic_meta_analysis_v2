# rm(list = ls(pattern = 'temp.*|test.*'))
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
source('functions_pre_preparation.R')
source('opts.R')
library(loggit)

dir.create(opts$cons_da$folder)
set_logfile(paste0(opts$cons_da$folder, "/construct_dataset.log"))


input_data <- list(
  'pre-input' = purrr::map(
    .x = list.files(path = opts$dir_r_downloaded_data, pattern = opts$cons_dat$select_prepared_data_regex), 
    .f = function(x){
      readr::read_tsv(
        file = paste0(opts$dir_r_downloaded_data, '/', x), 
        col_types = opts$cons_dat$input_col_types, 
        locale = readr::locale(decimal_mark = opts$decimal))
    })
)

are_vectors_the_same(chr_vec_list = purrr::map(
  .x = input_data[['pre-input']],
  .f =  colnames))

input_data$input <- rlist::list.rbind(.data = input_data[['pre-input']])

input_data$QA_checklist_input <- get_input_QA_checklist(
  input_table_ = input_data[['input']], 
  Pub._col_ = 'Pub.')

#Remove columns without single value
input_data$input$NP_ <- NULL
input_data$input$XP_ <- NULL
input_data$input$Protein <- NULL
input_data$input$Unknown <- NULL

# Checked - ID, Probe, NR_, Nucleotide, Unigene, ENSG_, ENST_, Gene_ID, Accession
input_data$input$NM_ <- stringr::str_remove_all(string = input_data$input$NM_, pattern = '\\.[0-9]{1,}')
input_data$input$XM_ <- stringr::str_remove_all(string = input_data$input$XM_, pattern = '\\.[0-9]{1,}')
input_data$input$XR_ <- stringr::str_remove_all(string = input_data$input$XR_, pattern = '\\.[0-9]{1,}')

input_data$input$Symbol <- input_data_symbol_uncorruption_wrapper(manual_input_symbol = input_data$input$Symbol)

input_data$input_reformated <- purrr::map_dfc(
  .x = input_data$input, 
  .f = function(x){
    if (class(x) == 'character') {
      remove_corrupting_symbols_from_chrvec(
        chr_vec = x, 
        repeated_spaces = T, 
        trailing_spaces = T, 
        character_NAs = T, 
        change_to_lower = F)
    } else x
})

input_data$input_reformated <- input_data$input_reformated[order(input_data$input_reformated[['Pub.']], input_data$input_reformated[['Exp.']]),]

write.table(
  x = input_data$input_reformated, 
  file = 'input_data_for_analysis.tsv', 
  sep = '\t', 
  dec = ',', 
  row.names = F)







# test <- input_data$input_reformated
# 
# test <- test[order(test[['Pub.']], test[['Exp.']]),]
# 
# test_probes_yes <- subset(x = test, subset = !is.na(test$Probe), select = c('Probe', 'Pub.', 'Exp.'))
# 
# test_probes_no <- subset(x = test, subset = is.na(test$Probe), select = c('Probe', 'Pub.', 'Exp.'))
# 
# length(test_probes_yes[[1]]) + length(test_probes_no[[1]]) == length(test[[1]])
# 
# test_any_probe <- purrr::map_lgl(.x = unique(test_probes_yes$Pub.), .f = function(x){ x %in% unique(test_probes_no$Pub.)})





