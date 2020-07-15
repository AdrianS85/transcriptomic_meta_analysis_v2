source('functions_pre_preparation.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')



descriptions <- list('input' = openxlsx::read.xlsx(xlsxFile = 'Antidepressants_metadata_140520.xlsx', sheet = 'pubs'))

descriptions$input_qa <- verify_df(df_ = descriptions$input, sort_by_col = 'Pub.')



### Prepare duration columns ###
descriptions$input_qa$df$Duration[descriptions$input_qa$df$Duration == 'single dose'] <- '1 day'

descriptions$input_qa$df$Duration_d <- cleanup_differing_units(charvec_ = descriptions$input_qa$df$Duration, unit_identifier_regexList = list('h' = '[0-9]{1,} h.*', 'd' = '[0-9]{1,} d.*', 'w' = '[0-9]{1,} w.*', 'mo' = '[0-9]{1,} mo.*'), unit_multiplier_list = list('h' = 1/24, 'd' = 1, 'w' = 7, 'm' = 30))[[1]]

descriptions$input_qa$df$Duration_cat <- dplyr::case_when(
  descriptions$input_qa$df$Duration_d <= 1 ~ 'acute',
  descriptions$input_qa$df$Duration_d > 1 & descriptions$input_qa$df$Duration_d <= 7 ~ 'medium',
  descriptions$input_qa$df$Duration_d > 7 ~ 'prolonged'
)
### Prepare duration columns ###



### Prepare latency columns ###
descriptions$input_qa$df$Latency[descriptions$input_qa$df$Latency == 'immediately'] <- '0 h'

descriptions$input_qa$df$Latency_h <- cleanup_differing_units(charvec_ = descriptions$input_qa$df$Latency, unit_identifier_regexList = list('mn' = '[0-9]{1,} mi.*', 'h' = '[0-9]{1,} h.*', 'd' = '[0-9]{1,} d.*', 'w' = '[0-9]{1,} w.*'), unit_multiplier_list = list('mn' = 1/60, 'h' = 1, 'd' = 24, 'w' = 168))[[1]]

descriptions$input_qa$df$Latency_cat <- dplyr::case_when(
  descriptions$input_qa$df$Latency_h <= 1 ~ 'immediate',
  descriptions$input_qa$df$Latency_h > 1 & descriptions$input_qa$df$Latency_h <= 48 ~ 'medium',
  descriptions$input_qa$df$Latency_h > 49 ~ 'prolonged'
)
### Prepare latency columns ###




### Prepare age columns ###
descriptions$input_qa$df$Age_d <- cleanup_differing_units(charvec_ = descriptions$input_qa$df$Age, unit_identifier_regexList = list('ed' = '[0-9]{1,} ed$', 'd' = '[0-9]{1,} d.*', 'w' = '[0-9]{1,} w.*', 'mo' = '[0-9]{1,} mo.*'), unit_multiplier_list = list('ed' = -1, 'd' = 1, 'w' = 7, 'm' = 30))[[1]]

descriptions$input_qa$df$Age_cat <- age_category_wrapper(
  age_col = descriptions$input_qa$df$Age, 
  age_days_col = descriptions$input_qa$df$Age_d, 
  species_col = descriptions$input_qa$df$Species,
  return_test = F)

# Single mouse studied gave grams insteead of age
temp_mus_weight_age_40 <- descriptions$input_qa$df$Species == 'mouse' & descriptions$input_qa$df$Age == '40 g' & !is.na(descriptions$input_qa$df$Age)

descriptions$input_qa$df$Age_cat[temp_mus_weight_age_40] <- 'adult'
### Prepare age columns ###



descriptions$output <- dplyr::select(descriptions$input_qa$df, 1:6, 27:28, 7, 29:30, 8:10, 31:32, 11:14, 16:23, 25)

descriptions$output_qa <- verify_df(df_ = descriptions$output, sort_by_col = 'Pub.', only_qa = T)

write.table(x = descriptions$output, file = 'descriptions_for_analysis.tsv', sep = '\t', dec = ',', row.names = F)

# 1 um/l in dose
