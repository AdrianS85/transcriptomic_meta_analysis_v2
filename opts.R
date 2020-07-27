set_col_check <- function()
{
  col_check_ <- list('Exp.' = 'Exp.', 'Comp.' = 'Comp',	'Pub.' = 'Pub.',	'ID' = 'ID',	'adj,P,Val' = 'adj,P,Val',	'P,Value' = 'P,Value',	't' = 't',	'B' = 'B',	'logFC' = 'logFC', 'Probe' = 'Probe',	'Symbol' = 'Symbol',	'Description' = 'Description',	'NM_' = 'NM_',	'Accession' = 'Accession',	'Gene_ID' = 'Gene_ID',	'ENST_' = 'ENST_',	'ENSG_' = 'ENSG_', 'Unigene' = 'Unigene',	'XM_' = 'XM_',	'XR_' = 'XR_',	'NR_' = 'NR_',	'NP_' = 'NP_',	'XP_' = 'XP_',	'Nucleotide' = 'Nucleotide',	'Protein' = 'Protein',	'Unknown' = 'Unknown')
  
  return(col_check_)
}


opts <- list(
  'decimal' = ',',
  'check_cols_for_input' = set_col_check(),
  'dir_r_downloaded_data' = 'geo_ma',
  'dir_seq_input' = 'geo_seq',
  'geo_download_p' = 0.05,
  'pre_prep_reform' = NA,
  'cons_dat' = NA,
  'ann' = NA)

opts$pre_prep_reform <- list(
  'range_of_cols_with_IDs_in_output_start' = 9,
  'range_of_cols_with_IDs_in_output_end' = 25,
  'input_col_types' = 'cnnnnncccccccccccccccccccccccccccccccc',
  'folder' = 'pre_preparation_reformatting'
)

opts$cons_dat <- list(
  'select_prepared_data_regex' = '.*_prepared.*',
  'input_col_types' = 'nncnnnnnccccccccccccccccc',
  'folder' = 'construct_dataset'
)

opts$ann <- list('des_exp_id_col' = 'Exp.',
                'des_species_col' = 'Species',
                'input_exp_id_col' = 'Exp.',
                'des_platform_col' = 'Assay',
                'folder' = 'annotation',
                'splitting_message_ok' = 'Splitting of %s into %s and %s was ok.',
                'splitting_message_bad' = 'Annotated data do not have the same number of rows as finalized and leftover data. Splitting failed!'
)


dir.create(opts[['dir_r_downloaded_data']])
dir.create(opts[['dir_seq_input']])


source('functions_for_annotation.R')
source('functions_GEO_download.R')
source('functions_pre_preparation.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/bioinfo_little_helpers.R')
