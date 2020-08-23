library(loggit)
default::default(message) <- list(echo = F)
# rm(list = ls(pattern = 'temp.*|test.*'))

`%>%` <- dplyr::`%>%`
`%nin%` <- Hmisc::`%nin%`

source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')
source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/bioinfo_little_helpers.R')
source('functions_GEO_download.R')
source('functions_pre_preparation.R')
source('functions_for_annotation.R')
source('functions_for_data_prep.R')



set_col_check <- function()
{
  col_check_ <- list('Exp.' = 'Exp.', 'Comp.' = 'Comp',	'ID' = 'ID',	'adj,P,Val' = 'adj,P,Val',	'P,Value' = 'P,Value',	't' = 't',	'B' = 'B',	'logFC' = 'logFC', 'Probe' = 'Probe',	'Symbol' = 'Symbol',	'Description' = 'Description',	'NM_' = 'NM_',	'Accession' = 'Accession',	'Gene_ID' = 'Gene_ID',	'ENST_' = 'ENST_',	'ENSG_' = 'ENSG_', 'Unigene' = 'Unigene',	'XM_' = 'XM_',	'XR_' = 'XR_',	'NR_' = 'NR_',	'NP_' = 'NP_',	'XP_' = 'XP_',	'Nucleotide' = 'Nucleotide',	'Protein' = 'Protein',	'Unknown' = 'Unknown')
  
  return(col_check_)
}





opts <- list(
  'decimal' = ',',
  'check_cols_for_input' = set_col_check(),
  'dir_r_downloaded_data' = 'geo_ma',
  'dir_seq_input' = 'geo_seq',
  'dir_descriptions' = 'desc',
  'geo_download_p' = 0.05,
  'pre_prep_reform' = NA,
  'cons_dat' = NA,
  'ann' = NA,
  'str_sep' = ', ',
  'ext_gene_name_col' = 'external_gene_name')





opts$pre_prep_reform <- list(
  'range_of_cols_with_IDs_in_output_start' = 9,
  'range_of_cols_with_IDs_in_output_end' = 25,
  'input_col_types' = 'cnnnnncccccccccccccccccccccccccccccccc',
  'folder' = 'pre_preparation_reformatting',
  'bind_sep' = '___',
  'which_platform_to_get_exp_and_comp_nb_for_regex' = '^GPL.*',
  'exp_names_extraction_regex' = '^[[0-9]]{1,3}',
  'comp_names_extraction_regex' = '^[[0-9]]{1,3}_[[0-9]]{1,3}',
  'exp_comp_divider_symbol' = '_',
  'search_bys' = set_search_bys(search_by_n_ = 'identifer', check_cols_ = opts$check_cols_for_input),
  'GPL1261' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.*GPL1261.tsv')),
  'GPL5425' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.*GPL5425.tsv')),
  'GPL6246' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.GPL6246.tsv')),
  'GPL6885' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.*GPL6885.tsv')),
  'GPL6887' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.*GPL6887.tsv')),
  'GPL8160' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.*GPL8160.tsv')),
  'GPL10427' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.*GPL10427.tsv')),
  'GPL13912' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.*GPL13912.tsv')),
  'GPL16570' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.*GPL16570.tsv')),
  'GPL17223' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.*GPL17223.tsv')),
  'GPL25480' = list('file_names'  = list.files(path = opts$dir_r_downloaded_data, pattern = '.*GPL25480.tsv'))
)

opts$pre_prep_reform[['platform_to_get_exp_and_comp_nb_for']] <- names(opts$pre_prep_reform[stringr::str_detect(string = names(opts$pre_prep_reform), pattern = opts$pre_prep_reform[['which_platform_to_get_exp_and_comp_nb_for_regex']])]) # This is only for downloaded experiments

for (platform in opts$pre_prep_reform[['platform_to_get_exp_and_comp_nb_for']]) 
{
  opts$pre_prep_reform[[platform]][['exp_names']] <- stringr::str_extract(
    string = opts$pre_prep_reform[[platform]][['file_names']], 
    pattern = opts$pre_prep_reform[['exp_names_extraction_regex']])
  
  opts$pre_prep_reform[[platform]][['comp_names']] <- stringr::str_remove(
    string = stringr::str_extract(
      string = opts$pre_prep_reform[[platform]][['file_names']],
      pattern = opts$pre_prep_reform[['comp_names_extraction_regex']]),
    pattern = paste0('.*', opts$pre_prep_reform[['exp_comp_divider_symbol']]))
}





opts$ann <- list('exp_id_col' = opts$check_cols_for_input[['Exp.']],
                'des_species_col' = 'Species',
                'des_platform_col' = 'Assay',
                'folder' = 'annotation',
                'splitting_message_ok' = 'Splitting of %s into %s and %s was ok.',
                'splitting_message_bad' = 'Annotated data do not have the same number of rows as finalized and leftover data. Splitting failed!',
                'cutoff_for_logfc' = 0.05
)





opts$data <- list('folder' = 'data',
                  'folder_gene_groups' = 'gene_groups',
                  'folder_clusters' = 'clusters',
                  'exp_and_comp_col' = 'EC_ID',
                  'logfc_median_col' = 'logFC_median'
)




dir.create(opts[['dir_r_downloaded_data']])
dir.create(opts[['dir_seq_input']])
dir.create(opts$pre_prep_reform[['folder']])
dir.create(opts$ann[['folder']])
dir.create(opts$data[['folder']])
dir.create(paste0(
  opts$data[['folder']], 
  '/', 
  opts$data[['folder_gene_groups']]))








