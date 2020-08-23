### Read in the table with description of experiments. They need to have the same identifiers as data in PRE_DATA, that is: Paper(int), Experiment(chr). Also needs Species(chr) column, if probe_id identification function or ncbi query tool is to be used
id_platf <- list(
  'desc' = readr::read_tsv('descriptions_for_analysis.tsv', col_types = 'nnccccnccnccccncccccccnnccncc', locale = readr::locale(decimal_mark = ',')),
  'input' = readr::read_tsv('input_data_for_analysis.tsv', col_types = 'nncnnnnnccccccccccccc', locale = readr::locale(decimal_mark = ',')))

id_platf$sample_input <- dplyr::sample_n(tbl = id_platf$input, size = 2500)


### Prepare data for probe analysis, and for not  probe analysis
id_platf$input_pubs_with_probes <- unique(subset(x = id_platf$sample_input, subset = !is.na(id_platf$sample_input$Probe))[['Pub.']])

id_platf$pubs_for_probe_analysis$input <- subset(x = id_platf$sample_input, subset = id_platf$sample_input$`Pub.` %in% id_platf$input_pubs_with_probes)

id_platf$pubs_NOT_for_probe_analysis$input <- subset(x = id_platf$sample_input, subset = id_platf$sample_input$`Pub.` %nin% id_platf$input_pubs_with_probes)

length(id_platf$sample_input[[1]]) == length(id_platf$pubs_for_probe_analysis$input[[1]]) + length(id_platf$pubs_NOT_for_probe_analysis$input[[1]])
### Prepare data for probe analysis, and not for probe analysis

id_platf$pubs_for_probe_analysis$highest_hit_analysis <- master_annotator(
  descriptions_df = id_platf$desc, 
  des_paper_id_col = 'Pub.', 
  des_species_col = 'Species', 
  input_df = id_platf$pubs_for_probe_analysis$input, 
  input_paper_col = 'Pub.', 
  input_id_col = 'Probe',
  PERFORM_B_get_the_highest_hit_returning_id_type = T,
  B_highest_hit_int_Probe_IDs_to_test = 50)

### !!! Add saving all the highest_hit dataset?


data_prepared_for_annotation <- list(
  'data_for_probe_annotation' = id_platf$pubs_for_probe_analysis$input,
  'platform_ids_for_probe_annotation' = id_platf$pubs_for_probe_analysis$highest_hit_analysis$ids_to_be_used_for_annotation,
  'data_NOT_for_probe_annotation' = id_platf$pubs_NOT_for_probe_analysis$input,
  'descriptions' = id_platf$desc
  )

save(data_prepared_for_annotation, file = 'data_prepared_for_annotation')

