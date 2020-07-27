source('opts.R')
# rm(list = ls(pattern = 'temp.*|test.*'))

# opts_comparisons_76 <- prepare_comparisons_for_exp_GSE119291_wrapper()
# checked - means 1) length of dataset and 2 random genes logFC vs geo2r script for given dataset selected during different batch check-up, 2) list values vs meta-data descriptions


#########################
### PREPARE PLATFORMS ###
#########################
opts_geo_dl <- list(
  'GPL1261' = list(
    'platform' = "GPL1261",
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
    'AnnotGPL' = T),
  'GPL5425' = list(
    'platform' = 'GPL5425',
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
    'AnnotGPL' = T),
  'GPL6246' = list(
    'platform' = "GPL6246",
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
    'AnnotGPL' = T),
  'GPL6885' = list(
    'platform' = "GPL6885",
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","GI","GenBank.Accession"),
    'AnnotGPL' = T),
  'GPL6887' = list(
    'platform' = "GPL6887",
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","UniGene.title","UniGene.symbol","UniGene.ID","Nucleotide.Title","GI","GenBank.Accession"),
    'AnnotGPL' = T),
  'GPL8160' = list(
    'platform' = "GPL8160",
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC"),
    'AnnotGPL' = F),
  'GPL10427' = list(
    'platform' = "GPL10427",
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","accessions","ProbeUID","ProbeName","GeneName","SystematicName","Description","GB_LIST"),
    'AnnotGPL' = F),
  'GPL13912' = list(
    'platform' = "GPL13912",
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","NAME","GB_ACC","GENE_ID","GENE_SYMBOL","GENE_NAME","UNIGENE_ID","ENSEMBL_ID","ACCESSION_STRING","DESCRIPTION"),
    'AnnotGPL' = F),
  'GPL16570' = list(
    'platform' = "GPL16570",
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","probeset_id","gene_assignment","mrna_assignment","swissprot","unigene"),
    'AnnotGPL' = F),
  "GPL17223" = list(
    'platform' = "GPL17223",
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","Search_key","ProbeId","Gid","Transcript","GB_ACC","Symbol","Definition"),
    'AnnotGPL' = F),
  'GPL25480' = list(
    'platform' = 'GPL25480',
    'col_names' = c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID"),
    'AnnotGPL' = F)
)
#########################
### PREPARE PLATFORMS ###
#########################



###########################
### PREPARE COMPARISONS ###
###########################
opts_geo_dl[['GPL1261']][['comps']] <- list(
  '1' = list(
    'Pub.' = 70,
    'Exp.' = 37,
    'Comp.' = 1,
    'sample_comp' = 'S100a10_bacTRAP',
    'series' = "GSE35761",
    'group_design' = paste0('101010')), #checked
  '2' = list(
    'Pub.' = 70,
    'Exp.' = 37,
    'Comp.' = 2,
    'sample_comp' = 'Glt25d2_bacTRAP',
    'series' = "GSE35763",
    'group_design' = paste0('101010')), #checked
  '3' = list(
    'Pub.' = 70,
    'Exp.' = 37,
    'Comp.' = 3,
    'sample_comp' = 'p11_KO_S100a10_bacTRAP',
    'series' = "GSE35765",
    'group_design' = paste0('110010')), #checked
  '4' = list(
    'Pub.' = 77,
    'Exp.' = 43,
    'Comp.' = 1,
    'sample_comp' = '30_mg',
    'series' = "GSE63469",
    'group_design' = paste0('011X0X')), #checked
  '5' = list(
    'Pub.' = 77,
    'Exp.' = 43,
    'Comp.' = 2,
    'sample_comp' = '100_mg',
    'series' = "GSE63469",
    'group_design' = paste0('0XX101')), #checked
  '6' = list(
    'Pub.' = 78,
    'Exp.' = 44,
    'Comp.' = 1,
    'sample_comp' = 'cortex',
    'series' = "GSE118668",
    'group_design' = paste0('1111111100000000')), #checked
  '7' = list(
    'Pub.' = 78,
    'Exp.' = 44,
    'Comp.' = 2,
    'sample_comp' = 'hippocampus',
    'series' = "GSE118669",
    'group_design' = paste0('1111111100000000')), #checked
  '8' = list(
    'Pub.' = 66,
    'Exp.' = 34,
    'Comp.' = 1,
    'sample_comp' = 'none',
    'series' = "GSE6476",
    'group_design' = paste0('0011')) #checked
)





opts_geo_dl[['GPL5425']][['comps']] <- list(
  '1' = list(
    'Pub.' = 71,
    'Exp.' = 38,
    'Comp.' = 2,
    'sample_comp' = 'paroxetine_1d',
    'series' = "GSE59895",
    'group_design' = GSE59895_group_design()[['paroxetine_1d']]), #checked
  '2' = list(
    'Pub.' = 71,
    'Exp.' = 38,
    'Comp.' = 3,
    'sample_comp' = 'paroxetine_3d',
    'series' = "GSE59895",
    'group_design' = GSE59895_group_design()[['paroxetine_3d']]), #checked
  '3' = list(
    'Pub.' = 71,
    'Exp.' = 38,
    'Comp.' = 4,
    'sample_comp' = 'paroxetine_5d',
    'series' = "GSE59895",
    'group_design' = GSE59895_group_design()[['paroxetine_5d']]), #checked
  '4' = list(
    'Pub.' = 71,
    'Exp.' = 38,
    'Comp.' = 6,
    'sample_comp' = 'phenelzine_1d',
    'series' = "GSE59895",
    'group_design' = GSE59895_group_design()[['phenelzine_1d']]), #checked
  '5' = list(
    'Pub.' = 71,
    'Exp.' = 38,
    'Comp.' = 7,
    'sample_comp' = 'phenelzine_3d',
    'series' = "GSE59895",
    'group_design' = GSE59895_group_design()[['phenelzine_3d']]), #checked
  '6' = list(
    'Pub.' = 71,
    'Exp.' = 38,
    'Comp.' = 8,
    'sample_comp' = 'phenelzine_5d',
    'series' = "GSE59895",
    'group_design' = GSE59895_group_design()[['phenelzine_5d']]), #checked
  '7' = list(
    'Pub.' = 71,
    'Exp.' = 38,
    'Comp.' = 14,
    'sample_comp' = 'tramadol_1d',
    'series' = "GSE59895",
    'group_design' = GSE59895_group_design()[['tramadol_1d']]), #checked
  '8' = list(
    'Pub.' = 71,
    'Exp.' = 38,
    'Comp.' = 15,
    'sample_comp' = 'tramadol_3d',
    'series' = "GSE59895",
    'group_design' = GSE59895_group_design()[['tramadol_3d']]), #checked
  '9' = list(
    'Pub.' = 71,
    'Exp.' = 38,
    'Comp.' = 16,
    'sample_comp' = 'tramadol_5d',
    'series' = "GSE59895",
    'group_design' = GSE59895_group_design()[['tramadol_5d']]) #checked
)





opts_geo_dl[['GPL6246']][['comps']] <- list(
  '1' = list(
    'Pub.' = 67,
    'Exp.' = 35,
    'Comp.' = 1,
    'sample_comp' = 'none',
    'series' = "GSE26364",
    'group_design' = paste0('0011')) #checked
)





opts_geo_dl[['GPL6885']][['comps']] <- list(
  '1' = list(
    'Pub.' = 68,
    'Exp.' = 36,
    'Comp.' = 1,
    'sample_comp' = 'cortex',
    'series' = "GSE26836",
    'group_design' = paste0('XXXXXX111000')), #checked
  '2' = list(
    'Pub.' = 68,
    'Exp.' = 36,
    'Comp.' = 2,
    'sample_comp' = 'hippocampus',
    'series' = "GSE26836",
    'group_design' = paste0('000111XXXXXX')), #checked
  '3' = list(
    'Pub.' = 28,
    'Exp.' = 18,
    'Comp.' = 1,
    'sample_comp' = 'HA',
    'series' = "GSE27532",
    'group_design' = paste0('XXXXXXXX11110000')), #checked
  '4' = list(
    'Pub.' = 28,
    'Exp.' = 18,
    'Comp.' = 2,
    'sample_comp' = 'LA',
    'series' = "GSE27532",
    'group_design' = paste0('11110000XXXXXXXX')) #checked
)




opts_geo_dl[['GPL6887']][['comps']] <- list(
  '1' = list(
    'Pub.' = 72,
    'Exp.' = 39,
    'Comp.' = 1,
    'sample_comp' = 'hippocampus_1h',
    'series' = "GSE73798",
    'group_design' = paste0('1XXXX0XXXXXX1XXXXXXXXXXXXXXX0XXXXXXXXXXXXXXXXX0X1XXXXXXXXXXX')), #checked
  '2' = list(
    'Pub.' = 72,
    'Exp.' = 39,
    'Comp.' = 2,
    'sample_comp' = 'hippocampus_2h',
    'series' = "GSE73798",
    'group_design' = paste0('XXXXXXXXXXXXXXXXX01XXXX0XXXXXXXXXXXX1XXXXXXXXXXXXXXXX01XXXXX')), #checked
  '3' = list(
    'Pub.' = 72,
    'Exp.' = 39,
    'Comp.' = 3,
    'sample_comp' = 'hippocampus_4h',
    'series' = "GSE73798",
    'group_design' = paste0('XXXXXX1XXXXXXXXXXXXXXXXXXXXXXX1XXXX0XXXXX01XXXXXXXXXXXXXXXX0')), #checked
  '4' = list(
    'Pub.' = 72,
    'Exp.' = 39,
    'Comp.' = 4,
    'sample_comp' = 'hippocampus_8h',
    'series' = "GSE73798",
    'group_design' = paste0('X1XXXXXXXXX0XXXXXXXXXXXX1XXXX0XXXXXXXXXXXXXXXXX0X1XXXXXXXXXX')), #checked
  '5' = list(
    'Pub.' = 72,
    'Exp.' = 39,
    'Comp.' = 5,
    'sample_comp' = 'striatum_1h',
    'series' = "GSE73799",
    'group_design' = paste0('XXXXXX1XXXX0XXXX0X1XXXXXXXXXXXXXXXXXXXXX0X1XXXXXXXXXXXXXXXXX')), #checked
  '6' = list(
    'Pub.' = 72,
    'Exp.' = 39,
    'Comp.' = 6,
    'sample_comp' = 'striatum_2h',
    'series' = "GSE73799",
    'group_design' = paste0('XXXXXXXXXXXXXXXXXXXXXXX01XXXX01XXXXXXXXXXXXXXXX01XXXXXXXXXXX')), #checked
  '7' = list(
    'Pub.' = 72,
    'Exp.' = 39,
    'Comp.' = 7,
    'sample_comp' = 'striatum_4h',
    'series' = "GSE73799",
    'group_design' = paste0('1XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX01XXXXXXXXXXXXXXXX01XXXX0')), #checked
  '8' = list(
    'Pub.' = 72,
    'Exp.' = 39,
    'Comp.' = 8,
    'sample_comp' = 'striatum_8h',
    'series' = "GSE73799",
    'group_design' = paste0('XXXXX0X1XXXX1XXXX0XXXXXXXXXXXXXXXXXXXXXXX0X1XXXXXXXXXXXXXXXX')), #checked
  '9' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 1,
    'sample_comp' = 'mianserin_1h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['mianserin_1h']]),  #checked
  '10' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 2,
    'sample_comp' = 'imipramine_1h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['imipramine_1h']]), #checked
  '11' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 3,
    'sample_comp' = 'fluoxetine_1h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['fluoxetine_1h']]), #checked
  '12' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 4,
    'sample_comp' = 'bupropion_1h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['bupropion_1h']]), #checked
  '13' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 5,
    'sample_comp' = 'tianeptine_1h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['tianeptine_1h']]), #checked
  '14' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 6,
    'sample_comp' = 'tranylcypromine_1h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['tranylcypromine_1h']]), #checked
  '15' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 7,
    'sample_comp' = 'mianserin_2h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['mianserin_2h']]), #checked
  '16' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 8,
    'sample_comp' = 'imipramine_2h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['imipramine_2h']]), #checked
  '17' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 9,
    'sample_comp' = 'fluoxetine_2h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['fluoxetine_2h']]), #checked
  '18' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 10,
    'sample_comp' = 'bupropion_2h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['bupropion_2h']]), #checked
  '19' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 11,
    'sample_comp' = 'tianeptine_2h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['tianeptine_2h']]), #checked
  '20' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 12,
    'sample_comp' = 'tranylcypromine_2h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['tranylcypromine_2h']]), #checked
  '21' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 13,
    'sample_comp' = 'mianserin_4h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['mianserin_4h']]), #checked
  '22' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 14,
    'sample_comp' = 'imipramine_4h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['imipramine_4h']]), #checked
  '23' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 15,
    'sample_comp' = 'fluoxetine_4h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['fluoxetine_4h']]), #checked
  '24' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 16,
    'sample_comp' = 'bupropion_4h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['bupropion_4h']]), #checked
  '25' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 17,
    'sample_comp' = 'tianeptine_4h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['tianeptine_4h']]), #checked
  '26' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 18,
    'sample_comp' = 'tranylcypromine_4h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['tranylcypromine_4h']]), #checked
  '27' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 19,
    'sample_comp' = 'mianserin_8h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['mianserin_8h']]), #checked
  '28' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 20,
    'sample_comp' = 'imipramine_8h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['imipramine_8h']]), #checked
  '29' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 21,
    'sample_comp' = 'fluoxetine_8h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['fluoxetine_8h']]), #checked
  '30' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 22,
    'sample_comp' = 'bupropion_8h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['bupropion_8h']]), #checked
  '31' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 23,
    'sample_comp' = 'tianeptine_8h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['tianeptine_8h']]), #checked
  '32' = list(
    'Pub.' = 33,
    'Exp.' = 20,
    'Comp.' = 24,
    'sample_comp' = 'tranylcypromine_8h',
    'series' = "GSE48955",
    'group_design' = GSE48955_group_design()[['tranylcypromine_8h']]) #checked
)





opts_geo_dl[['GPL8160']][['comps']] <- list(
  '1' = list(
    'Pub.' = 34,
    'Exp.' = 21,
    'Comp.' = 1,
    'sample_comp' = 'none',
    'series' = "GSE42940",
    'group_design' = paste0('00001111')) #checked
)





opts_geo_dl[['GPL10427']][['comps']] <- list(
  '1' = list(
    'Pub.' = 60,
    'Exp.' = 32,
    'Comp.' = 1,
    'sample_comp' = 'hippocampus',
    'series' = "GSE97916",
    'group_design' = paste0('0001111XXXXXXXXX')), #checked
  '2' = list(
    'Pub.' = 60,
    'Exp.' = 32,
    'Comp.' = 2,
    'sample_comp' = 'hypothalamus',
    'series' = "GSE82016",
    'group_design' = paste0('0001111XXXXXXXXX')),#checked
  '3' = list(
    'Pub.' = 45,
    'Exp.' = 26,
    'Comp.' = 1,
    'sample_comp' = 'none',
    'series' = "GSE76110",
    'group_design' = paste0('00011111XXXXXXXX')) #checked
)





opts_geo_dl[['GPL13912']][['comps']] <- list(
  '1' = list(
    'Pub.' = 52,
    'Exp.' = 29,
    'Comp.' = 1,
    'sample_comp' = 'hippocampus',
    'series' = "GSE84183",
    'group_design' = paste0('00000000XXXXXXXXXXXXXXXX11111111XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')), #checked
  '2' = list(
    'Pub.' = 52,
    'Exp.' = 29,
    'Comp.' = 2,
    'sample_comp' = 'cortex',
    'series' = "GSE84183",
    'group_design' = paste0('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX00000000XXXXXXXXXXXXXXXX11111111')) #checked
)





opts_geo_dl[['GPL16570']][['comps']] <- list(
  '1' = list(
    'Pub.' = 56,
    'Exp.' = 30,
    'Comp.' = 1,
    'sample_comp' = 'none',
    'series' = "GSE93041",
    'group_design' = paste0('000XXX111')), #checked
  '2' = list(
    'Pub.' = 50,
    'Exp.' = 28,
    'Comp.' = 1,
    'sample_comp' = 'amitriptyline',
    'series' = "GSE63005",
    'group_design' = paste0('11XX00XXXXXX')), #checked
  '3' = list(
    'Pub.' = 50,
    'Exp.' = 28,
    'Comp.' = 2,
    'sample_comp' = 'citalopram',
    'series' = "GSE63005",
    'group_design' = paste0('XX1100XXXXXX')), #checked
  '4' = list(
    'Pub.' = 50,
    'Exp.' = 28,
    'Comp.' = 3,
    'sample_comp' = 'duloxetine',
    'series' = "GSE63005",
    'group_design' = paste0('XXXX0011XXXX')), #checked
  '5' = list(
    'Pub.' = 50,
    'Exp.' = 28,
    'Comp.' = 4,
    'sample_comp' = 'imipramine',
    'series' = "GSE63005",
    'group_design' = paste0('XXXX00XX11XX')), #checked
  '6' = list(
    'Pub.' = 50,
    'Exp.' = 28,
    'Comp.' = 5,
    'sample_comp' = 'mirtazapine',
    'series' = "GSE63005",
    'group_design' = paste0('XXXX00XXXX11')) #checked
)





opts_geo_dl[['GPL17223']][['comps']] <- list(
  '1' = list(
    'Pub.' = 36,
    'Exp.' = 22,
    'Comp.' = 1,
    'sample_comp' = 'cortex',
    'series' = "GSE47541",
    'group_design' = paste0('XXXXXXX11110000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')), #checked
  '2' = list(
    'Pub.' = 36,
    'Exp.' = 22,
    'Comp.' = 2,
    'sample_comp' = 'hippocampus',
    'series' = "GSE47541",
    'group_design' = paste0('XXXXXXXXXXXXXXXXXXXXXXX11110000XXXXXXXXXXXXXXX')), #checked
  '3' = list(
    'Pub.' = 36,
    'Exp.' = 22,
    'Comp.' = 3,
    'sample_comp' = 'dorsal_raphe',
    'series' = "GSE47541",
    'group_design' = paste0('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX11110000'))
) #checked




opts_geo_dl[['GPL25480']][['comps']] <- list(
  '1' = list(
    'Pub.' = 76,
    'Exp.' = 42,
    'Comp.' = 1,
    'sample_comp' = paste0('bupropion_control'),
    'series' = "GSE119291",
    'group_design' = GSE119291_group_design()[['one']]), #checked
  '2' = list(
    'Pub.' = 76,
    'Exp.' = 42,
    'Comp.' = 2,
    'sample_comp' = paste0('phenelzine_control'),
    'series' = "GSE119291",
    'group_design' = GSE119291_group_design()[['two']]), #checked
  '3' = list(
    'Pub.' = 76,
    'Exp.' = 42,
    'Comp.' = 3,
    'sample_comp' = paste0('trazodone_control'),
    'series' = "GSE119291",
    'group_design' = GSE119291_group_design()[['three']]), #checked
  '4' = list(
    'Pub.' = 76,
    'Exp.' = 42,
    'Comp.' = 4,
    'sample_comp' = paste0('tranylcypromine_control'),
    'series' = "GSE119291",
    'group_design' = GSE119291_group_design()[['four']]), #checked
  '5' = list(
    'Pub.' = 76,
    'Exp.' = 42,
    'Comp.' = 5,
    'sample_comp' = paste0('bupropion_SZ'),
    'series' = "GSE119291",
    'group_design' = GSE119291_group_design()[['five']]), #checked
  '6' = list(
    'Pub.' = 76,
    'Exp.' = 42,
    'Comp.' = 6,
    'sample_comp' = paste0('phenelzine_SZ'),
    'series' = "GSE119291",
    'group_design' = GSE119291_group_design()[['six']]), #checked
  '7' = list(
    'Pub.' = 76,
    'Exp.' = 42,
    'Comp.' = 7,
    'sample_comp' = paste0('trazodone_SZ'),
    'series' = "GSE119291",
    'group_design' = GSE119291_group_design()[['seven']]), #checked
  '8' = list(
    'Pub.' = 76,
    'Exp.' = 42,
    'Comp.' = 8,
    'sample_comp' = paste0('tranylcypromine_SZ'),
    'series' = "GSE119291",
    'group_design' = GSE119291_group_design()[['eight']]) #checked
)
###########################
### PREPARE COMPARISONS ###
###########################



#####################
### DOWNLOAD DATA ###
#####################
geo2r_download <- purrr::map(.x = opts_geo_dl, .f = function(platform_)
{
  comparison_download <- purrr::map(
    .x = platform_[['comps']],
    .f = function(comp){
      
      name <- paste(comp[['Exp.']], comp[['Comp.']], comp[['sample_comp']], comp[['series']], platform_[['platform']], sep = '_')
      
      table <- get_full_topTable(
        series_ = comp[['series']],
        platform_ = platform_[['platform']],
        group_names_ = comp[['group_design']],
        col_names_ = platform_[['col_names']],
        p_cutoff_ = 0.05,
        AnnotGPL_ = platform_[['AnnotGPL']]) #checked

      write.table(
        x = table,
        file = paste0(opts$dir_r_downloaded_data, '/', name, '.tsv'),
        sep = '\t',
        row.names = F,
        dec = ',')

      named_table <- list()
      named_table[[name]] <- table
      
      return(named_table)
    })
  
  return(comparison_download)
})

geo2r_download <- rlist::list.ungroup(rlist::list.ungroup(geo2r_download))

save(geo2r_download, file = paste0(opts$dir_r_downloaded_data, '/geo2r_download'))
#####################
### DOWNLOAD DATA ###
#####################










###############
### TESTING ###
###############






tT <- topTable(fit2, adjust="fdr", sort.by="B", number=length(fit2$genes[[1]]))

tT <- subset(tT, tT$P.Value < 0.05)

rm(cont.matrix, design, ex, fit, fit2, gset, tT, fl, gsms, i, idx, lofC, qx, sml)
###############
### TESTING ###
###############