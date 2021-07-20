library(annotatr)

annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic','hg19_genes_intronexonboundaries')
annotations = build_annotations(genome = 'hg19', annotations = annots)

dm_regions <- read_regions(con = "ml_top50.bed", genome = 'hg19', format = 'bed')
dm_annotated = annotate_regions(regions = dm_regions, annotations = annotations, ignore.strand = TRUE, quiet = FALSE)
df_dm_annotated <- data.frame(dm_annotated)
dm_random_regions <- randomize_regions(regions = dm_regions, allow.overlaps = TRUE, per.chromosome = TRUE)
dm_random_annotated <- annotate_regions(regions = dm_random_regions, annotations = annotations, ignore.strand = TRUE, quiet = TRUE)
dm_annsum <- summarize_annotations(annotated_regions = dm_annotated,quiet = TRUE)
dm_annsum_rnd <- summarize_annotations(annotated_regions = dm_annotated, annotated_random = dm_random_annotated, quiet = TRUE)
dm_numsum <- summarize_numerical(annotated_regions = dm_annotated, by = c('annot.type', 'annot.id'), over = c('diff_meth'), quiet = TRUE)

annots_order = c( 
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_genes_1to5kb', 
  'hg19_genes_promoters', 
  'hg19_genes_5UTRs', 
  'hg19_genes_exons', 
  'hg19_genes_intronexonboundaries', 
  'hg19_genes_introns', 
  'hg19_genes_intergenic') 
dm_vs_kg_annotations = plot_annotation( 
  annotated_regions = dm_annotated, 
  annotation_order = annots_order, 
  plot_title = '# of Sites Tested for DM annotated on chr9', 
  x_label = 'knownGene Annotations', 
  y_label = 'Count') 
print(dm_vs_kg_annotations)

annots_order = c( 
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_cpg_inter',
  'hg19_cpg_shelves',
  'hg19_genes_1to5kb', 
  'hg19_genes_promoters', 
  'hg19_genes_5UTRs', 
  'hg19_genes_exons', 
  'hg19_genes_intronexonboundaries', 
  'hg19_genes_introns', 
  'hg19_genes_intergenic') 

dm_vs_kg_annotations_wrandom = plot_annotation( 
  annotated_regions = dm_annotated, 
  annotated_random = dm_random_annotated, 
  annotation_order = annots_order, 
  plot_title = 'Dist. of Sites Tested for DM (with rndm.)', 
  x_label = 'Annotations', 
  y_label = 'Count') 
print(dm_vs_kg_annotations_wrandom)

annots_order = c( 
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_genes_1to5kb', 
  'hg19_genes_promoters', 
  'hg19_genes_5UTRs', 
  'hg19_genes_exons', 
  'hg19_genes_intronexonboundaries', 
  'hg19_genes_introns', 
  'hg19_genes_intergenic') 
dm_vs_coannotations = plot_coannotations( 
  annotated_regions = dm_annotated, 
  annotation_order = annots_order, 
  axes_label = 'Annotations', 
  plot_title = 'Regions in Pairs of Annotations') 
print(dm_vs_coannotations)

dm_vs_regions_annot = plot_numerical( 
  annotated_regions = dm_annotated, 
  x = 'mu0', 
  facet = 'annot.type', 
  facet_order = c( 
    'hg19_cpg_islands',
    'hg19_cpg_shores',
    'hg19_genes_1to5kb', 
    'hg19_genes_promoters', 
    'hg19_genes_5UTRs', 
    'hg19_genes_exons', 
    'hg19_genes_intronexonboundaries', 
    'hg19_genes_introns', 
    'hg19_genes_intergenic'), 
  bin_width = 5, 
  plot_title = 'Group 0 Region Methylation In Genes', 
  x_label = 'Group 0') 
print(dm_vs_regions_annot)

dm_vs_regions_name = plot_numerical( 
  annotated_regions = dm_annotated, 
  x = 'mu0', 
  y = 'mu1', 
  facet = 'annot.type', 
  facet_order = c('hg19_genes_1to5kb','hg19_genes_promoters', 
                  'hg19_genes_5UTRs','hg19_genes_3UTRs', 'hg19_custom_ezh2', 
                  'hg19_genes_intergenic', 'hg19_cpg_islands', 'hg19_cpg_shores'), 
  plot_title = 'Region Methylation: Group 0 vs Group 1', 
  x_label = 'Group 0', 
  y_label = 'Group 1') 
print(dm_vs_regions_name)

dm_vs_num_co = plot_numerical_coannotations( 
  annotated_regions = dm_annotated, 
  x = 'mu0', 
  annot1 = 'hg19_cpg_islands', 
  annot2 = 'hg19_genes_promoters', 
  bin_width = 5, 
  plot_title = 'Group 0 Perc. Meth. in CpG Islands and Promoters', 
  x_label = 'Percent Methylation') 
print(dm_vs_num_co)

# The orders for the x-axis labels. This is also a subset
# of the labels (hyper, hypo, none). 
x_order = c( 
  'hyper', 
  'hypo') 
# The orders for the fill labels. Can also use this
# parameter to subset annotation types to fill. 
fill_order = c( 
  'hg19_cpg_islands', 
  'hg19_cpg_shores', 
  'hg19_cpg_shelves', 
  'hg19_cpg_inter') 
# Make a barplot of the data class where each bar
# is composed of the counts of CpG annotations. 
dm_vs_cpg_cat1 = plot_categorical( 
  annotated_regions = dm_annotated, x='DM_status', fill='annot.type', 
  x_order = x_order, fill_order = fill_order, position='stack', 
  plot_title = 'DM Status by CpG Annotation Counts', 
  legend_title = 'Annotations', 
  x_label = 'DM status', 
  y_label = 'Count') 
print(dm_vs_cpg_cat1)

# Use the same order vectors as the previous code block,
# but use proportional fill instead of counts. 

# Make a barplot of the data class where each bar
# is composed of the *proportion* of CpG annotations. 
dm_vs_cpg_cat2 = plot_categorical( 
  annotated_regions = dm_annotated, x='DM_status', fill='annot.type', 
  x_order = x_order, fill_order = fill_order, position='fill', 
  plot_title = 'DM Status by CpG Annotation Proportions', 
  legend_title = 'Annotations', 
  x_label = 'DM status', 
  y_label = 'Proportion') 
print(dm_vs_cpg_cat2)

# Add in the randomized annotations for "Random Regions" bar 

# Make a barplot of the data class where each bar
# is composed of the *proportion* of CpG annotations, and
# includes "All" regions tested for DM and "Random Regions"
# regions consisting of randomized regions. 
dm_vs_cpg_cat_random = plot_categorical( 
  annotated_regions = dm_annotated, annotated_random = dm_random_annotated, 
  x='DM_status', fill='annot.type', 
  x_order = x_order, fill_order = fill_order, position='fill', 
  plot_title = 'DM Status by CpG Annotation Proportions', 
  legend_title = 'Annotations', 
  x_label = 'DM status', 
  y_label = 'Proportion') 
print(dm_vs_cpg_cat_random)

x_order = c( 
  'hg19_custom_ezh2', 
  'hg19_genes_1to5kb', 
  'hg19_genes_promoters', 
  'hg19_genes_5UTRs', 
  'hg19_genes_exons', 
  'hg19_genes_introns', 
  'hg19_genes_3UTRs', 
  'hg19_genes_intergenic') 
# The orders for the fill labels. 
fill_order = c( 
  'hyper', 
  'hypo', 
  'none') 
dm_vs_kg_cat = plot_categorical( 
  annotated_regions = dm_annotated, x='annot.type', fill='DM_status', 
  x_order = x_order, fill_order = fill_order, position='fill', 
  legend_title = 'DM Status', 
  x_label = 'knownGene Annotations', 
  y_label = 'Proportion') 
print(dm_vs_kg_cat)