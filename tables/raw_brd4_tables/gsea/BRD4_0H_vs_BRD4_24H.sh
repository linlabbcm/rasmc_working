#!/usr/bin/bash 

echo "Running GSEA for BRD4 0H vs BRD4 24H"

#java -cp /storage/cylin/home/rhirsch/gsea-3.0.jar -Xmx4096m xtools.gsea.Gsea -res /storage/cylin/grail/projects/rasmc_all/tables/raw_brd4_tables/BRD4_0H_vs_BRD4_24H.gct -cls /storage/cylin/grail/projects/rasmc_all/tables/raw_brd4_tables/BRD4_0H_vs_BRD4_24H.cls#BRD4_0H_versus_BRD4_24H -gmx /storage/cylin/grail/projects/rasmc_all/tables/raw_brd4_tables/CLUSTER_GENE_SETS.gmt -chip gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip -collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label BRD4_0H_versus_BRD4_24H_clusters -metric log2_Ratio_of_Classes -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out /storage/cylin/grail/projects/rasmc_all/tables/raw_brd4_tables/gsea/ -gui false

java -Xmx4000m -cp /storage/cylin/home/cl6/gsea2-3.0_beta_2.jar xtools.gsea.Gsea -res /storage/cylin/grail/projects/rasmc_all/tables/raw_brd4_tables/BRD4_0H_vs_BRD4_24H.gct -cls /storage/cylin/grail/projects/rasmc_all/tables/raw_brd4_tables/BRD4_0H_vs_BRD4_24H.cls#BRD4_0H_versus_BRD4_24H -gmx /storage/cylin/grail/projects/rasmc_all/tables/raw_brd4_tables/CLUSTER_GENE_SETS.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label BRD4_0H_versus_BRD4_24H_clusters -metric log2_Ratio_of_Classes -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out /storage/cylin/grail/projects/rasmc_all/gsea/ -gui false
