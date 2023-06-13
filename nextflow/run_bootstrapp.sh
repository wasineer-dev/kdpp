nextflow run vlmc_bootstrap.nf --ndup 5 --sampletarget=320000000 --output=distance_10x.tsv
nextflow clean -f
nextflow run vlmc_bootstrap.nf --ndup 5 --sampletarget=160000000 --output=distance_5x.tsv
nextflow clean -f
nextflow run vlmc_bootstrap.nf --ndup 5 --sampletarget=96000000 --output=distance_3x.tsv
nextflow clean -f
