#这里进行单位点的FST计算，但是单位点的结果并不是特别准确，需要结合windows的结果结合来看

vcftools --vcf 825_filter_maf005_geno01_mind01.vcf --weir-fst-pop high_height_ind.txt --weir-fst-pop low_height_ind.txt --out horse_high_low_height_fst
