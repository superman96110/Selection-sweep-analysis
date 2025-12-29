
vcftools --vcf 825_filter_maf005_geno01_mind01.vcf --keep low_height_ind.txt --TajimaD 50000 --out low

vcftools --vcf 825_filter_maf005_geno01_mind01.vcf --keep high_height_ind.txt --TajimaD 50000 --out high
