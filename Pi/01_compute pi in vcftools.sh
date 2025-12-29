vcftools --vcf 825_filter_maf005_geno01_mind01.vcf --keep high_height_ind.txt --window-pi 50000 --window-pi-step 10000 --out high.pi

vcftools --vcf 825_filter_maf005_geno01_mind01.vcf --keep low_height_ind.txt --window-pi 50000 --window-pi-step 10000 --out low.pi
