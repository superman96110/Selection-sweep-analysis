#计算FST(windows size 50k, wndows step 10k)
nohup vcftools --vcf /workspace/home/goat/supeng/bodysize/horse/height_male_horse/fst/filter_horse_gtex.recode.vcf --weir-fst-pop /workspace/home/goat/supeng/bodysize/horse/height_male_horse/fst/low.txt --weir-fst-pop /workspace/home/goat/supeng/bodysize/horse/height_male_horse/fst/high.txt --fst-window-size 50000 --fst-window-step 10000 --out horse_fst_height_male.fst &

#计算FST snp
nohup vcftools --vcf /workspace/home/goat/supeng/bodysize/horse/height_male_horse/fst/filter_horse_gtex.recode.vcf --weir-fst-pop /workspace/home/goat/supeng/bodysize/horse/height_male_horse/fst/low.txt --weir-fst-pop /workspace/home/goat/supeng/bodysize/horse/height_male_horse/fst/high.txt --out horse_fst_height_male.fst &
