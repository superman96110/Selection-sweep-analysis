#根据bed的注释文件，对top1%和top5%的fst位点进行功能注释
#bed文件时根据gff文件生成的，用于bedtools的注释

#!/bin/bash

input_file="/data/supeng/bodysize/horse/fst/height/horse_high_low_height_fst.windowed.weir.fst"
output_dir="/data/supeng/bodysize/horse/fst/height/"
annotation_file="/data/supeng/bodysize/horse/fst/height/horse3-gene_nochr.bed"

# 创建输出目录（如果不存在）
mkdir -p $output_dir

# 1. 对文件进行排序，按FST值从高到低排序
sort -k5,5gr $input_file > ${output_dir}sorted_all.fst.50kb

# 2. 计算总行数
total_lines=$(wc -l < ${output_dir}sorted_all.fst.50kb)

# 检查文件是否为空
if [ "$total_lines" -eq 0 ]; then
            echo "输入文件为空，终止运行。"
                exit 1
fi

# 3. 计算前1%和前5%的行数
top_1_percent_lines=$(echo "($total_lines * 0.01)/1" | bc)
top_5_percent_lines=$(echo "($total_lines * 0.05)/1" | bc)

# 检查行数是否有效
if [ "$top_1_percent_lines" -le 0 ] || [ "$top_5_percent_lines" -le 0 ]; then
            echo "行数计算出错，请检查输入文件。"
                exit 1
fi

# 4. 提取最高的前1%和5%的行，并保留所需的列
head -n $top_1_percent_lines ${output_dir}sorted_all.fst.50kb | awk '{print $1"\t"$2"\t"$3"\t"$5}' > ${output_dir}top_1_percent.bed
head -n $top_5_percent_lines ${output_dir}sorted_all.fst.50kb | awk '{print $1"\t"$2"\t"$3"\t"$5}' > ${output_dir}top_5_percent.bed

# 5. 使用 bedtools 进行注释
bedtools intersect -a ${output_dir}top_1_percent.bed -b $annotation_file -wa -wb > ${output_dir}top_1_percent.ann
bedtools intersect -a ${output_dir}top_5_percent.bed -b $annotation_file -wa -wb > ${output_dir}top_5_percent.ann

# 6. 提取基因列表并去重
awk '{print $NF}' ${output_dir}top_1_percent.ann | sort | uniq > ${output_dir}top1gene.txt
awk '{print $NF}' ${output_dir}top_5_percent.ann | sort | uniq > ${output_dir}top5gene.txt

echo "Processing complete. Files saved in $output_dir"


