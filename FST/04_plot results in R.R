# ===== 1) 颜色：31条染色体 =====
colorset <- c(
    "#3E0A52", "#423D77", "#3F678B", "#468C8D", "#5FB47F", "#9FD55C",
    "#F9E956", "#A3D2CA", "#84A9AC", "#627D98", "#5A5C6B", "#C7D2FE",
    "#A9C0F5", "#7BAAF7", "#536DFE", "#3A3A90", "#002147", "#0047AB",
    "#1E90FF", "#00BFFF", "#87CEFA", "#ADD8E6", "#B0E0E6", "#00CED1",
    "#48D1CC", "#FF6347", "#FFD700", "#ADFF2F", "#40E0D0", "#FF69B4",
    "#BA55D3"
)
# 确保正好31种
colorset <- colorset[1:31]

# 染色体标签 1-31
chr_labels <- as.character(1:31)

# ===== 2) 读 windowed fst 结果（马）=====
# 你的文件名：horse_high_low_height_fst.windowed.weir.fst
fst_raw <- fread("horse_high_low_height_fst.windowed.weir.fst", header = TRUE)

# 统一CHR格式（去chr前缀）并限制1-31
fst_raw[, CHR := as.numeric(gsub("^chr", "", as.character(CHROM), ignore.case = TRUE))]
fst_raw <- fst_raw[!is.na(CHR) & CHR >= 1 & CHR <= 31]

# CMplot输入：SNP, CHR, BP, Fst
# BP建议用 BIN_START（或者用MID都行，但你现在习惯BIN_START就保持）
fst_data <- fst_raw[, .(
    SNP = paste0("chr", CHR, ":", as.numeric(BIN_START)),
    CHR = CHR,
    BP  = as.numeric(BIN_START),
    Fst = as.numeric(WEIGHTED_FST)
)]
fst_data <- fst_data[is.finite(Fst)]
fst_data <- na.omit(fst_data)

# 阈值（1% / 5%）
threshold_1percent_fst <- as.numeric(quantile(fst_data$Fst, 0.99, na.rm = TRUE))
threshold_5percent_fst <- as.numeric(quantile(fst_data$Fst, 0.95, na.rm = TRUE))

# ===== 3) 读 top1 注释（马）=====
# 你的文件名：top_1_percent.ann
# 你之前定义过：CHROM START END SOME_VALUE SOME_ID GENE_START GENE_END GENE_NAME
ann <- fread("top_1_percent.ann", header = FALSE)
setnames(ann, c("CHROM","START","END","SOME_VALUE","SOME_ID","GENE_START","GENE_END","GENE_NAME"))

ann[, CHR := as.numeric(gsub("^chr", "", as.character(CHROM), ignore.case = TRUE))]
ann <- ann[!is.na(CHR) & CHR >= 1 & CHR <= 31]

# 取Top基因（你要标几个自己调，下面默认30）
top_n_gene <- 10
ann_sorted <- ann[order(-SOME_VALUE)]
top_genes_unique <- unique(ann_sorted$GENE_NAME)[1:top_n_gene]
top_ann <- ann_sorted[GENE_NAME %in% top_genes_unique][order(-SOME_VALUE), .SD[1], by = GENE_NAME]

# 把基因对应到“窗口点”(用该基因区间的START作为窗口定位；也可改成MID)
top_fst_SNPs <- paste0("chr", top_ann$CHR, ":", as.numeric(top_ann$START))
top_fst_gene_names <- top_ann$GENE_NAME

# 关键：确保这些highlight点在fst_data里真的存在（否则CMplot不显示/报错）
hit <- top_fst_SNPs %in% fst_data$SNP
top_fst_SNPs <- top_fst_SNPs[hit]
top_fst_gene_names <- top_fst_gene_names[hit]

# ===== 4) 画 Fst 曼哈顿图（保持你原来参数风格）=====
CMplot(
    fst_data[, .(SNP, CHR, BP, Fst)],
    type = "p",
    plot.type = "m",
    LOG10 = FALSE,
    
    ylim = c(0, 1),          # 你要求最大1（如果你的Fst>1会被截断）
    ylab = "FST",
    col = colorset,
    
    # 注意你原来 threshold 顺序是(1%,5%)；我保持一致
    threshold = c(threshold_1percent_fst, threshold_5percent_fst),
    threshold.col = c("blue", "red"),
    
    chr.labels = chr_labels,
    main = "Fst Manhattan Plot (Horse chr1-31)",
    
    highlight = top_fst_SNPs,
    highlight.col = "green",
    highlight.text = top_fst_gene_names,
    
    amplify = FALSE,
    file = NULL
)
