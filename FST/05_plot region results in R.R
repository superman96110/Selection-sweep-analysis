#这个是绘制区间FST折线图的代码
library(ggplot2)
library(dplyr)
library(scales)

# =========================
# 1) 读取数据
# =========================
fst_data <- read.table(
  "horse_high_low_height_fst.windowed.weir.fst",
  header = TRUE, stringsAsFactors = FALSE
)

# （可选）检查需要的列是否存在
need_cols <- c("CHROM", "BIN_START", "WEIGHTED_FST")
miss <- setdiff(need_cols, names(fst_data))
if (length(miss) > 0) stop("缺少列：", paste(miss, collapse = ", "))

# 如果没有 BIN_END，就用 BIN_START + window_size 近似一个 BIN_END（不推荐，但可用）
if (!"BIN_END" %in% names(fst_data)) {
  warning("数据中没有 BIN_END 列，将用相邻 BIN_START 推断窗口大小来近似 BIN_END。")
  fst_data <- fst_data %>%
    arrange(CHROM, BIN_START) %>%
    group_by(CHROM) %>%
    mutate(win = lead(BIN_START) - BIN_START,
           win = ifelse(is.na(win), median(win, na.rm = TRUE), win),
           BIN_END = BIN_START + win) %>%
    ungroup()
}

# =========================
# 2) 设定要绘制的区域
# =========================
chrom_id  <- "6"
start_pos <- 82502985
end_pos   <- 82731843

# 灰色高亮区间（可选）
highlight_start <- 82610000
highlight_end   <- 82650000

# =========================
# 3) 取区间数据（关键修复：用“重叠窗口”）
#    只要窗口与 [start_pos, end_pos] 有重叠就保留
# =========================
interval_data <- fst_data %>%
  filter(CHROM == chrom_id) %>%
  filter(BIN_START <= end_pos, BIN_END >= start_pos) %>%   # ✅ 修复缺点
  mutate(POS = (BIN_START + BIN_END) / 2)

# =========================
# 4) 作图
# =========================
ymax <- max(interval_data$WEIGHTED_FST, na.rm = TRUE)

p <- ggplot(interval_data, aes(x = POS, y = WEIGHTED_FST)) +
  # 灰色高亮块（不需要就把 highlight_* 设为 NA）
  { if (!anyNA(c(highlight_start, highlight_end)))
      geom_rect(
        aes(xmin = highlight_start, xmax = highlight_end,
            ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "grey70", alpha = 0.6
      )
  } +
  geom_line(color = "black", linewidth = 0.4) +
  coord_cartesian(ylim = c(0, ymax * 1.05), expand = FALSE) +
  scale_x_continuous(
    limits = c(start_pos, end_pos),
    expand = c(0, 0),  
    breaks = pretty_breaks(n = 6),
    labels = label_number(accuracy = 1, big.mark = "")  # 不加逗号
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Fst") +
  theme_classic(base_size = 9) +
  theme(
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.3, color = "grey30"),
    axis.text = element_text(size = 7, color = "grey10"),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(2, 2, 2, 2)
  )

print(p)
