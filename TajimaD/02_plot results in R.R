library(ggplot2)
library(dplyr)
library(scales)

# =========================
# 1) 参数：染色体、基因区间、扩展长度
# =========================
chrom_id   <- "6"

gene_start <- 82552985
gene_end   <- 82681843

flank_bp   <- 400000
plot_start <- max(0, gene_start - flank_bp)
plot_end   <- gene_end + flank_bp

# TajimaD 的窗口大小（vcftools --TajimaD 后面的数字）
tajima_win <- 50000

# =========================
# 2) 读取 TajimaD 结果
# =========================
high_td <- read.table("high.Tajima.D", header = TRUE, stringsAsFactors = FALSE)
low_td  <- read.table("low.Tajima.D",  header = TRUE, stringsAsFactors = FALSE)

# 统一 TajimaD 列名为 TAJIMA_D（兼容不同 vcftools 输出列名）
rename_tajima_col <- function(df) {
  cand <- names(df)[grepl("tajima", names(df), ignore.case = TRUE)]
  if (length(cand) == 0) stop("找不到 TajimaD 列，请检查文件列名：", paste(names(df), collapse = ", "))
  df %>% rename(TAJIMA_D = all_of(cand[1]))
}

high_td <- rename_tajima_col(high_td)
low_td  <- rename_tajima_col(low_td)

# 必要列检查
stopifnot(all(c("CHROM", "BIN_START", "TAJIMA_D") %in% names(high_td)))
stopifnot(all(c("CHROM", "BIN_START", "TAJIMA_D") %in% names(low_td)))

# 如果没有 BIN_END，则用窗口大小补一个 BIN_END
if (!"BIN_END" %in% names(high_td)) high_td$BIN_END <- high_td$BIN_START + tajima_win
if (!"BIN_END" %in% names(low_td))  low_td$BIN_END  <- low_td$BIN_START + tajima_win

high_td <- high_td %>% mutate(GROUP = "High")
low_td  <- low_td  %>% mutate(GROUP = "Low")

td_data <- bind_rows(high_td, low_td) %>%
  filter(CHROM == chrom_id) %>%
  # 重叠筛选（避免左边界缺点）
  filter(BIN_START <= plot_end, BIN_END >= plot_start) %>%
  mutate(POS = (BIN_START + BIN_END) / 2)

# y轴范围：不对称，按数据范围自动，并留边距
ymin <- min(td_data$TAJIMA_D, na.rm = TRUE)
ymax <- max(td_data$TAJIMA_D, na.rm = TRUE)
ypad <- (ymax - ymin) * 0.05
if (!is.finite(ypad) || ypad == 0) ypad <- 0.1  # 极端情况下防止 0 范围

# =========================
# 3) 作图（无0参考线 + y轴不对称）
# =========================
p <- ggplot(td_data, aes(x = POS, y = TAJIMA_D, group = GROUP, color = GROUP)) +

  # 基因区间灰色高亮
  annotate("rect",
           xmin = gene_start, xmax = gene_end,
           ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.6) +

  geom_line(linewidth = 0.45, na.rm = TRUE) +

  coord_cartesian(
    xlim = c(plot_start, plot_end),
    ylim = c(ymin - ypad, ymax + ypad),
    expand = FALSE
  ) +

  scale_x_continuous(
    breaks = pretty_breaks(n = 6),
    labels = label_number(accuracy = 1, big.mark = "")
  ) +

  labs(x = NULL, y = "Tajima's D", color = NULL) +

  scale_color_manual(values = c(High = "#1f77b4", Low = "#d62728")) +

  theme_classic(base_size = 9) +
  theme(
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.3, color = "grey30"),
    axis.text = element_text(size = 7, color = "grey10"),
    axis.title.y = element_text(size = 8),

    # 图例右上角（panel 内）
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.key = element_rect(fill = NA, color = NA),
    legend.text = element_text(size = 7),

    plot.margin = margin(2, 2, 2, 2)
  )

print(p)
