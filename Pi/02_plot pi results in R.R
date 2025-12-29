setwd("F:/caas/毕业课题/第三章_选择信号/horse data/selection sweep/pi/258_ind/")
library(ggplot2)
library(dplyr)
library(scales)

# =========================
# 参数：基因区间 + 扩展
# =========================
chrom_id   <- "6"

gene_start <- 82552985
gene_end   <- 82681843

flank_bp   <- 400000
plot_start <- max(0, gene_start - flank_bp)
plot_end   <- gene_end + flank_bp

# =========================
# 读取 pi 结果（vcftools 输出一般是 *.windowed.pi）
# =========================
high_pi <- read.table("high.pi.windowed.pi", header = TRUE, stringsAsFactors = FALSE)
low_pi  <- read.table("low.pi.windowed.pi",  header = TRUE, stringsAsFactors = FALSE)

# 必要列检查
stopifnot(all(c("CHROM","BIN_START","BIN_END","PI") %in% names(high_pi)))
stopifnot(all(c("CHROM","BIN_START","BIN_END","PI") %in% names(low_pi)))

high_pi <- high_pi %>% mutate(GROUP = "High")
low_pi  <- low_pi  %>% mutate(GROUP = "Low")

pi_data <- bind_rows(high_pi, low_pi) %>%
    filter(CHROM == chrom_id) %>%
    # ✅ 重叠筛选（避免左边界缺点）
    filter(BIN_START <= plot_end, BIN_END >= plot_start) %>%
    mutate(POS = (BIN_START + BIN_END) / 2)

ymax <- max(pi_data$PI, na.rm = TRUE)

# =========================
# 作图（颜色区分 + 右上角图例）
# =========================
p <- ggplot(pi_data, aes(x = POS, y = PI, group = GROUP, color = GROUP)) +
    # 基因区间灰色高亮
    annotate("rect",
             xmin = gene_start, xmax = gene_end,
             ymin = -Inf, ymax = Inf,
             fill = "grey70", alpha = 0.6) +
    
    geom_line(linewidth = 0.45, na.rm = TRUE) +
    
    coord_cartesian(
        xlim = c(plot_start, plot_end),
        ylim = c(0, ymax * 1.05),
        expand = FALSE
    ) +
    
    scale_x_continuous(
        breaks = pretty_breaks(n = 6),
        labels = label_number(accuracy = 1, big.mark = "")
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    
    labs(x = NULL, y = expression(pi), color = NULL) +
    
    # 颜色（你可以按图里习惯换成自己想要的两种颜色）
    scale_color_manual(values = c(High = "#1f77b4", Low = "#d62728")) +
    
    theme_classic(base_size = 9) +
    theme(
        panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.3, color = "grey30"),
        axis.text = element_text(size = 7, color = "grey10"),
        axis.title.y = element_text(size = 8),
        
        # ✅ 图例右上角（panel 内）
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.key = element_rect(fill = NA, color = NA),
        legend.text = element_text(size = 7),
        
        plot.margin = margin(2, 2, 2, 2)
    )

print(p)
