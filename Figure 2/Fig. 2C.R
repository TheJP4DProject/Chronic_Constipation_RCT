library(ggplot2)
library(dplyr)

df <- read.csv("tableS2_moutsp.csv") %>%
  filter(qval < 0.1, abs(coef) > 0.15)

plot_data <- df %>%
  mutate(
    coef_abs = abs(coef),
    direction = ifelse(coef > 0, "Positive", "Negative"),
    display_name = gsub("\\[ID:(\\d+)\\]", "\\1", GTDB),
    signif_label = case_when(
      qval < 0.01 ~ "**",
      qval < 0.05 ~ "*",
      TRUE ~ ""
    ),
    signif_y = coef_abs + 0.02
  ) %>%
  arrange(desc(coef)) %>%
  mutate(display_name = factor(display_name, levels = display_name))

panel_theme <- theme_minimal(base_size = 5, base_family = "Helvetica") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 5, face = "bold"),
    axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none"
  )

p <- ggplot(plot_data, aes(display_name, coef_abs, fill = direction)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.1) +
  geom_text(
    aes(y = signif_y, label = signif_label),
    size = 5 / 2.835
  ) +
  scale_fill_manual(
    values = c("Positive" = "#D55E00", "Negative" = "#0072B2")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "CC maa (FDR < 0.1, |coef| > 0.15)",
    y = "Effect size",
    x = NULL
  ) +
  panel_theme

ggsave(
  "Fig2C.pdf",
  p,
  width = 5,
  height = 3
)
