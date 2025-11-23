library(ggplot2)
library(cowplot)
bands <- data.frame(
  ybottom = c(-15, log(1/3), log(3)),
  ytop    = c(log(1/3), log(3), 9),
  color   = c("green", "grey", "red")
)

plot_data_summary <- function(data, extra_title = "") {
  # Process the input data
  data_summary2 <- data %>%
    mutate(
      logJAB = log(JAB),
      N_shape = case_when(
        N < 30 ~ "n < 30",
        N >= 30 & N < 900 ~ "30 \u2264 n < 900",
        N >= 900 & N < 30000 ~ "900 \u2264 n < 30k",
        N >= 30000 & N < 900000 ~ "30k \u2264 n < 900k",
        N >= 900000 ~ "n \u2265 900k"
      ),
      N_shape = factor(
        N_shape,
        levels = c(
          "n < 30",
          "30 \u2264 n < 900",
          "900 \u2264 n < 30k",
          "30k \u2264 n < 900k",
          "n \u2265 900k"
        )
      )
    )

  total_analyses <- data_summary2 %>%
    pull(analysisId) %>%
    unique() %>%
    length()

  test_counts <- data_summary2 %>%
    group_by(statisticalMethod) %>%
    summarise(n = n_distinct(analysisId), .groups = "drop")

  test_labels <- setNames(
    paste0(test_counts$statisticalMethod, " (n=", test_counts$n, ")"),
    test_counts$statisticalMethod
  )

  data_summary2 <- data_summary2 %>%
    mutate(statisticalMethod = factor(statisticalMethod, levels = names(test_labels)))

  # Main plot
  p <- ggplot(data_summary2, aes(x = pValue, y = logJAB, shape = N_shape)) +
    geom_rect(data = bands, aes(xmin = 0, xmax = 1, ymin = ybottom, ymax = ytop, fill = color),
              alpha = 0.2, inherit.aes = FALSE) +
    geom_point(aes(color = statisticalMethod), size = 1) +
    ylim(-15, 9) +
    scale_x_continuous(name = "pValue", limits = c(0, 1)) +
    scale_y_continuous(name = "ln(eJAB01)", breaks = seq(-15, 9, 2)) +
    scale_fill_identity() +
    scale_color_brewer(
      name = "Test Type",
      palette = "Set3", labels = test_labels
    ) +
    scale_shape_manual(
      name = "Sample Size",
      values = c(
        "n < 30" = 17,
        "30 \u2264 n < 900" = 15,
        "900 \u2264 n < 30k" = 16,
        "30k \u2264 n < 900k" = 1,
        "n \u2265 900k" = 8
      )
    ) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "black") +
    geom_hline(yintercept = log(1/3), linetype = "dashed", color = "black") +
    geom_hline(yintercept = log(3), linetype = "dashed", color = "black") +
    theme_minimal() +
    ggtitle(paste("ln(eJAB01) vs pValue - Analyses:", total_analyses, extra_title))

  # Inset plot
  inset_plot <- p +
    coord_cartesian(xlim = c(0, 0.05), ylim = c(-3, 3)) +
    scale_x_continuous(name = "", breaks = c(0, 0.025, 0.05)) +
    scale_y_continuous(name = "", breaks = c(-3, 0, 3)) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.background = element_rect(fill = "white", color = "black", size = 1.5),
      axis.text = element_text(size = 8),
      axis.ticks = element_line(size = 0.5)
    )

  print(ggdraw(p) +
          draw_plot(inset_plot, x = 0.2, y = 0.1, width = 0.4, height = 0.4))
}
