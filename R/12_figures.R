#-----------------------------------------------------------------------------------------------------------------------

# Title: Construct plots and visualizations for methods paper
# Author: Marcus Becker

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# Set path to data folders; note that this will be updated once data has been stored.
root <- ""

# Set ggplot theme to stay consistent
ggplot2::theme_set(theme_light())

#-----------------------------------------------------------------------------------------------------------------------

# 1. WMU plots

# Load WMU model and raw data
load(paste0(root, "results/methods-paper/WMU-analysis/model_raw_data_2021-07-22.rdata"))

# Make predictions for initial plotting
x <- seq(0, 0.77, 0.01) # 0.77 is the highest aerial density
pred <- predict(m.lm, newdata = data.frame(aerial_avg = x), se.fit = TRUE)
pred.line <- data.frame(x, pred$fit)
pred.line.u <- data.frame(x, u = pred$fit + 1.65*pred$se.fit)
pred.line.l <- data.frame(x, l = pred$fit - 1.65*pred$se.fit)

all.pred <- pred.line %>% left_join(pred.line.u) %>% left_join(pred.line.l) %>% mutate(pass = "initial")

# Make initial aerial/camera WMU density comparison plot:
df_sum_dens %>%
  ggplot(aes(x = aerial_avg, y = density_avg)) +
  geom_linerange(aes(ymin = density_lci_0.9, ymax = density_uci_0.9),
                 color = "gray70", size = 0.7) +
  geom_linerange(aes(xmin = aerial_lci_0.9, xmax = aerial_uci_0.9),
                 color = "gray70", size = 0.7) +
  geom_ribbon(data = all.pred, inherit.aes = FALSE, aes(x = x, ymax = u-0.005, ymin = l+0.005), fill = "grey", alpha = 0.3) +
  geom_line(data = all.pred, aes(x = x, y = pred.fit), color = "gray30", size = 1.3) +
  geom_line(data = all.pred, aes(x = x, y = l), color = "grey40", linetype = "dotted", size = 0.5) +
  geom_line(data = all.pred, aes(x = x, y = u), color = "grey60", linetype = "dotted") +
  geom_abline(intercept = c(0, 0), linetype = 2) +
  geom_point(aes(size = n_deployments), color = "grey20", fill = "grey60", shape = 21, stroke = 1.15) +
  scale_size_continuous(labels = c("10-20", "21-35", "36-50", "51-80", "81+"),
                        breaks = c(9, 22, 35, 50, 85),
                        limits = c(8, 200),
                        range = c(2, 10)) +
  scale_y_continuous(breaks = seq(0, 4, 0.25), limits = c(0, 4)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  guides(size = guide_legend(title = "Number of Cameras:")) +
  labs(y = expression(Camera~density~(animals~per~km^2)),
       x = expression(Aerial~density~(animals~per~km^2))) +
  theme(legend.position = c(0.88, 0.80),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0), size = 13),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0), size = 13),
        panel.grid = element_blank(),
        legend.spacing.x = unit(0.15, 'cm'))

ggsave(filename = paste0(root, "results/methods-paper/figures/density-comp.png"), dpi = 500, width = 7, height = 5)

# Make predictions for corrected plotting
x <- seq(0, 0.77, 0.01) # 0.77 is the highest aerial density
pred <- predict(m.lm.corr, newdata = data.frame(aerial_avg = x), se.fit = TRUE)
pred.line <- data.frame(x, pred$fit)
pred.line.u <- data.frame(x, u = pred$fit + 1.65*pred$se.fit)
pred.line.l <- data.frame(x, l = pred$fit - 1.65*pred$se.fit)

all.pred <- pred.line %>% left_join(pred.line.u) %>% left_join(pred.line.l) %>% mutate(pass = "corrected")

# Make corrected aerial/camera WMU density comparison plot:
df_sum_dens_corr %>%
  ggplot(aes(x = aerial_avg, y = density_avg)) +
  geom_linerange(aes(ymin = density_lci_0.9, ymax = density_uci_0.9),
                 color = "gray70", size = 0.7) +
  geom_linerange(aes(xmin = aerial_lci_0.9, xmax = aerial_uci_0.9),
                 color = "gray70", size = 0.7) +
  geom_ribbon(data = all.pred, inherit.aes = FALSE, aes(x = x, ymax = u-0.005, ymin = l+0.005), fill = "grey", alpha = 0.3) +
  geom_line(data = all.pred, aes(x = x, y = pred.fit), color = "gray30", size = 1.3) +
  geom_line(data = all.pred, aes(x = x, y = l), color = "grey40", linetype = "dotted", size = 0.5) +
  geom_line(data = all.pred, aes(x = x, y = u), color = "grey60", linetype = "dotted") +
  geom_abline(intercept = c(0, 0), linetype = 2) +
  geom_point(aes(size = n_deployments), color = "grey20", fill = "grey60", shape = 21, stroke = 1.15) +
  scale_size_continuous(labels = c("10-20", "21-35", "36-50", "51-80", "81+"),
                        breaks = c(9, 22, 35, 50, 85),
                        limits = c(8, 200),
                        range = c(2, 10)) +
  scale_y_continuous(breaks = seq(0, 4, 0.25), limits = c(0, 4)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  guides(size = guide_legend(title = "Number of Cameras:")) +
  labs(y = expression(Camera~density~(animals~per~km^2)),
       x = expression(Aerial~density~(animals~per~km^2))) +
  theme(legend.position = c(0.88, 0.80),
        legend.title = element_text(size = 10),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0), size = 13),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0), size = 13),
        panel.grid = element_blank(),
        legend.spacing.x = unit(0.15, 'cm'))

ggsave(filename = paste0(root, "results/methods-paper/figures/density-comp-corr.png"), dpi = 500, width = 7, height = 5)

#-----------------------------------------------------------------------------------------------------------------------

d <- df_sum_dens %>%
  ungroup() %>%
  mutate(Camera = (density_uci_0.9 - density_lci_0.9) / density_avg,
         Aerial = (aerial_uci_0.9 - aerial_lci_0.9) / aerial_avg) %>%
  select(WMUNIT_COD, n_deployments, Aerial, Camera) %>%
  pivot_longer(cols = c(Aerial, Camera), names_to = "Method:", values_to = "width")

d1 <- d %>% filter(`Method:` == "Aerial")
max.aerial <- max(d1$width, na.rm = TRUE)
min.aerial <- min(d1$width, na.rm = TRUE)

d %>%
  filter(`Method:` == "Camera") %>%
  ggplot(aes(x = n_deployments, y = width)) +
  geom_hline(yintercept = max.aerial + 0.01, color = "grey20") +
  geom_hline(yintercept = min.aerial, color = "grey20") +
  scale_y_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3)) +
  scale_x_continuous(breaks = seq(0, 220, 20), limits = c(0, 220)) +
  annotate(geom = "rect", ymax = max.aerial, ymin = min.aerial, xmax = Inf, xmin = -Inf, fill = "grey", alpha = 0.3) +
  annotate(geom = "text", label = "Range of aerial survey precision", x = 30, y = 0.41, size = 3.25) +
  geom_point(size = 4, color = "grey20", fill = "grey60", shape = 21, stroke = 1.15) +
  coord_cartesian(xlim = c(0, 220)) +
  labs(y = "90% Confidence Interval Width / Mean",
       x = "Number of Cameras in the WMU") +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = c(0.86, 0.84),
        axis.title.x = element_text(margin = margin(15, 0, 0, 0), size = 13),
        axis.title.y = element_text(margin = margin(0, 15, 0, 0), size = 13),
        panel.grid.minor = element_blank(),
        legend.spacing.x = unit(0.15, 'cm'))

ggsave(filename = paste0(root, "results/methods-paper/figures/precision-comp.png"), dpi = 500, width = 7, height = 5)

#-----------------------------------------------------------------------------------------------------------------------
