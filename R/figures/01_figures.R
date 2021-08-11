#-----------------------------------------------------------------------------------------------------------------------

# Title: Construct plots and visualizations for methods paper
# Author: Marcus Becker

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)

# Set preferred ggplot2 theme
ggplot2::theme_set(theme_light())

# Load previously processed data:
# Use `data/processed` - user can change if they wish. This folder is in the .gitignore.
processed <- "./data/processed/"

# Image parameters
height <- 5
width <- 6
dpi <- 500

#-----------------------------------------------------------------------------------------------------------------------

# Figure 1
# Estimated probability of moose leaving the camera field-of-view based on length of interval between sequential images

d <- read_csv("./data/processed/gap-leave-prob_raw-data.csv")

# Moose predictions
d_pred_moose <- d %>%
  filter(common_name == "Moose") %>%
  group_by(common_name) %>%
  nest() %>%
  # Define model for each gap group, and then make predictions
  mutate(model = map(.x = data, ~ smooth.spline(x = .$diff_time, y = .$left, df = 3)),
         pred = map(.x = model, ~ predict(., x = 20:120))) %>%
  select(common_name, pred) %>%
  unnest_wider(pred) %>%
  unnest(cols = c(x, y)) %>%
  rename(diff_time = x, pred = y) %>%
  ungroup()

# Figure
d %>%
  filter(common_name == "Moose") %>%
  ggplot(mapping = aes(x = diff_time, y = left)) +
  geom_jitter(height = 0.05, size = 2, alpha = 0.4, color = "gray20") +
  geom_line(data = d_pred_moose, aes(x = diff_time, y = pred), color = "gray40", size = 1.3) +
  scale_y_continuous(labels = c("Stayed", 0.2, 0.4, 0.6, 0.8, "Left"), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(20, 120, by = 20), limits = c(20, 120)) +
  labs(y = "Probability of leaving field-of-view",
       x = "Length of gap between images (seconds)") +
  theme(axis.title.x = element_text(margin = margin(10, 0, 0, 0), size = 12),
        axis.title.y = element_text(margin = margin(0, -15, 0, 0)),
        axis.text.y = element_text(face = c("bold", "plain", "plain", "plain", "plain", "bold"),
                                   size = c(13, 10, 10, 10, 10, 13)))

# Save figure
ggsave(filename = "./results/figures/Figure 1.png", height = height, width = width, dpi = dpi)

#-----------------------------------------------------------------------------------------------------------------------

# Figure 2
# Distribution of calculated moose densities across the 2,990 cameras.

d <- read_csv("./data/processed/abmi-cmu_all-years_density.csv")

# Only moose, only at ABMI core sites (no CMU, no off-grid)
d_moose_dens <- d %>%
  filter(common_name == "Moose",
         str_detect(project, "ABMI"),
         !str_detect(location, "OG")) %>%
  # Summarise density across the two seasons (summer and winter)
  group_by(location, common_name) %>%
  summarise(density = mean(density_km2, na.rm = TRUE))

# Figure
d_moose_dens %>%
  filter(density < 10) %>%
  ggplot(mapping = aes(x = density)) +
  geom_histogram(fill = "gray60", color = "gray20", bins = 40) +
  coord_cartesian(ylim = c(0, 500)) +
  scale_x_continuous(breaks = seq(-1, 10, 1), expand = c(0.02, 0)) +
  labs(x = expression(Camera~density~(animals~per~km^2)),
       y = "Number of Cameras") +
  theme(legend.position = "none",
        strip.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 13, margin = margin(0, 10, 0, 0)),
        axis.title.x = element_text(size = 13, margin = margin(10, 0, 0, 0)))

# Save figure
ggsave(filename = "./results/figures/Figure 2.png", height = height, width = width, dpi = dpi)

#-----------------------------------------------------------------------------------------------------------------------

# Figure 3
# Moose densities at deployments located in typical treed vs open microhabitats in both deciduous and coniferous stands

# Make a simple dataframe of the relevant data for Moose
d <- data.frame(common_name = c("Moose", "Moose", "Moose", "Moose", "Moose", "Moose"),
                stand_type = c("Deciduous", "Deciduous", "Deciduous", "Conifer", "Conifer", "Conifer"),
                opening = c("Treed", "Open\n(Low)", "Open\n(Productive)", "Treed", "Open\n(Low)", "Open\n(Productive)"),
                density = c(0.581, 0.279, 0.94, 0.21, 0.162, 0.755),
                lci = c(0.283, 0.072, 0.445, 0.087, 0.076, 0.351),
                uci = c(0.948, 0.592, 1.553, 0.36, 0.288, 1.257))

# Labels
labels <- as_labeller(c(`Deciduous` = "(A) Deciduous Forest",
                        `Conifer` = "(B) Coniferous Forest"))

# Figure
d %>%
  mutate(opening = factor(opening, levels = c("Treed", "Open\n(Low)", "Open\n(Productive)")),
         stand_type = factor(stand_type, levels = c("Deciduous", "Conifer"))) %>%
  ggplot(aes(y = density, x = opening)) +
  geom_point(size = 5, shape = 18) +
  geom_errorbar(aes(ymax = uci, ymin = lci), width = 0.21, size = 0.7) +
  facet_wrap(~ stand_type, labeller = labels) +
  scale_y_continuous(breaks = seq(0, 1.6, 0.2), limits = c(0, 1.6)) +
  labs(x = "", y = expression(Moose~Density~(animals~per~km^2))) +
  theme(strip.text = element_text(color = "black", size = 12),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0), size = 12),
        axis.text.x = element_text(size = 11))

# Save figure
ggsave(filename = "./results/figures/Figure 3.png", height = height, width = width, dpi = dpi)

#-----------------------------------------------------------------------------------------------------------------------

# Figure 4
# Proportion of time that moose spend investigating the camera and pole, and time spent investigating plus associated

d <- read_csv("./data/lookup/behaviour-moose-results.csv")

labels <- as_labeller(c(`invest` = " (A) Investigating Time Only", `investassociate` = "(B) Investigating and Associated Time"))

d %>%
  mutate(estimate = mean,
         mult = 1 / (1 - estimate)) %>%
  ggplot(aes(y = estimate, x = veghf)) +
  geom_point(size = 4, shape = 18) +
  geom_errorbar(aes(ymax = uci, ymin = lci), width = 0.3) +
  facet_wrap(~ analysis, labeller = labels) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1),
                     sec.axis = sec_axis(~ 1 / (1 - .), breaks = c(1, 2, 3, 5, 10, 100), name = "Density Factor")) +
  labs(y = "Proportion of Time",
       x = "Habitat Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.text = element_text(color = "black", size = 10),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0), size = 12),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0), size = 12))

# Save figure
ggsave(filename = "./results/figures/Figure 4.png", height = height, width = width, dpi = dpi)

#-----------------------------------------------------------------------------------------------------------------------

# Figures 6, 7, & 8 - WMU application

# Load WMU model and raw data
load("./results/models/wmu_application_model_output.rdata")

# Make predictions for initial plotting
x <- seq(0, 0.77, 0.01) # 0.77 is the highest aerial density
pred <- predict(m.lm, newdata = data.frame(aerial_avg = x), se.fit = TRUE)
pred.line <- data.frame(x, pred$fit)
pred.line.u <- data.frame(x, u = pred$fit + 1.65*pred$se.fit)
pred.line.l <- data.frame(x, l = pred$fit - 1.65*pred$se.fit)

all.pred <- pred.line %>% left_join(pred.line.u) %>% left_join(pred.line.l) %>% mutate(pass = "initial")

# Figure 6
# Relationship between density estimated with cameras and aerial surveys

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
  scale_y_continuous(breaks = seq(0, 3.5, 0.25), limits = c(0, 3.5)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  guides(size = guide_legend(title = "Number of Cameras:")) +
  labs(y = expression(Camera~density~(animals~per~km^2)),
       x = expression(Aerial~density~(animals~per~km^2))) +
  theme(legend.position = c(0.87, 0.81),
        legend.title = element_text(size = 8.5),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0), size = 13),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0), size = 13),
        panel.grid = element_blank(),
        legend.spacing.x = unit(0.1, 'cm'))
  # Warning is due to one aerial estimate with unreported CI.

# Save figure
ggsave(filename = "./results/figures/Figure 6.png", height = height, width = width, dpi = dpi)

# Make predictions for corrected plotting
pred <- predict(m.lm.corr, newdata = data.frame(aerial_avg = x), se.fit = TRUE)
pred.line <- data.frame(x, pred$fit)
pred.line.u <- data.frame(x, u = pred$fit + 1.65*pred$se.fit)
pred.line.l <- data.frame(x, l = pred$fit - 1.65*pred$se.fit)

all.pred <- pred.line %>% left_join(pred.line.u) %>% left_join(pred.line.l) %>% mutate(pass = "corrected")

# Figure 7
# Relationship between density estimated with cameras and aerial surveys, with adjustment made for time investigating

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
  # Using same limits to help readers compare between the two figures
  scale_y_continuous(breaks = seq(0, 3.5, 0.25), limits = c(0, 3.5)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  guides(size = guide_legend(title = "Number of Cameras:")) +
  labs(y = expression(Camera~density~(animals~per~km^2)),
       x = expression(Aerial~density~(animals~per~km^2))) +
  theme(legend.position = c(0.87, 0.81),
        legend.title = element_text(size = 8.5),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0), size = 13),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0), size = 13),
        panel.grid = element_blank(),
        legend.spacing.x = unit(0.1, 'cm'))
  # Warning is due to one aerial estimate with unreported CI.

# Save figure
ggsave(filename = "./results/figures/Figure 7.png", height = height, width = width, dpi = dpi)

#-----------------------------------------------------------------------------------------------------------------------

# Figure 8
# Relationship between the number of cameras in each WMU and the precision of the TIFC moose density estimate

d <- df_sum_dens %>%
  ungroup() %>%
  mutate(Camera = (density_uci_0.9 - density_lci_0.9) / density_avg,
         Aerial = (aerial_uci_0.9 - aerial_lci_0.9) / aerial_avg) %>%
  select(WMUNIT_COD, n_deployments, Aerial, Camera) %>%
  pivot_longer(cols = c(Aerial, Camera), names_to = "Method:", values_to = "width")

d1 <- d %>% filter(`Method:` == "Aerial")
max.aerial <- max(d1$width, na.rm = TRUE)
min.aerial <- min(d1$width, na.rm = TRUE)

# Figure 8
d %>%
  filter(`Method:` == "Camera") %>%
  ggplot(aes(x = n_deployments, y = width)) +
  geom_hline(yintercept = max.aerial + 0.01, color = "grey20") +
  geom_hline(yintercept = min.aerial, color = "grey20") +
  scale_y_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3)) +
  scale_x_continuous(breaks = seq(0, 220, 20), limits = c(0, 220)) +
  annotate(geom = "rect", ymax = max.aerial, ymin = min.aerial, xmax = Inf, xmin = -Inf, fill = "grey", alpha = 0.3) +
  annotate(geom = "text", label = "Range of aerial survey precision", x = 35, y = 0.41, size = 3.25) +
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

ggsave(filename = "./results/figures/Figure 8.png", height = height, width = width, dpi = dpi)

#-----------------------------------------------------------------------------------------------------------------------
