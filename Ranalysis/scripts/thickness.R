library(tidyr)
library(ggplot2)
library(dplyr)

setwd("/home/gospozha/haifa/O.patagonica/microCT/")
# Read data
df <- read.csv("skeletal_thickness.csv", check.names = FALSE) %>%
  rename(sample = `Sample number`, thickness = Avg_size_um, sd = sdv_um) %>%
  mutate(site = factor(str_to_title(str_trim(Site)), levels = c("Shallow", "Deep")))

# Summary statistics
stats <- df %>%
  group_by(site) %>%
  summarise(n = n(), mean = mean(thickness), sd = sd(thickness),
            median = median(thickness), min = min(thickness), max = max(thickness),
            .groups = "drop")

stats

# Statistical tests
welch <- t.test(thickness ~ site, data = df)
welch
# t = 0.26614, df = 4.5097, p-value = 0.8019


# Plot: individual colony means ± within-colony SD;
# diamonds show group means ± SD among colonies

p <- ggplot(df, aes(x = site, y = thickness)) +
  # Individual colonies: mean ± within-colony SD
  geom_errorbar(aes(ymin = pmax(0, thickness - sd), ymax = thickness + sd,
                    group = sample),
                width = 0.08, linewidth = 0.5, color = "grey55",
                position = position_dodge(width = 0.55)) +
  geom_point(aes(group = sample), shape = 21, size = 2.8,
             fill = "grey75", color = "grey35",
             position = position_dodge(width = 0.55)) +
  
  # Group mean ± SD among colonies
  geom_errorbar(data = stats, aes(x = site, y = mean,
                                  ymin = pmax(0, mean - sd), ymax = mean + sd),
                inherit.aes = FALSE, width = 0.16,
                linewidth = 0.8, color = "black") +
  geom_point(data = stats, aes(x = site, y = mean),
             inherit.aes = FALSE, shape = 23, size = 3.8,
             fill = "white", color = "black", stroke = 0.9) +
  
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Skeletal thickness (µm)") +
  theme_classic(base_size = 13) +
  theme(axis.text = element_text(color = "black"),
        axis.line = element_line(linewidth = 0.5))

print(p)

# Save 
ggsave("skeletal_thickness.pdf", p, width = 4, height = 4)
ggsave("skeletal_thickness.png", p, width = 4, height = 4, dpi = 600)