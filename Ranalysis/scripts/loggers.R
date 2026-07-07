library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)
library(tidyr)
library(patchwork)
library(ggbreak)
library(emmeans)
library(ggpattern)
library(anytime)

Sys.setlocale("LC_TIME", "en_US.UTF-8")

setwd("/home/gospozha/haifa/O.patagonica/loggers/data/")

# prepare each dataset from different depths

m10=read.csv("./Licor general code/10m_Cat_MO.csv", header=TRUE, stringsAsFactors = TRUE)
m10$Depth="10"
m10$Depth=as.factor(m10$Depth)
m10$Local_timestamp=as.POSIXct(m10$Unix_Timestamp)
str(m10)

m25=read.csv("./Licor general code/25m_Cat_MO.csv", header=TRUE, stringsAsFactors = TRUE)
m25$Depth="25"
m25$Depth=as.factor(m25$Depth)
m25$Local_timestamp=as.POSIXct(m25$Unix_Timestamp)
str(m25)

m45=read.csv("./Licor general code/45m_Cat_MO.csv", header=TRUE, stringsAsFactors = TRUE)
m45$Depth="45"
m45$Depth=as.factor(m45$Depth)
m45$Local_timestamp=as.POSIXct(m45$Unix_Timestamp)
str(m45)

# merge all depths

m_all=rbind(m10, m25, m45)
str(m_all)

# cat datasets so they start and end at the same date
df <- m_all %>%
  filter(Local_timestamp >= "2024-06-20 11:19:00",
         Local_timestamp <="2024-09-29 13:39:00") %>%
  dplyr::select(Local_timestamp, Depth, Temperature, PAR) %>%
  mutate(
    datetime = ymd_hms(Local_timestamp),
    date = as.Date(datetime)
  )


# Temperature - as one line

# Temperature: daily min and max per depth
daily_temp <- df %>%
  group_by(date, Depth) %>%
  summarise(
    min_temp = min(Temperature, na.rm = TRUE),
    max_temp = max(Temperature, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # reshape to long format for plotting as one zig-zag line
  tidyr::pivot_longer(cols = c(min_temp, max_temp),
                      names_to = "stat",
                      values_to = "temp") %>%
  arrange(Depth, date, stat)

# PAR: midday average per depth, only until July 20
# Make sure datetime and date are correctly formatted
df <- df %>%
  mutate(
    datetime = as.POSIXct(datetime),
    date = as.Date(date),
    time = format(datetime, "%H:%M:%S")
  )

# daily_par <- df %>%
#   filter(date <= as.Date("2024-07-20")) %>%
#   group_by(date, Depth) %>%
#   summarise(
#     max_PAR = max(PAR, na.rm = TRUE),
#     .groups = "drop"
#   )
# PAR: daily midday average per depth, only until July 20
daily_par_midday <- df %>%
  filter(
    date <= as.Date("2024-07-20"),
    time >= "12:00:00",
    time <= "14:00:00"
  ) %>%
  group_by(date, Depth) %>%
  summarise(
    mean_midday_PAR = mean(PAR, na.rm = TRUE),
    n_measurements = sum(!is.na(PAR)),
    .groups = "drop"
  )

# Temp plot
temp_plot <- ggplot(daily_temp, aes(x = date, y = temp, color = Depth, group = Depth)) +
  geom_line(size = 0.4) +
  scale_x_date(date_labels = "%d/%m", date_breaks = "14 days") +
  labs(
    #title = "Daily temperature",
    x = "Date (DD/MM)",
    y = "°C",
    color = "Depth"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# PAR plot
# par_plot <- ggplot(daily_par, aes(x = date, y = max_PAR, color = Depth)) +
#   geom_line(size = 0.4) +
#   scale_x_date(date_labels = "%d/%m", date_breaks = "14 days") +
#   labs(
#     title = "Daily maximum PAR",
#     x = "Date (DD/MM)",
#     y = "PAR (µmol/(s·m²))",
#     color = "Depth (m)"
#   ) +
#   theme_minimal(base_size = 10) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
depth_cols <- c(
  "10" = "#DC267F",
  "25" = "#785EF0",
  "45" = "#648FFF"
)

par_plot <- ggplot(daily_par_midday, aes(x = date, y = mean_midday_PAR, color = Depth)) +
  geom_line(linewidth = 0.6) +
  scale_x_date(date_labels = "%d/%m", date_breaks = "14 days") +
  scale_color_manual(values = depth_cols) +
  labs(
    #title = "Midday average PAR",
    x = "Date (DD/MM)",
    y = "µmol m⁻² s⁻¹",
    color = "Depth"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

par_plot
combined <- (par_plot / temp_plot) +
  plot_annotation(tag_levels = "A")

combined

ggsave("combined_loggers.jpg", combined, width = 8, height = 7)


# Temperature: as ribbon
daily_temp <- df %>%
  group_by(date, Depth) %>%
  summarise(
    min_temp = min(Temperature, na.rm = TRUE),
    max_temp = max(Temperature, na.rm = TRUE),
    .groups = "drop"
  )

# Temperature plot (min and max as ribbons )
temp_plot <- ggplot(daily_temp, aes(x = date)) +
  geom_ribbon(aes(ymax = max_temp, ymin = min_temp,fill = Depth), alpha = 0.7) +
  scale_x_date(date_labels = "%d/%m", date_breaks = "14 days") +
  scale_color_manual(values = depth_cols) +
  labs(
    #title = "Daily temperature",
    x = "Date (DD/MM)",
    y = "°C",
    color = "Depth",
    linetype = "Statistic"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ribbons and lines
daily_temp <- df %>%
  group_by(date, Depth) %>%
  summarise(
    min_temp = min(Temperature, na.rm = TRUE),
    max_temp = max(Temperature, na.rm = TRUE),
    .groups = "drop"
  )

temp_plot <- ggplot(daily_temp, aes(x = date)) +
  geom_ribbon(
    aes(ymax = max_temp, ymin = min_temp, fill = Depth),
    alpha = 0.35
  ) +
  geom_line(
    aes(y = max_temp, color = Depth),
    linewidth = 0.4
  ) +
  geom_line(
    aes(y = min_temp, color = Depth),
    linewidth = 0.4
  ) +
  scale_fill_manual(values = depth_cols) +
  scale_color_manual(values = depth_cols) +
  scale_x_date(date_labels = "%d/%m", date_breaks = "14 days") +
  labs(
    x = "Date (DD/MM)",
    y = "°C",
    fill = "Depth",
    color = "Depth"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

temp_plot

combined <- (par_plot / temp_plot) +
  plot_annotation(tag_levels = "A")

combined

ggsave("combined_loggers_midday.jpg", combined, width = 8, height = 7)


#### statistics ####
library(lme4)
library(lmerTest)
library(multcomp)
library(multcompView)
library(broom.mixed)

# Prepare data

df <- df %>%
  mutate(
    datetime = as.POSIXct(datetime),
    date = as.Date(datetime),
    Depth = factor(Depth, levels = c("10", "25", "45")),
    time = format(datetime, "%H:%M:%S")
  )

# Daily temperature summaries

daily_temp_stats <- df %>%
  group_by(date, Depth) %>%
  summarise(
    temp_mean = mean(Temperature, na.rm = TRUE),
    temp_min = min(Temperature, na.rm = TRUE),
    temp_max = max(Temperature, na.rm = TRUE),
    temp_delta = temp_max - temp_min,
    n_measurements = sum(!is.na(Temperature)),
    .groups = "drop"
  )



# Mean ± SE across daily summaries

temp_summary_article <- daily_temp_stats %>%
  group_by(Depth) %>%
  summarise(
    mean_temp = mean(temp_mean, na.rm = TRUE),
    se_temp = sd(temp_mean, na.rm = TRUE) / sqrt(n()),
    
    min_temp = mean(temp_min, na.rm = TRUE),
    se_min_temp = sd(temp_min, na.rm = TRUE) / sqrt(n()),
    
    max_temp = mean(temp_max, na.rm = TRUE),
    se_max_temp = sd(temp_max, na.rm = TRUE) / sqrt(n()),
    
    delta_temp = mean(temp_delta, na.rm = TRUE),
    se_delta_temp = sd(temp_delta, na.rm = TRUE) / sqrt(n()),
    
    n_days = n(),
    .groups = "drop"
  )%>%
  mutate(
    across(
      .cols = where(is.numeric) & !n_days,
      .fns = ~ round(.x, 1)
    )
  )

temp_summary_article


# Export temperature tables

write_csv(temp_summary_article, "temperature_summary_article.csv")


# Temperature mixed models
# Model: daily metric ~ Depth + (1 | date)


m_temp_mean <- lmer(temp_mean ~ Depth + (1 | date), data = daily_temp_stats)
m_temp_delta <- lmer(temp_delta ~ Depth + (1 | date), data = daily_temp_stats)

# ANOVA tables

anova_temp_mean <- anova(m_temp_mean)
anova_temp_delta <- anova(m_temp_delta)

anova_temp_mean
# < 2.2e-16
anova_temp_delta
# < 2.2e-16

# Pairwise comparisons

emm_temp_mean <- emmeans(m_temp_mean, pairwise ~ Depth, adjust = "BH")
emm_temp_delta <- emmeans(m_temp_delta, pairwise ~ Depth, adjust = "BH")

emm_temp_mean
# all signif
emm_temp_delta
# 10-45, 25-45 - 45 is much smoother


# Daily PAR summaries
# Midday PAR only, from 12:00 to 14:00
# Only until 2024-07-20


daily_par_stats <- df %>%
  filter(
    date <= as.Date("2024-07-20"),
    time >= "12:00:00",
    time <= "14:00:00"
  ) %>%
  group_by(date, Depth) %>%
  summarise(
    midday_PAR_mean = mean(PAR, na.rm = TRUE),
    midday_PAR_min = min(PAR, na.rm = TRUE),
    midday_PAR_max = max(PAR, na.rm = TRUE),
    midday_PAR_delta = midday_PAR_max - midday_PAR_min,
    n_measurements = sum(!is.na(PAR)),
    .groups = "drop"
  )



# Mean ± SE across daily midday summaries


par_summary_article <- daily_par_stats %>%
  group_by(Depth) %>%
  summarise(
    mean_midday_PAR = mean(midday_PAR_mean, na.rm = TRUE),
    se_midday_PAR = sd(midday_PAR_mean, na.rm = TRUE) / sqrt(n()),
    
    min_midday_PAR = mean(midday_PAR_min, na.rm = TRUE),
    se_min_midday_PAR = sd(midday_PAR_min, na.rm = TRUE) / sqrt(n()),
    
    max_midday_PAR = mean(midday_PAR_max, na.rm = TRUE),
    se_max_midday_PAR = sd(midday_PAR_max, na.rm = TRUE) / sqrt(n()),
    
    delta_midday_PAR = mean(midday_PAR_delta, na.rm = TRUE),
    se_delta_midday_PAR = sd(midday_PAR_delta, na.rm = TRUE) / sqrt(n()),
    
    n_days = n(),
    .groups = "drop"
  )%>%
  mutate(
    across(
      .cols = where(is.numeric) & !n_days,
      .fns = ~ round(.x, 1)
    )
  )

par_summary_article


# Export PAR tables


write_csv(par_summary_article, "PAR_summary_article.csv")


# Log-transform PAR variables


daily_par_stats <- daily_par_stats %>%
  mutate(
    log_mean_midday_PAR = log10(midday_PAR_mean + 1),
    log_max_midday_PAR = log10(midday_PAR_max + 1),
    log_delta_midday_PAR = log10(midday_PAR_delta + 1)
  )


# PAR mixed models
# Model: daily PAR metric ~ Depth + (1 | date)


m_par_mean <- lmer(log_mean_midday_PAR ~ Depth + (1 | date), data = daily_par_stats)
m_par_delta <- lmer(log_delta_midday_PAR ~ Depth + (1 | date), data = daily_par_stats)

# ANOVA tables

anova_par_mean <- anova(m_par_mean)
anova_par_delta <- anova(m_par_delta)

anova_par_mean
# < 2.2e-16 
anova_par_delta
# < 2.2e-16 
# Pairwise comparisons

emm_par_mean <- emmeans(m_par_mean, pairwise ~ Depth, adjust = "BH")
emm_par_delta <- emmeans(m_par_delta, pairwise ~ Depth, adjust = "BH")

emm_par_mean
# all signif
emm_par_delta
# all signif
