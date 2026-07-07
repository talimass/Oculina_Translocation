library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(car)
library(cowplot)
library(gridExtra)
library(multcompView)
library(tidyverse)
library(patchwork)
library(effsize)
library(lubridate)
library(lme4)
library(lmerTest)
library(fitdistrplus)
library(emmeans)
setwd("/home/gospozha/haifa/O.patagonica/physio/")
load(("physio.RData"))

#### physiology ####
# opening file, checking, replacing characters to vectors when needed
physio <- read.csv('physiological_data.csv', stringsAsFactors = F)
str(physio)
#View(physio)
physio$Depth = as.factor(physio$Depth)
physio$Depth <- factor(
  physio$Depth,
  levels = c("10", "25", "45"),
  labels = c("T10", "T25", "T45")
)
physio <- na.omit(physio)
physio$Colony = as.factor(physio$Colony)
# adding the same theme to each plot
mytheme = theme_bw()+
  theme(plot.title = element_text(size = 10), legend.position = "none",
        axis.text = element_text(colour = "black", size = 8), axis.title = element_text(size = 9))


#### monthly growth rate ####

physio <- physio %>%
  mutate(
    start_date = case_when(
      Depth == "T10" ~ dmy("08.05.2023"),
      Depth == "T25" ~ dmy("11.02.2023"),
      Depth == "T45" ~ dmy("01.06.2023")
    ),
    end_date = case_when(
      Depth == "T10" ~ dmy("28.09.2023"),
      Depth == "T25" ~ dmy("17.10.2023"),
      Depth == "T45" ~ dmy("18.12.2023")
    ),
    experiment_days = as.numeric(end_date - start_date),
    experiment_months = experiment_days / 30.4375,
    growth_cm2_month = growrth_cm.2 / experiment_months
  )

physio %>%
  distinct(Depth, start_date, end_date, experiment_days, experiment_months)

physio.growth <- physio %>%
  filter(growth_cm2_month > 0)

shapiro.test(physio.growth$growth_cm2_month)
leveneTest(growth_cm2_month ~ Depth, data = physio.growth)

shapiro.test(sqrt(physio.growth$growth_cm2_month))
leveneTest(sqrt(growth_cm2_month) ~ Depth, data = physio.growth)

model <- lmer(
  sqrt(growth_cm2_month) ~ Depth + (1 | Colony),
  data = physio.growth
)

summary(model)
anova(model)
nobs(model)

plot(model)
qqnorm(resid(model))
qqline(resid(model))


growth <- ggplot(physio.growth, aes(y = growth_cm2_month, x = Depth)) +
  geom_boxplot(
    outlier.shape = NA,
    aes(fill = Depth),
    fatten = 0.5,
    alpha = 0.7,
    lwd = 0.4
  ) +
  geom_jitter(width = 0.2, size = 0.8) +
  labs(
    x = "Depth (m)",
    y = expression(cm^2~month^-1),
    title = "Monthly growth rate"
  ) +
  mytheme

growth
#### protein ####

shapiro.test(physio$protein_ug_ml) #  normal  
leveneTest(protein_ug_ml~Depth,d=physio) # variance is ok

ggplot(physio, aes(x = protein_ug_ml, fill=Depth, alpha=0.1)) +
  geom_density() +
  mytheme

# anova + Tukey post-hoc 
#model <- lm(data = physio, protein_ug_ml ~ Depth)
model <- lmer(protein_ug_ml ~ Depth + (1 | Colony), data = physio)
summary(model)  
#plot(model) # Q-Q plot is ok
anova(model) 
nobs(model)
qqnorm(resid(model))
qqline(resid(model))




# posthoc test
posthoc <- emmeans(model, pairwise ~ Depth, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
Names <- gsub("Depth", "", Names)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Depth <- rownames(letters.df) 
placement <- physio %>% 
  group_by(Depth) %>%
  summarise(quantile(protein_ug_ml, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 
# boxplot with letters

protein <- ggplot(physio, aes(y = (protein_ug_ml), x = Depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  labs(x = "Depth (m)", y= "µg mL⁻¹", title = "Coral host protein concentration")+
  mytheme +
  geom_text(data = letters.df, aes(x = Depth, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

protein


#### cell/cm2 ####
# normality, homoscedasticity
shapiro.test(physio$zoox_cm.2) #  normal  
leveneTest(zoox_cm.2~Depth,d=physio) # variance is ok

ggplot(physio, aes(x = (zoox_cm.2), fill=Depth, alpha=0.1)) +
  geom_density() +
  mytheme

# anova  
#model <- lm(data = physio, (zoox_cm.2) ~ Depth)
model <- lmer(zoox_cm.2 ~ Depth + (1 | Colony), data = physio)
summary(model)  
#plot(model) # Q-Q plot is ok
anova(model) 
nobs(model)
qqnorm(resid(model))
qqline(resid(model))


# boxplot 

cellcm <- ggplot(physio, aes(y = (zoox_cm.2), x = Depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  labs(x = "Depth (m)", y= ~cells ~x~cm^-2, title = "Symbiont cell count per surface area")+
  mytheme 

cellcm


#### chl/ug ####
# normality, homoscedasticity
shapiro.test((physio$chlorophyl_ug_algea)) #  not normal  
leveneTest(chlorophyl_ug_algea~Depth,d=physio) # variance is ok

shapiro.test(sqrt(physio$chlorophyl_ug_algea)) #  not normal  

ggplot(physio, aes(x = (chlorophyl_ug_algea), fill=Depth, alpha=0.1)) +
  geom_density() +
  mytheme
max(physio$chlorophyl_ug_algea)

physio_no_outlier <- physio %>%
  filter(chlorophyl_ug_algea != max(chlorophyl_ug_algea, na.rm = TRUE))
shapiro.test((physio_no_outlier$chlorophyl_ug_algea)) #  not normal  
leveneTest(chlorophyl_ug_algea~Depth,d=physio_no_outlier) # variance is ok

model <- lmer(sqrt(chlorophyl_ug_algea) ~ Depth + (1 | Colony), data = physio_no_outlier)
summary(model)  
#plot(model) # Q-Q plot is ok
anova(model) 
nobs(model)
qqnorm(resid(model))
qqline(resid(model))
# boxplot

chl <- ggplot(physio_no_outlier, aes(y = chlorophyl_ug_algea, x = Depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  labs(x = "Depth (m)", title="Chlorophyll per ug algae", y="~chlorophyll[a] ~x ~µg⁻¹") +
  mytheme 
 
chl


#### photophysiology ####

fire.data <- read.csv("Annotated_FIRe_data.csv")
fire.data$Depth = factor(fire.data$depth, levels = c("4", "4_October", "10", "25"),
                         labels = c("I3", "C3", "T10", "T25"))

fire.data$Colony <- as.factor(fire.data$Colony)
fire.data <- na.omit(fire.data)

mycolors <- c('#C77CFF',  '#F8766D', '#00BA38',"#00BFC4"  )

# reduce dataset 
fire <- fire.data %>%
  drop_na() %>%
  group_by(Depth, Colony) %>%
  summarise(
    `Fv.Fm` = mean(`fv_fm`, na.rm = TRUE),
    Sigma = mean(Sigma, na.rm = TRUE),
    `Pmax.e.s` = mean(`Pmax.e.s`, na.rm = TRUE),
    p = mean(p, na.rm = TRUE),
    .groups = "drop"  # optional: ungroups the result
  )

#### Functional absorption cross-section of PSII ####

# checking for normality and homoscedasticity
shapiro.test((fire$Sigma))  # data is normal
leveneTest(Sigma~Depth, d=fire) # no heteroscedasticity of variance

ggplot(fire, aes(x = (Sigma), fill=Depth, alpha=0.1)) +
  geom_density() +
  mytheme

model <- lmer((Sigma) ~ Depth + (1 | Colony), data = fire)
summary(model)  
#plot(model) # Q-Q plot is ok
anova(model) 
nobs(model)
qqnorm(resid(model))
qqline(resid(model))


# boxplot

sigma <- ggplot(fire, aes(y = (Sigma), x = Depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  scale_fill_manual(values = mycolors )+
  labs(title="Functional absorption cross-section of PSII", y= "σPSII’(A2)", x= "Depth (m)") +
  mytheme
sigma


#### Quantum yield of photochemistry in PSII ####

# checking for normality and homoscedasticity

shapiro.test((fire$Fv.Fm))  # data is not normal (p.value < 0.05)
leveneTest(Fv.Fm~Depth, d=fire) # heteroscedasticity is ok


model <- lmer((Fv.Fm) ~ Depth + (1 | Colony), data = fire)
summary(model)  
#plot(model) # Q-Q plot is ok
anova(model) 
nobs(model)
qqnorm(resid(model))
qqline(resid(model))

# posthoc test
posthoc <- emmeans(model, pairwise ~ Depth, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
Names <- gsub("Depth", "", Names)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Depth <- rownames(letters.df) 
placement <- fire %>% 
  group_by(Depth) %>%
  summarise(quantile(Fv.Fm, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters
# groups with the same letter are the same
# groups that are significantly different get different letters

FvFm = ggplot(fire, aes(y = Fv.Fm, x = Depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  scale_fill_manual(values = mycolors )+
  labs(title="Algal quantum yield of photochemistry in PSII", y="Fv’/Fm’", x= "Depth (m)") +
  mytheme +
  geom_text(data = letters.df, aes(x = Depth, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

FvFm

#### Maximum photosynthetic rate ####
# checking for normality and homoscedasticity
shapiro.test(log(fire$Pmax.e.s))  # data is not normal (p.value < 0.05)
leveneTest(log(Pmax.e.s)~Depth, d=fire) # heteroscedasticity is ok


model <- lmer(log(Pmax.e.s) ~ Depth + (1 | Colony), data = fire)
summary(model)  
#plot(model) # Q-Q plot is ok
anova(model) 
nobs(model)
qqnorm(resid(model))
qqline(resid(model))
shapiro.test(resid(model))

# posthoc test
posthoc <- emmeans(model, pairwise ~ Depth, adjust = "tukey")
p_values <- as.data.frame(posthoc$contrasts)
d <- p_values$p.value < 0.05
Names <- gsub(" ", "", p_values$contrast)
Names <- gsub("Depth", "", Names)
names(d) <- Names
# compact letter display
letters <- multcompLetters(d)
letters.df <- data.frame(letters$Letters)
colnames(letters.df)[1] <- "Letter"
letters.df$Depth <- rownames(letters.df) 
placement <- fire %>% 
  group_by(Depth) %>%
  summarise(quantile(Pmax.e.s, na.rm = TRUE)[4])
colnames(placement)[2] <- "Placement.Value"
letters.df <- left_join(letters.df, placement) 

# boxplot with letters
# groups with the same letter are the same
# groups that are significantly different get different letters

Pmax = ggplot(fire, aes(y = Pmax.e.s, x = Depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  scale_fill_manual(values = mycolors )+
  labs(title="Maximum photosynthetic rate", y="electron s⁻¹ PSII⁻¹", x= "Depth (m)") +
  mytheme +
  geom_text(data = letters.df, aes(x = Depth, y = Placement.Value, label = Letter),
            size = 3, color = "black", hjust = -0.5, vjust = -0.8, fontface = "italic")

Pmax


#### connectivity parameter ####
# checking for normality and homoscedasticity
shapiro.test(log(fire$p))  # ok
leveneTest(log(p)~Depth, d=fire) # ok


model <- lmer(log(p) ~ Depth + (1 | Colony), data = fire)
summary(model)  
#plot(model) # Q-Q plot is ok
anova(model) 
nobs(model)
qqnorm(resid(model))
qqline(resid(model))
shapiro.test(resid(model))


p <- ggplot(fire, aes(y = (p), x = Depth)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Depth), fatten = 0.5, alpha = 0.7, lwd = 0.4)+
  geom_jitter(width=0.2, size = 0.8) +
  scale_fill_manual(values = mycolors )+
  labs(title="Connectivity parameter", y= "p") +
  mytheme

p

#### combined plots ####

condition_colors <- c(
  "I3" = "#FFB000",
  "C3"        = "#FE6100",
  "T10"        = "#DC267F",
  "T25"        = "#785EF0",
  "T45" = "#648FFF"
)


combined_phys <- (
  (growth + protein) / (cellcm + chl) / (FvFm + sigma) / (Pmax + p) ) +
  plot_annotation(tag_levels = "A")& 
  theme(plot.tag = element_text(face = "bold", size = 9)) 

combined_phys <- 
  combined_phys &
  theme(plot.title = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        text = element_text(size = 7))

combined_phys <- combined_phys &
  scale_fill_manual(values = condition_colors, name = "Depth")

combined_phys

ggsave(
  "boxplots_combined.pdf",
  plot = combined_phys,
  width = 18,
  height = 19,
  units = "cm",
  device = "pdf",
  useDingbats = FALSE
)
ggsave(
  "boxplots_combined.jpg",
  plot = combined_phys,
  width = 18,
  height = 19,
  units = "cm",
  device = "jpg"
)

save.image("physio.RData")
