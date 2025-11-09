# Episodic thinking: endpoint-focused vs. present-focused
# Behavioral data analysis & viso
# This script requires one file: "behavioral_rawdata.xlsx"
# Programmed by Feng XIAO (updated on 2025.11.5)
############################################################################################################

### Preparation
## Load required packages for analysis
package_list <- c('car','tidyr','dplyr','readxl','effsize','afex','emmeans',
                  'e1071','lmtest','broom',
                  'ggplot2','patchwork')
lapply(package_list, require, character.only = TRUE)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


### Data input
rd_pretest <- read_excel('behavioral_rawdata.xlsx', sheet = 'pretest', na = "---")
rd_posttest <- read_excel('behavioral_rawdata.xlsx', sheet = 'posttest', na = "---")
rd_k <- read_excel('behavioral_rawdata.xlsx', sheet = 'k_value', na = "---")


### Demographic data analysis
## Original sample size
# Age
mean(rd_pretest$Age) #20.93
sd(rd_pretest$Age) #2.09
# Gender
dim(filter(rd_pretest, Gender == 1))[1] #9 males
dim(filter(rd_pretest, Gender == 2))[1] #31 females

## Final sample size (excluding 2 with MDD history)
rd_pretest_filtered <- rd_pretest %>%
  filter(`Mental disease/Long time medicine` == 0)
# Age
mean(rd_pretest_filtered$Age) #20.95
sd(rd_pretest_filtered$Age) #2.14
mean((filter(rd_pretest_filtered, Gender == 1))$Age) #male: 21.89
sd((filter(rd_pretest_filtered, Gender == 1))$Age) #male: 2.26
mean((filter(rd_pretest_filtered, Gender == 2))$Age) #female: 20.66
sd((filter(rd_pretest_filtered, Gender == 2))$Age) #female: 2.06
t.test((filter(rd_pretest_filtered, Gender == 1))$Age, (filter(rd_pretest_filtered, Gender == 2))$Age,
       paired =FALSE, alternative = "two.sided", var.equal=FALSE,
       conf.level=0.95) #ages did not differ within males and females
# Gender
dim(filter(rd_pretest_filtered, Gender == 1))[1] #9 males
dim(filter(rd_pretest_filtered, Gender == 2))[1] #29 females
# Handedness
dim(filter(rd_pretest_filtered, Handedness == 1))[1] #left-hand: 0
dim(filter(rd_pretest_filtered, Handedness == 2))[1] #right-hand: 38
# Meditation experience
dim(filter(rd_pretest_filtered, meditation == 1))[1] #Yes: 10
dim(filter(rd_pretest_filtered, meditation == 2))[1] #No: 28

## Posttest
# Headsize
rd_posttest_filtered <- rd_posttest %>%
  filter(`Mental disease/Long time medicine` == 0)
mean(rd_posttest_filtered$`Headsize (cm)`) #56.35
sd(rd_posttest_filtered$`Headsize (cm)`) #1.70

### Group comparisons
## Involvement level
inv <- c(rd_posttest_filtered$Involvement_P, rd_posttest_filtered$Involvement_E)
mean(inv) #overall: 76.05
sd(inv) #overall: 14.97
mean(rd_posttest_filtered$Involvement_P) #Present: 76.92
sd(rd_posttest_filtered$Involvement_P) #Present: 15.55
mean(rd_posttest_filtered$Involvement_E) #Endpoint: 75.18
sd(rd_posttest_filtered$Involvement_E) #Endpoint: 14.52
t.test(rd_posttest_filtered$Involvement_E,
       rd_posttest_filtered$Involvement_P,
       paired = TRUE, alternative = "two.sided") #NS


## Duration estimation
## Analysis
# Before vs. After (Present) vs. After (Endpoint)
df_time_diff1 <- data.frame(rd_pretest_filtered$Time_diff1,
                            rd_posttest_filtered$TimeP_diff1,
                            rd_posttest_filtered$TimeE_diff1)
colnames(df_time_diff1) <- c('Pretest','Present','Endpoint')
df_time_diff2 <- data.frame(rd_pretest_filtered$Time_diff2,
                            rd_posttest_filtered$TimeP_diff2,
                            rd_posttest_filtered$TimeE_diff2)
colnames(df_time_diff2) <- c('Pretest','Present','Endpoint')
df_time_diff1$id <- 1:nrow(df_time_diff1)
df_time_diff2$id <- 1:nrow(df_time_diff2)
df_time1 <- rbind(df_time_diff1, df_time_diff2)
df_time1 <- na.omit(df_time1) #delete the missing data
df_time2 <- pivot_longer(df_time1,
                         cols = c("Pretest", "Present", "Endpoint"),
                         names_to = "group",
                         values_to = "timing_speed")
# Skewness (less than 3 is acceptable)
skewness(df_time1$Pretest) #-0.08
skewness(df_time1$Present) #-0.16
skewness(df_time1$Endpoint) #-0.66
# Lavene's tests
leveneTest(df_time2$timing_speed, factor(df_time2$group)) #equal sds
# ANOVA (within-group factor: pretest, present, endpoint)
aov_time <- aov_ez(
  id = "id",                          
  dv = "timing_speed",
  within = "group",
  data = df_time2,
  type = 3
)
summary(aov_time)
  #F(2,74)=6.63, p=.002;but Mauchly's test is significant (p = .014)
  # Greenhouse-Geisser results for corection: p=.004 
afex::nice(aov_time, es = "pes") #F(1.65,61.92)=6.63, p=.004, partial etasq=0.15
pairs(emmeans(aov_time, ~ group), adjust = "bonferroni")
cohen.d(df_time1$Pretest, df_time1$Endpoint,
        paired = TRUE,
        hedges.correction = FALSE)
  #Endpoint(before>after): t(37)=3.11, p=.011, d=0.51

##Plotting
dat_mean <- df_time1 %>%
  group_by(id) %>%
  summarise(
    Pretest  = mean(Pretest,  na.rm = TRUE),
    Present  = mean(Present,  na.rm = TRUE),
    Endpoint = mean(Endpoint, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(Pretest, Present, Endpoint),
               names_to = "group", values_to = "timing_speed")
lines_present <- dat_mean %>%
  filter(group %in% c("Pretest", "Present")) %>%
  mutate(Intervention = "Present",
         Time = if_else(group == "Pretest", "Before", "After"))
lines_endpoint <- dat_mean %>%
  filter(group %in% c("Pretest", "Endpoint")) %>%
  mutate(Intervention = "Endpoint",
         Time = if_else(group == "Pretest", "Before", "After"))
df_lines <- bind_rows(lines_present, lines_endpoint) %>%
  mutate(
    Time = factor(Time, levels = c("Before", "After")),
    Intervention = factor(Intervention, levels = c("Endpoint", "Present"))
  )
df_sum <- df_lines %>%
  group_by(Intervention, Time) %>%
  summarise(
    mean = mean(timing_speed, na.rm = TRUE),
    se   = sd(timing_speed, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
pd <- position_dodge(width = 0.3)

p_time<-ggplot(df_sum, aes(x = Time, y = mean,
                           color = Intervention, group = Intervention)) +
  geom_line(linewidth = 0.45, position = pd) +
  geom_point(aes(shape = Intervention),
             size = 1.0, stroke = 0.7, position = pd) +   
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.10, linewidth = 0.45, position = pd) +
  scale_x_discrete(limits = c("Before", "After"),
                   labels = c("Before", "After")) +
  scale_y_continuous(limits = c(-4, 0),
                     breaks = seq(-4, 0, by = 1),
                     expand = c(0, 0)) +
  scale_color_manual(values = c("Endpoint" = "#B22222", "Present" = "#4169E1")) +
  scale_shape_manual(values = c("Endpoint" = 1, "Present" = 2)) + 
  labs(x = NULL, y = "Duration estimation difference (s)", title = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.35),
    axis.title       = element_text(size = 7, colour = "black"),
    axis.text        = element_text(size = 7, colour = "black"),
    panel.background = element_rect(fill = "transparent"),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.position  = "none",
    plot.margin      = ggplot2::margin(3, 3, 2, 2)
  ) +
  patchwork::plot_annotation(
    title = "(A) Interval-timing speed",
    theme = theme(plot.title = element_text(size = 7, colour = "black",
                                            face = "bold",
                                            margin = ggplot2::margin(b = 2)))
  )
ggsave("pic_time.pdf", plot = p_time, width = 2, height = 2, units = "in", device = cairo_pdf)


## Discounting rates (k-values)
## Analysis
df_k <- data.frame(rd_k$Pretest, rd_k$Present, rd_k$Endpoint)
colnames(df_k) <- c('Pretest','Present','Endpoint')
df_logk <- log(df_k) #log k values
df_logk$id <- 1:nrow(df_logk)
df_logk1 <- pivot_longer(df_logk, cols = c("Pretest", "Present", "Endpoint"),
                         names_to = "group",
                         values_to = "discounting_rate")
# Skewness (less than 3 is acceptable)
skewness(df_logk$Pretest) #-0.62
skewness(df_logk$Present) #0.93
skewness(df_logk$Endpoint) #0.90
# Lavene's tests
leveneTest(df_logk1$discounting_rate, factor(df_logk1$group)) #equal sds
# ANOVA (within-group factor: pretest, present, endpoint)
aov_logk <- aov_ez(
  id = "id",                          
  dv = "discounting_rate",
  within = "group",
  data = df_logk1,
  type = 3
)
summary(aov_logk)
  #F(2,74)=6.84, p=.002
afex::nice(aov_logk, es = "pes") #partial etasq=0.16.
pairs(emmeans(aov_logk, ~ group), adjust = "bonferroni")
cohen.d(df_logk$Pretest, df_logk$Endpoint,
        paired = TRUE,
        hedges.correction = FALSE)
  #Endpoint(before<after): t(37)=3.29, p=.007, d=0.45
cohen.d(df_logk$Pretest, df_logk$Present,
        paired = TRUE,
        hedges.correction = FALSE)
  #Present(before<after): t(37)=3.27, p=.007, d=0.50

##Plotting
dat_mean <- df_logk %>%
  group_by(id) %>%
  summarise(
    Pretest  = mean(Pretest,  na.rm = TRUE),
    Present  = mean(Present,  na.rm = TRUE),
    Endpoint = mean(Endpoint, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(Pretest, Present, Endpoint),
               names_to = "group", values_to = "discounting_rate")
lines_present <- dat_mean %>%
  filter(group %in% c("Pretest", "Present")) %>%
  mutate(Intervention = "Present",
         Time = if_else(group == "Pretest", "Before", "After"))
lines_endpoint <- dat_mean %>%
  filter(group %in% c("Pretest", "Endpoint")) %>%
  mutate(Intervention = "Endpoint",
         Time = if_else(group == "Pretest", "Before", "After"))
df_lines <- bind_rows(lines_present, lines_endpoint) %>%
  mutate(
    Time = factor(Time, levels = c("Before", "After")),
    Intervention = factor(Intervention, levels = c("Endpoint", "Present"))
  )
df_sum <- df_lines %>%
  group_by(Intervention, Time) %>%
  summarise(
    mean = mean(discounting_rate, na.rm = TRUE),
    se   = sd(discounting_rate, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_logk<-ggplot(df_sum, aes(x = Time, y = mean,
                           color = Intervention, group = Intervention)) +
  geom_line(linewidth = 0.45, position = pd) +
  geom_point(aes(shape = Intervention),
             size = 1.0, stroke = 0.7, position = pd) +   
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.10, linewidth = 0.45, position = pd) +
  scale_x_discrete(limits = c("Before", "After"),
                   labels = c("Before", "After")) +
  scale_y_continuous(limits = c(-8, 0),
                     breaks = seq(-8, 0, by = 2),
                     expand = c(0, 0)) +
  scale_color_manual(values = c("Endpoint" = "#B22222", "Present" = "#4169E1")) +
  scale_shape_manual(values = c("Endpoint" = 1, "Present" = 2)) + 
  labs(x = NULL, y = "Log k-value", title = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.35),
    axis.title       = element_text(size = 7, colour = "black"),
    axis.text        = element_text(size = 7, colour = "black"),
    panel.background = element_rect(fill = "transparent"),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.position  = "none",
    plot.margin      = ggplot2::margin(3, 3, 2, 2)
  ) +
  patchwork::plot_annotation(
    title = "(B) Delay discounting rate",
    theme = theme(plot.title = element_text(size = 7, colour = "black",
                                            face = "bold",
                                            margin = ggplot2::margin(b = 2)))
  )
ggsave("pic_logk.pdf", plot = p_logk, width = 2, height = 2, units = "in", device = cairo_pdf)


## Monetary allocation task
## Analysis
if (!"id" %in% names(rd_pretest_filtered))  rd_pretest_filtered$id  <- seq_len(nrow(rd_pretest_filtered))
if (!"id" %in% names(rd_posttest_filtered)) rd_posttest_filtered$id <- rd_pretest_filtered$id
# Pretest
pre_DS  <- rd_pretest_filtered$Prop_ShortSaving - rd_pretest_filtered$Prop_ShortSpend
pre_DL  <- rd_pretest_filtered$Prop_LongSaving  - rd_pretest_filtered$Prop_ShortSpend
pre_DSL <- rd_pretest_filtered$Prop_LongSaving  - rd_pretest_filtered$Prop_ShortSaving
# Present
pres_DS  <- rd_posttest_filtered$PropP_ShortSaving - rd_posttest_filtered$PropP_ShortSpend
pres_DL  <- rd_posttest_filtered$PropP_LongSaving  - rd_posttest_filtered$PropP_ShortSpend
pres_DSL <- rd_posttest_filtered$PropP_LongSaving  - rd_posttest_filtered$PropP_ShortSaving
# Endpoint
end_DS  <- rd_posttest_filtered$PropE_ShortSaving - rd_posttest_filtered$PropE_ShortSpend
end_DL  <- rd_posttest_filtered$PropE_LongSaving  - rd_posttest_filtered$PropE_ShortSpend
end_DSL <- rd_posttest_filtered$PropE_LongSaving  - rd_posttest_filtered$PropE_ShortSaving

mk_long <- function(id, pre, pres, end, metric_name){
  tibble(
    id = id,
    Pretest  = pre,
    Present  = pres,
    Endpoint = end
  ) |>
    drop_na() |>
    pivot_longer(c(Pretest, Present, Endpoint),
                 names_to = "Session", values_to = "Score") |>
    mutate(Session = factor(Session, levels = c("Pretest","Present","Endpoint")),
           Metric  = metric_name)
}

df_DS_long  <- mk_long(rd_pretest_filtered$id, pre_DS,  pres_DS,  end_DS,  "DS")
df_DL_long  <- mk_long(rd_pretest_filtered$id, pre_DL,  pres_DL,  end_DL,  "DL")
df_DSL_long <- mk_long(rd_pretest_filtered$id, pre_DSL, pres_DSL, end_DSL, "DSL")

## Monetary allocation task ¡ª DS
## DS = short-term saving ??? short-term spending
df_DS_wide <- tidyr::pivot_wider(
  data  = as.data.frame(df_DS_long),
  id_cols = id,
  names_from  = Session,
  values_from = Score
)
# Skewness (less than 3 is acceptable)
skewness(df_DS_wide$Pretest,  na.rm = TRUE) #-2.13
skewness(df_DS_wide$Present,  na.rm = TRUE) #-0.97
skewness(df_DS_wide$Endpoint, na.rm = TRUE) #-0.44
# Levene's tests
leveneTest(df_DS_long$Score, df_DS_long$Session) #equal sds
# ANOVA (within-group factor: pretest, present, endpoint)
aov_DS <- aov_ez(
  id     = "id",
  dv     = "Score",
  within = "Session",
  data   = df_DS_long,
  type   = 3
)
summary(aov_DS) #NS; but short-term saving>spending: F(1,37)=17.27, p<.001
# partial etasq for the intercept effect
F_intercept <- 17.27
df_error <- 37
eta_sq <- F_intercept / (F_intercept + df_error)
eta_sq #0.32

## Monetary allocation task ¡ª DL
## DL = long-term saving ??? short-term spending
df_DL_wide <- tidyr::pivot_wider(
  data  = as.data.frame(df_DL_long),
  id_cols = id,
  names_from  = Session,
  values_from = Score
)
# Skewness
skewness(df_DL_wide$Pretest,  na.rm = TRUE) #-1.50
skewness(df_DL_wide$Present,  na.rm = TRUE) #-1.16
skewness(df_DL_wide$Endpoint, na.rm = TRUE) #-0.81
# Levene's tests
leveneTest(df_DL_long$Score, df_DL_long$Session) #equal sds
# ANOVA (within-group factor: pretest, present, endpoint)
aov_DL <- aov_ez(
  id     = "id",
  dv     = "Score",
  within = "Session",
  data   = df_DL_long,
  type   = 3
)
summary(aov_DL) #NS,but long-term saving>spending: F(1,37)=13.38, p=.001
# partial etasq for the intercept effect
F_intercept <- 13.38
df_error <- 37
eta_sq <- F_intercept / (F_intercept + df_error)
eta_sq #0.27

## Monetary allocation task ¡ª DSL
## DSL = long-term saving ??? short-term saving
df_DSL_wide <-  tidyr::pivot_wider(
  data  = as.data.frame(df_DSL_long),
  id_cols = id,
  names_from  = Session,
  values_from = Score
)
# Skewness
skewness(df_DSL_wide$Pretest,  na.rm = TRUE) #-0.22
skewness(df_DSL_wide$Present,  na.rm = TRUE) #-0.34
skewness(df_DSL_wide$Endpoint, na.rm = TRUE) #-0.25
# Levene's tests
leveneTest(df_DSL_long$Score, df_DSL_long$Session) #equal sds
# ANOVA (within-group factor: pretest, present, endpoint)
aov_DSL <- aov_ez(
  id     = "id",
  dv     = "Score",
  within = "Session",
  data   = df_DSL_long,
  type   = 3
)
summary(aov_DSL) #NS

## Plotting
if (!"id" %in% names(rd_pretest_filtered))  rd_pretest_filtered$id  <- seq_len(nrow(rd_pretest_filtered))
if (!"id" %in% names(rd_posttest_filtered)) rd_posttest_filtered$id <- rd_pretest_filtered$id
pre_long <- rd_pretest_filtered %>%
  transmute(
    id,
    Session  = "Before",
    ShortSpend  = Prop_ShortSpend,
    ShortSaving = Prop_ShortSaving,
    LongSaving  = Prop_LongSaving
  ) %>%
  pivot_longer(cols = c(ShortSpend, ShortSaving, LongSaving),
               names_to = "Category", values_to = "Prop")
end_long <- rd_posttest_filtered %>%
  transmute(
    id,
    Session  = "Endpoint-focused",
    ShortSpend  = PropE_ShortSpend,
    ShortSaving = PropE_ShortSaving,
    LongSaving  = PropE_LongSaving
  ) %>%
  pivot_longer(cols = c(ShortSpend, ShortSaving, LongSaving),
               names_to = "Category", values_to = "Prop")
pres_long <- rd_posttest_filtered %>%
  transmute(
    id,
    Session  = "Present-focused",
    ShortSpend  = PropP_ShortSpend,
    ShortSaving = PropP_ShortSaving,
    LongSaving  = PropP_LongSaving
  ) %>%
  pivot_longer(cols = c(ShortSpend, ShortSaving, LongSaving),
               names_to = "Category", values_to = "Prop")
alloc_long <- bind_rows(pre_long, end_long, pres_long) %>%
  mutate(
    Session  = factor(Session, levels = c("Before","Endpoint-focused","Present-focused")),
    Category = factor(Category, levels = c("ShortSpend","ShortSaving","LongSaving"),
                      labels = c("Short-term spending","Short-term saving","Long-term saving"))
  )
alloc_sum <- alloc_long %>%
  group_by(Session, Category) %>%
  summarise(
    mean_pct = mean(Prop, na.rm = TRUE) * 100,
    se_pct   = (sd(Prop, na.rm = TRUE) / sqrt(sum(!is.na(Prop)))) * 100,
    .groups  = "drop"
  )
pd <- position_dodge(width = 0.15)

p_ma<-ggplot(alloc_sum, aes(x = Session, y = mean_pct,
                      color = Category, shape = Category,
                      group = Category)) +
  geom_line(linewidth = 0.5, position = pd) +
  geom_point(size = 1.0, stroke = 0.7, position = pd) +
  geom_errorbar(aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct),
                width = 0.10, linewidth = 0.45, position = pd) +
  scale_color_manual(values = c(
    "Short-term spending" = "#B22222",
    "Short-term saving"   = "#4169E1",
    "Long-term saving"    = "#6F4F28"
  )) +
  scale_shape_manual(values = c(
    "Short-term spending" = 21,  
    "Short-term saving"   = 22,  
    "Long-term saving"    = 24   
  )) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 25),
                     labels = function(x) paste0(x, "%"),
                     expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Allocation proportion (%)",
    color = NULL,
    title = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.35),
    axis.title       = element_text(size = 7, colour = "black"),
    axis.text        = element_text(size = 7, colour = "black"),
    panel.background = element_rect(fill = "transparent"),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.position  = "none",
    plot.margin      = ggplot2::margin(3, 3, 2, 2)
  ) +
  patchwork::plot_annotation(
    title = "(C) Monetary allocation strategy",
    theme = theme(plot.title = element_text(size = 7, colour = "black",
                                            face = "bold",
                                            margin = ggplot2::margin(b = 2)))
  )
ggsave("pic_ma.pdf", plot = p_ma, width = 2, height = 2, units = "in", device = cairo_pdf)


## Regression: examine whether logk changes predicts allocation changes
if (!"id" %in% names(rd_pretest_filtered))  rd_pretest_filtered$id  <- seq_len(nrow(rd_pretest_filtered))
if (!"id" %in% names(rd_posttest_filtered)) rd_posttest_filtered$id <- rd_pretest_filtered$id
if (!"id" %in% names(rd_k))                 rd_k$id                 <- rd_pretest_filtered$id

df_k <- rd_k %>%
  transmute(id,
            Pretest  = Pretest,
            Present  = Present,
            Endpoint = Endpoint) %>%
  mutate(across(c(Pretest, Present, Endpoint), log, .names = "log_{.col}"))

# Allocation changes calculation: endpoint and present
delta_present <- rd_pretest_filtered %>%
  transmute(id,
            dP_ShortSpend  = rd_posttest_filtered$PropP_ShortSpend  - Prop_ShortSpend,
            dP_ShortSaving = rd_posttest_filtered$PropP_ShortSaving - Prop_ShortSaving,
            dP_LongSaving  = rd_posttest_filtered$PropP_LongSaving  - Prop_LongSaving) %>%
  left_join(df_k %>% transmute(id, dP_logk = log_Present - log_Pretest),
            by = "id")
delta_endpoint <- rd_pretest_filtered %>%
  transmute(id,
            dE_ShortSpend  = rd_posttest_filtered$PropE_ShortSpend  - Prop_ShortSpend,
            dE_ShortSaving = rd_posttest_filtered$PropE_ShortSaving - Prop_ShortSaving,
            dE_LongSaving  = rd_posttest_filtered$PropE_LongSaving  - Prop_LongSaving) %>%
  left_join(df_k %>% transmute(id, dE_logk = log_Endpoint - log_Pretest),
            by = "id")

# Standardization for regression
std_lm <- function(y, x) {
  dat <- na.omit(data.frame(y = scale(y), x = scale(x)))
  lm(y ~ x, data = dat)
}
fit_P_SS_std  <- std_lm(delta_present$dP_ShortSpend,  delta_present$dP_logk)
fit_P_SSh_std <- std_lm(delta_present$dP_ShortSaving, delta_present$dP_logk)
fit_P_LS_std  <- std_lm(delta_present$dP_LongSaving,  delta_present$dP_logk)
fit_E_SS_std  <- std_lm(delta_endpoint$dE_ShortSpend,  delta_endpoint$dE_logk)
fit_E_SSh_std <- std_lm(delta_endpoint$dE_ShortSaving, delta_endpoint$dE_logk)
fit_E_LS_std  <- std_lm(delta_endpoint$dE_LongSaving,  delta_endpoint$dE_logk)

summarize_fit <- function(fit, outcome, group){
  broom::tidy(fit) %>%
    filter(term == "x") %>%
    transmute(Group = group,
              Outcome = outcome,
              Beta = estimate,
              SE   = std.error,
              t    = statistic,
              p    = p.value)
}

# Benjamini¨CHochberg corrections
tab_std <- dplyr::bind_rows(
  summarize_fit(fit_P_SS_std,  "Short-term spending", "Present"),
  summarize_fit(fit_P_SSh_std, "Short-term saving",   "Present"),
  summarize_fit(fit_P_LS_std,  "Long-term saving",    "Present"),
  summarize_fit(fit_E_SS_std,  "Short-term spending", "Endpoint"),
  summarize_fit(fit_E_SSh_std, "Short-term saving",   "Endpoint"),
  summarize_fit(fit_E_LS_std,  "Long-term saving",    "Endpoint")
) |>
  dplyr::group_by(Group) |>
  dplyr::mutate(
    p_bh   = p.adjust(p, method = "BH")   
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(Group, Outcome)

# Calculation of Rsq and fsq
r2_table <- tibble::tibble(
  Group = c(rep("Present", 3), rep("Endpoint", 3)),
  Outcome = rep(c("Short-term spending","Short-term saving","Long-term saving"), 2),
  R2 = c(summary(fit_P_SS_std)$r.squared,
         summary(fit_P_SSh_std)$r.squared,
         summary(fit_P_LS_std)$r.squared,
         summary(fit_E_SS_std)$r.squared,
         summary(fit_E_SSh_std)$r.squared,
         summary(fit_E_LS_std)$r.squared)
) |>
  dplyr::mutate(f2 = R2/(1 - R2))

tab_std
r2_table

## Plotting
extract_ci <- function(fit, outcome, group){
  tt <- broom::tidy(fit, conf.int = TRUE) |>
    dplyr::filter(term == "x") |>
    dplyr::transmute(
      Group   = group,
      Outcome = outcome,
      Beta    = estimate,
      CI_low  = conf.low,
      CI_high = conf.high,
      p       = p.value
    )
  tt
}

coef_df <- dplyr::bind_rows(
  extract_ci(fit_P_SS_std,  "Short-term spending", "Present"),
  extract_ci(fit_P_SSh_std, "Short-term saving",   "Present"),
  extract_ci(fit_P_LS_std,  "Long-term saving",    "Present"),
  extract_ci(fit_E_SS_std,  "Short-term spending", "Endpoint"),
  extract_ci(fit_E_SSh_std, "Short-term saving",   "Endpoint"),
  extract_ci(fit_E_LS_std,  "Long-term saving",    "Endpoint")
) |>
  dplyr::left_join(tab_std |>
                     dplyr::select(Group, Outcome, p_bh),
                   by = c("Group","Outcome")) |>
  dplyr::mutate(
    Outcome = factor(Outcome,
                     levels = c("Short-term spending","Short-term saving","Long-term saving"))
  )

pd <- position_dodge(width = 0.45)

p_reg<-pd <- position_dodge(width = 0.45)

p_reg <- ggplot(coef_df,
                aes(x = Outcome, y = Beta, color = Group, shape = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.18, linewidth = 0.55, position = pd) +
  geom_point(size = 1.2, stroke = 0.7, position = pd) +
  scale_color_manual(values = c("Endpoint" = "#B22222", "Present" = "#4169E1")) +
  scale_shape_manual(values = c("Endpoint" = 1, "Present" = 2)) +
  scale_y_continuous(
    limits = c(-2, 2),                
    breaks = seq(-2, 2, by = 1),       
    labels = scales::number_format(accuracy = 1),  
    expand = c(0, 0)
  ) +
  labs(
    y = "Standardized ¦Â",
    x = NULL,
    title = NULL,
    color = NULL, shape = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.35),
    axis.title       = element_text(size = 7, colour = "black"),
    axis.text        = element_text(size = 7, colour = "black"),
    panel.background = element_rect(fill = "transparent"),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.position  = "none",
    plot.margin      = ggplot2::margin(3, 3, 2, 2)
  ) +
  patchwork::plot_annotation(
    title = "(D) Allocation changes predicted by ¦¤log k",
    theme = theme(plot.title = element_text(size = 7, colour = "black",
                                            face = "bold",
                                            margin = ggplot2::margin(b = 2)))
  )
ggsave("pic_reg.pdf", plot = p_reg, width = 2, height = 2, units = "in", device = cairo_pdf)