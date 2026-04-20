# Episodic thinking: endpoint-focused vs. present-focused
# Behavioral data analysis & viso
# This script requires one file: "behavioral_rawdata.xlsx"
# Programmed by Feng XIAO (updated on 2026.4.19)
############################################################################################################

### Preparation
## Load required packages for analysis
package_list <- c('car','tidyr','dplyr','readxl','effsize','afex','emmeans',
                  'e1071','lmtest','broom','boot','lme4','lmerTest','effectsize',
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

## Add order (EP vs PE) and create Post1/Post2 mapping
## Condition column is in rd_pretest sheet: "EP" or "PE"
## EP = Endpoint first, then Present
## PE = Present first, then Endpoint
rd_pretest_filtered <- rd_pretest_filtered %>%
  mutate(
    order = factor(Condition, levels = c("PE","EP"))  # PE: Present-first; EP: Endpoint-first
  )
# Carry order into posttest (same row order / same participants)
rd_posttest_filtered <- rd_posttest_filtered %>%
  mutate(order = rd_pretest_filtered$order)
# Helper: map (Present, Endpoint) -> (Post1, Post2) by order
map_post12 <- function(present, endpoint, order){
  post1 <- ifelse(order == "EP", endpoint, present)
  post2 <- ifelse(order == "EP", present, endpoint)
  list(post1 = post1, post2 = post2)
}

## Mediation analysis
# M = difference in involvement level (Endpoint - Present)
# Y1 = difference in interval-timing speed (Endpoint - Present)
# Y2 = difference in logk (Endpoint - Present)
id_vec <- 1:nrow(rd_posttest_filtered) #id creation
# M
M_vec <- rd_posttest_filtered$Involvement_E - rd_posttest_filtered$Involvement_P
# Y1
Y1_vec <- (rd_posttest_filtered$TimeE_diff1 + rd_posttest_filtered$TimeE_diff2)/2 -
  (rd_posttest_filtered$TimeP_diff1 + rd_posttest_filtered$TimeP_diff2)/2
# Y2
df_k <- data.frame(
  id       = id_vec,
  Pretest  = rd_k$Pretest,
  Present  = rd_k$Present,
  Endpoint = rd_k$Endpoint
)
df_logk <- df_k %>%
  mutate(
    logk_Pretest  = log(Pretest),
    logk_Present  = log(Present),
    logk_Endpoint = log(Endpoint)
  )
Y2_vec <- df_logk$logk_Endpoint - df_logk$logk_Present
df_med_Y1 <- data.frame(id = id_vec, M = M_vec, Y1 = Y1_vec) %>% na.omit()
df_med_Y2 <- data.frame(id = id_vec, M = M_vec, Y2 = Y2_vec) %>% na.omit()
# Function
run_ws_mediation <- function(df, y_name, R = 5000, seed = 666) {
  set.seed(seed)
  y <- df[[y_name]]
  cat("\n====================================================\n")
  cat(sprintf("Mediation outcome: %s\n", y_name))
  cat(sprintf("N (complete cases) = %d\n", nrow(df)))
  cat("----------------------------------------------------\n")
  # Path a: is mean(M) different from 0?
  cat("Path a (mean(M) vs 0):\n")
  print(t.test(df$M))
  # Path b: does M predict Y?
  cat("\nPath b (Y ~ M):\n")
  lm_b <- lm(y ~ M, data = df)
  print(summary(lm_b))
  # Total effect: is mean(Y) different from 0?
  cat("\nTotal effect (mean(Y) vs 0):\n")
  print(t.test(y))
  # Bootstrap indirect effect a*b
  med_fun <- function(data, i) {
    d <- data[i, ]
    a <- mean(d$M)
    b <- coef(lm(d[[y_name]] ~ d$M))[2]
    a * b
  }
  boot_res <- boot(df, med_fun, R = R)
  cat("\nBootstrap indirect effect (a*b), percentile CI:\n")
  print(boot.ci(boot_res, type = "perc"))
  invisible(list(lm_b = lm_b, boot_res = boot_res))
}
# Analysis
res_Y1 <- run_ws_mediation(df_med_Y1, "Y1", R = 5000)
res_Y2 <- run_ws_mediation(df_med_Y2, "Y2", R = 5000)


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
# Distribution checks (bounded/ceiling + shape)
invP <- rd_posttest_filtered$Involvement_P
invE <- rd_posttest_filtered$Involvement_E
invD <- invE - invP
# Skewness (less than 3 is acceptable)
cat(sprintf("Skewness P=%.2f, E=%.2f, Diff=%.2f\n",
            skewness(invP, na.rm=TRUE),
            skewness(invE, na.rm=TRUE),
            skewness(invD, na.rm=TRUE)))
 #Skewness P=-0.45, E=-0.55, Diff=0.56
# Ceiling / floor proportions (helps justify bounded-scale concerns)
cat(sprintf("Ceiling(>=95) P=%.1f%%, E=%.1f%%\n",
            100*mean(invP >= 95, na.rm=TRUE),
            100*mean(invE >= 95, na.rm=TRUE)))
 #Ceiling(>=95) P=13.2%, E=10.5%
cat(sprintf("Floor(<=5)  P=%.1f%%, E=%.1f%%\n",
            100*mean(invP <= 5,  na.rm=TRUE),
            100*mean(invE <= 5,  na.rm=TRUE)))
 #Floor(<=5)  P=0.0%, E=0.0%
wilcox.test(invE,
            invP,
            paired = TRUE,
            exact = FALSE,
            alternative = "two.sided") #NS

## Duration estimation
## Analysis
# Average 2 trials within each condition per subj
time_pre <- rowMeans(cbind(rd_pretest_filtered$Time_diff1, rd_pretest_filtered$Time_diff2), na.rm=TRUE)
time_P <- rowMeans(cbind(rd_posttest_filtered$TimeP_diff1, rd_posttest_filtered$TimeP_diff2), na.rm=TRUE)
time_E <- rowMeans(cbind(rd_posttest_filtered$TimeE_diff1, rd_posttest_filtered$TimeE_diff2), na.rm=TRUE)
df_time_subj <- data.frame(
  id    = 1:nrow(rd_pretest_filtered),
  order = rd_pretest_filtered$order,
  Pretest  = time_pre,
  Present  = time_P,
  Endpoint = time_E
)
df_time_long <- df_time_subj %>%
  tidyr::pivot_longer(cols = c("Pretest","Present","Endpoint"),
                      names_to = "condition",
                      values_to = "timing_speed") %>%
  mutate(condition = factor(condition, levels = c("Pretest","Present","Endpoint")))
df_time_long <- df_time_long %>%
  left_join(
    data.frame(
      id = 1:nrow(rd_posttest_filtered),
      invP = rd_posttest_filtered$Involvement_P,
      invE = rd_posttest_filtered$Involvement_E
    ),
    by = "id"
  ) %>%
  mutate(
    involvement = dplyr::case_when(
      condition == "Present"  ~ invP,
      condition == "Endpoint" ~ invE,
      TRUE ~ NA_real_
    )
  )
# Skewness checks
skewness(df_time_subj$Pretest,  na.rm = TRUE) #-0.22
skewness(df_time_subj$Present,  na.rm = TRUE) #-0.14
skewness(df_time_subj$Endpoint, na.rm = TRUE) #-0.55
# Main model: condition + order
m_time0 <- lmer(timing_speed ~ condition + order + (1|id), data = df_time_long, REML = FALSE)
anova(m_time0)
anova_res <- anova(m_time0)
eta_squared(anova_res, partial = TRUE) #condition: F(2,76)=6.58, p=.002,etas=0.15
emm_time <- emmeans(m_time0, ~ condition)
pairs(emm_time, adjust = "bonferroni") 
eff_size(emm_time, sigma = sigma(m_time0), edf = df.residual(m_time0))
 #endpoint(pre-post): t(78)=3.44,p=.003,d=0.80
 #present(pre-post): t(78)=2.57,p=.036,d=0.60
# Add involvement as covariate (post conditions only)
df_time_post <- df_time_long %>% filter(condition %in% c("Present","Endpoint"))
m_time_inv <- lmer(timing_speed ~ condition + order + scale(involvement) + (1|id),
                   data = df_time_post, REML = FALSE)
anova(m_time_inv) #NS

## Plotting
df_time_bar_sum <- dplyr::bind_rows(
  tibble::tibble(
    timepoint = "Pretest",
    condition = "Pretest",
    mean = mean(time_pre, na.rm = TRUE),
    se   = sd(time_pre, na.rm = TRUE) / sqrt(sum(!is.na(time_pre))),
    n    = sum(!is.na(time_pre))
  ),
  tibble::tibble(
    timepoint = "Post1",
    condition = "Endpoint",
    mean = mean(time_E[rd_pretest_filtered$order == "EP"], na.rm = TRUE),
    se   = sd(time_E[rd_pretest_filtered$order == "EP"], na.rm = TRUE) /
      sqrt(sum(!is.na(time_E[rd_pretest_filtered$order == "EP"]))),
    n    = sum(!is.na(time_E[rd_pretest_filtered$order == "EP"]))
  ),
  tibble::tibble(
    timepoint = "Post1",
    condition = "Present",
    mean = mean(time_P[rd_pretest_filtered$order == "PE"], na.rm = TRUE),
    se   = sd(time_P[rd_pretest_filtered$order == "PE"], na.rm = TRUE) /
      sqrt(sum(!is.na(time_P[rd_pretest_filtered$order == "PE"]))),
    n    = sum(!is.na(time_P[rd_pretest_filtered$order == "PE"]))
  ),
  tibble::tibble(
    timepoint = "Post2",
    condition = "Endpoint",
    mean = mean(time_E[rd_pretest_filtered$order == "PE"], na.rm = TRUE),
    se   = sd(time_E[rd_pretest_filtered$order == "PE"], na.rm = TRUE) /
      sqrt(sum(!is.na(time_E[rd_pretest_filtered$order == "PE"]))),
    n    = sum(!is.na(time_E[rd_pretest_filtered$order == "PE"]))
  ),
  tibble::tibble(
    timepoint = "Post2",
    condition = "Present",
    mean = mean(time_P[rd_pretest_filtered$order == "EP"], na.rm = TRUE),
    se   = sd(time_P[rd_pretest_filtered$order == "EP"], na.rm = TRUE) /
      sqrt(sum(!is.na(time_P[rd_pretest_filtered$order == "EP"]))),
    n    = sum(!is.na(time_P[rd_pretest_filtered$order == "EP"]))
  )
) %>%
  dplyr::mutate(
    timepoint = factor(timepoint, levels = c("Pretest", "Post1", "Post2")),
    condition = factor(condition, levels = c("Pretest", "Endpoint", "Present"))
  )

p_time <- ggplot(
  df_time_bar_sum,
  aes(x = timepoint, y = mean, fill = condition)
) +
  geom_col(
    width = 0.62,
    colour = "black",
    linewidth = 0.35,
    position = position_dodge(width = 0.70)
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.10,
    linewidth = 0.45,
    position = position_dodge(width = 0.70)
  ) +
  scale_x_discrete(labels = c(
    "Pretest" = "Pretest",
    "Post1"   = "Posttest 1",
    "Post2"   = "Posttest 2"
  )) +
  scale_y_continuous(
    limits = c(-6, 0),
    breaks = seq(-6, 0, by = 1.5),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Duration estimation difference (s)",
    title = NULL
  ) +
  scale_fill_manual(values = c(
    "Pretest"  = "grey70",
    "Endpoint" = "#B22222",
    "Present"  = "#4169E1"
  )) +
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
    theme = theme(
      plot.title = element_text(
        size = 7,
        colour = "black",
        face = "bold",
        margin = ggplot2::margin(b = 2)
      )
    )
  )
ggsave("pic_time.pdf", plot = p_time, width = 2, height = 2, units = "in", device = cairo_pdf)

## Delay discounting
## Analysis
df_k <- data.frame(
  id = 1:nrow(rd_k),
  order = rd_pretest_filtered$order,
  Pretest  = rd_k$Pretest,
  Present  = rd_k$Present,
  Endpoint = rd_k$Endpoint
)
# log-transform
df_k$Pretest  <- log(df_k$Pretest)
df_k$Present  <- log(df_k$Present)
df_k$Endpoint <- log(df_k$Endpoint)
df_logk_long <- df_k %>%
  pivot_longer(cols = c("Pretest","Present","Endpoint"),
               names_to = "condition",
               values_to = "discounting_rate") %>%
  mutate(condition = factor(condition, levels=c("Pretest","Present","Endpoint")))
# Skewness checks
skewness(df_k$Pretest,  na.rm = TRUE) #-0.62
skewness(df_k$Present,  na.rm = TRUE) #0.93
skewness(df_k$Endpoint, na.rm = TRUE) #0.90
# Main model: condition + order
m_logk0 <- lmer(discounting_rate ~ condition + order + (1|id), data=df_logk_long, REML=FALSE)
anova(m_logk0)
anova_res <- anova(m_logk0)
eta_squared(anova_res, partial = TRUE) #condition: F(2,76)=7.02, p=.002, etas=0.16
emm_logk <- emmeans(m_logk0, ~ condition)
pairs(emm_logk, adjust = "bonferroni")
eff_size(emm_logk, sigma = sigma(m_logk0), edf = df.residual(m_logk0))
 #endpoint(pre-post): t(78)=-3.06,p=.009,d=0.71
 #present(pre-post): t(78)=-3.33,p=.004,d=0.77
# Add involvement as covariate (post conditions only)
inv_map <- data.frame(
  id   = 1:nrow(rd_posttest_filtered),
  invP = rd_posttest_filtered$Involvement_P,
  invE = rd_posttest_filtered$Involvement_E
)
df_logk_long <- df_logk_long %>%
  left_join(inv_map, by = "id")
df_logk_post <- df_logk_long %>%
  filter(condition %in% c("Present","Endpoint")) %>%
  mutate(
    involvement = case_when(
      condition == "Present"  ~ invP,
      condition == "Endpoint" ~ invE,
      TRUE ~ NA_real_
    )
  )
df_logk_post$condition <- droplevels(df_logk_post$condition)
m_logk_inv <- lmer(discounting_rate ~ condition + order +
                     scale(involvement) + (1|id),
                   data = df_logk_post,
                   REML = FALSE)

anova(m_logk_inv) #NS

## Plotting
tmp_post_k <- map_post12(df_k$Present, df_k$Endpoint, df_k$order)

df_logk_bar_sum <- dplyr::bind_rows(
  tibble::tibble(
    timepoint = "Pretest",
    condition = "Pretest",
    mean = mean(df_k$Pretest, na.rm = TRUE),
    se   = sd(df_k$Pretest, na.rm = TRUE) / sqrt(sum(!is.na(df_k$Pretest))),
    n    = sum(!is.na(df_k$Pretest))
  ),
  tibble::tibble(
    timepoint = "Post1",
    condition = "Endpoint",
    mean = mean(tmp_post_k$post1[df_k$order == "EP"], na.rm = TRUE),
    se   = sd(tmp_post_k$post1[df_k$order == "EP"], na.rm = TRUE) /
      sqrt(sum(!is.na(tmp_post_k$post1[df_k$order == "EP"]))),
    n    = sum(!is.na(tmp_post_k$post1[df_k$order == "EP"]))
  ),
  tibble::tibble(
    timepoint = "Post1",
    condition = "Present",
    mean = mean(tmp_post_k$post1[df_k$order == "PE"], na.rm = TRUE),
    se   = sd(tmp_post_k$post1[df_k$order == "PE"], na.rm = TRUE) /
      sqrt(sum(!is.na(tmp_post_k$post1[df_k$order == "PE"]))),
    n    = sum(!is.na(tmp_post_k$post1[df_k$order == "PE"]))
  ),
  tibble::tibble(
    timepoint = "Post2",
    condition = "Endpoint",
    mean = mean(tmp_post_k$post2[df_k$order == "PE"], na.rm = TRUE),
    se   = sd(tmp_post_k$post2[df_k$order == "PE"], na.rm = TRUE) /
      sqrt(sum(!is.na(tmp_post_k$post2[df_k$order == "PE"]))),
    n    = sum(!is.na(tmp_post_k$post2[df_k$order == "PE"]))
  ),
  tibble::tibble(
    timepoint = "Post2",
    condition = "Present",
    mean = mean(tmp_post_k$post2[df_k$order == "EP"], na.rm = TRUE),
    se   = sd(tmp_post_k$post2[df_k$order == "EP"], na.rm = TRUE) /
      sqrt(sum(!is.na(tmp_post_k$post2[df_k$order == "EP"]))),
    n    = sum(!is.na(tmp_post_k$post2[df_k$order == "EP"]))
  )
) %>%
  dplyr::mutate(
    timepoint = factor(timepoint, levels = c("Pretest", "Post1", "Post2")),
    condition = factor(condition, levels = c("Pretest", "Endpoint", "Present"))
  )

p_logk <- ggplot2::ggplot(
  df_logk_bar_sum,
  ggplot2::aes(x = timepoint, y = mean, fill = condition)
) +
  ggplot2::geom_col(
    width = 0.62,
    colour = "black",
    linewidth = 0.35,
    position = ggplot2::position_dodge(width = 0.70)
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean - se, ymax = mean + se),
    width = 0.10,
    linewidth = 0.45,
    position = ggplot2::position_dodge(width = 0.70)
  ) +
  ggplot2::scale_x_discrete(labels = c(
    "Pretest" = "Pretest",
    "Post1"   = "Posttest 1",
    "Post2"   = "Posttest 2"
  )) +
  ggplot2::scale_y_continuous(
    limits = c(-12, 0),
    breaks = seq(-12, 0, by = 3),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Log k-value",
    title = NULL
  ) +
  ggplot2::scale_fill_manual(values = c(
    "Pretest"  = "grey70",
    "Endpoint" = "#B22222",
    "Present"  = "#4169E1"
  )) +
  ggplot2::theme_classic(base_size = 8) +
  ggplot2::theme(
    axis.line        = ggplot2::element_line(colour = "black", linewidth = 0.35),
    axis.title       = ggplot2::element_text(size = 7, colour = "black"),
    axis.text        = ggplot2::element_text(size = 7, colour = "black"),
    panel.background = ggplot2::element_rect(fill = "transparent"),
    plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.position  = "none",
    plot.margin      = ggplot2::margin(3, 3, 2, 2)
  ) +
  patchwork::plot_annotation(
    title = "(B) Delay discounting rate",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 7,
        colour = "black",
        face = "bold",
        margin = ggplot2::margin(b = 2)
      )
    )
  )
ggsave("pic_logk.pdf", plot = p_logk, width = 2, height = 2, units = "in", device = cairo_pdf)

## Monetary allocation task
## Analysis
# ID + order (use SubjectNumber as the only id key)
if (!"SubjectNumber" %in% names(rd_pretest_filtered)) {
  stop("rd_pretest_filtered must contain SubjectNumber.")
}
if (!"SubjectNumber" %in% names(rd_posttest_filtered)) {
  rd_posttest_filtered$SubjectNumber <- rd_pretest_filtered$SubjectNumber
}

rd_pretest_filtered$id  <- rd_pretest_filtered$SubjectNumber
rd_posttest_filtered$id <- rd_posttest_filtered$SubjectNumber

if ("Condition" %in% names(rd_pretest_filtered) && !"order" %in% names(rd_pretest_filtered)) {
  rd_pretest_filtered$order <- rd_pretest_filtered$Condition
}
rd_pretest_filtered$order <- factor(rd_pretest_filtered$order, levels = c("EP","PE"))

# Pretest (proportion diffs)
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
# involvement map keyed by SubjectNumber
inv_map <- tibble::tibble(
  id   = rd_posttest_filtered$SubjectNumber,
  invP = rd_posttest_filtered$Involvement_P,
  invE = rd_posttest_filtered$Involvement_E
) %>% dplyr::distinct(id, .keep_all = TRUE)

# Helper: build long dataframe with order + involvement (post only)
mk_long_alloc <- function(id, order, pre, pres, end, metric_name, inv_map) {
  
  df_wide <- tibble::tibble(
    id       = id,
    order    = order,
    Pretest  = pre,
    Present  = pres,
    Endpoint = end
  ) %>%
    tidyr::drop_na()
  
  df_long <- df_wide %>%
    tidyr::pivot_longer(
      cols = c("Pretest","Present","Endpoint"),
      names_to  = "condition",
      values_to = "Score"
    ) %>%
    dplyr::mutate(
      condition = factor(condition, levels = c("Pretest","Present","Endpoint")),
      Metric    = metric_name
    ) %>%
    dplyr::left_join(inv_map, by = "id") %>%
    dplyr::mutate(
      involvement = dplyr::case_when(
        condition == "Present"  ~ invP,
        condition == "Endpoint" ~ invE,
        TRUE ~ NA_real_
      )
    )
  
  list(wide = df_wide, long = df_long)
}

# Build data
alloc_DS  <- mk_long_alloc(
  id    = rd_pretest_filtered$SubjectNumber,
  order = rd_pretest_filtered$order,
  pre   = pre_DS,  pres = pres_DS,  end = end_DS,
  metric_name = "DS",
  inv_map = inv_map
)
alloc_DL  <- mk_long_alloc(
  id    = rd_pretest_filtered$SubjectNumber,
  order = rd_pretest_filtered$order,
  pre   = pre_DL,  pres = pres_DL,  end = end_DL,
  metric_name = "DL",
  inv_map = inv_map
)
alloc_DSL <- mk_long_alloc(
  id    = rd_pretest_filtered$SubjectNumber,
  order = rd_pretest_filtered$order,
  pre   = pre_DSL, pres = pres_DSL, end = end_DSL,
  metric_name = "DSL",
  inv_map = inv_map
)

# Use sum-to-zero contrasts for Type III tests (restore later if needed)
op_contr <- options(contrasts = c("contr.sum","contr.poly"))

# DS = short-term saving - short-term spending
df_DS_wide <- alloc_DS$wide
df_DS_long <- alloc_DS$long
# Skewness checks (per condition)
skewness(df_DS_wide$Pretest,  na.rm = TRUE)  # -2.13
skewness(df_DS_wide$Present,  na.rm = TRUE)  # -0.97
skewness(df_DS_wide$Endpoint, na.rm = TRUE)  # -0.44
# Main model: condition + order
df_DS_long$condition <- factor(df_DS_long$condition, levels = c("Pretest","Present","Endpoint"))
df_DS_long$order     <- factor(df_DS_long$order, levels = c("EP","PE"))
m_DS0 <- lmer(Score ~ condition + order + (1|id), data = df_DS_long, REML = FALSE)

# (A) Condition / order effects
a_DS0 <- anova(m_DS0)
print(a_DS0)
cat(sprintf("DS main effects p: condition=%.4f, order=%.4f\n",
            a_DS0["condition","Pr(>F)"], a_DS0["order","Pr(>F)"]))
 #DS main effects p: condition=0.9700, order=0.1104

# (B) Intercept (grand mean) test: is overall DS different from 0?
summ_fix_DS <- summary(m_DS0)$coefficients
print(summ_fix_DS["(Intercept)", ])
grand_mean_DS <- emmeans(m_DS0, ~ 1)
print(summary(grand_mean_DS, infer = c(TRUE, TRUE)))
mean_est_DS <- summary(grand_mean_DS)$emmean
sd_res_DS   <- sigma(m_DS0)
d_DS_intercept <- mean_est_DS / sd_res_DS
cat(sprintf("DS intercept d (mean/sigma)=%.3f\n", d_DS_intercept))
 #short-term saving>short-term spending: t(38)=3.91,p<.001,d=0.53

# Add involvement as covariate (post only)
df_DS_post <- df_DS_long %>%
  dplyr::filter(condition %in% c("Present","Endpoint")) %>%
  dplyr::mutate(condition = droplevels(condition))
m_DS_inv <- lmer(Score ~ condition + order + scale(involvement) + (1|id),
                 data = df_DS_post, REML = FALSE)
a_DS_inv <- anova(m_DS_inv)
print(a_DS_inv)
cat(sprintf("DS post-only p: condition=%.4f, order=%.4f, involvement=%.4f\n",
            a_DS_inv["condition","Pr(>F)"],
            a_DS_inv["order","Pr(>F)"],
            a_DS_inv["scale(involvement)","Pr(>F)"]))
print(emmeans(m_DS_inv, ~ order))
eta_p2_from_anova <- function(aov_tab, term){
  if (!term %in% rownames(aov_tab)) stop("term not found in anova table.")
  Fv  <- as.numeric(aov_tab[term, "F value"])
  df1 <- as.numeric(aov_tab[term, "NumDF"])
  df2 <- as.numeric(aov_tab[term, "DenDF"])
  eta_p2 <- (Fv * df1) / (Fv * df1 + df2)
  data.frame(term = term, F = Fv, df1 = df1, df2 = df2, eta_p2 = eta_p2)
}
a_DS_inv <- anova(m_DS_inv)
eta_p2_from_anova(a_DS_inv, "order")
 #Endpoint-focused first order preferred more to short-term saving,p=.043,etas=0.11

# Interaction check (post only)
m_DS_int <- lmer(Score ~ condition * order + scale(involvement) + (1|id),
                 data = df_DS_post, REML = FALSE)
a_DS_int <- anova(m_DS_int)
print(a_DS_int)
cat(sprintf("DS post-only interaction p: condition:order=%.4f\n",
            a_DS_int["condition:order","Pr(>F)"]))
 #No interaction effects between condition and order,p=.278

# DL = long-term saving - short-term spending
df_DL_wide <- alloc_DL$wide
df_DL_long <- alloc_DL$long
# Skewness checks (per condition)
skewness(df_DL_wide$Pretest,  na.rm = TRUE)  # -1.49
skewness(df_DL_wide$Present,  na.rm = TRUE)  # -1.16
skewness(df_DL_wide$Endpoint, na.rm = TRUE)  # -0.81

# Main model: condition + order
df_DL_long$condition <- factor(df_DL_long$condition, levels = c("Pretest","Present","Endpoint"))
df_DL_long$order     <- factor(df_DL_long$order, levels = c("EP","PE"))
m_DL0 <- lmer(Score ~ condition + order + (1|id), data = df_DL_long, REML = FALSE)

# (A) Condition / order effects
a_DL0 <- anova(m_DL0)
print(a_DL0)
cat(sprintf("DL main effects p: condition=%.4f, order=%.4f\n",
            a_DL0["condition","Pr(>F)"], a_DL0["order","Pr(>F)"]))
 #DL main effects p: condition=0.8768, order=0.7286

# (B) Intercept (grand mean) test: is overall DL different from 0?
summ_fix_DL <- summary(m_DL0)$coefficients
print(summ_fix_DL["(Intercept)", ])
grand_mean_DL <- emmeans(m_DL0, ~ 1)
print(summary(grand_mean_DL, infer = c(TRUE, TRUE)))
mean_est_DL <- summary(grand_mean_DL)$emmean
sd_res_DL   <- sigma(m_DL0)
d_DL_intercept <- mean_est_DL / sd_res_DL
cat(sprintf("DL intercept d (mean/sigma)=%.3f\n", d_DL_intercept))
 #long-term saving>short-term spending: t(38)=3.70,p=.001,d=1.03

# Add involvement as covariate (post only)
df_DL_post <- df_DL_long %>%
  dplyr::filter(condition %in% c("Present","Endpoint")) %>%
  dplyr::mutate(condition = droplevels(condition))
m_DL_inv <- lmer(Score ~ condition + order + scale(involvement) + (1|id),
                 data = df_DL_post, REML = FALSE)
a_DL_inv <- anova(m_DL_inv)
print(a_DL_inv)
cat(sprintf("DL post-only p: condition=%.4f, order=%.4f, involvement=%.4f\n",
            a_DL_inv["condition","Pr(>F)"],
            a_DL_inv["order","Pr(>F)"],
            a_DL_inv["scale(involvement)","Pr(>F)"]))
 #DL post-only p: condition=0.9222, order=0.8727, involvement=0.7143

# DSL = long-term saving - short-term saving
df_DSL_wide <- alloc_DSL$wide
df_DSL_long <- alloc_DSL$long
# Skewness checks (per condition)
skewness(df_DSL_wide$Pretest,  na.rm = TRUE)  # -0.22
skewness(df_DSL_wide$Present,  na.rm = TRUE)  # -0.34
skewness(df_DSL_wide$Endpoint, na.rm = TRUE)  # -0.25

# Main model: condition + order
df_DSL_long$condition <- factor(df_DSL_long$condition, levels = c("Pretest","Present","Endpoint"))
df_DSL_long$order     <- factor(df_DSL_long$order, levels = c("EP","PE"))
m_DSL0 <- lmer(Score ~ condition + order + (1|id), data = df_DSL_long, REML = FALSE)

# (A) Condition / order effects
a_DSL0 <- anova(m_DSL0)
print(a_DSL0)
cat(sprintf("DSL main effects p: condition=%.4f, order=%.4f\n",
            a_DSL0["condition","Pr(>F)"], a_DSL0["order","Pr(>F)"]))
 #DSL main effects p: condition=0.9814, order=0.1845

# (B) Intercept (grand mean) test: is overall DSL different from 0?
summ_fix_DSL <- summary(m_DSL0)$coefficients
print(summ_fix_DSL["(Intercept)", ]) #NS

# Add involvement as covariate (post only)
df_DSL_post <- df_DSL_long %>%
  dplyr::filter(condition %in% c("Present","Endpoint")) %>%
  dplyr::mutate(condition = droplevels(condition))
m_DSL_inv <- lmer(Score ~ condition + order + scale(involvement) + (1|id),
                  data = df_DSL_post, REML = FALSE)
a_DSL_inv <- anova(m_DSL_inv)
print(a_DSL_inv)
cat(sprintf("DSL post-only p: condition=%.4f, order=%.4f, involvement=%.4f\n",
            a_DSL_inv["condition","Pr(>F)"],
            a_DSL_inv["order","Pr(>F)"],
            a_DSL_inv["scale(involvement)","Pr(>F)"]))
print(emmeans(m_DSL_inv, ~ order)) #NS

# restore contrasts option
options(op_contr)

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