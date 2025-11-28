Cox Regression Analysis Project
# Cox Proportional Hazards Regression: Lung Cancer Survival

This project demonstrates basic **survival analysis** and **Cox proportional hazards regression** in R using lung cancer patient data.

We will:

1. Estimate and visualise **Kaplan–Meier survival curves**.
2. Fit **univariable** and **multivariable Cox models** for symptom status.
3. Briefly explain the **back-door path** concept and its relation to multiple regression.
4. Analyse **weight loss** as a potential predictor of survival and interpret hazard ratios.
5. Export model results to CSV for reporting.

---

## 1. Data and Setup

### 1.1 Required R packages

```r
# Install if needed:
# install.packages("survival")

library(survival)
# Set your working directory (edit to your own path if needed)
# setwd("C:/Users/YourName/Path/To/Project")

# Read data (update path if you store it in e.g. "data/lung.txt")
d_lung <- read.csv("lung.txt")

# Take a quick look at the data
d_lung
str(d_lung)

'data.frame': 228 obs. of 7 variables:
 $ id      : int
 $ time    : int  # follow-up time in days
 $ status  : int  # event indicator (1 = death, 0 = censored)
 $ age     : int
 $ sex     : int  # 1 = male, 2 = female
 $ symptoms: int  # 1 = symptomatic, 0 = asymptomatic
 $ wt.loss : int  # weight loss in pounds

2. Exercise 1 – Kaplan–Meier Survival Curves
2.1 Overall survival function

We first estimate the overall survival function using survfit():

m0 <- survfit(Surv(time, status) ~ 1, data = d_lung)

plot(
  m0,
  xlab = "Days",
  ylab = "Survival probability",
  mark.time = TRUE
)

m0

Call: survfit(formula = Surv(time, status) ~ 1, data = d_lung)

      n events median 0.95LCL 0.95UCL
[1,] 228    165    310     285     363

Included plot
![PLOT_1_m0](images/PLOT_1_m0.png)

Overall Kaplan–Meier curve (saved separately as PLOT_1_m0):

2.2 One-year survival probability

Because follow-up time is recorded in days, we use times = 365.25 for 1-year survival:

summary(m0, times = 365.25)

Call: survfit(formula = Surv(time, status) ~ 1, data = d_lung)

 time n.risk n.event survival std.err lower 95% CI upper 95% CI
  365     65     121    0.409  0.0358        0.345        0.486

Estimated 1-year survival ≈ 0.41 (41%).

Survival by symptom status
Now we compare survival between symptomatic and asymptomatic patients:
m1 <- survfit(Surv(time, status) ~ factor(symptoms), data = d_lung)

plot(
  m1,
  conf.int = FALSE,
  mark.time = TRUE,
  col = c("black", "orange"),
  xlab = "Days",
  ylab = "Survival probability"
)

legend(
  "topright",
  lty = c(1, 1), lwd = 2,
  col = c("orange", "black"),
  legend = c("Symptomatic", "Asymptomatic")
)

Included plot
![PLOT_2_m1](images/PLOT_2_m1.png)

3. Exercise 2 – Cox Regression for Symptom Status
3.1 Univariable Cox model (symptoms only)

We fit a simple Cox model with symptom status as the only explanatory variable:

mcox1 <- coxph(Surv(time, status) ~ factor(symptoms), data = d_lung)
summary(mcox1)

Call:
coxph(formula = Surv(time, status) ~ factor(symptoms), data = d_lung)

n = 228, number of events = 165 

                    coef exp(coef) se(coef)      z Pr(>|z|)   
factor(symptoms)1 0.5372   1.7112  0.1874  2.866  0.00416 **

exp(coef) exp(-coef) lower .95 upper .95
1.711       0.5844     1.185     2.471

Concordance= 0.562
Likelihood ratio test= 8.99  on 1 df, p=0.003
Score (logrank) test = 8.41  on 1 df, p=0.004

Log-rank test p-value

hazard_ratio     <- coef(summary(mcox1))[1, 2]
hazard_ratio_95ci <- summary(mcox1)$conf.int[1, 3:4]
hazard_ratio_p   <- coef(summary(mcox1))[1, 5]

hazard_ratio
hazard_ratio_95ci
hazard_ratio_p

HR ≈ 1.71

95% CI ≈ 1.19–2.47

p ≈ 0.004

Interpretation: Symptomatic patients have about 71% higher hazard of death compared with asymptomatic patients, and this association is statistically significant.

3.2 Multivariable Cox model (symptoms + age + sex)

We now adjust for age and sex:

mcox2 <- coxph(
  Surv(time, status) ~ factor(symptoms) + age + factor(sex),
  data = d_lung
)
summary(mcox2)


Example (abridged) output:

Call:
coxph(formula = Surv(time, status) ~ factor(symptoms) + age + 
         factor(sex), data = d_lung)

n = 228, number of events = 165 

                    coef exp(coef) se(coef)      z Pr(>|z|)    
factor(symptoms)1  0.559982  1.750640 0.188951  2.964 0.003040 ** 
age                0.014337  1.014440 0.009122  1.572 0.116022    
factor(sex)2      -0.559854  0.571292 0.168069 -3.331 0.000865 ***

Extract adjusted hazard ratio for symptom status
hazard_ratio2      <- coef(summary(mcox2))[1, 2]
hazard_ratio2_95ci <- summary(mcox2)$conf.int[1, 3:4]
hazard_ratio2_p    <- coef(summary(mcox2))[1, 5]

hazard_ratio2
hazard_ratio2_95ci
hazard_ratio2_p


Adjusted HR (symptoms) ≈ 1.75

95% CI ≈ 1.21–2.54

p ≈ 0.003

Interpretation: Even after adjusting for age and sex, symptomatic patients still have a substantially higher hazard of death compared with asymptomatic patients.

3.3 Summary table of symptom hazard ratios
results_symptoms <- data.frame(
  Model    = c("Unadjusted", "Adjusted"),
  HR       = c(hazard_ratio, hazard_ratio2),
  CI_lower = c(hazard_ratio_95ci[1], hazard_ratio2_95ci[1]),
  CI_upper = c(hazard_ratio_95ci[2], hazard_ratio2_95ci[2]),
  Pvalue   = c(hazard_ratio_p, hazard_ratio2_p)
)

results_symptoms


Example:

      Model       HR CI_lower CI_upper      Pvalue
1 Unadjusted 1.711208 1.185073 2.470929 0.004159160
2   Adjusted 1.750640 1.208819 2.535319 0.003040366


Save to CSV:

write.csv2(results_symptoms, file = "results_symptoms_survival.csv", row.names = FALSE)

4. Exercise 3 – Back-door Path Concept

A back-door path is any non-causal pathway between an exposure 
𝑋
X and an outcome 
𝑌
Y that:

Starts with an arrow pointing into 
𝑋
X, and

Can create spurious associations or confounding.

Multiple regression can help block these back-door paths by conditioning on appropriate confounders (variables that are common causes of both 
𝑋
X and 
𝑌
Y). By including these confounders in the model (but avoiding colliders and mediators), the regression aims to better isolate the causal effect of 
𝑋
X on 
𝑌
Y.

5. Exercise 4 – Weight Loss and Survival
5.1 Convert weight loss from pounds to kilograms
# 1 pound = 0.453592 kg
d_lung$wtloss_kg <- d_lung$wt.loss * 0.453592

5.2 Univariable Cox model for weight loss (per 1 kg)
mcox3 <- coxph(Surv(time, status) ~ wtloss_kg, data = d_lung)
summary(mcox3)


Example (abridged) output:

Call:
coxph(formula = Surv(time, status) ~ wtloss_kg, data = d_lung)

n = 214, number of events = 152 
(14 observations deleted due to missingness)

            coef exp(coef) se(coef)     z Pr(>|z|)
wtloss_kg 0.002908  1.002913 0.013402 0.217    0.828

exp(coef) exp(-coef) lower .95 upper .95
1.003     0.9971    0.9769      1.03

Concordance= 0.525
Likelihood ratio test= 0.05  on 1 df, p=0.8

Extract HR, CI, and p-value
summary(mcox3)$sctest["pvalue"]
summary(mcox3)$conf.int

hr_wtloss1      <- coef(summary(mcox3))[1, 2]
hr_wtloss1_95ci <- summary(mcox3)$conf.int[1, 3:4]
hr_wtloss1_p    <- coef(summary(mcox3))[1, 5]

hr_wtloss1
hr_wtloss1_95ci
hr_wtloss1_p


HR per 1 kg ≈ 1.003

95% CI ≈ 0.98–1.03

p ≈ 0.83

Interpretation: There is no statistically significant association between weight loss and survival. A hazard ratio of 1.003 means each additional kilogram of weight loss is associated with only a 0.3% increase in the hazard of death, which is negligible.

5.3 Adjusted Cox model (per 5 kg, adjusted for age and sex)

We rescale weight loss to per 5 kg:

d_lung$wtloss_5kg <- d_lung$wtloss_kg / 5

mcox4 <- coxph(Surv(time, status) ~ wtloss_5kg + sex + age, data = d_lung)
summary(mcox4)


Example (abridged) output:

Call:
coxph(formula = Surv(time, status) ~ wtloss_5kg + sex + age, data = d_lung)

n = 214, number of events = 152 
(14 observations deleted due to missingness)

              coef exp(coef) se(coef)      z Pr(>|z|)
wtloss_5kg  0.008373  1.008408 0.068271  0.123   0.9024   
sex        -0.521032  0.593907 0.174354 -2.988   0.0028 **
age         0.020088  1.020291 0.009664  2.079   0.0377 *

Extract HR per 5 kg
hr_wtloss2      <- coef(summary(mcox4))[1, 2]
hr_wtloss2_95ci <- summary(mcox4)$conf.int[1, 3:4]
hr_wtloss2_p    <- coef(summary(mcox4))[1, 5]

hr_wtloss2
hr_wtloss2_95ci
hr_wtloss2_p


HR per 5 kg ≈ 1.008

95% CI ≈ 0.88–1.15

p ≈ 0.90

Interpretation: After adjusting for age and sex, there is still no evidence that weight loss is associated with survival. A hazard ratio of 1.008 indicates that for every 5 kg increase in weight loss, the hazard of death increases by only 0.8%, which is not clinically meaningful.

5.4 Summary table of weight-loss hazard ratios
results_wtloss <- data.frame(
  Model    = c("Unadjusted", "Adjusted"),
  HR       = c(hr_wtloss1, hr_wtloss2),
  CI_lower = c(hr_wtloss1_95ci[1], hr_wtloss2_95ci[1]),
  CI_upper = c(hr_wtloss1_95ci[2], hr_wtloss2_95ci[2]),
  Pvalue   = c(hr_wtloss1_p, hr_wtloss2_p)
)

results_wtloss


Example:

      Model       HR  CI_lower CI_upper    Pvalue
1 Unadjusted 1.002913 0.9769125 1.029605 0.8281974
2   Adjusted 1.008408 0.8821128 1.152786 0.9023884


Save to CSV:

write.csv2(results_wtloss, file = "results_wtloss.csv", row.names = FALSE)