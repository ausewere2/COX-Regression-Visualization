Cox Regression Survival Analysis — Lung Cancer Dataset

This project demonstrates Kaplan–Meier survival analysis, Cox proportional hazards regression, confounder adjustment, and weight-loss survival modeling using the survival package in R.

All analyses use the lung cancer dataset (lung.txt) containing 228 observations and 7 variables.

📁 Project Structure
/images               <- folder containing all your JPG plots
README.md             <- this file
lung.txt              <- dataset
results_symptoms_survival.csv
results_wtloss.csv

1. Load Libraries and Data
library(survival)

setwd()
d_lung = read.csv("C:\\Users\\PAT EDOCHE\\Downloads\\R_FOLDER\\period 2\\lung.txt")

str(d_lung)

2. Kaplan–Meier Survival Function
m0 <- survfit(Surv(time, status) ~ 1, data = d_lung)

plot(m0, xlab = "Days", ylab = "Survival probability", mark.time = TRUE)

📸 Figure 1. Overall Survival Curve
![PLOT_1_m0](PLOT_1_m0.JPG)

1-Year Survival Estimate
summary(m0, times = 365.25)

3. Kaplan–Meier Curves by Symptom Status
m1 <- survfit(Surv(time, status) ~ factor(symptoms), data = d_lung)

plot(m1, conf.int = FALSE, mark.time = TRUE, col = c("black", "orange"),
     xlab = "Days", ylab = "Survival probability")

legend("topright",
       lty = c(1,1), lwd = 2,
       col = c("orange","black"),
       c("Symptomatic","Asymptomatic"))

📸 Figure 2. Survival by Symptom Status
![PLOT_2_m1](PLOT_2_m1.JPG)

4. Cox Regression — Symptoms Only (Unadjusted)
mcox1 <- coxph(Surv(time, status) ~ factor(symptoms), data = d_lung)
summary(mcox1)

Hazard Ratio

HR = 1.71 (95% CI 1.19–2.47)

p = 0.004

Symptomatic patients had 71% higher hazard of death than asymptomatic patients.

CSV Export
write.csv2(results, "results_symptoms_survival.csv")

5. Cox Regression — Adjusting for Age and Sex
mcox2 <- coxph(Surv(time, status) ~ factor(symptoms) + age + factor(sex), data = d_lung)
summary(mcox2)

Adjusted Hazard Ratio

HR = 1.75 (95% CI 1.21–2.54)

p = 0.003

Adjustment slightly strengthened the effect.

6. What Is a Back–Door Path?

A back-door path is any non-causal pathway between exposure (X) and outcome (Y) that begins with an arrow pointing into X.
Multiple regression blocks these back-door paths by conditioning on true confounders, preventing biased estimation.

7. Weight Loss and Survival — Unadjusted & Adjusted Models
Convert to kg
d_lung$wtloss_kg <- d_lung$wt.loss * 0.453592

Unadjusted Model
mcox3 <- coxph(Surv(time,status) ~ wtloss_kg, data=d_lung)
summary(mcox3)


HR = 1.003

95% CI = 0.98–1.03

p = 0.83

➡ No evidence that weight loss predicts survival.

Adjusted (per 5 kg weight loss)
d_lung$wtloss_5kg <- d_lung$wtloss_kg / 5
mcox4 <- coxph(Surv(time,status) ~ wtloss_5kg + sex + age, data=d_lung)
summary(mcox4)


HR = 1.008

95% CI = 0.88–1.15

p = 0.90

➡ Still no meaningful association.

8. Final Results Table
results = data.frame(
  Model     = c("Unadjusted", "Adjusted"),
  HR        = c(hr_wtloss1, hr_wtloss2),
  CI_lower  = c(hr_wtloss1_95ci[1], hr_wtloss2_95ci[1]),
  CI_upper  = c(hr_wtloss1_95ci[2], hr_wtloss2_95ci[2]),
  Pvalue    = c(hr_wtloss1_p, hr_wtloss2_p)
)


Saved as:

results_wtloss.csv

9. Required Plot Modifications (Exercise Requirement)

Your final figure includes:

✔ New x-axis tick labels
✔ New line type (lty=2)
✔ Title + annotation added

Example Code
plot(m0, xlab="Days", ylab="Survival Probability", lty=2)

axis(1, at=seq(0,1000,200), labels=paste0(seq(0,1000,200)," d"))

title("Modified Kaplan–Meier Curve")
text(600,0.2,"n = 228 patients", cex=0.8)

📸 Figure 3. Modified Final Plot

10. Interpretation Summary

Symptomatic patients have significantly lower survival.

Age is a weak predictor; sex shows a strong protective effect (females).

Weight loss does not predict survival.

Adjustment for confounders increases accuracy but does not change the weight-loss result.
