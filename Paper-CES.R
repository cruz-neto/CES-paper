## =============================================================================
## RECREATIONAL DEMAND AND ACCESS COSTS IN NATIONAL FOREST OF SILVANIA
## Individual Travel Cost Method (ITCM) - REVISED pipeline
## Addresses comments from Reviewer 1 and Reviewer 2 (Dossie Geoturismo e o
## Cerrado Brasileiro, avaliacao de 26/07/2026 e 16/08/2026)
##

options(scipen = 999)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(MASS); library(AER)
  library(lmtest); library(sandwich); library(performance); library(car)
})

## -----------------------------------------------------------------------
## 1. LOAD RAW DATA
## -----------------------------------------------------------------------
database <- read_csv("database.csv", show_col_types = FALSE)
names(database)[1:2] <- c("id", "idx")

## -----------------------------------------------------------------------
## 2. DATA CLEANING (documented, auditable)
## -----------------------------------------------------------------------

## 2.1 Remove rows with a corrupted fuel_consumption value (Unix-timestamp
##     leaked into a numeric field during export -> travel_cost ~ US$6.2bn)
corrupted_ids <- database$id[database$fuel_consumption > 1e6]
database <- database %>% filter(!id %in% corrupted_ids)

## 2.2 Fix distance = 0 caused by a failed lookup on accented municipality
##     names (Brasilia, Morrinhos, Goiania, Silvania all affected).
##     Reference distances: in-sample validated values are used whenever an
##     un-corrupted respondent from the same origin exists; otherwise the
##     official road distance (source: road-distance calculators, see
##     README) is used.
database$start_clean <- database$start
database$start_clean[grepl("Silv", database$start_clean)]  <- "Silvania"
database$start_clean[grepl("Goi", database$start_clean) & database$start_clean != "Gameleira de Goias"] <- "Goiania"

database <- database %>%
  mutate(
    distance_fixed = case_when(
      distance > 0                ~ distance,
      start_clean == "Silvania"   ~ 8,      # in-sample validated (n=76 non-zero rows, all = 8 km)
      start_clean == "Goiania"    ~ 73.1,   # in-sample validated (n=10 non-zero rows, all = 73.1 km)
      start_clean == "Brasilia"   ~ 180,    # road distance, external source
      start_clean == "Morrinhos"  ~ 187,    # road distance, external source
      TRUE ~ distance
    )
  )

## 2.3 Income midpoints (US$/month); open top bracket uses the conventional
##     1.5x-lower-bound rule for open-ended income brackets.
database <- database %>%
  mutate(
    income_mid = case_when(
      income == "ate 1"            ~ 121,
      income == "acima de 1 ate 3" ~ 485,
      income == "acima de 3 ate 5" ~ 970,
      income == "acima de 5"       ~ 1818,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(income_mid))   # drops 1 row with invalid/missing income code

## -----------------------------------------------------------------------
## 3. VARIABLE CONSTRUCTION
## -----------------------------------------------------------------------

## 3.1 Generalized access cost (synthetic construction):
##   TC_i = 2*d_i*c_km + k*(w_i/h_m)*(2*d_i/v)
##   d_i   = one-way distance (km)
##   c_km  = US$0.26/km   (fuel + vehicle wear)
##   k     = 0.30         (opportunity-cost-of-time coefficient)
##   w_i   = monthly income midpoint (US$)
##   h_m   = 220           (hours/month)
##   v     = 60 km/h       (assumed average travel speed)
c_km <- 0.26; k <- 0.30; h_m <- 220; v <- 60

database <- database %>%
  mutate(
    travel_cost_new = 2*distance_fixed*c_km + k*(income_mid/h_m)*(2*distance_fixed/v)
  )

## 3.2 Socio-economic dummies, rebuilt from raw ordinal variables using the
##     reference categories stated in Section 2.4 of the manuscript.
database <- database %>%
  mutate(
    age_old_dummy         = ifelse(age %in% c("36 a 59 anos","60 anos ou mais"), 1, 0),        # ref = 18-35
    schooling_below_dummy = ifelse(schooling %in% c("Ensino medio incompleto","Ensino completo","Graduacao incompleta"), 1, 0), # ref = undergrad+
    income_dummy          = ifelse(income %in% c("ate 1","acima de 1 ate 3"), 1, 0),           # ref = >3 MW
    sex_dummy              = ifelse(sex == "Feminino", 1, 0)
  )

cat("Final analytic sample size:", nrow(database), "\n")
write_csv(database, "database_final_with_dummies.csv")

## -----------------------------------------------------------------------
## 4. DESCRIPTIVE STATISTICS
## -----------------------------------------------------------------------
summary(database$visitation); var(database$visitation); mean(database$visitation)
summary(database$distance_fixed)
summary(database$travel_cost_new)

## -----------------------------------------------------------------------
## 5. COUNT-DATA MODELS
## -----------------------------------------------------------------------

## 5.1 Baseline Poisson vs Negative Binomial
pois <- glm(visitation ~ log(travel_cost_new), family = poisson(link = "log"), data = database)
summary(pois)
dispersiontest(pois)

nb <- glm.nb(visitation ~ log(travel_cost_new), data = database)
summary(nb)

lrtest(pois, nb)
AIC(pois, nb); BIC(pois, nb)

## 5.2 Full model with socio-economic controls (log specification, final model)
nb_full <- glm.nb(
  visitation ~ log(travel_cost_new) + income_dummy + age_old_dummy + schooling_below_dummy + sex_dummy,
  data = database
)
summary(nb_full)
cat("Theta:", nb_full$theta, " SE:", nb_full$SE.theta, "\n")

pois_full <- glm(
  visitation ~ log(travel_cost_new) + income_dummy + age_old_dummy + schooling_below_dummy + sex_dummy,
  family = poisson(link = "log"), data = database
)
dispersiontest(pois_full)
lrtest(pois_full, nb_full)
AIC(pois_full, nb_full); BIC(pois_full, nb_full)

## 5.3 Robustness: linear travel-cost specification
nb_full_lin <- glm.nb(
  visitation ~ travel_cost_new + income_dummy + age_old_dummy + schooling_below_dummy + sex_dummy,
  data = database
)
summary(nb_full_lin)
AIC(nb_full, nb_full_lin)   # log specification preferred

## -----------------------------------------------------------------------
## 6. DIAGNOSTICS REQUESTED BY REVIEWER 1
## -----------------------------------------------------------------------

## 6.1 Multicollinearity
vif(nb_full)

## 6.2 Heteroskedasticity-robust (sandwich) standard errors
coeftest(nb_full, vcov = sandwich)

## 6.3 Overdispersion / zero-inflation (performance package)
check_overdispersion(nb_full)
check_zeroinflation(nb_full)

## 6.4 Influential observations
cd <- cooks.distance(nb_full)
cat("Max Cook's distance:", max(cd), " | # obs above 4/n:", sum(cd > 4/nrow(database)), "\n")

## 6.5 Pseudo-R2
r2_mcfadden(nb_full)

## 6.6 Residual diagnostics (plots)
png("fig_residuals_vs_fitted.png", width = 1000, height = 700, res = 130)
plot(fitted(nb_full), residuals(nb_full, type = "deviance"),
     xlab = "Fitted values", ylab = "Deviance residuals",
     main = "Residuals vs Fitted - Negative Binomial (final model)")
abline(h = 0, lty = 2, col = "red")
dev.off()

png("fig_qq_residuals.png", width = 1000, height = 700, res = 130)
qqnorm(residuals(nb_full, type = "deviance"), main = "Normal Q-Q - Deviance residuals")
qqline(residuals(nb_full, type = "deviance"), col = "red")
dev.off()

png("fig_cooks_distance.png", width = 1000, height = 700, res = 130)
plot(cd, type = "h", ylab = "Cook's distance", main = "Influence diagnostics")
abline(h = 4/nrow(database), col = "red", lty = 2)
dev.off()

## 6.7 On-site (endogenous) truncation note
## visitation is bounded below at 1 because only actual visitors were
## interviewed on-site; this is a well-known source of bias in intercept
## surveys (Shaw, 1988; Englin & Shonkwiler, 1995) and is reported as a
## limitation in the manuscript rather than corrected via a truncated-NB
## re-estimation, to keep the model comparable across both reviewers'
## comments; a truncated-NB robustness check is straightforward future work.
table(database$visitation)

## 6.8 Spatial dependence proxy (no point coordinates available; only
## municipality of origin). Residual clustering by municipality is tested
## as a practical proxy; a formal Moran's I would require point coordinates.
database$resid_dev <- residuals(nb_full, type = "deviance")
summary(aov(resid_dev ~ factor(start_clean), data = database))

## -----------------------------------------------------------------------
## 7. PEARSON CORRELATIONS AMONG ACCESS-COST COMPONENTS (Figure 2, revised)
## -----------------------------------------------------------------------
library(GGally)
ggpairs(database[, c("travel_cost_new","distance_fixed","deslocation_time")],
        columnLabels = c("Access Cost (US$)","Distance (km)","Travel Time (min)"))

## -----------------------------------------------------------------------
## 8. FINAL RESULTS TABLE
## -----------------------------------------------------------------------
sink("table_models_final.txt")
cat("=== Table 1: Negative Binomial model (revised) ===\n")
print(summary(nb_full))
cat("\n=== Table 2: Model diagnostics ===\n")
cat("AIC (Poisson):", AIC(pois_full), " | AIC (NB):", AIC(nb_full), "\n")
cat("BIC (Poisson):", BIC(pois_full), " | BIC (NB):", BIC(nb_full), "\n")
cat("Dispersion test (Poisson): see dispersiontest(pois_full)\n")
cat("LR test Poisson vs NB: see lrtest(pois_full, nb_full)\n")
cat("Theta (dispersion parameter):", nb_full$theta, " SE:", nb_full$SE.theta, "\n")
cat("VIF range:", paste(round(range(vif(nb_full)),3), collapse=" - "), "\n")
cat("McFadden pseudo-R2:", r2_mcfadden(nb_full)$R2, "\n")
sink()
