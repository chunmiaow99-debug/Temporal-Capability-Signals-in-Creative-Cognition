# ============================================================
# JCB Study 1 — Main Analysis
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---------- Data path ----------
data_dir <- "/Users/wcm/Library/CloudStorage/OneDrive-个人/JCB-2/R1/JCB_R1/数据"

# ============================================================
# 1. Load data and construct analysis variables
# ============================================================

master <- read_csv(
  file.path(data_dir, "master_s1_long.csv"),
  show_col_types = FALSE
)

person_master <- master %>%
  group_by(case) %>%
  summarise(
    cse_m  = mean(csem),
    cse_sd = sd(csem),
    .groups = "drop"
  ) %>%
  left_join(
    master %>%
      filter(wave == 1) %>%
      select(
        case, cpem, cpe1m, cpe2m, cpe3m,
        bcpem, gender, age, tenure
      ),
    by = "case"
  )

# ============================================================
# 2. Descriptive statistics and correlations
# ============================================================

tab1_vars <- c(
  "cse_m", "cse_sd", "cpem", "cpe1m", "cpe2m", "cpe3m",
  "bcpem", "gender", "age", "tenure"
)

tab1_lbl <- c(
  "Average CSE level",
  "CSE instability",
  "Creative process engagement",
  "Identifying problems",
  "Information searching",
  "Idea generating",
  "CPE tendency",
  "Gender",
  "Age",
  "Work experience"
)

tab1_msd <- person_master %>%
  summarise(
    across(
      all_of(tab1_vars),
      list(M = ~mean(.x), SD = ~sd(.x))
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = c("var", "stat"),
    names_sep = "_(?=[MSD]+$)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(
    Variable = tab1_lbl[match(var, tab1_vars)]
  ) %>%
  select(Variable, M, SD) %>%
  mutate(
    across(c(M, SD), ~round(.x, 2))
  )

print(tab1_msd, n = 10)

R <- person_master %>%
  select(all_of(tab1_vars)) %>%
  cor() %>%
  round(3)

print(R)

# ============================================================
# 3. Main regression analyses
# ============================================================

d <- person_master

run_models <- function(dv, d) {

  f1 <- as.formula(
    paste(dv, "~ bcpem + gender + age + tenure")
  )

  f2 <- as.formula(
    paste(dv, "~ bcpem + gender + age + tenure + cse_m")
  )

  f3 <- as.formula(
    paste(dv, "~ bcpem + gender + age + tenure + cse_m + cse_sd")
  )

  list(
    lm(f1, data = d),
    lm(f2, data = d),
    lm(f3, data = d)
  )
}

dv_set <- c(
  cpem  = "cpem",
  ident = "cpe1m",
  info  = "cpe2m",
  idea  = "cpe3m"
)

all_models <- list()

for (nm in names(dv_set)) {
  all_models[[nm]] <- run_models(
    dv_set[[nm]],
    d
  )
}

# ============================================================
# 4. Main results
# ============================================================

# Full regression coefficients and model fit
for (nm in names(dv_set)) {

  cat("\n---", nm, "---\n")

  for (i in 1:3) {
    model <- all_models[[nm]][[i]]
    print(summary(model)$coefficients)
    cat("R2 =", round(summary(model)$r.squared, 3), "\n")
  }
}

# CSE instability coefficients from the final model for each outcome
cat("\n--- CSE instability coefficients ---\n")

for (nm in names(dv_set)) {

  s <- summary(
    all_models[[nm]][[3]]
  )$coefficients

  b  <- s["cse_sd", 1]
  se <- s["cse_sd", 2]
  p  <- s["cse_sd", 4]

  p_txt <- ifelse(
    p < .001,
    "< .001",
    paste0("= ", sprintf("%.3f", p))
  )

  cat(sprintf(
    "%-6s b=%.3f, SE=%.3f, p %s, R2=%.3f\n",
    nm,
    b,
    se,
    p_txt,
    summary(all_models[[nm]][[3]])$r.squared
  ))
}
