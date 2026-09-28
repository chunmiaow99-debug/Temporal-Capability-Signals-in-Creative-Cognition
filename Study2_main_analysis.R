# ============================================================
# JCB Study 2 — Main Analysis
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MplusAutomation)
})

# ---------- Paths ----------
data_dir <- "/Users/wcm/Library/CloudStorage/OneDrive-个人/JCB-2/R1/JCB_R1/数据"
out_dir  <- file.path(data_dir, "Study2_RSA_output")
mplus_cmd <- "mplus"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load data
# ============================================================

d <- read_csv(
  file.path(data_dir, "data12_full.csv"),
  show_col_types = FALSE
)

# ============================================================
# 2. Descriptive statistics and correlations
# ============================================================

core_vars <- c("cpem", "csem", "gainm", "lossm")
ctrl_vars <- c("gender", "exp", "edu", "opp")
all_vars  <- c(core_vars, ctrl_vars)

lbl <- c(
  "CPE", "CSE", "Gain", "Loss",
  "Gender", "Exp", "Edu", "Opp"
)

msd_core <- d %>%
  summarise(
    across(
      all_of(core_vars),
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
    Variable = lbl[match(var, core_vars)]
  )

msd_ctrl <- d %>%
  distinct(case, .keep_all = TRUE) %>%
  summarise(
    across(
      all_of(ctrl_vars),
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
    Variable = lbl[match(var, ctrl_vars) + 4]
  )

tab_msd <- bind_rows(msd_core, msd_ctrl) %>%
  select(Variable, M, SD) %>%
  mutate(
    across(c(M, SD), ~round(.x, 2))
  )

print(tab_msd)

pm <- d %>%
  group_by(case) %>%
  summarise(
    across(all_of(core_vars), mean),
    .groups = "drop"
  )

d_w <- d %>%
  left_join(
    pm,
    by = "case",
    suffix = c("", ".pm")
  ) %>%
  mutate(
    across(
      all_of(core_vars),
      ~ .x - get(paste0(cur_column(), ".pm")),
      .names = "w_{.col}"
    )
  )

Rw <- d_w %>%
  select(paste0("w_", core_vars)) %>%
  cor(use = "pairwise.complete.obs")

colnames(Rw) <- rownames(Rw) <- core_vars

plevel <- d %>%
  group_by(case) %>%
  summarise(
    across(all_of(core_vars), mean),
    .groups = "drop"
  ) %>%
  left_join(
    d %>%
      distinct(case, .keep_all = TRUE) %>%
      select(case, all_of(ctrl_vars)),
    by = "case"
  )

Rb <- plevel %>%
  select(all_of(all_vars)) %>%
  cor(use = "pairwise.complete.obs")

M <- matrix(
  NA,
  length(all_vars),
  length(all_vars),
  dimnames = list(lbl, lbl)
)

diag(M) <- 1

for (i in 2:length(all_vars)) {
  for (j in 1:(i - 1)) {

    ci <- all_vars[i]
    cj <- all_vars[j]

    M[i, j] <- if (i <= length(core_vars) &&
                    j <= length(core_vars)) {
      Rw[ci, cj]
    } else {
      Rb[ci, cj]
    }
  }
}

print(round(M, 3))

# ============================================================
# 3. Multilevel Bayesian RSA
# ============================================================

run_rsa <- function(core_var, med_var = "", tag) {

  use_vars <- paste(
    compact(
      c(
        "case", core_var,
        "xc", "x2c", "lagyc", "lagy2c", "xyc",
        med_var,
        "xm", "gender", "exp", "edu", "opp"
      )
    ),
    collapse = " "
  )

  d_model <- d %>%
    group_by(case) %>%
    mutate(
      xm = mean(csem, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(
      any_of(
        c(
          "case", core_var,
          "xc", "x2c", "lagyc", "lagy2c", "xyc",
          med_var,
          "xm", "gender", "exp", "edu", "opp"
        )
      )
    )

  med_within <- if (med_var == "") "" else paste0(
    "p71 | cpem ON xc;\n",
    "p72 | cpem ON lagyc;\n",
    "p73 | cpem ON x2c;\n",
    "p74 | cpem ON xyc;\n",
    "p75 | cpem ON lagy2c;\n",
    "p14 | cpem ON ", med_var, ";\n",
    "cpem WITH ", med_var, " @0;\n"
  )

  med_between <- if (med_var == "") "" else paste0(
    med_var, " ON xm gender exp edu opp;\n",
    "[p71] (p71);\n",
    "[p72] (p72);\n",
    "[p73] (p73);\n",
    "[p74] (p74);\n",
    "[p75] (p75);\n",
    "[p14] (p14);\n",
    "cpem WITH ", med_var, " @0;\n"
  )

  med_new <- if (med_var == "") "" else
    "inSA inSB inSC inSD\n"

  med_cons <- if (med_var == "") "" else paste0(
    "inSA = effA * p14;\n",
    "inSB = effB * p14;\n",
    "inSC = effC * p14;\n",
    "inSD = effD * p14;\n"
  )

  model <- mplusObject(

    TITLE = paste0("Study 2 RSA: ", tag),

    VARIABLE = paste0(
      "CLUSTER = case;\n",
      "USEVARIABLES = ", use_vars, ";\n",
      "WITHIN = xc x2c lagyc lagy2c xyc;\n",
      "BETWEEN = xm gender exp edu opp;\n"
    ),

    ANALYSIS = paste0(
      "TYPE = TWOLEVEL RANDOM;\n",
      "ESTIMATOR = BAYES;\n",
      "PROCESSORS = 2;\n",
      "BITERATIONS = (20000);\n",
      "BSEED = 20250704;\n"
    ),

    MODEL = paste0(
      "%WITHIN%\n",
      "b1 | ", core_var, " ON xc;\n",
      "b2 | ", core_var, " ON lagyc;\n",
      "b3 | ", core_var, " ON x2c;\n",
      "b4 | ", core_var, " ON xyc;\n",
      "b5 | ", core_var, " ON lagy2c;\n",
      med_within,
      "%BETWEEN%\n",
      if (med_var != "")
        paste0(med_var, " ON xm gender exp edu opp;\n"),
      "cpem ON xm gender exp edu opp;\n",
      "[b1] (gamma10);\n",
      "[b2] (gamma20);\n",
      "[b3] (gamma30);\n",
      "[b4] (gamma40);\n",
      "[b5] (gamma50);\n",
      med_between
    ),

    MODELCONSTRAINT = paste0(
      "NEW(effA effB effC effD ",
      med_new,
      ");\n",
      "effA = gamma10 + gamma20;\n",
      "effB = gamma30 + gamma40 + gamma50;\n",
      "effC = gamma10 - gamma20;\n",
      "effD = gamma30 - gamma40 + gamma50;\n",
      med_cons
    ),

    OUTPUT = "CINTERVAL;",

    rdata = d_model
  )

  inp <- file.path(
    out_dir,
    paste0("rsa_", tag, ".inp")
  )

  mplusModeler(
    model,
    modelout = inp,
    run = 0L,
    writeData = "always",
    hashfilename = FALSE,
    quiet = TRUE
  )

  runModels(
    target = out_dir,
    filefilter = paste0("rsa_", tag, ".inp"),
    Mplus_command = mplus_cmd,
    showOutput = FALSE,
    replaceOutfile = "always"
  )

  out <- sub("\\.inp$", ".out", inp)
  res <- readModels(out)

  pars <- res$parameters$unstandardized

  key <- pars %>%
    filter(
      paramHeader == "New.Additional.Parameters"
    ) %>%
    select(
      param,
      est,
      posterior_sd,
      pval,
      lower_2.5ci,
      upper_2.5ci
    ) %>%
    mutate(
      across(
        where(is.numeric),
        ~round(.x, 3)
      )
    )

  print(key)

  invisible(key)
}

# ============================================================
# 4. Main models
# ============================================================

direct_cpe <- run_rsa(
  core_var = "cpem",
  med_var  = "",
  tag      = "direct_cpe"
)

gain_model <- run_rsa(
  core_var = "gainm",
  med_var  = "gainm",
  tag      = "gain_mediation"
)

loss_model <- run_rsa(
  core_var = "lossm",
  med_var  = "lossm",
  tag      = "loss_mediation"
)
