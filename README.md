# metabolic-syndrome-meta-analysis
r code for data analaysis of metabolic sx in nepal
############################################################
# Meta-analysis: Metabolic Syndrome in Nepal (k ~ 8, I2 high)
# PRIMARY: GLMM (binomial-normal) with ML tau^2
# SENSITIVITY: Logit (Inverse-variance) with REML tau^2 + HK
# No funnel/Egger. Manual LOO and cumulative series.
############################################################

## 0) Packages --------------------------------------------------------------
need <- c("readxl","writexl","dplyr","tidyr","stringr","janitor",
          "meta","metafor","ggplot2","binom","readr")
new  <- need[!(need %in% installed.packages()[,"Package"])]
if(length(new)) install.packages(new, dependencies = TRUE)

suppressPackageStartupMessages({
  library(readxl); library(writexl)
  library(dplyr);  library(tidyr);  library(stringr); library(janitor)
  library(meta);   library(metafor)
  library(ggplot2);library(binom);  library(readr)
})

## 1) Paths & Input ---------------------------------------------------------
base_dir   <- "C:\\Users\\VYAS IT\\Desktop\\meta analysis\\metabolic syndrome"
in_xlsx    <- file.path(base_dir, "for_mets.xlsx")
in_xls     <- file.path(base_dir, "for_mets.xls")
in_csv     <- file.path(base_dir, "for_mets.csv")
infile     <- if (file.exists(in_xlsx)) in_xlsx else if (file.exists(in_xls)) in_xls else in_csv

if (!file.exists(infile)) stop("Input file not found. Check the path and filename.")

out_dir  <- file.path(base_dir, "outputs")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

dat_raw <- if (grepl("\\.csv$", infile, ignore.case = TRUE)) {
  readr::read_csv(infile, show_col_types = FALSE)
} else {
  readxl::read_excel(infile)
}

## 2) Clean & harmonize columns --------------------------------------------
# Expected columns (original names):
# Aut, year, M/U/R, ms, Male_Num, Num, prev_NCEP, wc/whr, htn, dm, high_tg, low_hdl, NOS_Score

dat <- dat_raw %>%
  janitor::clean_names() %>%   # -> aut, year, m_u_r, ms, male_num, num, prev_ncep, wc_whr, ...
  rename(
    author         = aut,
    site_code      = m_u_r,
    total_n        = num,
    prev_ncep_pct  = prev_ncep,
    wc_or_whr_pct  = wc_whr
  ) %>%
  mutate(
    author         = as.character(author),
    year           = suppressWarnings(as.integer(year)),
    site           = factor(site_code, levels = c(0,1,2),
                            labels = c("Mixed","Urban","Rural")),
    total_n        = suppressWarnings(as.integer(total_n)),
    ms             = suppressWarnings(readr::parse_number(as.character(ms))),
    prev_ncep_pct  = suppressWarnings(readr::parse_number(as.character(prev_ncep_pct))),
    wc_or_whr_pct  = suppressWarnings(readr::parse_number(as.character(wc_or_whr_pct))),
    htn            = suppressWarnings(readr::parse_number(as.character(htn))),
    dm             = suppressWarnings(readr::parse_number(as.character(dm))),
    high_tg        = suppressWarnings(readr::parse_number(as.character(high_tg))),
    low_hdl        = suppressWarnings(readr::parse_number(as.character(low_hdl))),
    nos_score      = suppressWarnings(readr::parse_number(as.character(nos_score)))
  )

## 3) Derive xi/ni (MetS cases and totals) ----------------------------------
derive_cases <- function(ms, prev_pct, n){
  if (is.na(n) || n <= 0) return(NA_integer_)
  if (!is.na(ms)) {
    if (ms <= 1)                 return(round(ms * n))        # proportion
    if (ms > 1 && ms <= 100 &&
        (is.na(prev_pct) || abs(ms - prev_pct) < 1e-6)) return(round(ms/100 * n)) # percent
    if (ms > 0 && ms <= n)       return(round(ms))            # count
  }
  if (is.na(ms) && !is.na(prev_pct)) return(round(prev_pct/100 * n))
  return(NA_integer_)
}

dat <- dat %>%
  rowwise() %>%
  mutate(
    xi = derive_cases(ms, prev_ncep_pct, total_n),
    ni = total_n
  ) %>%
  ungroup()

# Sanity checks
if (any(is.na(dat$xi) | is.na(dat$ni))) stop("Missing xi/ni after derivation.")
if (any(dat$xi > dat$ni)) stop("Found xi > ni; check inputs.")
if (any(dat$ni <= 0)) stop("Found ni <= 0; check inputs.")

## 4) Per-study prevalence table (exact binomial CI) ------------------------
stud_prev <- binom::binom.exact(dat$xi, dat$ni) %>%
  as_tibble() %>%
  transmute(prev = mean, ci_l = lower, ci_u = upper)

study_table <- dat %>%
  bind_cols(stud_prev) %>%
  mutate(prev_pct = 100*prev, ci_l_pct = 100*ci_l, ci_u_pct = 100*ci_u) %>%
  select(author, year, site, xi, ni, prev_pct, ci_l_pct, ci_u_pct,
         prev_ncep_pct, wc_or_whr_pct, high_tg, low_hdl, htn, dm, nos_score)

write_csv(study_table, file.path(out_dir, "01_study_characteristics_and_prevalence.csv"))

## 5) Helper for back-transformation & summary rows -------------------------
bt <- function(x, sm) meta::backtransf(x, sm = sm)
summarise_meta <- function(mdl, label){
  tibble(
    model     = label,
    k         = mdl$k,
    pooled    = bt(mdl$TE.random,    mdl$sm),
    ci_l      = bt(mdl$lower.random, mdl$sm),
    ci_u      = bt(mdl$upper.random, mdl$sm),
    pi_l      = bt(mdl$lower.predict,mdl$sm),
    pi_u      = bt(mdl$upper.predict,mdl$sm),
    I2        = mdl$I2,
    tau2      = mdl$tau2
  )
}

## 6) PRIMARY model: GLMM (ML) ---------------------------------------------
m_glmm <- metaprop(
  event = xi, n = ni, studlab = paste(dat$author, dat$year),
  data = dat,
  method = "GLMM",
  method.tau = "ML",          # REQUIRED for GLMM
  prediction = TRUE,
  title = "MetS prevalence (GLMM, ML)"
)
print(summary(m_glmm))

png(file.path(plot_dir, "forest_overall_GLMM_ML.png"), width = 1600, height = 1200, res = 200)
forest(m_glmm, prediction = TRUE, print.tau2 = TRUE, print.i2 = TRUE,
       leftcols = c("studlab","event","n"),
       leftlabs  = c("Study","MetS","N"),
       rightlabs = c("Prev [95% CI]"), smlab = "Pooled prevalence")
dev.off()

## 7) SENSITIVITY: Logit (Inverse-variance), REML + HK ----------------------
m_logit <- metaprop(
  event = xi, n = ni, studlab = paste(dat$author, dat$year),
  data = dat,
  sm = "PLOGIT",
  method = "Inverse",         # force inverse-variance (NOT GLMM)
  method.tau = "REML",
  method.random.ci = "HK",    # Hartung–Knapp
  prediction = TRUE,
  incr = 0.5, method.incr = "all",
  title = "MetS prevalence (logit, REML + HK)"
)
print(summary(m_logit))

png(file.path(plot_dir, "forest_overall_Logit_REML_HK.png"), width = 1600, height = 1200, res = 200)
forest(m_logit, prediction = TRUE, print.tau2 = TRUE, print.i2 = TRUE,
       leftcols = c("studlab","event","n"),
       leftlabs  = c("Study","MetS","N"),
       rightlabs = c("Prev [95% CI]"), smlab = "Pooled prevalence (logit)")
dev.off()

## 8) Subgroup (exploratory): by Site ---------------------------------------
m_glmm_site <- metaprop(
  event = xi, n = ni, studlab = paste(dat$author, dat$year),
  data = dat,
  method = "GLMM", method.tau = "ML", prediction = TRUE,
  byvar = site, print.byvar = TRUE,
  title = "MetS prevalence by site (GLMM, ML)"
)
print(summary(m_glmm_site))

# Export subgroup summary table
subgroup_tbl <- try({
  data.frame(
    Group      = m_glmm_site$bylevs,
    k          = m_glmm_site$k.w,
    pooled_pct = 100*bt(m_glmm_site$TE.random.w,    m_glmm_site$sm),
    ci_l_pct   = 100*bt(m_glmm_site$lower.random.w, m_glmm_site$sm),
    ci_u_pct   = 100*bt(m_glmm_site$upper.random.w, m_glmm_site$sm),
    I2         = m_glmm_site$I2.w,
    tau2       = m_glmm_site$tau2.w
  )
}, silent = TRUE)
if (!inherits(subgroup_tbl, "try-error")) {
  write_csv(subgroup_tbl, file.path(out_dir, "03_subgroup_by_site_GLMM.csv"))
}

png(file.path(plot_dir, "forest_by_site_GLMM_ML.png"), width = 1800, height = 1400, res = 200)
forest(m_glmm_site, prediction = TRUE, print.byvar = TRUE,
       print.i2 = TRUE, print.tau2 = TRUE,
       leftcols = c("studlab","event","n","byvar"),
       leftlabs  = c("Study","MetS","N","Site"),
       smlab = "Pooled prevalence by site")
dev.off()

## 9) Leave-one-out (GLMM) - manual & Logit metainf -------------------------
# Manual LOO for GLMM (robust)
loo_glmm_list <- lapply(seq_len(nrow(dat)), function(i) {
  fit <- metaprop(
    event = dat$xi[-i], n = dat$ni[-i],
    studlab = paste(dat$author[-i], dat$year[-i]),
    method = "GLMM", method.tau = "ML", prediction = TRUE
  )
  tibble(
    study_omitted = paste(dat$author[i], dat$year[i]),
    pooled        = bt(fit$TE.random,        fit$sm),
    ci_l          = bt(fit$lower.random,     fit$sm),
    ci_u          = bt(fit$upper.random,     fit$sm),
    pi_l          = bt(fit$lower.predict,    fit$sm),
    pi_u          = bt(fit$upper.predict,    fit$sm),
    I2            = fit$I2,
    tau2          = fit$tau2
  )
})
loo_glmm <- bind_rows(loo_glmm_list)
write_csv(loo_glmm, file.path(out_dir, "04_leave_one_out_GLMM_ML.csv"))

## Leave-one-out forest (no sorting)
png(file.path(plot_dir, "leave_one_out_Logit_REML_HK.png"),
    width = 1600, height = 1200, res = 200)
meta::forest(loo_logit,
             pooled.totals = FALSE,
             leftlabs = "Study omitted",
             smlab = "Pooled prevalence (logit model)")
dev.off()


png(file.path(plot_dir, "baujat_Logit_REML_HK.png"), width = 1400, height = 1200, res = 190)
baujat(m_logit)
dev.off()

## 10) Cumulative series (manual, GLMM & Logit) -----------------------------
# GLMM manual cumulative
dat_ord <- dat %>% arrange(year)
cum_glmm_list <- lapply(seq_len(nrow(dat_ord)), function(k){
  fit <- metaprop(
    event = dat_ord$xi[1:k], n = dat_ord$ni[1:k],
    studlab = paste(dat_ord$author[1:k], dat_ord$year[1:k]),
    method = "GLMM", method.tau = "ML", prediction = TRUE
  )
  tibble(
    k_included = k,
    last_study = paste(dat_ord$author[k], dat_ord$year[k]),
    pooled     = bt(fit$TE.random,      fit$sm),
    ci_l       = bt(fit$lower.random,   fit$sm),
    ci_u       = bt(fit$upper.random,   fit$sm),
    pi_l       = bt(fit$lower.predict,  fit$sm),
    pi_u       = bt(fit$upper.predict,  fit$sm),
    I2         = fit$I2,
    tau2       = fit$tau2
  )
})
cum_glmm <- bind_rows(cum_glmm_list)
write_csv(cum_glmm, file.path(out_dir, "05_cumulative_meta_GLMM_ML.csv"))

png(file.path(plot_dir, "cumulative_meta_GLMM_ML.png"), width = 1600, height = 1200, res = 200)
plot(cum_glmm$k_included, 100*cum_glmm$pooled, type = "b", pch = 16,
     xlab = "Studies included (by year order)",
     ylab = "Cumulative pooled prevalence (%)",
     main = "Cumulative meta-analysis (GLMM, ML)")
arrows(cum_glmm$k_included, 100*cum_glmm$ci_l,
       cum_glmm$k_included, 100*cum_glmm$ci_u,
       angle = 90, code = 3, length = 0.03)
dev.off()

# Logit manual cumulative
cum_logit_list <- lapply(seq_len(nrow(dat_ord)), function(k){
  fit <- metaprop(
    event = dat_ord$xi[1:k], n = dat_ord$ni[1:k],
    studlab = paste(dat_ord$author[1:k], dat_ord$year[1:k]),
    sm = "PLOGIT", method = "Inverse", method.tau = "REML",
    method.random.ci = "HK", prediction = TRUE,
    incr = 0.5, method.incr = "all"
  )
  tibble(
    k_included = k,
    last_study = paste(dat_ord$author[k], dat_ord$year[k]),
    pooled     = bt(fit$TE.random,      fit$sm),
    ci_l       = bt(fit$lower.random,   fit$sm),
    ci_u       = bt(fit$upper.random,   fit$sm),
    pi_l       = bt(fit$lower.predict,  fit$sm),
    pi_u       = bt(fit$upper.predict,  fit$sm),
    I2         = fit$I2,
    tau2       = fit$tau2
  )
})
cum_logit <- bind_rows(cum_logit_list)
write_csv(cum_logit, file.path(out_dir, "05_cumulative_meta_Logit_REML_HK.csv"))

png(file.path(plot_dir, "cumulative_meta_Logit_REML_HK.png"), width = 1600, height = 1200, res = 200)
plot(cum_logit$k_included, 100*cum_logit$pooled, type = "b", pch = 16,
     xlab = "Studies included (by year order)",
     ylab = "Cumulative pooled prevalence (%)",
     main = "Cumulative meta-analysis (Logit, REML + HK)")
arrows(cum_logit$k_included, 100*cum_logit$ci_l,
       cum_logit$k_included, 100*cum_logit$ci_u,
       angle = 90, code = 3, length = 0.03)
dev.off()

## 11) Component meta-analyses (GLMM, ML) -----------------------------------
comp_map <- c(
  wc_or_whr_pct = "High waist (WC/WHR)",
  high_tg       = "High triglycerides",
  low_hdl       = "Low HDL",
  htn           = "Hypertension",
  dm            = "Diabetes / high glucose"
)

comp_rows   <- list()
comp_models <- list()

for (v in names(comp_map)) {
  if (!v %in% names(dat)) next
  tmp <- dat %>%
    mutate(comp_pct = suppressWarnings(readr::parse_number(as.character(.data[[v]]))),
           comp_evt = ifelse(!is.na(comp_pct) & !is.na(ni), round(comp_pct/100 * ni), NA_integer_)) %>%
    filter(!is.na(comp_evt), !is.na(ni), comp_evt <= ni)
  
  if (nrow(tmp) >= 2) {
    m_comp <- metaprop(
      event = comp_evt, n = ni, studlab = paste(tmp$author, tmp$year),
      data = tmp, method = "GLMM", method.tau = "ML", prediction = TRUE,
      title = paste0(comp_map[[v]], " (GLMM, ML)")
    )
    comp_models[[v]] <- m_comp
    srow <- summarise_meta(m_comp, comp_map[[v]])
    comp_rows[[v]] <- srow
    
    png(file.path(plot_dir, paste0("forest_component_", gsub("[^A-Za-z0-9]+","_", v), "_GLMM_ML.png")),
        width = 1600, height = 1200, res = 200)
    forest(m_comp, prediction = TRUE, print.i2 = TRUE, print.tau2 = TRUE,
           leftcols = c("studlab","event","n"),
           leftlabs  = c("Study","Events","N"),
           smlab = comp_map[[v]])
    dev.off()
  }
}

comp_tbl <- if (length(comp_rows)) bind_rows(comp_rows) else tibble()
if (nrow(comp_tbl)) {
  comp_tbl_pct <- comp_tbl %>%
    mutate(across(c(pooled, ci_l, ci_u, pi_l, pi_u), ~ . * 100)) %>%
    rename(pooled_pct = pooled, ci_l_pct = ci_l, ci_u_pct = ci_u,
           pi_l_pct = pi_l, pi_u_pct = pi_u)
  write_csv(comp_tbl_pct, file.path(out_dir, "06_components_pooled_GLMM_ML.csv"))
  
  # Simple summary bar
  pdat <- comp_tbl_pct %>% arrange(desc(pooled_pct)) %>% mutate(component = factor(model, levels = model))
  png(file.path(plot_dir, "components_summary_bar_GLMM_ML.png"), width = 1600, height = 1100, res = 200)
  ggplot(pdat, aes(x = component, y = pooled_pct)) +
    geom_col() +
    geom_errorbar(aes(ymin = ci_l_pct, ymax = ci_u_pct), width = 0.2) +
    coord_flip() +
    labs(x = NULL, y = "Pooled prevalence (%)",
         title = "Metabolic syndrome components (pooled, GLMM ML)") +
    theme_minimal(base_size = 12)
  dev.off()
}

## 12) Export consolidated workbook & summary text --------------------------
main_summ <- bind_rows(
  summarise_meta(m_glmm,  "GLMM ML"),
  summarise_meta(m_logit, "Logit REML + HK")
) %>%
  mutate(across(c(pooled, ci_l, ci_u, pi_l, pi_u), ~ . * 100)) %>%
  rename(pooled_pct = pooled, ci_l_pct = ci_l, ci_u_pct = ci_u,
         pi_l_pct = pi_l, pi_u_pct = pi_u)

write_csv(main_summ, file.path(out_dir, "02_pooled_summaries_main.csv"))

wb <- list(
  "Study_Characteristics" = study_table,
  "Pooled_Summaries"      = main_summ,
  "Subgroup_Site"         = if (exists("subgroup_tbl")) subgroup_tbl else NULL,
  "LeaveOneOut_GLMM"      = loo_glmm,
  "Cumulative_GLMM"       = cum_glmm,
  "Cumulative_Logit"      = cum_logit,
  "Components_Pooled"     = if (exists("comp_tbl_pct")) comp_tbl_pct else NULL
)
# Remove NULL sheets
wb <- wb[!vapply(wb, is.null, logical(1))]
writexl::write_xlsx(wb, path = file.path(out_dir, "ALL_RESULTS_MetS_Nepal.xlsx"))

sink(file.path(out_dir, "SUMMARY_README.txt"))
cat("Meta-analysis of Metabolic Syndrome in Nepal (k =", m_glmm$k, ")\n",
    "PRIMARY: GLMM (ML). SENSITIVITY: Logit (REML + HK).\n",
    "Very high heterogeneity: emphasize prediction intervals and ranges.\n\n", sep="")
cat("== GLMM (ML) ==\n")
cat(sprintf("Pooled: %.1f%% [%.1f–%.1f], PI [%.1f–%.1f]; I2=%.1f%%; tau^2=%.4f\n\n",
            100*bt(m_glmm$TE.random, m_glmm$sm),
            100*bt(m_glmm$lower.random, m_glmm$sm),
            100*bt(m_glmm$upper.random, m_glmm$sm),
            100*bt(m_glmm$lower.predict, m_glmm$sm),
            100*bt(m_glmm$upper.predict, m_glmm$sm),
            m_glmm$I2, m_glmm$tau2))
cat("== Logit (REML + HK) ==\n")
cat(sprintf("Pooled: %.1f%% [%.1f–%.1f], PI [%.1f–%.1f]; I2=%.1f%%; tau^2=%.4f\n\n",
            100*bt(m_logit$TE.random, m_logit$sm),
            100*bt(m_logit$lower.random, m_logit$sm),
            100*bt(m_logit$upper.random, m_logit$sm),
            100*bt(m_logit$lower.predict, m_logit$sm),
            100*bt(m_logit$upper.predict, m_logit$sm),
            m_logit$I2, m_logit$tau2))
cat("Caveat: k small & I2 high; pooled mean is fragile; report prediction intervals prominently.\n")
sink()

message("Done. CSVs and plots saved in: ", out_dir, " / ", plot_dir)
## Overall forest: Logit REML + HK
pdf(file.path(plot_dir, "forest_overall_Logit_REML_HK.pdf"),
    width = 12, height = 10)   # bigger canvas
forest(m_logit, prediction = TRUE, print.tau2 = TRUE, print.i2 = TRUE,
       leftcols = c("studlab","event","n"),
       leftlabs  = c("Study","MetS","N"),
       rightlabs = c("Prev [95% CI]"),
       smlab = "Pooled prevalence (logit, REML + HK)")
dev.off()

tiff(file.path(plot_dir, "forest_overall_Logit_REML_HK.tiff"),
     width = 3200, height = 2400, res = 400, compression = "lzw")
forest(m_logit, prediction = TRUE, print.tau2 = TRUE, print.i2 = TRUE,
       leftcols = c("studlab","event","n"),
       leftlabs  = c("Study","MetS","N"),
       rightlabs = c("Prev [95% CI]"),
       smlab = "Pooled prevalence (logit, REML + HK)")
dev.off()
## Leave-one-out: Logit REML + HK
loo_logit <- metainf(m_logit)

pdf(file.path(plot_dir, "leave_one_out_Logit_REML_HK.pdf"),
    width = 12, height = 10)
forest(loo_logit,
       pooled.totals = FALSE,
       leftlabs = "Study omitted",
       smlab = "Pooled prevalence (logit model)")
dev.off()

tiff(file.path(plot_dir, "leave_one_out_Logit_REML_HK.tiff"),
     width = 3200, height = 2400, res = 400, compression = "lzw")
forest(loo_logit,
       pooled.totals = FALSE,
       leftlabs = "Study omitted",
       smlab = "Pooled prevalence (logit model)")
dev.off()

