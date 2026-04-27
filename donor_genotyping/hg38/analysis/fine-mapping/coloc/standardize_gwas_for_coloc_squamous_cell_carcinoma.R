library(data.table)

# ---------------------------
# User inputs (this trait)
# ---------------------------
trait_name <- "squamous_cell_carcinoma"
type <- "cc"              # "cc" or "quant"
prefer_snp <- "rsid"      # "rsid" or "variant_id"

infile  <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/squamous_cell_carcinoma/GCST90475583.h.tsv.gz"
outfile <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/squamous_cell_carcinoma/standardized_GCST90475583.h.tsv.gz"

# Optional overrides (use if needed)
N_override <- NA_real_
cases_total <- NA_real_     # e.g., 19217
controls_total <- NA_real_  # e.g., 419674
sdY_override <- NA_real_    # quant only (e.g., 1 for invnorm)

# ---------------------------
# Standardizer
# ---------------------------
standardize_gwas_for_coloc <- function(gwas_path, trait_name,
                                       type = c("cc","quant"),
                                       prefer_snp = c("rsid","variant_id"),
                                       N_override = NA_real_,
                                       cases_total = NA_real_,
                                       controls_total = NA_real_,
                                       sdY_override = NA_real_) {
  type <- match.arg(type)
  prefer_snp <- match.arg(prefer_snp)
  
  dt <- fread(
    gwas_path,
    sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE,
    na.strings = c("", "NA", "NaN", "nan", ".")
  )
  cn <- names(dt)
  pick <- function(...) { c(...)[c(... ) %in% cn][1] }
  
  # ---- core fields ----
  chr_col <- pick("chromosome","hm_chrom","chr","CHROM")
  pos_col <- pick("base_pair_location","hm_pos","pos","POS")
  
  ea_col  <- pick("effect_allele","hm_effect_allele","A1","EA","alt","ALT")
  oa_col  <- pick("other_allele","hm_other_allele","A2","OA","ref","REF")
  
  p_col   <- pick("p_value","p","P","pval")
  
  # EAF: prefer effect_allele_frequency; do NOT fall back to case/control AF for eaf
  eaf_col <- pick("effect_allele_frequency","hm_effect_allele_frequency","eaf","EAF","af")
  
  n_col   <- pick("cum_eff_sample_size","n","N","neff","effective_n","effective_sample_size")
  
  ncase_col <- pick("num_cases","cases","ncase","N_cases")
  nctrl_col <- pick("num_controls","controls","ncontrol","N_controls")
  
  rsid_col <- pick("rsid","hm_rsid")
  vid_col  <- pick("variant_id","hm_variant_id")
  
  snp_col <- if (prefer_snp == "variant_id") {
    pick("variant_id","hm_variant_id","rsid","hm_rsid")
  } else {
    pick("rsid","hm_rsid","variant_id","hm_variant_id")
  }
  
  beta_col <- pick("hm_beta","beta","BETA","effect","estimate")
  or_col   <- pick("odds_ratio","hm_odds_ratio","OR","or")
  
  se_col    <- pick("standard_error","se","SE","stderr","standarderror")
  ci_lo_col <- pick("ci_lower","hm_ci_lower","lower_ci","ci_l","lci","LCI")
  ci_hi_col <- pick("ci_upper","hm_ci_upper","upper_ci","ci_u","uci","UCI")
  
  # required columns
  req <- c(chr_col, pos_col, ea_col, oa_col, p_col, snp_col)
  if (any(is.na(req))) {
    stop("Missing required cols in ", gwas_path, "\nHave: ", paste(cn, collapse = ", "))
  }
  if (is.na(beta_col) && is.na(or_col)) {
    stop("No beta or OR column found in ", gwas_path)
  }
  
  # ---- start standardized table (build from dt, keep only standardized columns later) ----
  out <- dt[, .(
    trait = trait_name,
    type  = type,
    snp   = as.character(get(snp_col)),
    rsid  = if (!is.na(rsid_col)) as.character(get(rsid_col)) else NA_character_,
    variant_id = if (!is.na(vid_col)) as.character(get(vid_col)) else NA_character_,
    chr   = suppressWarnings(as.integer(get(chr_col))),
    pos   = suppressWarnings(as.integer(get(pos_col))),
    ea    = as.character(get(ea_col)),
    oa    = as.character(get(oa_col)),
    p     = suppressWarnings(as.numeric(get(p_col))),
    eaf   = if (!is.na(eaf_col)) suppressWarnings(as.numeric(get(eaf_col))) else NA_real_,
    num_cases    = if (!is.na(ncase_col)) suppressWarnings(as.numeric(get(ncase_col))) else NA_real_,
    num_controls = if (!is.na(nctrl_col)) suppressWarnings(as.numeric(get(nctrl_col))) else NA_real_,
    N_raw = if (!is.na(n_col)) suppressWarnings(as.numeric(get(n_col))) else NA_real_
  )]
  
  # ---- beta + effect_scale ----
  if (!is.na(beta_col)) {
    out[, beta := suppressWarnings(as.numeric(dt[[beta_col]]))]
    out[, effect_scale := "beta"]
  } else {
    OR <- suppressWarnings(as.numeric(dt[[or_col]]))
    out[, beta := suppressWarnings(log(OR))]
    out[, effect_scale := "logOR"]
  }
  
  # ---- se + se_source (priority: SE -> CI -> p+beta) ----
  out[, se := NA_real_]
  out[, se_source := NA_character_]
  
  # 1) explicit SE
  if (!is.na(se_col)) {
    se0 <- suppressWarnings(as.numeric(dt[[se_col]]))
    idx <- which(is.finite(se0) & se0 > 0)
    if (length(idx) > 0) {
      out[idx, se := se0[idx]]
      out[idx, se_source := "SE"]
    }
  }
  
  # 2) CI (assume 95%, infer on log scale; valid for OR-type effects)
  if (!is.na(ci_lo_col) && !is.na(ci_hi_col)) {
    lo <- suppressWarnings(as.numeric(dt[[ci_lo_col]]))
    hi <- suppressWarnings(as.numeric(dt[[ci_hi_col]]))
    ok_ci <- is.finite(lo) & is.finite(hi) & lo > 0 & hi > 0
    se_ci <- (log(hi) - log(lo)) / (2 * 1.96)
    
    idx <- which(is.na(out$se) & ok_ci & is.finite(se_ci) & se_ci > 0)
    if (length(idx) > 0) {
      out[idx, se := se_ci[idx]]
      out[idx, se_source := "CI"]
    }
  }
  
  # 3) p + beta fallback
  if (any(is.na(out$se))) {
    p <- out$p
    b <- out$beta
    z <- suppressWarnings(qnorm(p / 2, lower.tail = FALSE))
    se_pb <- abs(b) / z
    idx <- which(is.na(out$se) & is.finite(p) & p > 0 & p <= 1 &
                   is.finite(b) & is.finite(z) & z > 0 &
                   is.finite(se_pb) & se_pb > 0)
    if (length(idx) > 0) {
      out[idx, se := se_pb[idx]]
      out[idx, se_source := "p+beta"]
    }
  }
  
  # varbeta
  out[, varbeta := se^2]
  
  # ---- MAF ----
  out[, MAF := ifelse(is.na(eaf), NA_real_, pmin(eaf, 1 - eaf))]
  
  # ---- N + N_source ----
  out[, N := NA_real_]
  out[, N_source := NA_character_]
  
  # Prefer per-SNP effective N if provided
  if (!is.na(n_col)) {
    N0 <- out$N_raw
    idx <- which(is.finite(N0) & N0 > 0)
    if (length(idx) > 0) {
      out[idx, N := N0[idx]]
      out[idx, N_source := n_col]
    }
  }
  
  # If N missing/bad and per-SNP num_cases/num_controls exist: derive N = cases+controls
  idx_needN <- which(!(is.finite(out$N) & out$N > 0))
  if (length(idx_needN) > 0) {
    nc <- out$num_cases
    nt <- out$num_controls
    N_der <- nc + nt
    idx <- idx_needN[is.finite(N_der[idx_needN]) & N_der[idx_needN] > 0]
    if (length(idx) > 0) {
      out[idx, N := N_der[idx]]
      out[idx, N_source := "derived_cases_controls"]
    }
  }
  
  # If still missing/bad: N_override
  idx_needN <- which(!(is.finite(out$N) & out$N > 0))
  if (length(idx_needN) > 0 && is.finite(N_override) && N_override > 0) {
    out[idx_needN, N := N_override]
    out[idx_needN, N_source := "N_override"]
  }
  
  # ---- s + s_source ----
  out[, s := NA_real_]
  out[, s_source := "NA"]
  
  if (type == "cc") {
    uniq_nonmiss <- function(x) unique(x[is.finite(x)])
    u_cases <- uniq_nonmiss(out$num_cases)
    u_ctrls <- uniq_nonmiss(out$num_controls)
    
    has_counts_cols <- !is.na(ncase_col) && !is.na(nctrl_col)
    
    if (has_counts_cols && (length(u_cases) > 1 || length(u_ctrls) > 1)) {
      ok <- is.finite(out$num_cases) & is.finite(out$num_controls) & (out$num_cases + out$num_controls) > 0
      out[ok, s := num_cases / (num_cases + num_controls)]
      out[ok, s_source := "per_snp"]
    } else if (has_counts_cols && length(u_cases) == 1 && length(u_ctrls) == 1) {
      ok <- is.finite(out$num_cases) & is.finite(out$num_controls) & (out$num_cases + out$num_controls) > 0
      out[ok, s := num_cases / (num_cases + num_controls)]
      out[ok, s_source := "study_total"]
    } else if (is.finite(cases_total) && is.finite(controls_total) && (cases_total + controls_total) > 0) {
      out[, num_cases := cases_total]
      out[, num_controls := controls_total]
      out[, s := cases_total / (cases_total + controls_total)]
      out[, s_source := "study_total"]
    } else {
      out[, s := NA_real_]
      out[, s_source := "NA"]
    }
  } else {
    out[, s := NA_real_]
    out[, s_source := "NA"]
  }
  
  # ---- sdY + sdY_source (sdY only for quant) ----
  out[, sdY := NA_real_]
  out[, sdY_source := "NA"]
  if (type == "quant") {
    if (is.finite(sdY_override)) {
      out[, sdY := sdY_override]
      out[, sdY_source := "override"]
    }
  }
  
  # ---- replace blanks with NA (character columns) ----
  for (ccol in c("snp","rsid","variant_id","ea","oa","effect_scale","se_source","N_source","s_source","sdY_source")) {
    if (ccol %in% names(out)) {
      out[get(ccol) == "", (ccol) := NA_character_]
    }
  }
  
  # ---- cleaning ----
  out <- out[
    is.finite(beta) &
      is.finite(se) & se > 0 &
      is.finite(varbeta) & varbeta > 0 &
      is.finite(p) & p > 0 & p <= 1 &
      is.finite(chr) & is.finite(pos) &
      !is.na(ea) & ea != "" &
      !is.na(oa) & oa != "" &
      is.finite(N) & N > 0
  ]
  
  # do NOT require eaf/MAF; but if MAF present require (0,0.5)
  out <- out[is.na(MAF) | (is.finite(MAF) & MAF > 0 & MAF < 0.5)]
  
  # for cc: require s in (0,1) if provided; allow NA only if no study totals were given
  if (type == "cc") {
    have_any_totals <- FALSE
    if (!is.na(ncase_col) && !is.na(nctrl_col)) {
      u_cases2 <- unique(out$num_cases[is.finite(out$num_cases)])
      u_ctrls2 <- unique(out$num_controls[is.finite(out$num_controls)])
      if (length(u_cases2) >= 1 && length(u_ctrls2) >= 1) have_any_totals <- TRUE
    }
    if (is.finite(cases_total) && is.finite(controls_total) && (cases_total + controls_total) > 0) {
      have_any_totals <- TRUE
    }
    if (have_any_totals) {
      out <- out[is.finite(s) & s > 0 & s < 1]
    } else {
      out <- out[is.na(s) | (is.finite(s) & s > 0 & s < 1)]
    }
  }
  
  # ---- enforce final standardized column order (DO NOT CHANGE) ----
  out_final <- out[, .(
    trait,
    type,
    snp,
    rsid,
    variant_id,
    chr,
    pos,
    ea,
    oa,
    beta,
    se,
    varbeta,
    p,
    effect_scale,
    se_source,
    eaf,
    MAF,
    N,
    N_source,
    num_cases,
    num_controls,
    s,
    s_source,
    sdY,
    sdY_source
  )]
  
  out_final[]
}

# ---------------------------
# Run + write
# ---------------------------
out <- standardize_gwas_for_coloc(
  gwas_path   = infile,
  trait_name  = trait_name,
  type        = type,
  prefer_snp  = prefer_snp,
  N_override  = N_override,
  cases_total    = cases_total,
  controls_total = controls_total,
  sdY_override   = sdY_override
)

fwrite(out, outfile, sep = "\t", na = "NA", quote = FALSE)
cat(sprintf("[%s] wrote %d rows -> %s\n", trait_name, nrow(out), outfile))
