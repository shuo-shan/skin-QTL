library(data.table)

# ----------------------------
# User inputs / metadata
# ----------------------------
trait   <- "systemic_lupus_erythematosus"
infile  <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/systemic_lupus_erythematosus/GCST90476183.h.tsv.gz"
outfile <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/systemic_lupus_erythematosus/standardized_GCST90476183.h.tsv.gz"

type <- "cc"          # "cc" or "quant"
prefer_snp <- "rsid"  # "rsid" or "variant_id"

# Optional overrides (set if needed)
N_override <- NA_real_

# (cc only; optional) If provided, used only when file lacks usable cases/controls
cases_total <- NA_real_
controls_total <- NA_real_

# (quant only)
sdY_override <- NA_real_

# ----------------------------
# Standardization function
# ----------------------------
standardize_gwas_for_coloc <- function(gwas_path, trait_name,
                                       type = c("cc","quant"),
                                       prefer_snp = c("rsid","variant_id"),
                                       sdY_override = NA_real_,
                                       N_override = NA_real_,
                                       cases_total = NA_real_,
                                       controls_total = NA_real_) {
  
  type <- match.arg(type)
  prefer_snp <- match.arg(prefer_snp)
  
  dt <- fread(gwas_path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
  cn <- names(dt)
  
  pick <- function(...) {
    cand <- c(...)
    hit <- cand[cand %in% cn]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }
  
  # ---- Core columns ----
  chr_col <- pick("chromosome","chr","CHROM","hm_chrom")
  pos_col <- pick("base_pair_location","pos","POS","hm_pos")
  ea_col  <- pick("effect_allele","ea","EA","A1","ALT","alt","hm_effect_allele")
  oa_col  <- pick("other_allele","oa","OA","A2","REF","ref","hm_other_allele")
  
  p_col   <- pick("p_value","p","P","pval","pval_nominal")
  eaf_col <- pick("effect_allele_frequency","eaf","EAF","af","hm_effect_allele_frequency")
  
  rsid_col <- pick("rsid","hm_rsid")
  vid_col  <- pick("variant_id","hm_variant_id")
  
  snp_col <- if (prefer_snp == "variant_id") {
    pick("variant_id","hm_variant_id","rsid","hm_rsid")
  } else {
    pick("rsid","hm_rsid","variant_id","hm_variant_id")
  }
  
  # Effect sizes
  beta_col <- pick("beta","hm_beta")
  or_col   <- pick("odds_ratio","OR","or","hm_odds_ratio")
  
  # SE sources: explicit SE > CI > p+beta
  se_col    <- pick("standard_error","se","SE","stderr","hm_standard_error","hm_se")
  ci_lo_col <- pick("ci_lower","lower_ci","lci","ci_l","hm_ci_lower")
  ci_hi_col <- pick("ci_upper","upper_ci","uci","ci_u","hm_ci_upper")
  
  # N and case/control
  n_col     <- pick("cum_eff_sample_size","effective_N","n_eff","n","N")
  ncase_col <- pick("num_cases","cases","ncase","N_cases")
  nctrl_col <- pick("num_controls","controls","ncontrol","N_controls")
  
  # ---- Required checks ----
  req <- c(chr_col, pos_col, ea_col, oa_col, p_col)
  if (any(is.na(req))) {
    stop("Missing required cols in ", gwas_path, "\nHave: ", paste(cn, collapse = ", "))
  }
  if (is.na(beta_col) && is.na(or_col)) {
    stop("No beta or OR column found in ", gwas_path, "\nHave: ", paste(cn, collapse = ", "))
  }
  
  # ---- Build base output (all standardized columns exist) ----
  out <- dt[, .(
    trait        = trait_name,
    type         = type,
    snp          = if (!is.na(snp_col))  as.character(get(snp_col))  else NA_character_,
    rsid         = if (!is.na(rsid_col)) as.character(get(rsid_col)) else NA_character_,
    variant_id   = if (!is.na(vid_col))  as.character(get(vid_col))  else NA_character_,
    chr          = suppressWarnings(as.integer(get(chr_col))),
    pos          = suppressWarnings(as.integer(get(pos_col))),
    ea           = as.character(get(ea_col)),
    oa           = as.character(get(oa_col)),
    beta         = NA_real_,
    se           = NA_real_,
    varbeta      = NA_real_,
    p            = suppressWarnings(as.numeric(get(p_col))),
    effect_scale = NA_character_,
    se_source    = NA_character_,
    eaf          = if (!is.na(eaf_col)) suppressWarnings(as.numeric(get(eaf_col))) else NA_real_,
    MAF          = NA_real_,
    N            = NA_real_,
    N_source     = NA_character_,
    num_cases    = if (!is.na(ncase_col)) suppressWarnings(as.numeric(get(ncase_col))) else NA_real_,
    num_controls = if (!is.na(nctrl_col)) suppressWarnings(as.numeric(get(nctrl_col))) else NA_real_,
    s            = NA_real_,
    s_source     = NA_character_,
    sdY          = NA_real_,
    sdY_source   = NA_character_
  )]
  
  # Normalize blanks -> NA (true NA; later fwrite(na="NA") prints NA)
  for (ccol in c("snp","rsid","variant_id","ea","oa")) {
    out[get(ccol) == "", (ccol) := NA_character_]
  }
  
  # ---- beta + effect_scale ----
  if (!is.na(beta_col)) {
    out[, beta := suppressWarnings(as.numeric(dt[[beta_col]]))]
    out[, effect_scale := "beta"]
  } else {
    OR <- suppressWarnings(as.numeric(dt[[or_col]]))
    OR[!is.finite(OR) | OR <= 0] <- NA_real_
    out[, beta := log(OR)]
    out[, effect_scale := "logOR"]
  }
  
  # ---- se + se_source (priority: SE > CI > p+beta) ----
  out[, se := NA_real_]
  out[, se_source := NA_character_]
  
  # 1) explicit SE column (only where finite >0)
  if (!is.na(se_col)) {
    se0 <- suppressWarnings(as.numeric(dt[[se_col]]))
    idx <- which(is.finite(se0) & se0 > 0)
    if (length(idx) > 0) {
      out[idx, se := se0[idx]]
      out[idx, se_source := "SE"]
    }
  }
  
  # 2) CI (assume 95%; infer on log scale)
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
  
  # 3) p + beta: se = |beta| / z, z from two-sided p
  z <- suppressWarnings(qnorm(out$p / 2, lower.tail = FALSE))
  se_pb <- abs(out$beta) / z
  idx <- which(is.na(out$se) & is.finite(se_pb) & se_pb > 0)
  if (length(idx) > 0) {
    out[idx, se := se_pb[idx]]
    out[idx, se_source := "p+beta"]
  }
  
  # varbeta
  out[, varbeta := se^2]
  
  # ---- eaf/MAF ----
  out[, MAF := ifelse(is.na(eaf), NA_real_, pmin(eaf, 1 - eaf))]
  
  # ---- N + N_source (priority: per-SNP N col > derived cases/controls > override) ----
  if (!is.na(n_col)) {
    N0 <- suppressWarnings(as.numeric(dt[[n_col]]))
    idx <- which(is.finite(N0) & N0 > 0)
    if (length(idx) > 0) {
      out[idx, N := N0[idx]]
      out[idx, N_source := n_col]
    }
  }
  
  idxN <- which(!(is.finite(out$N) & out$N > 0) &
                  is.finite(out$num_cases) & out$num_cases > 0 &
                  is.finite(out$num_controls) & out$num_controls > 0)
  if (length(idxN) > 0) {
    out[idxN, N := num_cases + num_controls]
    out[idxN, N_source := "derived_cases_controls"]
  }
  
  if (is.finite(N_override) && N_override > 0) {
    idxO <- which(!(is.finite(out$N) & out$N > 0))
    if (length(idxO) > 0) {
      out[idxO, N := N_override]
      out[idxO, N_source := "N_override"]
    }
  }
  
  # ---- s + s_source (cc only; with per-SNP vs repeated-total detection) ----
  if (type == "cc") {
    out[, s := NA_real_]
    out[, s_source := "NA"]
    
    nc_vals <- unique(out[is.finite(num_cases), num_cases])
    nt_vals <- unique(out[is.finite(num_controls), num_controls])
    have_cases_ctrl_cols <- (length(nc_vals) > 0 && length(nt_vals) > 0)
    
    if (have_cases_ctrl_cols && length(nc_vals) == 1 && length(nt_vals) == 1) {
      denom <- out$num_cases + out$num_controls
      out[is.finite(denom) & denom > 0, s := num_cases / denom]
      out[is.finite(s), s_source := "study_total"]
    } else if (have_cases_ctrl_cols && (length(nc_vals) > 1 || length(nt_vals) > 1)) {
      denom <- out$num_cases + out$num_controls
      out[is.finite(denom) & denom > 0, s := num_cases / denom]
      out[is.finite(s), s_source := "per_snp"]
    } else if (is.finite(cases_total) && is.finite(controls_total) && cases_total > 0 && controls_total > 0) {
      out[, s := cases_total / (cases_total + controls_total)]
      out[, s_source := "study_total"]
    } else {
      out[, s := NA_real_]
      out[, s_source := "NA"]
    }
    
    out[, sdY := NA_real_]
    out[, sdY_source := "NA"]
  } else {
    out[, s := NA_real_]
    out[, s_source := "NA"]
    
    if (is.finite(sdY_override)) {
      out[, sdY := sdY_override]
      out[, sdY_source := "override"]
    } else {
      out[, sdY := NA_real_]
      out[, sdY_source := "NA"]
    }
  }
  
  # ---- Cleaning ----
  out <- out[is.finite(beta)]
  out <- out[is.finite(se) & se > 0]
  out <- out[is.finite(varbeta) & varbeta > 0]
  out <- out[is.finite(p) & p > 0 & p <= 1]
  out <- out[is.finite(chr) & is.finite(pos)]
  out <- out[!is.na(ea) & ea != "" & !is.na(oa) & oa != ""]
  out <- out[is.finite(N) & N > 0]
  out <- out[is.na(MAF) | (is.finite(MAF) & MAF > 0 & MAF < 0.5)]
  
  if (type == "cc") {
    had_totals <- FALSE
    if ((is.finite(cases_total) && is.finite(controls_total) && cases_total > 0 && controls_total > 0)) had_totals <- TRUE
    if (length(unique(out[is.finite(num_cases), num_cases])) == 1 &&
        length(unique(out[is.finite(num_controls), num_controls])) == 1 &&
        nrow(out[is.finite(num_cases) & is.finite(num_controls)]) > 0) {
      had_totals <- TRUE
    }
    if (had_totals) {
      out <- out[is.finite(s) & s > 0 & s < 1]
    } else {
      out <- out[is.na(s) | (is.finite(s) & s > 0 & s < 1)]
    }
  }
  
  # ---- Final: replace blank strings with "NA" (requested) ----
  # NOTE: fwrite(na="NA") prints true NA as "NA". Here we also eliminate empty strings.
  for (ccol in c("trait","type","snp","rsid","variant_id","ea","oa",
                 "effect_scale","se_source","N_source","s_source","sdY_source")) {
    out[is.na(get(ccol)) | get(ccol) == "", (ccol) := "NA"]
  }
  
  # ---- Enforce exact column order (sdY + sdY_source at end) ----
  std_cols <- c(
    "trait","type","snp","rsid","variant_id","chr","pos","ea","oa",
    "beta","se","varbeta","p","effect_scale","se_source",
    "eaf","MAF","N","N_source",
    "num_cases","num_controls","s","s_source",
    "sdY","sdY_source"
  )
  out <- out[, ..std_cols]
  out[]
}

# ----------------------------
# Run + write
# ----------------------------
out <- standardize_gwas_for_coloc(
  gwas_path      = infile,
  trait_name     = trait,
  type           = type,
  prefer_snp     = prefer_snp,
  sdY_override   = sdY_override,
  N_override     = N_override,
  cases_total    = cases_total,
  controls_total = controls_total
)

fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("[%s] wrote %d rows -> %s\n", trait, nrow(out), outfile))
