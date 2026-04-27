library(data.table)

# -----------------------------
# User inputs / metadata
# -----------------------------
trait_name  <- "crohns_disease"
type        <- "cc"                 # "cc" or "quant"
prefer_snp  <- "rsid"               # "rsid" or "variant_id"

infile  <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/crohns_disease/GCST90475318.h.tsv.gz"
outfile <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/crohns_disease/standardized_GCST90475318.h.tsv.gz"

# Optional overrides (set to NA if not used)
N_override      <- NA_real_
cases_total     <- NA_real_   # only used if cc and file lacks usable case/control cols
controls_total  <- NA_real_
sdY_override    <- NA_real_   # quant only (set e.g. 1 for invnorm/z-scored traits)

# -----------------------------
# Standardizer
# -----------------------------
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
    hits <- c(...)
    hits <- hits[hits %in% cn]
    if (length(hits) == 0) return(NA_character_)
    hits[1]
  }
  
  # --- Map columns from this header ---
  chr_col <- pick("chromosome","hm_chrom","chr","CHROM")
  pos_col <- pick("base_pair_location","hm_pos","pos","POS")
  ea_col  <- pick("effect_allele","hm_effect_allele","A1","EA","alt","ALT")
  oa_col  <- pick("other_allele","hm_other_allele","A2","OA","ref","REF")
  
  p_col   <- pick("p_value","p","P","pval")
  eaf_col <- pick("effect_allele_frequency","hm_effect_allele_frequency","eaf","EAF","af")
  
  # N preference (per-SNP effective N if provided)
  n_col   <- pick("cum_eff_sample_size","n","N")
  
  # case/control counts (may be per-SNP or study totals repeated)
  ncase_col <- pick("num_cases","cases","ncase","N_cases")
  nctrl_col <- pick("num_controls","controls","ncontrol","N_controls")
  
  rsid_col <- pick("rsid","hm_rsid")
  vid_col  <- pick("variant_id","hm_variant_id")
  
  snp_col <- if (prefer_snp == "variant_id") {
    pick("variant_id","hm_variant_id","rsid","hm_rsid")
  } else {
    pick("rsid","hm_rsid","variant_id","hm_variant_id")
  }
  
  beta_col <- pick("hm_beta","beta")
  or_col   <- pick("odds_ratio","hm_odds_ratio","OR","or")
  
  se_col    <- pick("standard_error","se","SE")
  ci_lo_col <- pick("ci_lower","hm_ci_lower","lower_ci","ci_l","lci")
  ci_hi_col <- pick("ci_upper","hm_ci_upper","upper_ci","ci_u","uci")
  
  # Required basics
  req <- c(chr_col, pos_col, ea_col, oa_col, p_col)
  if (any(is.na(req))) {
    stop("Missing required cols in ", gwas_path, "\nHave: ", paste(cn, collapse = ", "))
  }
  if (is.na(beta_col) && is.na(or_col)) {
    stop("No beta or OR column found in ", gwas_path, "\nHave: ", paste(cn, collapse = ", "))
  }
  
  # Helper: replace blank/NA characters with literal "NA"
  fix_char <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- "NA"
    x
  }
  
  out <- dt[, .(
    trait      = trait_name,
    type       = type,
    snp        = if (!is.na(snp_col)) as.character(get(snp_col)) else NA_character_,
    rsid       = if (!is.na(rsid_col)) as.character(get(rsid_col)) else NA_character_,
    variant_id = if (!is.na(vid_col))  as.character(get(vid_col))  else NA_character_,
    chr        = suppressWarnings(as.integer(get(chr_col))),
    pos        = suppressWarnings(as.integer(get(pos_col))),
    ea         = as.character(get(ea_col)),
    oa         = as.character(get(oa_col)),
    p          = suppressWarnings(as.numeric(get(p_col))),
    eaf        = if (!is.na(eaf_col)) suppressWarnings(as.numeric(get(eaf_col))) else NA_real_,
    N_raw      = if (!is.na(n_col))   suppressWarnings(as.numeric(get(n_col)))   else NA_real_,
    num_cases  = if (!is.na(ncase_col)) suppressWarnings(as.numeric(get(ncase_col))) else NA_real_,
    num_controls = if (!is.na(nctrl_col)) suppressWarnings(as.numeric(get(nctrl_col))) else NA_real_
  )]
  
  # ---- beta + effect_scale ----
  out[, beta := NA_real_]
  out[, effect_scale := NA_character_]
  if (!is.na(beta_col)) {
    out[, beta := suppressWarnings(as.numeric(dt[[beta_col]]))]
    out[, effect_scale := "beta"]
  } else {
    OR <- suppressWarnings(as.numeric(dt[[or_col]]))
    out[, beta := suppressWarnings(log(OR))]
    out[, effect_scale := "logOR"]
  }
  
  # ---- se + se_source priority: SE -> CI -> p+beta ----
  out[, se := NA_real_]
  out[, se_source := NA_character_]
  
  # (1) explicit SE
  if (!is.na(se_col)) {
    se0 <- suppressWarnings(as.numeric(dt[[se_col]]))
    ok <- is.finite(se0) & se0 > 0
    out[ok, se := se0[ok]]
    out[ok, se_source := "SE"]
  }
  
  # (2) CI (assume 95% on OR/log scale)
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
  
  # (3) infer from p + beta (two-sided)
  p_safe <- out$p
  p_safe <- ifelse(is.finite(p_safe) & p_safe > 0, p_safe, NA_real_)
  # avoid p=0 underflow
  p_safe <- ifelse(!is.na(p_safe), pmax(p_safe, .Machine$double.xmin), NA_real_)
  z <- suppressWarnings(qnorm(p_safe / 2, lower.tail = FALSE))
  se_pb <- suppressWarnings(abs(out$beta) / z)
  
  idx <- which(is.na(out$se) & is.finite(se_pb) & se_pb > 0)
  if (length(idx) > 0) {
    out[idx, se := se_pb[idx]]
    out[idx, se_source := "p+beta"]
  }
  
  out[, varbeta := se^2]
  
  # ---- eaf -> MAF ----
  out[, MAF := ifelse(is.na(eaf), NA_real_, pmin(eaf, 1 - eaf))]
  
  # ---- N + N_source ----
  out[, N := NA_real_]
  out[, N_source := NA_character_]
  
  # prefer per-SNP N column if usable
  okN <- is.finite(out$N_raw) & out$N_raw > 0
  if (!is.na(n_col)) {
    out[okN, N := N_raw]
    out[okN, N_source := n_col]  # column name
  }
  
  # if missing/bad N, derive from cases+controls if present
  ok_ccN <- is.finite(out$num_cases) & is.finite(out$num_controls) & out$num_cases >= 0 & out$num_controls >= 0
  idx <- which((!is.finite(out$N) | out$N <= 0) & ok_ccN)
  if (length(idx) > 0) {
    out[idx, N := num_cases + num_controls]
    out[idx, N_source := "derived_cases_controls"]
  }
  
  # if still missing/bad, use N_override
  idx <- which(!is.finite(out$N) | out$N <= 0)
  if (length(idx) > 0 && is.finite(N_override) && N_override > 0) {
    out[idx, N := N_override]
    out[idx, N_source := "N_override"]
  }
  
  # ---- s + s_source (cc only) with per-SNP vs study-total detection ----
  out[, s := NA_real_]
  out[, s_source := "NA"]
  
  if (type == "cc") {
    # detect whether num_cases/num_controls vary across SNPs
    ncase_vals <- unique(out[is.finite(num_cases), num_cases])
    nctrl_vals <- unique(out[is.finite(num_controls), num_controls])
    
    has_file_counts <- length(ncase_vals) >= 1 && length(nctrl_vals) >= 1
    
    if (has_file_counts && (length(ncase_vals) > 1 || length(nctrl_vals) > 1)) {
      # per-SNP
      idx <- which(is.finite(out$num_cases) & is.finite(out$num_controls) & (out$num_cases + out$num_controls) > 0)
      out[idx, s := num_cases / (num_cases + num_controls)]
      out[idx, s_source := "per_snp"]
    } else if (has_file_counts && length(ncase_vals) == 1 && length(nctrl_vals) == 1) {
      # study totals repeated per row
      denom <- ncase_vals[1] + nctrl_vals[1]
      if (is.finite(denom) && denom > 0) {
        out[, s := ncase_vals[1] / denom]
        out[, s_source := "study_total"]
      }
    } else if (is.finite(cases_total) && is.finite(controls_total) && (cases_total + controls_total) > 0) {
      out[, s := cases_total / (cases_total + controls_total)]
      out[, s_source := "study_total"]
    } else {
      out[, s := NA_real_]
      out[, s_source := "NA"]
    }
  }
  
  # ---- sdY + sdY_source (quant only) ----
  out[, sdY := NA_real_]
  out[, sdY_source := "NA"]
  if (type == "quant") {
    if (is.finite(sdY_override)) {
      out[, sdY := sdY_override]
      out[, sdY_source := "override"]
    }
  }
  
  # ---- Cleaning rules ----
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
  
  # allow MAF NA; if present require (0,0.5)
  out <- out[is.na(MAF) | (is.finite(MAF) & MAF > 0 & MAF < 0.5)]
  
  # cc: require s in (0,1) if provided; allow NA only if no study totals were given
  if (type == "cc") {
    if (is.finite(cases_total) && is.finite(controls_total)) {
      out <- out[is.finite(s) & s > 0 & s < 1]
    } else {
      # if s_source indicates we had totals or per_snp, enforce; else allow NA
      out <- out[
        (is.na(s) & s_source == "NA") |
          (is.finite(s) & s > 0 & s < 1)
      ]
    }
  } else {
    out[, s := NA_real_]
    out[, s_source := "NA"]
  }
  
  # ---- Replace blanks/NA in character fields with literal "NA" ----
  char_cols <- c("trait","type","snp","rsid","variant_id","ea","oa","effect_scale",
                 "se_source","N_source","s_source","sdY_source")
  for (cc in char_cols) out[, (cc) := fix_char(get(cc))]
  
  # ---- Final standardized column order (DO NOT CHANGE) ----
  out_final <- out[, .(
    trait, type, snp, rsid, variant_id, chr, pos, ea, oa,
    beta, se, varbeta, p, effect_scale, se_source,
    eaf, MAF, N, N_source, num_cases, num_controls, s, s_source,
    sdY, sdY_source
  )]
  
  out_final[]
}

# -----------------------------
# Run + write
# -----------------------------
out <- standardize_gwas_for_coloc(
  gwas_path = infile,
  trait_name = trait_name,
  type = type,
  prefer_snp = prefer_snp,
  sdY_override = sdY_override,
  N_override = N_override,
  cases_total = cases_total,
  controls_total = controls_total
)

fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("[%s] wrote %d rows -> %s\n", trait_name, nrow(out), outfile))
