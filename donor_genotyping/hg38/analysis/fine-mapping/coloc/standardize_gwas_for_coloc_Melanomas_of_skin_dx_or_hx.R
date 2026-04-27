library(data.table)

# ---------------- User inputs ----------------
trait_name   <- "Melanomas_of_skin_dx_or_hx"
type         <- "cc"          # "cc" or "quant"
prefer_snp   <- "rsid"        # "rsid" or "variant_id"

infile       <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/Melanomas_of_skin_dx_or_hx/GCST90475577.h.tsv.gz"
outfile      <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/Melanomas_of_skin_dx_or_hx/standardized_GCST90475577.h.tsv.gz"

N_override       <- NA_real_  # optional constant; NA = don't use
cases_total      <- NA_real_  # cc only (optional)
controls_total   <- NA_real_  # cc only (optional)
sdY_override     <- NA_real_  # quant only (optional)

# ---------------- Standardizer ----------------
standardize_gwas_for_coloc <- function(gwas_path,
                                       trait_name,
                                       type = c("cc","quant"),
                                       prefer_snp = c("rsid","variant_id"),
                                       N_override = NA_real_,
                                       cases_total = NA_real_,
                                       controls_total = NA_real_,
                                       sdY_override = NA_real_) {
  
  type <- match.arg(type)
  prefer_snp <- match.arg(prefer_snp)
  
  dt <- fread(gwas_path, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
  cn <- names(dt)
  
  # Treat blank/whitespace-only as NA (character cols)
  for (j in which(vapply(dt, is.character, logical(1)))) {
    x <- dt[[j]]
    x <- trimws(x)
    x[x == ""] <- NA_character_
    dt[[j]] <- x
  }
  
  # robust column picker
  pick <- function(...) {
    hits <- c(...)
    hits[hits %in% cn][1]
  }
  
  # --- map columns for this header (and common alternates) ---
  chr_col <- pick("chromosome","chr","CHROM","hm_chrom")
  pos_col <- pick("base_pair_location","pos","POS","bp","hm_pos")
  
  ea_col  <- pick("effect_allele","EA","A1","alt","ALT","hm_effect_allele")
  oa_col  <- pick("other_allele","OA","A2","ref","REF","hm_other_allele")
  
  p_col   <- pick("p_value","p","P","pval","p_value_hm","hm_p_value")
  
  rsid_col <- pick("rsid","hm_rsid")
  vid_col  <- pick("variant_id","hm_variant_id")
  
  # effect size
  beta_col <- pick("beta","BETA","hm_beta")
  or_col   <- pick("odds_ratio","OR","or","hm_odds_ratio")
  
  # uncertainty (priority: SE -> CI -> p+beta)
  se_col    <- pick("standard_error","se","SE","stderr","hm_standard_error","hm_se")
  ci_lo_col <- pick("ci_lower","lower_ci","lci","ci_l","hm_ci_lower")
  ci_hi_col <- pick("ci_upper","upper_ci","uci","ci_u","hm_ci_upper")
  
  # allele freq / sample size
  eaf_col <- pick("effect_allele_frequency","eaf","EAF","af","hm_effect_allele_frequency")
  n_col   <- pick("cum_eff_sample_size","n","N","Neff","effective_n","effN","N_eff")
  
  # case/control counts
  ncase_col <- pick("num_cases","ncase","cases","N_cases")
  nctrl_col <- pick("num_controls","ncontrol","controls","N_controls")
  
  # required minimal fields
  req <- c(chr_col, pos_col, ea_col, oa_col, p_col)
  if (any(is.na(req))) {
    stop("Missing required cols in ", gwas_path, "\nHave: ", paste(cn, collapse = ", "))
  }
  if (is.na(beta_col) && is.na(or_col)) {
    stop("No beta or odds_ratio column found in ", gwas_path)
  }
  if (prefer_snp == "rsid" && is.na(rsid_col) && is.na(vid_col)) {
    stop("prefer_snp='rsid' but rsid/variant_id not found in header.")
  }
  if (prefer_snp == "variant_id" && is.na(vid_col) && is.na(rsid_col)) {
    stop("prefer_snp='variant_id' but variant_id/rsid not found in header.")
  }
  
  # --- base output skeleton (DO NOT change order later) ---
  out <- dt[, .(
    trait      = trait_name,
    type       = type,
    
    snp        = NA_character_,  # filled below
    rsid       = if (!is.na(rsid_col)) as.character(get(rsid_col)) else NA_character_,
    variant_id = if (!is.na(vid_col))  as.character(get(vid_col))  else NA_character_,
    
    chr        = suppressWarnings(as.integer(get(chr_col))),
    pos        = suppressWarnings(as.integer(get(pos_col))),
    
    ea         = toupper(as.character(get(ea_col))),
    oa         = toupper(as.character(get(oa_col))),
    
    beta       = NA_real_,
    se         = NA_real_,
    varbeta    = NA_real_,
    p          = suppressWarnings(as.numeric(get(p_col))),
    
    effect_scale = NA_character_,
    se_source    = NA_character_,
    
    eaf       = if (!is.na(eaf_col)) suppressWarnings(as.numeric(get(eaf_col))) else NA_real_,
    MAF       = NA_real_,
    
    N         = if (!is.na(n_col)) suppressWarnings(as.numeric(get(n_col))) else NA_real_,
    N_source  = NA_character_,
    
    num_cases    = if (!is.na(ncase_col)) suppressWarnings(as.numeric(get(ncase_col))) else NA_real_,
    num_controls = if (!is.na(nctrl_col)) suppressWarnings(as.numeric(get(nctrl_col))) else NA_real_,
    
    s         = NA_real_,
    s_source  = NA_character_,
    
    sdY       = NA_real_,
    sdY_source= NA_character_
  )]
  
  # snp column (merge identifier)
  if (prefer_snp == "variant_id") {
    out[, snp := fifelse(!is.na(variant_id), variant_id, rsid)]
  } else {
    out[, snp := fifelse(!is.na(rsid), rsid, variant_id)]
  }
  
  # ---- beta + effect_scale ----
  if (!is.na(beta_col)) {
    out[, beta := suppressWarnings(as.numeric(dt[[beta_col]]))]
    out[, effect_scale := "beta"]
  } else {
    OR <- suppressWarnings(as.numeric(dt[[or_col]]))
    out[, beta := log(OR)]
    out[, effect_scale := "logOR"]
  }
  
  # ---- se with priority: SE -> CI(95%) -> p+beta ----
  out[, se := NA_real_]
  out[, se_source := NA_character_]
  
  # 1) explicit SE
  if (!is.na(se_col)) {
    se0 <- suppressWarnings(as.numeric(dt[[se_col]]))
    ok  <- is.finite(se0) & se0 > 0
    out[ok, se := se0[ok]]
    out[ok, se_source := "SE"]
  }
  
  # 2) CI bounds (95%)
  if (!is.na(ci_lo_col) && !is.na(ci_hi_col)) {
    lo <- suppressWarnings(as.numeric(dt[[ci_lo_col]]))
    hi <- suppressWarnings(as.numeric(dt[[ci_hi_col]]))
    
    if (out$effect_scale[1] == "logOR") {
      # CI provided on OR scale in this dataset -> use log(hi/lo)
      ok_ci <- is.finite(lo) & is.finite(hi) & lo > 0 & hi > 0
      se_ci <- (log(hi) - log(lo)) / (2 * 1.96)
    } else {
      # fallback: assume CI already on beta scale
      ok_ci <- is.finite(lo) & is.finite(hi)
      se_ci <- (hi - lo) / (2 * 1.96)
    }
    
    fill <- is.na(out$se) & ok_ci & is.finite(se_ci) & se_ci > 0
    out[fill, se := se_ci[fill]]
    out[fill, se_source := "CI"]
  }
  
  # 3) infer from p + beta (two-sided)
  fill <- is.na(out$se) & is.finite(out$beta) & is.finite(out$p) & out$p > 0 & out$p <= 1
  if (any(fill)) {
    z <- qnorm(out$p[fill] / 2, lower.tail = FALSE)
    se_pb <- abs(out$beta[fill]) / z
    ok_pb <- is.finite(se_pb) & se_pb > 0
    idx <- which(fill)[ok_pb]
    if (length(idx) > 0) {
      out[idx, se := se_pb[ok_pb]]
      out[idx, se_source := "p+beta"]
    }
  }
  
  out[, varbeta := se^2]
  
  # ---- eaf / MAF ----
  out[, MAF := ifelse(is.na(eaf), NA_real_, pmin(eaf, 1 - eaf))]
  
  # ---- N + N_source ----
  out[, N_source := NA_character_]
  if (!is.na(n_col)) {
    okN <- is.finite(out$N) & out$N > 0
    out[okN, N_source := n_col]
  }
  
  # derive from per-SNP cases/controls if N missing/bad
  badN <- !(is.finite(out$N) & out$N > 0)
  if (any(badN) && (!is.na(ncase_col) || !is.na(nctrl_col))) {
    ok_ccN <- badN &
      is.finite(out$num_cases) & is.finite(out$num_controls) &
      out$num_cases >= 0 & out$num_controls >= 0 &
      (out$num_cases + out$num_controls) > 0
    out[ok_ccN, N := num_cases + num_controls]
    out[ok_ccN, N_source := "derived_cases_controls"]
  }
  
  # fallback to N_override
  badN <- !(is.finite(out$N) & out$N > 0)
  if (any(badN) && is.finite(N_override) && N_override > 0) {
    out[badN, N := N_override]
    out[badN, N_source := "N_override"]
  }
  
  # ---- s + s_source ----
  out[, s := NA_real_]
  out[, s_source := "NA"]
  
  if (type == "cc") {
    ok_counts <- is.finite(out$num_cases) & is.finite(out$num_controls) &
      out$num_cases >= 0 & out$num_controls >= 0 &
      (out$num_cases + out$num_controls) > 0
    
    # detect per-SNP vs repeated totals by uniqueness over non-missing
    uniq_cases <- unique(out$num_cases[is.finite(out$num_cases)])
    uniq_ctrls <- unique(out$num_controls[is.finite(out$num_controls)])
    counts_look_study_total <- (length(uniq_cases) == 1L) && (length(uniq_ctrls) == 1L)
    
    if (any(ok_counts)) {
      if (counts_look_study_total) {
        cases_total_dt    <- uniq_cases[1]
        controls_total_dt <- uniq_ctrls[1]
        if (is.finite(cases_total_dt) && is.finite(controls_total_dt) &&
            (cases_total_dt + controls_total_dt) > 0) {
          out[, s := cases_total_dt / (cases_total_dt + controls_total_dt)]
          out[, s_source := "study_total"]
        }
      } else {
        out[ok_counts, s := num_cases / (num_cases + num_controls)]
        out[ok_counts, s_source := "per_snp"]
      }
    }
    
    # If still NA, fall back to user-provided totals (if given)
    if (all(is.na(out$s))) {
      if (is.finite(cases_total) && is.finite(controls_total) && (cases_total + controls_total) > 0) {
        out[, s := cases_total / (cases_total + controls_total)]
        out[, s_source := "study_total"]
      } else {
        out[, s := NA_real_]
        out[, s_source := "NA"]
      }
    }
  } else {
    out[, s := NA_real_]
    out[, s_source := "NA"]
  }
  
  # ---- sdY + sdY_source ----
  if (type == "quant") {
    if (is.finite(sdY_override) && sdY_override > 0) {
      out[, sdY := sdY_override]
      out[, sdY_source := "override"]
    } else {
      out[, sdY := NA_real_]
      out[, sdY_source := "NA"]
    }
  } else {
    out[, sdY := NA_real_]
    out[, sdY_source := "NA"]
  }
  
  # ---- Cleaning ----
  out <- out[
    is.finite(beta) &
      is.finite(se) & se > 0 &
      is.finite(varbeta) & varbeta > 0 &
      is.finite(p) & p > 0 & p <= 1 &
      is.finite(chr) & is.finite(pos) &
      !is.na(ea) &
      !is.na(oa) &
      is.finite(N) & N > 0
  ]
  
  # allow MAF NA; if present require 0<MAF<0.5
  out <- out[is.na(MAF) | (MAF > 0 & MAF < 0.5)]
  
  # cc-specific s rule
  if (type == "cc") {
    out <- out[is.na(s) | (s > 0 & s < 1)]
    totals_given <- is.finite(cases_total) && is.finite(controls_total) && (cases_total + controls_total) > 0
    if (totals_given) out <- out[!is.na(s)]
  }
  
  # enforce EXACT standardized column order (sdY + sdY_source at end)
  out <- out[, .(
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
  
  out[]
}

# ---------------- Run + write ----------------
out <- standardize_gwas_for_coloc(
  gwas_path      = infile,
  trait_name     = trait_name,
  type           = type,
  prefer_snp     = prefer_snp,
  N_override     = N_override,
  cases_total    = cases_total,
  controls_total = controls_total,
  sdY_override   = sdY_override
)

fwrite(out, outfile, sep = "\t", compress = "gzip", na = "NA", quote=F)
cat(sprintf("[%s] wrote %d rows -> %s\n", trait_name, nrow(out), outfile))
