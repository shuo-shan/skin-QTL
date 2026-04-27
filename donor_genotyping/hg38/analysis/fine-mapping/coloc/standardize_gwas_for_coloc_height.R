library(data.table)

# ----------------------------
# User inputs / metadata
# ----------------------------
trait_name  <- "height"
type        <- "quant"          # "cc" or "quant"
prefer_snp  <- "rsid"           # "rsid" or "variant_id"

infile  <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/height/GCST90475362.h.tsv.gz"
outfile <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/height/standardized_GCST90475362.h.tsv.gz"

# quant-only
sdY_override <- 1               # inv-normal / z-scored trait => ~ N(0,1)

# used only if per-SNP N missing/bad
N_override <- 424305            # optional (file already has n)

# cc-only (not used here)
cases_total    <- NA_real_
controls_total <- NA_real_


# ----------------------------
# Standardizer (general; works for this header)
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
    hits <- c(...)
    hits[hits %in% cn][1]
  }
  
  # Core locus / alleles
  chr_col <- pick("chromosome","chr","CHROM","hm_chrom")
  pos_col <- pick("base_pair_location","pos","POS","hm_pos")
  
  ea_col  <- pick("effect_allele","EA","A1","alt","ALT","hm_effect_allele")
  oa_col  <- pick("other_allele","OA","A2","ref","REF","hm_other_allele")
  
  # IDs
  rsid_col <- pick("rsid","RSID","hm_rsid")
  vid_col  <- pick("variant_id","variantID","hm_variant_id","hm_variantid","hm_codvariant_id") # tolerate weird header truncation
  
  snp_col <- if (prefer_snp == "variant_id") vid_col else rsid_col
  
  # Effects: beta OR OR/logOR
  beta_col <- pick("beta","hm_beta","BETA")
  or_col   <- pick("OR","or","odds_ratio","oddsratio")
  
  # SE priority: SE col else CI else p+beta
  se_col   <- pick("standard_error","se","SE","stderr","SE_beta","hm_standard_error","hm_se")
  ci_lo_col <- pick("ci_lower","CI_LOWER","lower_ci","lci","ci_l","lower")
  ci_hi_col <- pick("ci_upper","CI_UPPER","upper_ci","uci","ci_u","upper")
  
  # p-value
  p_col <- pick("p_value","p","P","pval","p_value_neg_log10","pvalue","hm_p_value","hm_p")
  
  # EAF
  eaf_col <- pick("effect_allele_frequency","eaf","EAF","af","ALT_AF","hm_effect_allele_frequency","hm_eaf")
  
  # Sample size (per-SNP effective N preferred)
  n_col <- pick("cum_eff_sample_size","effective_sample_size","effN","n_eff","n","N","samplesize")
  
  # Required minimal set for writing a usable coloc table
  req_min <- c(chr_col, pos_col, ea_col, oa_col, p_col, snp_col)
  if (any(is.na(req_min))) {
    stop("Missing required cols in ", gwas_path,
         "\nMissing: ", paste(req_min[is.na(req_min)], collapse = ", "),
         "\nHave: ", paste(cn, collapse = ", "))
  }
  
  # ---- build beta + effect_scale
  beta_vec <- rep(NA_real_, nrow(dt))
  effect_scale <- rep(NA_character_, nrow(dt))
  
  if (!is.na(beta_col)) {
    beta_vec <- suppressWarnings(as.numeric(dt[[beta_col]]))
    effect_scale <- "beta"
  } else if (!is.na(or_col)) {
    orv <- suppressWarnings(as.numeric(dt[[or_col]]))
    beta_vec <- suppressWarnings(log(orv))
    effect_scale <- "logOR"
  } else {
    stop("No beta or OR column found in: ", gwas_path)
  }
  
  # ---- build p
  p_vec <- suppressWarnings(as.numeric(dt[[p_col]]))
  # handle -log10(p) if provided (rare; safeguard)
  if (grepl("neg_log10", p_col, ignore.case = TRUE)) {
    p_vec <- suppressWarnings(10^(-as.numeric(dt[[p_col]])))
  }
  
  # ---- build se + se_source
  se_vec <- rep(NA_real_, nrow(dt))
  se_source <- rep(NA_character_, nrow(dt))
  
  if (!is.na(se_col)) {
    se_vec <- suppressWarnings(as.numeric(dt[[se_col]]))
    se_source <- "SE"
  } else if (!is.na(ci_lo_col) && !is.na(ci_hi_col)) {
    lo <- suppressWarnings(as.numeric(dt[[ci_lo_col]]))
    hi <- suppressWarnings(as.numeric(dt[[ci_hi_col]]))
    # assume 95% CI on the effect scale (beta or logOR)
    se_vec <- (hi - lo) / (2 * 1.96)
    se_source <- "CI"
  } else {
    # infer from p + beta
    z <- suppressWarnings(abs(qnorm(p_vec / 2, lower.tail = FALSE)))
    se_vec <- suppressWarnings(abs(beta_vec) / z)
    se_source <- "p+beta"
  }
  
  # ---- eaf / MAF
  eaf_vec <- if (!is.na(eaf_col)) suppressWarnings(as.numeric(dt[[eaf_col]])) else rep(NA_real_, nrow(dt))
  maf_vec <- if (!is.na(eaf_col)) pmin(eaf_vec, 1 - eaf_vec) else rep(NA_real_, nrow(dt))
  
  # ---- N + N_source
  N_vec <- rep(NA_real_, nrow(dt))
  N_source <- rep(NA_character_, nrow(dt))
  
  if (!is.na(n_col)) {
    N_vec <- suppressWarnings(as.numeric(dt[[n_col]]))
    N_source <- n_col
  }
  
  # cc: num_cases/num_controls and s detection (for quant: NA)
  num_cases <- rep(NA_real_, nrow(dt))
  num_controls <- rep(NA_real_, nrow(dt))
  s_vec <- rep(NA_real_, nrow(dt))
  s_source <- rep("NA", nrow(dt))
  
  if (type == "cc") {
    cases_col <- pick("num_cases","cases","N_cases","ncases","case_n","case_count")
    ctrls_col <- pick("num_controls","controls","N_controls","ncontrols","control_n","control_count")
    if (!is.na(cases_col) && !is.na(ctrls_col)) {
      num_cases <- suppressWarnings(as.numeric(dt[[cases_col]]))
      num_controls <- suppressWarnings(as.numeric(dt[[ctrls_col]]))
      
      u_cases <- unique(num_cases[is.finite(num_cases)])
      u_ctrls <- unique(num_controls[is.finite(num_controls)])
      
      if (length(u_cases) > 1 || length(u_ctrls) > 1) {
        s_vec <- num_cases / (num_cases + num_controls)
        s_source <- "per_snp"
      } else if (length(u_cases) == 1 && length(u_ctrls) == 1) {
        s_vec <- num_cases / (num_cases + num_controls)
        s_source <- "study_total"
      }
    } else if (is.finite(cases_total) && is.finite(controls_total)) {
      s_vec <- cases_total / (cases_total + controls_total)
      s_source <- "study_total"
    }
  }
  
  # If N missing/bad: try derive from cases/controls (cc only), else N_override
  badN <- !(is.finite(N_vec) & N_vec > 0)
  if (any(badN) && type == "cc") {
    can_derive <- is.finite(num_cases) & is.finite(num_controls) & (num_cases + num_controls) > 0
    idx <- which(badN & can_derive)
    if (length(idx) > 0) {
      N_vec[idx] <- num_cases[idx] + num_controls[idx]
      N_source[idx] <- "derived_cases_controls"
    }
    badN <- !(is.finite(N_vec) & N_vec > 0)
  }
  if (any(badN) && is.finite(N_override) && N_override > 0) {
    N_vec[badN] <- N_override
    N_source[badN] <- "N_override"
  }
  
  # ---- sdY (quant only)
  sdY_vec <- rep(NA_real_, nrow(dt))
  sdY_source <- rep("NA", nrow(dt))
  if (type == "quant") {
    if (is.finite(sdY_override)) {
      sdY_vec <- sdY_override
      sdY_source <- "override"
    }
  }
  
  # ---- assemble standardized output (EXACT column order)
  out <- data.table(
    trait = trait_name,
    type  = type,
    snp   = as.character(dt[[snp_col]]),
    rsid  = if (!is.na(rsid_col)) as.character(dt[[rsid_col]]) else NA_character_,
    variant_id = if (!is.na(vid_col)) as.character(dt[[vid_col]]) else NA_character_,
    chr   = suppressWarnings(as.integer(dt[[chr_col]])),
    pos   = suppressWarnings(as.integer(dt[[pos_col]])),
    ea    = as.character(dt[[ea_col]]),
    oa    = as.character(dt[[oa_col]]),
    beta  = beta_vec,
    se    = se_vec,
    varbeta = se_vec^2,
    p     = p_vec,
    effect_scale = effect_scale,
    se_source = se_source,
    eaf   = eaf_vec,
    MAF   = maf_vec,
    N     = N_vec,
    N_source = N_source,
    num_cases = num_cases,
    num_controls = num_controls,
    s     = if (type == "cc") s_vec else NA_real_,
    s_source = if (type == "cc") s_source else "NA",
    sdY   = sdY_vec,
    sdY_source = sdY_source
  )
  
  # Replace blank strings with NA (and later we fwrite(na="NA"))
  char_cols <- names(out)[vapply(out, is.character, logical(1))]
  for (cc in char_cols) {
    set(out, which(out[[cc]] == ""), cc, NA_character_)
  }
  
  # ----------------------------
  # Cleaning rules
  # ----------------------------
  out <- out[is.finite(beta) & is.finite(se) & se > 0 & is.finite(varbeta) & varbeta > 0]
  out <- out[is.finite(p) & p > 0 & p <= 1]
  out <- out[is.finite(chr) & is.finite(pos)]
  out <- out[!is.na(ea) & !is.na(oa) & ea != "" & oa != ""]
  
  out <- out[is.finite(N) & N > 0]
  
  # Do NOT require eaf/MAF; but if MAF is present, require (0<MAF<0.5)
  out <- out[is.na(MAF) | (is.finite(MAF) & MAF > 0 & MAF < 0.5)]
  
  # For cc: require s in (0,1) if provided; allow NA only if no study totals were given
  if (type == "cc") {
    have_totals <- is.finite(cases_total) && is.finite(controls_total)
    out <- out[is.na(s) | (is.finite(s) & s > 0 & s < 1)]
    if (have_totals) {
      # if totals given, we expect s to be defined (study_total)
      out <- out[!is.na(s)]
    }
  }
  
  # Ensure exact column order (safety)
  std_cols <- c(
    "trait","type","snp","rsid","variant_id","chr","pos","ea","oa",
    "beta","se","varbeta","p","effect_scale","se_source","eaf","MAF",
    "N","N_source","num_cases","num_controls","s","s_source","sdY","sdY_source"
  )
  out <- out[, ..std_cols]
  
  # estimate sdY
  estimate_sdY_from_out <- function(out, eaf_min = 0.05, eaf_max = 0.95) {
    vg <- 2 * out$eaf * (1 - out$eaf)
    sdY_sq <- (out$se^2) * out$N * vg
    
    ok <- is.finite(sdY_sq) &
      is.finite(out$se) & out$se > 0 & out$se < 1 &
      is.finite(out$N) & out$N > 0 &
      is.finite(out$eaf) & out$eaf > eaf_min & out$eaf < eaf_max
    
    sdY_vec <- sqrt(sdY_sq[ok])
    list(
      sdY_hat = sqrt(stats::median(sdY_sq[ok], na.rm = TRUE)),
      n_used  = length(sdY_vec),
      iqr     = stats::IQR(sdY_vec, na.rm = TRUE)
    )
  }
  
  if (unique(out$type) == "quant") {
    est <- estimate_sdY_from_out(out)
    message(sprintf("[sdY] %s: sdY=%.3f (n=%d, IQR=%.3f)",
                    unique(out$trait), est$sdY_hat, est$n_used, est$iqr))
  }
  
  out[]
}


# ----------------------------
# Run for this file
# ----------------------------
out <- standardize_gwas_for_coloc(
  gwas_path   = infile,
  trait_name  = trait_name,
  type        = type,
  prefer_snp  = prefer_snp,
  sdY_override = sdY_override,
  N_override  = N_override,
  cases_total = cases_total,
  controls_total = controls_total
)

fwrite(out, outfile, sep = "\t", na = "NA", quote=F)
cat(sprintf("[%s] wrote %d rows -> %s\n", trait_name, nrow(out), outfile))
