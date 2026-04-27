library(data.table)

trait <- "atopic_dermatitis"
infile  <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/atopic_dermatitis/GCST90503109.h.tsv.gz"
outfile <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/atopic_dermatitis/standardized_GCST90503109.h.tsv.gz"

# trait metadata
type <- "cc"
prefer_snp <- "rsid"
N_override <- 451435
cases_total <- 42963
controls_total <- 408472
sdY_override <- NA_real_  # cc -> ignored

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
    x <- c(...)
    x[x %in% cn][1]
  }
  
  # required-ish
  chr_col <- pick("chromosome","hm_chrom","chr","CHROM")
  pos_col <- pick("base_pair_location","hm_pos","pos","POS")
  ea_col  <- pick("effect_allele","hm_effect_allele","A1","EA","alt","ALT")
  oa_col  <- pick("other_allele","hm_other_allele","A2","OA","ref","REF")
  p_col   <- pick("p_value","p","P","pval","pval_meta")
  
  # identifiers
  rsid_col <- pick("rsid","hm_rsid","RSID")
  vid_col  <- pick("variant_id","hm_variant_id","VID","ID")
  
  snp_col <- if (prefer_snp == "variant_id") {
    pick("variant_id","hm_variant_id","rsid","hm_rsid")
  } else {
    pick("rsid","hm_rsid","variant_id","hm_variant_id")
  }
  
  # effect size sources
  beta_col <- pick("hm_beta","beta","BETA","effect")
  or_col   <- pick("odds_ratio","hm_odds_ratio","OR")
  
  # SE sources
  se_col    <- pick("standard_error","se","SE","stderr","standardError")
  ci_lo_col <- pick("ci_lower","hm_ci_lower","lower_ci","ci_l","lci","LCI","lower")
  ci_hi_col <- pick("ci_upper","hm_ci_upper","upper_ci","ci_u","uci","UCI","upper")
  
  # frequency + N
  eaf_col <- pick("effect_allele_frequency","hm_effect_allele_frequency","eaf","EAF","af","freq1","A1FREQ")
  n_col   <- pick("cum_eff_sample_size","n","N","Neff","effective_n","N_effective")
  
  # optional case/control columns (not present in your header, but keep robust)
  ncase_col <- pick("num_cases","cases","ncase","N_cases","N_CASES")
  nctrl_col <- pick("num_controls","controls","ncontrol","N_controls","N_CONTROLS")
  
  # sanity checks
  req <- c(chr_col, pos_col, ea_col, oa_col, p_col)
  if (any(is.na(req))) {
    stop("Missing required cols in ", gwas_path, "\nHave: ", paste(cn, collapse = ", "))
  }
  if (is.na(beta_col) && is.na(or_col)) {
    stop("No beta or OR column found in ", gwas_path)
  }
  
  # ---- build beta + effect_scale
  beta_vec <- if (!is.na(beta_col)) as.numeric(dt[[beta_col]]) else log(as.numeric(dt[[or_col]]))
  effect_scale <- if (!is.na(beta_col)) "beta" else "logOR"
  
  # ---- build p
  p_vec <- as.numeric(dt[[p_col]])
  
  # ---- build se + se_source (priority: SE > CI > p+beta)
  se_source <- NA_character_
  se_vec <- rep(NA_real_, nrow(dt))
  
  if (!is.na(se_col)) {
    se_vec <- as.numeric(dt[[se_col]])
    se_source <- "SE"
  } else if (!is.na(ci_lo_col) && !is.na(ci_hi_col)) {
    lo <- as.numeric(dt[[ci_lo_col]])
    hi <- as.numeric(dt[[ci_hi_col]])
    # assume 95% CI on the *same scale as beta_vec* (your rules: infer on log scale)
    se_vec <- (hi - lo) / (2 * 1.96)
    se_source <- "CI"
  } else {
    # infer from p + beta: se = |beta| / z, z = qnorm(1 - p/2)
    z <- suppressWarnings(qnorm(1 - p_vec / 2))
    se_vec <- abs(beta_vec) / z
    se_source <- "p+beta"
  }
  
  # ---- eaf / MAF
  eaf_vec <- if (!is.na(eaf_col)) as.numeric(dt[[eaf_col]]) else rep(NA_real_, nrow(dt))
  maf_vec <- if (!is.na(eaf_col)) pmin(eaf_vec, 1 - eaf_vec) else rep(NA_real_, nrow(dt))
  
  # ---- N + N_source rules
  N_vec <- if (!is.na(n_col)) as.numeric(dt[[n_col]]) else rep(NA_real_, nrow(dt))
  N_source <- if (!is.na(n_col)) n_col else NA_character_
  
  # if N missing/bad and per-SNP case/control exist, derive N
  have_cc_cols <- !is.na(ncase_col) && !is.na(nctrl_col)
  ncase_vec <- if (have_cc_cols) as.numeric(dt[[ncase_col]]) else rep(NA_real_, nrow(dt))
  nctrl_vec <- if (have_cc_cols) as.numeric(dt[[nctrl_col]]) else rep(NA_real_, nrow(dt))
  
  badN <- !is.finite(N_vec) | N_vec <= 0
  if (any(badN) && have_cc_cols) {
    derivedN <- ncase_vec + nctrl_vec
    idx <- badN & is.finite(derivedN) & derivedN > 0
    if (any(idx)) {
      N_vec[idx] <- derivedN[idx]
      N_source[idx] <- "derived_cases_controls"
    }
  }
  
  # if still missing/bad, use N_override
  if (any(!is.finite(N_vec) | N_vec <= 0) && is.finite(N_override) && N_override > 0) {
    idx <- !is.finite(N_vec) | N_vec <= 0
    N_vec[idx] <- N_override
    N_source[idx] <- "N_override"
  }
  
  # ---- s + s_source (cc only)
  s_vec <- rep(NA_real_, nrow(dt))
  s_source <- rep("NA", nrow(dt))
  
  if (type == "cc") {
    if (have_cc_cols) {
      # detect whether per-SNP or study totals repeated
      u_cases <- unique(ncase_vec[is.finite(ncase_vec)])
      u_ctrls <- unique(nctrl_vec[is.finite(nctrl_vec)])
      
      cases_vary <- length(u_cases) > 1
      ctrls_vary <- length(u_ctrls) > 1
      
      if (cases_vary || ctrls_vary) {
        denom <- ncase_vec + nctrl_vec
        s_vec <- ncase_vec / denom
        s_source <- "per_snp"
      } else if (length(u_cases) == 1 && length(u_ctrls) == 1) {
        denom <- ncase_vec + nctrl_vec
        s_vec <- ncase_vec / denom
        s_source <- "study_total"
      } else if (is.finite(cases_total) && is.finite(controls_total)) {
        s_vec <- rep(cases_total / (cases_total + controls_total), nrow(dt))
        s_source <- "study_total"
      } else {
        # leave NA
      }
    } else if (is.finite(cases_total) && is.finite(controls_total)) {
      s_vec <- rep(cases_total / (cases_total + controls_total), nrow(dt))
      s_source <- "study_total"
    }
  }
  
  # ---- sdY + sdY_source (quant only; must be last)
  sdY_vec <- rep(NA_real_, nrow(dt))
  sdY_source <- rep("NA", nrow(dt))
  if (type == "quant") {
    if (is.finite(sdY_override)) {
      sdY_vec <- rep(sdY_override, nrow(dt))
      sdY_source <- rep("override", nrow(dt))
    }
  }
  
  # ---- build output with EXACT required columns + order
  out <- data.table(
    trait = trait_name,
    type  = type,
    
    snp = as.character(if (!is.na(snp_col)) dt[[snp_col]] else NA_character_),
    rsid = as.character(if (!is.na(rsid_col)) dt[[rsid_col]] else NA_character_),
    variant_id = as.character(if (!is.na(vid_col)) dt[[vid_col]] else NA_character_),
    
    chr = as.integer(dt[[chr_col]]),
    pos = as.integer(dt[[pos_col]]),
    
    ea = as.character(dt[[ea_col]]),
    oa = as.character(dt[[oa_col]]),
    
    beta = as.numeric(beta_vec),
    se   = as.numeric(se_vec),
    varbeta = as.numeric(se_vec)^2,
    p = as.numeric(p_vec),
    
    effect_scale = effect_scale,
    se_source = se_source,
    
    eaf = as.numeric(eaf_vec),
    MAF = as.numeric(maf_vec),
    
    N = as.numeric(N_vec),
    N_source = as.character(N_source),
    
    num_cases = as.numeric(if (have_cc_cols) ncase_vec else NA_real_),
    num_controls = as.numeric(if (have_cc_cols) nctrl_vec else NA_real_),
    
    s = as.numeric(s_vec),
    s_source = as.character(s_source),
    
    sdY = as.numeric(sdY_vec),
    sdY_source = as.character(sdY_source)
  )
  
  # turn blank strings into NA (then fwrite(na="NA") prints "NA")
  blank_to_na <- function(x) {
    if (is.character(x)) {
      x[nchar(x) == 0] <- NA_character_
    }
    x
  }
  for (j in names(out)) set(out, j = j, value = blank_to_na(out[[j]]))
  
  # ---- cleaning rules
  out <- out[
    is.finite(beta) &
      is.finite(se) & se > 0 &
      is.finite(varbeta) & varbeta > 0 &
      is.finite(p) & p > 0 & p <= 1 &
      is.finite(chr) &
      is.finite(pos) &
      !is.na(ea) & !is.na(oa) &
      is.finite(N) & N > 0
  ]
  
  # do NOT require eaf/MAF; but if MAF is present (finite), require (0, 0.5)
  out <- out[is.na(MAF) | (is.finite(MAF) & MAF > 0 & MAF < 0.5)]
  
  # cc: require s in (0,1) if provided; allow NA only if no study totals were given
  if (type == "cc") {
    have_study_totals <- is.finite(cases_total) && is.finite(controls_total)
    if (have_study_totals) {
      out <- out[is.finite(s) & s > 0 & s < 1]
    } else {
      out <- out[is.na(s) | (is.finite(s) & s > 0 & s < 1)]
    }
  }
  
  # enforce exact column order (and only these columns)
  std_cols <- c(
    "trait","type","snp","rsid","variant_id","chr","pos","ea","oa",
    "beta","se","varbeta","p","effect_scale","se_source","eaf","MAF",
    "N","N_source","num_cases","num_controls","s","s_source","sdY","sdY_source"
  )
  out <- out[, ..std_cols]
  
  out[]
}

out <- standardize_gwas_for_coloc(
  gwas_path   = infile,
  trait_name  = trait,
  type        = type,
  prefer_snp  = prefer_snp,
  sdY_override = sdY_override,
  N_override  = N_override,
  cases_total = cases_total,
  controls_total = controls_total
)

fwrite(out, outfile, sep = "\t", na = "NA", quote=F)
cat(sprintf("[%s] wrote %d rows -> %s\n", trait, nrow(out), outfile))
