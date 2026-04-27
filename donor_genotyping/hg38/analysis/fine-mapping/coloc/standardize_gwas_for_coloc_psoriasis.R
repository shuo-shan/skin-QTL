library(data.table)

trait   <- "psoriasis"
type    <- "cc"  # "cc" or "quant"
infile  <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/psoriasis/GCST90472771.h.tsv.gz"
outfile <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/psoriasis/standardized_GCST90472771.h.tsv.gz"
prefer  <- "rsid"

# Optional: used only when per-SNP N is missing/bad
N_override     <- 494544
cases_total    <- 36466
controls_total <- 458078
sdY_override   <- NA_real_  # quant-only; e.g., 1 for z-scored / invnorm traits

gwas_path <- infile
trait_name <- trait
type <- "cc"
prefer_snp <- prefer

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
  cn_l <- tolower(cn)
  
  pick <- function(...) {
    cands <- c(...)
    cands_l <- tolower(cands)
    hit <- match(cands_l, cn_l)
    if (all(is.na(hit))) return(NA_character_)
    cn[hit[which(!is.na(hit))[1]]]
  }
  
  # ---- map columns (robust picking) ----
  chr_col <- pick("chromosome","chr","chrom","hm_chrom","CHROM")
  pos_col <- pick("base_pair_location","pos","bp","hm_pos","POS")
  ea_col  <- pick("effect_allele","ea","a1","alt","ALT","hm_effect_allele","EA")
  oa_col  <- pick("other_allele","oa","a2","ref","REF","hm_other_allele","OA")
  
  p_col   <- pick("p_value","p","pval","pvalue","P")
  eaf_col <- pick("effect_allele_frequency","eaf","af","EAF","hm_effect_allele_frequency")
  
  # per-SNP N (effective N etc.)
  n_col   <- pick("cum_eff_sample_size","n","N","neff","n_eff","effective_sample_size")
  
  # case/control counts (may be per-SNP OR repeated totals)
  ncase_col <- pick("num_cases","cases","ncase","N_cases","n_cases")
  nctrl_col <- pick("num_controls","controls","ncontrol","N_controls","n_controls")
  
  rsid_col <- pick("rsid","hm_rsid","RSID","rsID")
  vid_col  <- pick("variant_id","hm_variant_id","snpid","variant","VID")
  
  snp_col <- if (prefer_snp == "variant_id") {
    pick("variant_id","hm_variant_id","rsid","hm_rsid")
  } else {
    pick("rsid","hm_rsid","variant_id","hm_variant_id")
  }
  
  beta_col <- pick("beta","hm_beta","BETA")
  or_col   <- pick("odds_ratio","or","hm_odds_ratio","OR")
  
  se_col    <- pick("standard_error","se","SE","stderr","standarderror")
  ci_lo_col <- pick("ci_lower","lower_ci","lci","ci_l","hm_ci_lower")
  ci_hi_col <- pick("ci_upper","upper_ci","uci","ci_u","hm_ci_upper")
  
  # ---- required checks ----
  req <- c(chr_col, pos_col, ea_col, oa_col, p_col)
  if (any(is.na(req))) {
    stop("Missing required cols in ", gwas_path, "\nHave: ", paste(cn, collapse = ", "))
  }
  if (is.na(beta_col) && is.na(or_col)) {
    stop("No beta or OR column found in ", gwas_path, "\nHave: ", paste(cn, collapse = ", "))
  }
  
  # helper: numeric safely
  num <- function(x) suppressWarnings(as.numeric(x))
  
  # ---- build standardized table with ALL columns ALWAYS present (exact schema later) ----
  out <- dt[, .(
    trait      = trait_name,
    type       = type,
    
    snp        = if (!is.na(snp_col))  as.character(get(snp_col))  else NA_character_,
    rsid       = if (!is.na(rsid_col)) as.character(get(rsid_col)) else NA_character_,
    variant_id = if (!is.na(vid_col))  as.character(get(vid_col))  else NA_character_,
    
    chr        = suppressWarnings(as.integer(get(chr_col))),
    pos        = suppressWarnings(as.integer(get(pos_col))),
    
    ea         = as.character(get(ea_col)),
    oa         = as.character(get(oa_col)),
    
    beta       = NA_real_,
    se         = NA_real_,
    varbeta    = NA_real_,
    p          = num(get(p_col)),
    
    effect_scale = NA_character_,
    se_source    = NA_character_,
    
    eaf       = if (!is.na(eaf_col)) num(get(eaf_col)) else NA_real_,
    MAF       = NA_real_,
    
    N         = if (!is.na(n_col)) num(get(n_col)) else NA_real_,
    N_source  = if (!is.na(n_col)) n_col else NA_character_,
    
    num_cases    = if (!is.na(ncase_col)) num(get(ncase_col)) else NA_real_,
    num_controls = if (!is.na(nctrl_col)) num(get(nctrl_col)) else NA_real_,
    
    s        = NA_real_,
    s_source = NA_character_,
    
    sdY        = NA_real_,
    sdY_source = NA_character_
  )]
  
  # ---- beta policy: beta or log(OR) ----
  if (!is.na(beta_col)) {
    out[, beta := num(dt[[beta_col]])]
    out[, effect_scale := "beta"]
  } else {
    OR <- num(dt[[or_col]])
    out[, beta := log(OR)]
    out[, effect_scale := "logOR"]
  }
  
  # ---- se policy: SE > CI > p+beta ----
  # 1) SE column
  if (!is.na(se_col)) {
    se0 <- num(dt[[se_col]])
    out[, se := se0]
    out[is.finite(se) & se > 0, se_source := "SE"]
  }
  
  # 2) CI bounds (assume 95%, infer on log scale)
  if (!is.na(ci_lo_col) && !is.na(ci_hi_col)) {
    lo <- num(dt[[ci_lo_col]])
    hi <- num(dt[[ci_hi_col]])
    ok_ci <- is.finite(lo) & is.finite(hi) & lo > 0 & hi > 0
    se_ci <- (log(hi) - log(lo)) / (2 * 1.96)
    
    idx <- which((!is.finite(out$se) | out$se <= 0) & ok_ci & is.finite(se_ci) & se_ci > 0)
    if (length(idx) > 0) {
      out[idx, `:=`(se = se_ci[idx], se_source = "CI")]
    }
  }
  
  # 3) p + beta (two-sided)
  idx2 <- which(
    (!is.finite(out$se) | out$se <= 0) &
      is.finite(out$beta) &
      is.finite(out$p) & out$p > 0 & out$p < 1
  )
  if (length(idx2) > 0) {
    z <- qnorm(out$p[idx2] / 2, lower.tail = FALSE)
    se_p <- abs(out$beta[idx2]) / z
    ok <- is.finite(se_p) & se_p > 0
    if (any(ok)) {
      out[idx2[ok], `:=`(se = se_p[ok], se_source = "p+beta")]
    }
  }
  
  out[, varbeta := se^2]
  
  # ---- MAF policy ----
  out[, MAF := ifelse(is.na(eaf), NA_real_, pmin(eaf, 1 - eaf))]
  
  # ---- N policy ----
  # If per-SNP N missing/bad and per-SNP cases/controls exist -> derive
  out[( !is.finite(N) | N <= 0 ) &
        is.finite(num_cases) & is.finite(num_controls) &
        (num_cases + num_controls) > 0,
      `:=`(N = num_cases + num_controls, N_source = "derived_cases_controls")]
  
  # Fallback: N_override
  if (is.finite(N_override) && N_override > 0) {
    out[!is.finite(N) | N <= 0, `:=`(N = N_override, N_source = "N_override")]
  }
  
  # ---- s policy (case fraction; detect per-SNP vs repeated totals) ----
  totals_from_cols_available <- FALSE
  
  if (type == "cc") {
    # default
    out[, `:=`(s = NA_real_, s_source = "NA")]
    
    has_counts <- is.finite(out$num_cases) & is.finite(out$num_controls) & (out$num_cases + out$num_controls) > 0
    
    # detection based on uniqueness among non-missing values
    uniq_cases <- unique(out$num_cases[is.finite(out$num_cases)])
    uniq_ctrls <- unique(out$num_controls[is.finite(out$num_controls)])
    
    nuniq_cases <- length(uniq_cases)
    nuniq_ctrls <- length(uniq_ctrls)
    
    if (nuniq_cases > 1 || nuniq_ctrls > 1) {
      # per-SNP varying counts
      out[has_counts, `:=`(s = num_cases / (num_cases + num_controls), s_source = "per_snp")]
      totals_from_cols_available <- TRUE
    } else if (nuniq_cases == 1 && nuniq_ctrls == 1) {
      # study totals repeated on every row
      sc <- uniq_cases[1]
      sn <- uniq_ctrls[1]
      if (is.finite(sc) && is.finite(sn) && (sc + sn) > 0) {
        out[, `:=`(s = sc / (sc + sn), s_source = "study_total")]
        totals_from_cols_available <- TRUE
      }
    } else if (is.finite(cases_total) && is.finite(controls_total) && (cases_total + controls_total) > 0) {
      out[, `:=`(s = cases_total / (cases_total + controls_total), s_source = "study_total")]
      totals_from_cols_available <- TRUE
    }
  } else {
    out[, `:=`(s = NA_real_, s_source = "NA")]
  }
  
  # ---- sdY policy ----
  if (type == "quant") {
    if (is.finite(sdY_override) && sdY_override > 0) {
      out[, `:=`(sdY = sdY_override, sdY_source = "override")]
    } else {
      out[, `:=`(sdY = NA_real_, sdY_source = "NA")]
    }
  } else {
    out[, `:=`(sdY = NA_real_, sdY_source = "NA")]
  }
  
  # ---- normalize blanks for character fields (and we will also write na="NA") ----
  char_cols <- c("trait","type","snp","rsid","variant_id","ea","oa","effect_scale","se_source","N_source","s_source","sdY_source")
  for (cc in char_cols) {
    if (!cc %in% names(out)) next
    out[[cc]] <- as.character(out[[cc]])
    out[[cc]][is.na(out[[cc]]) | out[[cc]] == ""] <- "NA"
  }
  
  # ---- cleaning ----
  out <- out[is.finite(beta)]
  out <- out[is.finite(se) & se > 0]
  out <- out[is.finite(varbeta) & varbeta > 0]
  out <- out[is.finite(p) & p > 0 & p <= 1]
  out <- out[is.finite(chr) & is.finite(pos)]
  out <- out[!(is.na(ea) | ea == "" | ea == "NA") & !(is.na(oa) | oa == "" | oa == "NA")]
  out <- out[is.finite(N) & N > 0]
  
  # MAF optional, but if present must be (0,0.5)
  out <- out[is.na(MAF) | (is.finite(MAF) & MAF > 0 & MAF < 0.5)]
  
  # cc: require s in (0,1) IF provided; allow NA only if NO study totals were given anywhere
  if (type == "cc") {
    if (totals_from_cols_available || (is.finite(cases_total) && is.finite(controls_total) && (cases_total + controls_total) > 0)) {
      out <- out[is.finite(s) & s > 0 & s < 1]
    } else {
      out <- out[is.na(s) | (is.finite(s) & s > 0 & s < 1)]
    }
  }
  
  # ---- enforce EXACT final column order (DO NOT CHANGE) ----
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

out <- standardize_gwas_for_coloc(
  gwas_path      = infile,
  trait_name     = trait,
  type           = type,
  prefer_snp     = prefer,
  sdY_override   = sdY_override,
  N_override     = N_override,
  cases_total    = cases_total,
  controls_total = controls_total
)

fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("[%s] wrote %d rows -> %s\n", trait, nrow(out), outfile))