library(data.table)

# ---------------- User metadata / paths ----------------
trait_name  <- "alopecia_areata"
type        <- "cc"                 # alopecia areata is case-control
prefer_snp  <- "rsid"               # or "variant_id"
infile      <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/alopecia_areata/GCST90478824.h.tsv.gz"
outfile     <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/alopecia_areata/standardized_GCST90478824.h.tsv.gz"

# Optional overrides (leave as NA_real_ if not used)
N_override      <- NA_real_
cases_total     <- NA_real_
controls_total  <- NA_real_
sdY_override    <- NA_real_         # quant-only; ignored for cc

gwas_path <- infile

# ---------------- Standardizer ----------------
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
  pick <- function(...) { c(...)[c(... ) %in% cn][1] }
  
  # ---- map columns ----
  chr_col <- pick("chromosome","chr","CHROM","hm_chrom")
  pos_col <- pick("base_pair_location","pos","POS","hm_pos")
  
  ea_col  <- pick("effect_allele","ea","EA","A1","alt","ALT","hm_effect_allele")
  oa_col  <- pick("other_allele","oa","OA","A2","ref","REF","hm_other_allele")
  
  p_col   <- pick("p_value","p","P","pval","pval_nominal")
  
  eaf_col <- pick("effect_allele_frequency","eaf","EAF","af","freq","hm_effect_allele_frequency")
  
  rsid_col <- pick("rsid","hm_rsid")
  vid_col  <- pick("variant_id","hm_variant_id")
  
  snp_col <- if (prefer_snp == "variant_id") {
    pick("variant_id","hm_variant_id","rsid","hm_rsid")
  } else {
    pick("rsid","hm_rsid","variant_id","hm_variant_id")
  }
  
  beta_col <- pick("beta","hm_beta","effect","BETA")
  or_col   <- pick("odds_ratio","OR","hm_odds_ratio")
  
  se_col   <- pick("standard_error","se","SE","stderr","hm_standard_error","hm_se")
  ci_lo_col <- pick("ci_lower","lower_ci","lci","ci_l","hm_ci_lower")
  ci_hi_col <- pick("ci_upper","upper_ci","uci","ci_u","hm_ci_upper")
  
  n_col <- pick("cum_eff_sample_size","n","N","Neff","effective_N","effN","total_n")
  
  ncase_col <- pick("num_cases","cases","ncase","N_cases","case_n")
  nctrl_col <- pick("num_controls","controls","ncontrol","N_controls","control_n")
  
  # ---- required checks ----
  req <- c(chr_col, pos_col, ea_col, oa_col, p_col)
  if (any(is.na(req))) {
    stop("Missing required cols in ", gwas_path, "\nHave: ", paste(cn, collapse = ", "))
  }
  if (is.na(beta_col) && is.na(or_col)) {
    stop("No beta or odds_ratio column in ", gwas_path)
  }
  
  # ---- base output skeleton in FINAL REQUIRED ORDER ----
  out <- dt[, .(
    trait = as.character(trait_name),
    type  = as.character(type),
    
    snp   = if (!is.na(snp_col)) as.character(get(snp_col)) else NA_character_,
    rsid  = if (!is.na(rsid_col)) as.character(get(rsid_col)) else NA_character_,
    variant_id = if (!is.na(vid_col)) as.character(get(vid_col)) else NA_character_,
    
    chr = suppressWarnings(as.integer(get(chr_col))),
    pos = suppressWarnings(as.integer(get(pos_col))),
    
    ea = as.character(get(ea_col)),
    oa = as.character(get(oa_col)),
    
    beta = NA_real_,
    se   = NA_real_,
    varbeta = NA_real_,
    p    = suppressWarnings(as.numeric(get(p_col))),
    
    effect_scale = NA_character_,
    se_source    = NA_character_,
    
    eaf = if (!is.na(eaf_col)) suppressWarnings(as.numeric(get(eaf_col))) else NA_real_,
    MAF = NA_real_,
    
    N = NA_real_,
    N_source = NA_character_,
    
    num_cases    = if (!is.na(ncase_col)) suppressWarnings(as.numeric(get(ncase_col))) else NA_real_,
    num_controls = if (!is.na(nctrl_col)) suppressWarnings(as.numeric(get(nctrl_col))) else NA_real_,
    
    s = NA_real_,
    s_source = NA_character_,
    
    sdY = NA_real_,
    sdY_source = NA_character_
  )]
  
  # normalize blanks -> NA (you’ll write NA as literal "NA" via fwrite(na="NA"))
  for (ccol in c("snp","rsid","variant_id","ea","oa")) {
    out[get(ccol) %chin% c("", " "), (ccol) := NA_character_]
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
  
  # ---- se + se_source (priority: SE, then CI, then p+beta) ----
  out[, se := NA_real_]
  out[, se_source := NA_character_]
  
  # 1) explicit SE column
  if (!is.na(se_col)) {
    se0 <- suppressWarnings(as.numeric(dt[[se_col]]))
    out[, se := se0]
    out[is.finite(se) & se > 0, se_source := "SE"]
  }
  
  # 2) CI (assume 95%; infer on log scale)
  if (!is.na(ci_lo_col) && !is.na(ci_hi_col)) {
    lo <- suppressWarnings(as.numeric(dt[[ci_lo_col]]))
    hi <- suppressWarnings(as.numeric(dt[[ci_hi_col]]))
    ok_ci <- is.finite(lo) & is.finite(hi) & lo > 0 & hi > 0
    se_ci <- (log(hi) - log(lo)) / (2 * 1.96)
    
    idx <- which((!is.finite(out$se) | out$se <= 0) & ok_ci & is.finite(se_ci) & se_ci > 0)
    if (length(idx) > 0) {
      out[idx, se := se_ci[idx]]
      out[idx, se_source := "CI"]
    }
  }
  
  # 3) infer from p + beta (two-sided)
  idx_pb <- which((!is.finite(out$se) | out$se <= 0) & is.finite(out$p) & out$p > 0 & out$p <= 1 & is.finite(out$beta) & out$beta != 0)
  if (length(idx_pb) > 0) {
    z <- abs(qnorm(out$p[idx_pb] / 2, lower.tail = FALSE))
    se_pb <- abs(out$beta[idx_pb]) / z
    ok_pb <- is.finite(se_pb) & se_pb > 0
    if (any(ok_pb)) {
      ii <- idx_pb[ok_pb]
      out[ii, se := se_pb[ok_pb]]
      out[ii, se_source := "p+beta"]
    }
  }
  
  out[, varbeta := se^2]
  
  # ---- eaf -> MAF ----
  out[, MAF := ifelse(is.na(eaf), NA_real_, pmin(eaf, 1 - eaf))]
  
  # ---- N + N_source ----
  N0 <- if (!is.na(n_col)) suppressWarnings(as.numeric(dt[[n_col]])) else rep(NA_real_, nrow(out))
  out[, N := N0]
  out[is.finite(N) & N > 0, N_source := if (!is.na(n_col)) n_col else NA_character_]
  
  # if N missing/bad and per-SNP cases/controls exist: N = cases + controls
  needN <- which(!(is.finite(out$N) & out$N > 0) & is.finite(out$num_cases) & is.finite(out$num_controls))
  if (length(needN) > 0) {
    Ncc <- out$num_cases[needN] + out$num_controls[needN]
    okNcc <- is.finite(Ncc) & Ncc > 0
    if (any(okNcc)) {
      ii <- needN[okNcc]
      out[ii, N := Ncc[okNcc]]
      out[ii, N_source := "derived_cases_controls"]
    }
  }
  
  # if still missing/bad: N_override
  if (is.finite(N_override) && N_override > 0) {
    needN2 <- which(!(is.finite(out$N) & out$N > 0))
    if (length(needN2) > 0) {
      out[needN2, N := N_override]
      out[needN2, N_source := "N_override"]
    }
  } else {
    out[!(is.finite(N) & N > 0), N_source := NA_character_]
  }
  
  # ---- s + s_source ----
  if (type == "cc") {
    # detect whether num_cases/num_controls are per-SNP vs study totals repeated
    uc <- unique(out[is.finite(num_cases), num_cases])
    uC <- unique(out[is.finite(num_controls), num_controls])
    uc <- uc[is.finite(uc)]
    uC <- uC[is.finite(uC)]
    
    if (length(uc) == 1 && length(uC) == 1) {
      # study totals repeated per row
      denom <- uc + uC
      if (is.finite(denom) && denom > 0) {
        out[, s := uc / denom]
        out[, s_source := "study_total"]
      } else {
        out[, s := NA_real_]
        out[, s_source := "NA"]
      }
    } else if (length(uc) > 1 || length(uC) > 1) {
      # per-SNP
      denom <- out$num_cases + out$num_controls
      out[is.finite(denom) & denom > 0 & is.finite(num_cases),
          `:=`(s = num_cases / denom, s_source = "per_snp")]
      out[!(is.finite(s) & s > 0 & s < 1), `:=`(s = NA_real_, s_source = "NA")]
    } else if (is.finite(cases_total) && is.finite(controls_total) && (cases_total + controls_total) > 0) {
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
  
  # ---- sdY + sdY_source ----
  if (type == "quant") {
    if (is.finite(sdY_override)) {
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
  
  # ---- cleaning (per your rules) ----
  out <- out[
    is.finite(beta) &
      is.finite(se) & se > 0 &
      is.finite(varbeta) & varbeta > 0 &
      is.finite(p) & p > 0 & p <= 1 &
      is.finite(chr) & is.finite(pos) &
      !is.na(ea) & !is.na(oa) &
      is.finite(N) & N > 0
  ]
  
  # do NOT require eaf/MAF; but if MAF present, require (0 < MAF < 0.5)
  out <- out[is.na(MAF) | (MAF > 0 & MAF < 0.5)]
  
  # for cc: require s in (0,1) if provided; otherwise allow NA only if no study totals were given
  if (type == "cc") {
    # if we produced s from study totals or per_snp, it must be valid
    out <- out[is.na(s) | (s > 0 & s < 1)]
  }
  
  # ---- ensure exact column order (including sdY at end) ----
  std_cols <- c(
    "trait","type","snp","rsid","variant_id","chr","pos","ea","oa",
    "beta","se","varbeta","p","effect_scale","se_source",
    "eaf","MAF","N","N_source","num_cases","num_controls",
    "s","s_source","sdY","sdY_source"
  )
  out <- out[, ..std_cols]
  
  out[]
}

# ---------------- run ----------------
out <- standardize_gwas_for_coloc(
  gwas_path     = infile,
  trait_name    = trait_name,
  type          = type,
  prefer_snp    = prefer_snp,
  sdY_override  = sdY_override,
  N_override    = N_override,
  cases_total   = cases_total,
  controls_total= controls_total
)

fwrite(out, outfile, sep = "\t", na = "NA", quote=F)
cat(sprintf("[%s] wrote %d rows -> %s\n", trait_name, nrow(out), outfile))
