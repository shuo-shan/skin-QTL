library(rtracklayer)
library(tidyverse)

chain_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/data/hg19ToHg38.over.chain"
jin_dir    <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/vitiligo_jin2016"
out_dir    <- file.path(jin_dir, "hg38")
dir.create(out_dir, showWarnings = FALSE)

# load chain once
chain <- import.chain(chain_file)

chroms <- c(1:22, "X")

for (ch in chroms) {
  cat(sprintf("Processing chr%s...\n", ch))
  
  in_file <- file.path(jin_dir, sprintf("GWAS123chr%scmh.txt.gz", ch))
  if (!file.exists(in_file)) {
    cat(sprintf("  File not found, skipping\n"))
    next
  }
  
  jin_raw <- read_table(in_file, show_col_types = FALSE) %>%
    rename_with(str_to_lower) %>%
    mutate(
      snp     = str_to_lower(snp),
      beta    = if ("orx" %in% names(.)) log(orx) else log(or),  # ← 改这里
      varbeta = se^2
    )
  
  gr_hg19 <- GRanges(
    seqnames = paste0("chr", jin_raw$chr),
    ranges   = IRanges(start = jin_raw$bp, end = jin_raw$bp)
  )
  
  gr_hg38      <- liftOver(gr_hg19, chain)
  mapped       <- lengths(gr_hg38) == 1
  gr_hg38_flat <- unlist(gr_hg38[mapped])
  cat(sprintf("  Total: %d  |  Mapped: %d  |  Dropped: %d\n",
              nrow(jin_raw), sum(mapped), sum(!mapped)))
  
  pos_df <- tibble(
    chr = as.character(seqnames(gr_hg38_flat)),
    pos = start(gr_hg38_flat)
  )
  
  jin_hg38 <- bind_cols(jin_raw[mapped, ] %>% select(-chr, -bp), pos_df) %>%
    select(chr, pos, snp, everything())
  
  out_file <- file.path(out_dir, sprintf("jin2016_vitiligo_chr%s_hg38.tsv.gz", ch))
  write_tsv(jin_hg38, gzfile(out_file))
  cat(sprintf("  Saved: %s\n", out_file))
}

cat("\nAll done!\n")

# sanity check
jin19 <- read_tsv(
  file.path(out_dir, "jin2016_vitiligo_chr19_hg38.tsv.gz"),
  show_col_types = FALSE
)
jin19 %>% filter(snp == "rs6510827") %>% select(snp, chr, pos, p, beta)


# sanity check
jin19 <- read_tsv(file.path(out_dir, "jin2016_vitiligo_chr19_hg38.tsv.gz"), show_col_types = FALSE)
jin19 %>% filter(snp == "rs6510827") %>% select(snp, chr, pos, p, beta)