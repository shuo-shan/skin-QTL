## QTL Gene Discovery Visualization
## Upset plots (per celltype, eQTL & reQTL) + Venn diagrams (per condition x QTLtype)

library(tidyverse)
library(UpSetR)
library(ggvenn)   # install.packages("ggvenn") -- or use VennDiagram package

# ---- Load data ----
dat <- read.delim("all_QTL_genes_FDR05.txt", stringsAsFactors = FALSE)

celltypes  <- unique(dat$celltype)    # FRB, KRT, MEL
conditions <- unique(dat$condition)   # PBS, IFNB, IFNG, TNF
qtltypes   <- unique(dat$QTLtype)     # eQTL, reQTL

# ============================================================
# 1. UPSET PLOTS: per celltype, genes across conditions
#    Separately for eQTL and reQTL
# ============================================================

make_upset <- function(celltype_val, qtltype_val) {
  sub <- dat %>%
    filter(celltype == celltype_val, QTLtype == qtltype_val)
  
  if (nrow(sub) == 0) {
    message(sprintf("No data for %s %s -- skipping", celltype_val, qtltype_val))
    return(invisible(NULL))
  }
  
  # Build binary membership matrix: rows = genes, cols = conditions
  genes_all <- unique(sub$gene)
  conds_present <- unique(sub$condition)
  
  mat <- sapply(conds_present, function(cond) {
    as.integer(genes_all %in% sub$gene[sub$condition == cond])
  })
  rownames(mat) <- genes_all
  mat <- as.data.frame(mat)
  
  # Order conditions: PBS first, then interferons
  cond_order <- intersect(c("PBS", "IFNB", "IFNG", "TNF"), colnames(mat))
  mat <- mat[, cond_order, drop = FALSE]
  
  title_str <- sprintf("%s %s — Genes per Condition", celltype_val, qtltype_val)
  
  upset(mat,
        sets = cond_order,
        order.by = "freq",
        main.bar.color = ifelse(qtltype_val == "eQTL", "#2C7BB6", "#D7191C"),
        sets.bar.color = "grey40",
        text.scale = 1.4,
        mainbar.y.label = "Gene Intersection Size",
        sets.x.label = "Genes per Condition",
        mb.ratio = c(0.6, 0.4))
  
  title(main = title_str, line = 2.5)
}

# Save all upset plots
pdf("QTL_UpsetPlots.pdf", width = 10, height = 6)
for (ct in celltypes) {
  message("eQTL upset: ", ct)
  make_upset(ct, "eQTL")
  
  message("reQTL upset: ", ct)
  make_upset(ct, "reQTL")
}
dev.off()

message("Upset plots saved to QTL_UpsetPlots.pdf")

# ============================================================
# 2. VENN DIAGRAMS: per condition x QTLtype, genes per celltype
# ============================================================

# Color palette per celltype
ct_colors <- c(FRB = "#E69F00", KRT = "#009E73", MEL = "#CC79A7")

make_venn <- function(cond_val, qtltype_val) {
  sub <- dat %>%
    filter(condition == cond_val, QTLtype == qtltype_val)
  
  if (nrow(sub) == 0) {
    message(sprintf("No data for %s %s -- skipping", cond_val, qtltype_val))
    return(invisible(NULL))
  }
  
  cts_present <- intersect(celltypes, unique(sub$celltype))
  gene_sets <- setNames(
    lapply(cts_present, function(ct) unique(sub$gene[sub$celltype == ct])),
    cts_present
  )
  
  p <- ggvenn(gene_sets,
              fill_color = ct_colors[cts_present],
              stroke_size = 0.5,
              set_name_size = 5,
              text_size = 4) +
    labs(title = sprintf("%s %s — Gene Overlap Across Cell Types", cond_val, qtltype_val)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  return(p)
}

# Save all venn diagrams
pdf("QTL_VennDiagrams.pdf", width = 7, height = 7)
for (qt in qtltypes) {
  for (cond in c("PBS", "IFNB", "IFNG", "TNF")) {
    message("Venn: ", cond, " ", qt)
    p <- make_venn(cond, qt)
    if (!is.null(p)) print(p)
  }
}
dev.off()

message("Venn diagrams saved to QTL_VennDiagrams.pdf")
message("Done!")
