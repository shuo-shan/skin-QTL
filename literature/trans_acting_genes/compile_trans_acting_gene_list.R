library(msigdbr)
library(data.table)
library(tidyverse)
library(magrittr)
library(dplyr)


dir.literature = "/pi/manuel.garber-umw/human/skin/eQTLs/literature/trans_acting_genes/"

# ----- Kinases ----- 
# Original — catalytic subunits
kinase_go_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl("KINASE_ACTIVITY", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Add: kinase regulator activity (regulatory + scaffolding subunits)
kinase_regulator <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl("KINASE_REGULATOR_ACTIVITY|KINASE_ACTIVATOR_ACTIVITY|
                KINASE_INHIBITOR_ACTIVITY", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Add: GO:BP — kinase cascade membership catches context-specific components
kinase_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  filter(grepl("PROTEIN_PHOSPHORYLATION|KINASE_CASCADE|
                MAPK_CASCADE|JAK_STAT_CASCADE|
                PHOSPHORYLATION", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Add: KinHub / Manning kinome — gold standard curated human kinome (518 genes)
# Download: http://kinase.com/human/kinome/ -> "Kinase Gene List" CSV
# This is the most complete source for the canonical kinome
kinhub <- read.csv("/pi/manuel.garber-umw/human/skin/eQTLs/literature/trans_acting_genes/kinase_database_from_Kinome_Sudarsanam_2002.txt", header=F)
kinase_kinhub <- kinhub %>% pull(V1) %>% unique()

# Add: REACTOME — catches pathway-specific kinase complexes
kinase_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  filter(grepl("KINASE|PHOSPHORYLATION_OF|CDK|MAPK|AKT|
                JAK|SRC_FAMILY", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Pseudokinase seed (no clean GO term)
pseudokinase_seed <- c(
  "ERBB3","STK40","TRIB1","TRIB2","TRIB3",
  "STRADB","STRADA","HASPIN","ADCK1","ADCK2",
  "WNK1","WNK2","WNK3","WNK4",   # atypical Lys-less kinases
  "DSTYK","HSPB8","PTK7","ROR1","ROR2",
  "EPHB6","STRAD"
)

kinase_genes_full <- Reduce(union, list(
  kinase_go_mf,
  kinase_regulator,
  kinase_reactome,
  pseudokinase_seed,
  kinase_bp,
  kinase_kinhub  
))


cat("Original GO:MF kinases:", length(kinase_go_mf), "\n")
cat("Expanded kinase list:  ", length(kinase_genes_full), "\n")
rm(kinase_go_mf,kinase_regulator,kinase_reactome,pseudokinase_seed,kinase_bp,kinase_kinhub, kinhub)

# ----- RNA binding protein genes -----
# Original
rbp_go_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl("RNA_BINDING", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Add: more specific RBP-relevant GO:MF terms
rbp_specific <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl(paste(c(
    "MRNA_BINDING",
    "MRNA_3_UTR_BINDING",
    "MRNA_5_UTR_BINDING",
    "AU_RICH_ELEMENT_BINDING",
    "POLY_A_RNA_BINDING",
    "MRNA_STABILITY",
    "RRNA_BINDING",
    "NCRNA_BINDING",
    "TELOMERIC_RNA_BINDING"
  ), collapse = "|"), gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Add: GO:BP for post-transcriptional processes
rbp_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  filter(grepl(paste(c(
    "MRNA_STABILIZATION",
    "MRNA_DESTABILIZATION",
    "REGULATION_OF_MRNA_STABILITY",
    "MRNA_PROCESSING",
    "RNA_SPLICING",
    "NUCLEAR_MRNA_SPLICING",
    "MRNA_EXPORT",
    "REGULATION_OF_TRANSLATION",
    "NONSENSE_MEDIATED_MRNA_DECAY"
  ), collapse = "|"), gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Gold standard: ENCODE eCLIP-validated RBPs (356 proteins, hg38)
# This is the best available experimentally validated RBP list
# Download: https://www.encodeproject.org/eclip/ 
#           or pre-compiled at https://rbpdb.ccbr.utoronto.ca/
rbp_encode <- read.csv(paste0(dir.literature,"/RNA_binding_proteins_ENCODE_eCLIP_targets.txt"), header=F)
rbp_rbpencode <- rbp_encode %>% pull(V1) %>% unique()

rbpdb <- read.csv(paste0(dir.literature,"/RBPDB_v1.3.1_proteins_human_2012-11-21.csv"), header=F)
rbp_rbpdb <- rbpdb %>% pull(V5) %>% unique()
rbp_rbpdb <- rbp_rbpdb[rbp_rbpdb != ""]

# Also good: Gerstberger et al. 2014 (Nat Rev Genetics) curated ~1,500 human RBPs
# Available as Supplementary Table from the paper
rbp_gerstberger <- read.csv(paste0(dir.literature, "/RNA_binding_proteins_Gerstberger_2014_NatRevGenetics.txt"), header=F)
rbp_rbp.gerstberger <- rbp_gerstberger %>% pull(V1) %>% unique()

# Genes to EXCLUDE — RNA-binding but not post-transcriptional regulators
exclude_from_rbp <- c(
  # Ribosomal proteins (structural, not regulatory)
  grep("^RPL|^RPS", rbp_go_mf, value = TRUE),
  # Core RNA Pol subunits
  "POLR1A","POLR1B","POLR2A","POLR2B","POLR3A","POLR3B",
  # Aminoacyl-tRNA synthetases (bind tRNA but not mRNA regulators)
  "AARS1","CARS1","DARS1","EARS1","FARS1","GARS1",
  "HARS1","IARS1","KARS1","LARS1","MARS1","NARS1",
  "PARS1","QARS1","RARS1","SARS1","TARS1","VARS1","WARS1","YARS1"
)

rbp_genes_full <- Reduce(union, list(
  rbp_go_mf, rbp_specific, rbp_bp, rbp_rbpencode, rbp_rbpdb, rbp_rbp.gerstberger
)) %>%
  setdiff(exclude_from_rbp)

cat("Original GO:MF RBPs:  ", length(rbp_go_mf), "\n")
cat("After adding terms:   ", length(unique(c(rbp_go_mf, rbp_specific, rbp_bp, rbp_rbpencode, rbp_rbpdb, rbp_rbp.gerstberger))), "\n")
cat("After exclusions:     ", length(rbp_genes_full), "\n")
rm(rbp_go_mf, rbp_specific, rbp_bp, rbp_rbpencode, rbp_rbpdb, rbp_rbp.gerstberger, rbp_encode, rbp_gerstberger, rbpdb, exclude_from_rbp)

# ---- Transcription factors ----
tf_file="/pi/manuel.garber-umw/human/skin/eQTLs/literature/Lambert_2018_human_TFs.txt"
tf_genes <- unique(fread(tf_file, header=F)$V1)

cat("TF genes: ", length(tf_genes),"\n")
rm(tf_file)

# ---- Chromatin remodelers ----
# Chromatin remodeling complexes (Biological Process)
chromatin_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  filter(grepl("CHROMATIN_REMODELING|CHROMATIN_ORGANIZATION", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Chromatin remodeler enzymatic activity (Molecular Function)
chromatin_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl("CHROMATIN_BINDING|HISTONE_DEACETYLASE|HISTONE_ACETYLTRANSFERASE|
                HISTONE_METHYLTRANSFERASE|HISTONE_DEMETHYLASE", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# REACTOME curated complexes — often more precise than GO for named complexes
chromatin_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  filter(grepl("CHROMATIN|HDAC|HAT_MEDIATED|PRC2|SWI_SNF|NURF|NURD", gs_name)) %>%
  pull(gene_symbol) %>% unique()

chromatin_genes <- unique(c(chromatin_bp, chromatin_mf, chromatin_reactome))

cat("chromatin remodeler genes: ", length(chromatin_genes),"\n")
rm(chromatin_bp, chromatin_mf, chromatin_reactome)

# ---- Cytokine Receptors and signaling adaptors ----
# KEGG cytokine pathways
receptor_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY") %>%
  filter(grepl("CYTOKINE|JAK_STAT|TOLL|NF_KAPPA|TNF|IL|CHEMOKINE", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# REACTOME — more fine-grained
receptor_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  filter(grepl("CYTOKINE_SIGNALING|INTERLEUKIN|INTERFERON|JAK|TOLL_LIKE_RECEPTOR|
                SIGNALING_BY_INTERLEUKINS|MYD88|TRIF", gs_name)) %>%
  pull(gene_symbol) %>% unique()

receptor_genes <- union(receptor_kegg, receptor_reactome)

cat("cytokine receptor genes: ", length(receptor_genes),"\n")
rm(receptor_kegg, receptor_reactome)


# ---- Phosphatases  ----
phosphatase_genes <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl("PHOSPHATASE_ACTIVITY|PROTEIN_PHOSPHATASE|PHOSPHOPROTEIN_PHOSPHATASE", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# ---- Ubiquitin / proteasome / protein degradation machinery ----
ubiquitin_genes <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl("UBIQUITIN_LIGASE|UBIQUITIN_PROTEIN_LIGASE|DEUBIQUITINASE|
                PROTEASOME|CULLIN", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Supplement with REACTOME
ubiquitin_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  filter(grepl("UBIQUITIN|PROTEASOME|SCF_COMPLEX|APC_C", gs_name)) %>%
  pull(gene_symbol) %>% unique()

ubiquitin_genes <- union(ubiquitin_genes, ubiquitin_reactome)
rm(ubiquitin_reactome)

# ---- Metabolic Enzymes (The "Moonlighting" Regulators) ----
metabolic_genes <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  filter(grepl("TCA_CYCLE|GLYCOLYSIS|PENTOSE_PHOSPHATE|ONE_CARBON_POOL|
                FOLATE|ACETYL_COA|NAD_METABOLISM|SIRTUIN", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# ---- Coactivators and Mediator complex ----
coactivator_genes <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl("TRANSCRIPTION_COACTIVATOR|TRANSCRIPTION_COREGULATOR|
                MEDIATOR_COMPLEX|ENHANCER_BINDING", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# ---- miRNA/lncRNA pathway genes (RISC, DROSHA, DICER axis) ----
mirna_genes <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  filter(grepl("MIRNA_PROCESSING|MIRNA_BIOGENESIS|RNA_SILENCING|
                POSTTRANSCRIPTIONAL_GENE_SILENCING", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# ---- Nuclear transport (importins/exportins) ----
transport_genes <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl("NUCLEAR_IMPORT|NUCLEAR_EXPORT|IMPORTIN|EXPORTIN|
                NUCLEAR_TRANSPORT", gs_name)) %>%
  pull(gene_symbol) %>% unique()

# ---- Second messenger enzymes ----
second_messenger_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl(paste(c(
    "SECOND_MESSENGER",
    "PHOSPHOLIPASE_ACTIVITY",
    "ADENYLYL_CYCLASE_ACTIVITY",
    "GUANYLYL_CYCLASE_ACTIVITY",
    "PHOSPHOINOSITIDE_3_KINASE_ACTIVITY",
    "PHOSPHATIDYLINOSITOL_KINASE_ACTIVITY",
    "DIACYLGLYCEROL_KINASE_ACTIVITY",
    "CYCLIC_NUCLEOTIDE"
  ), collapse = "|"), gs_name)) %>%
  pull(gene_symbol) %>% unique()

second_messenger_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  filter(grepl(paste(c(
    "CAMP_MEDIATED_SIGNALING",
    "CGMP_MEDIATED_SIGNALING",
    "INOSITOL_PHOSPHATE_METABOLIC",
    "PHOSPHATIDYLINOSITOL_SIGNALING",
    "SECOND_MESSENGER_MEDIATED_SIGNALING",
    "REGULATION_OF_CAMP_LEVELS",
    "IP3_SIGNALING"
  ), collapse = "|"), gs_name)) %>%
  pull(gene_symbol) %>% unique()

# REACTOME is particularly good for PI3K isoforms
second_messenger_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  filter(grepl(paste(c(
    "PI3K",
    "PIP3",
    "ADENYLATE_CYCLASE",
    "PHOSPHOLIPASE_C",
    "DAG_AND_IP3_SIGNALING",
    "SECOND_MESSENGER"
  ), collapse = "|"), gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Hard-coded seed for key isoforms that GO/REACTOME terms sometimes miss
second_messenger_seed <- c(
  # Adenylyl cyclases
  "ADCY1","ADCY2","ADCY3","ADCY4","ADCY5","ADCY6","ADCY7","ADCY8","ADCY9","ADCY10",
  # PI3K catalytic + regulatory subunits
  "PIK3CA","PIK3CB","PIK3CD","PIK3CG",
  "PIK3R1","PIK3R2","PIK3R3","PIK3R5","PIK3R6",
  # Phospholipase C isoforms
  "PLCB1","PLCB2","PLCB3","PLCB4",
  "PLCG1","PLCG2",
  "PLCD1","PLCD3","PLCD4",
  "PLCE1","PLCZ1",
  # Phospholipase A2 isoforms
  "PLA2G1B","PLA2G2A","PLA2G4A","PLA2G4B","PLA2G4C",
  "PLA2G5","PLA2G6","PLA2G7",
  # Phosphodiesterases (degrade cAMP/cGMP — equally important for tuning signal)
  "PDE1A","PDE2A","PDE3A","PDE3B","PDE4A","PDE4B","PDE4C","PDE4D",
  "PDE5A","PDE7A","PDE7B","PDE8A","PDE8B"
)

second_messenger_genes <- Reduce(union, list(
  second_messenger_mf,
  second_messenger_bp,
  second_messenger_reactome,
  second_messenger_seed
))

cat("Second messenger enzymes:", length(second_messenger_genes), "\n")

rm(second_messenger_mf, second_messenger_bp, second_messenger_reactome,second_messenger_seed)

# ---- RNA Pol II transcriptional regulators ----
polii_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  filter(grepl(paste(c(
    "RNA_POLYMERASE_II",
    "TRANSCRIPTION_ELONGATION",
    "TRANSCRIPTION_INITIATION_FROM_RNA_POLYMERASE_II",
    "PAUSE_RELEASE",
    "PROMOTER_PROXIMAL_PAUSING",
    "TRANSCRIPTION_PREINITIATION_COMPLEX"
  ), collapse = "|"), gs_name)) %>%
  pull(gene_symbol) %>% unique()

polii_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  filter(grepl(paste(c(
    "RNA_POLYMERASE_II_CTD",
    "GENERAL_TRANSCRIPTION_INITIATION_FACTOR",
    "TRANSCRIPTION_FACTOR_ACTIVITY_RNA_POLYMERASE_II"
  ), collapse = "|"), gs_name)) %>%
  pull(gene_symbol) %>% unique()

# REACTOME has the most precise terms for pause-release machinery
polii_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  filter(grepl(paste(c(
    "RNA_POL_II",
    "RNAP_II",
    "CTD_PHOSPHORYLATION",
    "TRANSCRIPTION_COUPLED",
    "PTEFB",          # CDK9-containing complex
    "MEDIATOR",
    "GENERAL_TRANSCRIPTION_FACTORS",
    "FORMATION_OF_THE_TRANSCRIPTION_PREINITIATION_COMPLEX",
    "PAUSING_AND_RECOVERY"
  ), collapse = "|"), gs_name)) %>%
  pull(gene_symbol) %>% unique()

# Seed list — pause-release factors are poorly represented in MSigDB terms
polii_seed <- c(
  # CDKs that phosphorylate Pol II CTD
  "CDK7","CDK8","CDK9","CDK12","CDK13",
  # Cyclin partners
  "CCNH","CCNT1","CCNT2","CCNK",
  # TFIIH complex
  "GTF2H1","GTF2H2","GTF2H3","GTF2H4","GTF2H5",
  "ERCC2","ERCC3","MNAT1",
  # P-TEFb / SEC complex (pause release)
  "BRD4","HEXIM1","HEXIM2","LARP7","MePCE",
  # DSIF complex (DRB-sensitivity inducing factor)
  "SUPT4H1","SUPT5H",
  # NELF complex (negative elongation factor)
  "NELFA","NELFB","NELFC","NELFE",
  # SPT6 / FACT (elongation-coupled histone chaperones)
  "SUPT6H","SSRP1","SUPT16H",
  # Super elongation complex (SEC)
  "ELL","ELL2","ELL3","AFF1","AFF4","ENL","AF9",
  # General transcription factors
  "GTF2B","GTF2E1","GTF2E2","GTF2F1","GTF2F2",
  "TAF1","TAF4","TAF5","TAF6","TAF7","TAF9","TAF12",
  # XPB/XPD (TFIIH helicases)
  "ERCC2","ERCC3"
)

polii_genes <- Reduce(union, list(
  polii_bp,
  polii_mf,
  polii_reactome,
  polii_seed
))

cat("RNA Pol II regulators:", length(polii_genes), "\n")

rm(polii_bp, polii_mf, polii_reactome, polii_seed)

# ---- Phase separation / condensate proteins ----
# No clean MSigDB term exists — strategy is:
#   (a) curated seed from key papers
#   (b) IDR-based expansion using flDPnn or FuzDrop scores
#   (c) GO "nuclear speckle / P-body" as a proxy compartment

# (a) Curated seed from Boija et al. 2018 (Cell), Sabari et al. 2018 (Science),
#     Guo et al. 2019, and condensate reviews up to 2024
condensate_seed <- c(
  # FET family (archetypical IDR-containing)
  "FUS","EWSR1","TAF15",
  # Mediator / coactivator condensates
  "MED1","MED12","MED13","MED14","MED15","MED24",
  "EP300","CREBBP",
  # Master TFs known to form condensates
  "MYC","MYB","OCT4","SOX2","NANOG",    # pluripotency — relevant as benchmarks
  "RELA","IRF3","IRF7","STAT3",          # immune-relevant condensate TFs
  # Super-enhancer coactivators
  "BRD4","CDK8","CDK19",
  # Stress granule / P-body components with trans-regulatory roles
  "G3BP1","G3BP2","CAPRIN1",
  "DDX3X","DDX6","EDC4",
  "TIA1","TIAL1",
  # hnRNP proteins with IDRs (overlap with RBPs but condensate-specific role)
  "HNRNPA1","HNRNPA2B1","HNRNPD",
  # Nuclear speckle scaffold
  "SRSF1","SRSF2","SON","MALAT1",       # MALAT1 = lncRNA but organizes speckles
  # Paraspeckle
  "NEAT1","SFPQ","NONO","FUS",
  # Polycomb condensates
  "EZH2","EED","SUZ12","CBX2","CBX8",
  # Coactivator condensates in immune activation
  "RBBP5","ASH2L","WDR5",               # MLL/COMPASS complex
  "KAT6A","KAT6B","KAT7"               # MOZ/MORF acetyltransferases
)

# (b) GO proxy: proteins annotated to condensate-like compartments
condensate_go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>%
  filter(grepl(paste(c(
    "NUCLEAR_SPECKLE",
    "CAJAL_BODY",
    "STRESS_GRANULE",
    "P_BODY",
    "PROCESSING_BODY",
    "PARASPECKLE",
    "NUCLEAR_STRESS_BODY",
    "TRANSCRIPTION_REGULATOR_COMPLEX"   # enriched for condensate-forming coactivators
  ), collapse = "|"), gs_name)) %>%
  pull(gene_symbol) %>% unique()

condensate_genes <- union(condensate_seed, condensate_go)

cat("Phase separation / condensate proteins:", length(condensate_genes), "\n")

rm(condensate_seed, condensate_go)

# ----  Final union + expressed-gene filter ----
# Combine everything
all_candidate_regulators <- Reduce(union, list(
  chromatin_genes,
  coactivator_genes,
  condensate_genes,
  kinase_genes_full, 
  metabolic_genes,
  mirna_genes,
  phosphatase_genes,
  polii_genes,
  rbp_genes_full,
  receptor_genes,
  second_messenger_genes,
  tf_genes, 
  ubiquitin_genes,
  transport_genes
))


# Annotate with category membership (useful for hotspot interpretation later)
category_map <- bind_rows(
  data.frame(gene = chromatin_genes,        category = "chromatin"),
  data.frame(gene = coactivator_genes,      category = "coactivator"),
  data.frame(gene = condensate_genes,       category = "condensate_protein"),
  data.frame(gene = kinase_genes_full,      category = "kinase"),
  data.frame(gene = metabolic_genes,        category = "metabolic"),
  data.frame(gene = mirna_genes,            category = "miRNA_regulator"),
  data.frame(gene = phosphatase_genes,      category = "phosphatase"),
  data.frame(gene = polii_genes,            category = "polii_regulator"),
  data.frame(gene = rbp_genes_full,         category = "RBP"),
  data.frame(gene = receptor_genes,         category = "receptor"),
  data.frame(gene = second_messenger_genes, category = "second_messenger"),
  data.frame(gene = tf_genes,               category = "transcription_factor"),
  data.frame(gene = ubiquitin_genes,        category = "ubiquitin"),
  data.frame(gene = transport_genes,        category = "transport")
) %>%
  group_by(gene) %>%
  summarise(categories = paste(unique(category), collapse = ";"), .groups = "drop")

cat("\nCompiled trans-acting candidate genes:", nrow(category_map), "\n")
cat("Genes in multiple new categories (potential high-priority):\n")
category_map %>% filter(grepl(";", categories)) %>% print()

# write to file
fwrite(category_map, file = paste0(dir.literature, "/compiled_trans_acting_candidate_genes_and_category.txt"),
       quote=F, col.names=T, sep="\t")
