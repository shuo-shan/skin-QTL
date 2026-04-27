library(plinkQC)
library(ggplot2)
indir="/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata"
qcdir="/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata"
name <- "data" # # Because your files are test.bed, test.bim, test.fam
path2plink <- "/share/pkg/plink/1.90b6.27/plink"

# individual level QC
fail_individuals <- perIndividualQC(
  indir = indir,
  qcdir = qcdir,
  name = name,
  path2plink = path2plink,
  refSamplesFile = file.path(qcdir, "HapMap_ID2Pop.txt"),
  refColorsFile = file.path(qcdir, "HapMap_PopColors.txt"),
  prefixMergedDataset = file.path(qcdir, "data.HapMapIII"),  # only if you downloaded this .bed/.bim/.fam set (not just eigenvec)
  do.run_check_ancestry = FALSE,
  interactive = TRUE,
  verbose = TRUE,
  showPlinkOutput=TRUE
)

overview_individuals <- overviewPerIndividualQC(fail_individuals,
                                                interactive=TRUE)

saveRDS(overview_individuals, file = file.path(qcdir, "overview_individuals.rds"))


# Save all individual QC plots
for (plot_name in names(overview_individuals)) {
  ggsave(
    filename = file.path(qcdir, paste0("individualQC_", plot_name, ".png")),
    plot = overview_individuals[[plot_name]],
    width = 6, height = 4, dpi = 300
  )
}



# marker level QC:
fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
                            path2plink=path2plink,
                            verbose=TRUE, interactive=TRUE,
                            showPlinkOutput=FALSE)


overview_marker <- overviewPerMarkerQC(fail_markers, interactive=TRUE)

saveRDS(overview_marker, file = file.path(qcdir, "overview_marker.rds"))

# Save all marker QC plots
for (plot_name in names(overview_marker)) {
  ggsave(
    filename = file.path(qcdir, paste0("markerQC_", plot_name, ".png")),
    plot = overview_marker[[plot_name]],
    width = 6, height = 4, dpi = 300
  )
}


# Create cleanData of individuals and variants that passed QC
Ids  <- cleanData(indir=indir, qcdir=qcdir, name=name, path2plink=path2plink,
                  verbose=TRUE, showPlinkOutput=FALSE)



