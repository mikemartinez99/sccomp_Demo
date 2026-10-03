#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Title: sccomp.R
# Description: Demo sccomp on Pijuan-Sala Tal1-/- dataset.
#
# Author: Mike Martinez
# Project: Compute and Conquer Blog Demo
# Date created: 2026-10-03
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Install additional packages
BiocManager::install("MouseGastrulationData")


#----- Libraries
library(sessioninfo)
library(Seurat)
library(scuttle)
library(dplyr)
library(ggplot2)
library(sccomp)
library(MouseGastrulationData)

#----- Set project directory
wd <- "/Users/mike/Desktop/sccomp_demo/"

#----- Specify directories
outputDir <- paste0(wd, "outputs/")
figDir <- paste0(outputDir, "figures/")

#----- Create directories
dirs <- c(outputDir, figDir)
lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD IN THE KANG DATASET
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Read in the Tal1 chimera data
sce.chimera <-  Tal1ChimeraData()

#----- Drop QC labels (not real cell types)
sce.chimera <- sce.chimera[, !sce.chimera$celltype.mapped %in% c("Doublet", "Stripped")]
sce.chimera <- logNormCounts(sce.chimera)

#----- Convert singleCellObj to Seurat
seurat <- as.Seurat(sce.chimera, counts = "counts", data = "logcounts")
table(seurat$sample, seurat$tomato)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SCCOMP RUN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.libPaths(c(path.expand("~/Rlib-sccomp"), .libPaths()))
sccompCache <- file.path(path.expand("~"), ".sccomp_models",
                         as.character(packageVersion("sccomp")))
dir.create(sccompCache, recursive = TRUE, showWarnings = FALSE)
unlockBinding("sccomp_stan_models_cache_dir", asNamespace("sccomp"))
assignInNamespace("sccomp_stan_models_cache_dir", sccompCache, ns = "sccomp")

#----- Factor sample and tomato
seurat$sample <- factor(seurat$sample, levels = c("3", "4", "1", "2"))
seurat$tomato <- factor(seurat$tomato, levels = c("FALSE", "TRUE"))

#----- Run sccomp estimate
sccomp_estimated <-
    seurat |>
    sccomp_estimate(
        formula_composition = ~ tomato,
        .sample = sample,
        .cell_group = celltype.mapped,
        inference_method = "hmc",
        prior_mean = list(intercept = c(0, 1), coefficients = c(0, 3)),
        cores = 1,
        verbose = TRUE,
        output_directory = paste0(outputDir, "sccomp_draws_files")
    )

#----- Test the compostional effect
sccomp_result <- sccomp_test(sccomp_estimated)

#----- Flat table of effects, drop the list-columns so it can be written to csv
sccomp_flat <- sccomp_result |>
    select(where(~ !is.list(.x)))

write.csv(sccomp_flat,
          file = paste0(outputDir, "DA-results.csv"),
          row.names = FALSE)

#----- Extract the compositional effect for our comparison
sccomp_flat <- sccomp_result |>
    filter(parameter == "tomatoTRUE") |>
    select(where(~ !is.list(.x)))

#----- Significant if FDR < 0.05
sccomp_flat$sig <- factor(ifelse(sccomp_flat$c_FDR < 0.05, "Significant", "Not significant"),
                          levels = c("Significant", "Not significant"))


#----- Look at the results
p <- ggplot(sccomp_flat, aes(x = c_effect, y = reorder(celltype.mapped, c_effect))) +
            geom_vline(xintercept = 0, colour = "grey50") +
            geom_vline(xintercept = 0.1, linetype = "dashed", colour = "grey50") +
            geom_vline(xintercept = -0.1, linetype = "dashed", colour = "grey50") +
            geom_errorbarh(aes(xmin = c_lower, xmax = c_upper, colour = sig), height = 0.25) +
            geom_point(aes(fill = sig), shape = 21, colour = "black", size = 2.5, stroke = 0.6) +
            scale_colour_manual(values = c("Significant" = "#D7301F",
                                            "Not significant" = "grey60"),
                                name = "FDR < 0.05", drop = FALSE) +
            scale_fill_manual(values = c("Significant" = "#D7301F",
                                        "Not significant" = "grey60"),
                                name = "FDR < 0.05", drop = FALSE) +
            labs(x = "Composition effect (logit scale)",
                y = NULL) +
            theme_classic(base_size = 12) +
            theme(legend.position = "bottom",
                    panel.grid.minor = element_blank())
ggsave(paste0(figDir, "sccomp_output.png"), width = 8, height = 8)
