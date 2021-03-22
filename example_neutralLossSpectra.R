# load required libraries ------------------------------------------------------
library(Spectra)
library(MsBackendMassbank)
library(MsBackendMgf)
library(MsCoreUtils)
library(pheatmap)

# some general settings --------------------------------------------------------
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "red"))(paletteLength)
myBreaks <- c(seq(0,1, by = 1/50))

# neutral loss similiarity example ---------------------------------------------
filePath <- "examples/libraries/pos"

massbank_files <- list.files(filePath, pattern = ".txt",
                             full.names = TRUE)

lipid_spectra <- Spectra(massbank_files,
                         source = MsBackendMassbank(),
                         backend = MsBackendDataFrame())

#filter PE spectra -------------------------------------------------------------
pe_spectra <- lipid_spectra[which(grepl("PE", lipid_spectra$name))]

pe_spectra_comb <- combineSpectra(pe_spectra,
                                  f = pe_spectra$name,
                                  #p = ms2_spectra$CLUSTER_ID,
                                  intensityFun = base::sum,
                                  mzFun = base::mean,
                                  tolerance = 0.005,
                                  ppm = 0,
                                  minProp = 0.5,
                                  peaks = "intersect",
                                  weighted = TRUE)

# calculate NL spectrum (Spectra style) ----------------------------------------
pe_spectra_nl2 <- pe_spectra_comb
pe_spectra_nl2 <- applyProcessing(pe_spectra_nl2)
ms2 <- msLevel(pe_spectra_nl2) == 2L
mz(pe_spectra_nl2@backend)[ms2] <- mz(pe_spectra_nl2)[ms2] - precursorMz(pe_spectra_nl2)[ms2]

# create similarity matrix
similarity <- compareSpectra(pe_spectra_comb, tolerance = 0.005)
colnames(similarity) <-rownames(similarity) <- pe_spectra_comb$name

similarity_nl2 <- compareSpectra(pe_spectra_nl2, tolerance = 0.005)
colnames(similarity_nl2) <- rownames(similarity_nl2) <- pe_spectra_nl2$name

# Plot the heatmap
pheatmap(similarity, color=myColor, breaks=myBreaks, main = "Dotproduct")
pheatmap(similarity_nl2, color=myColor, breaks=myBreaks, main = "NL-Dotproduct ver 2")

plotSpectraMirror(pe_spectra_comb[1], pe_spectra_comb[2])
