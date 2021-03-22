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

# ceramide examples ------------------------------------------------------------
# read spectra
ceramide_file <- "examples/ceramides/AutoMSn_CelegansLipids.mgf"
ceramide_spectra <- Spectra(ceramide_file,
                            source = MsBackendMgf(),
                            backend = MsBackendDataFrame())

# read reference spectra
sphingo_files <- list.files("examples/ceramides/MS2_lib", full.names = TRUE)
sphingo_lib <- Spectra(sphingo_files,
                       source = MsBackendMgf(),
                       backend = MsBackendDataFrame())

# overview
plot(rtime(ceramide_spectra),
     precursorMz(ceramide_spectra),
     xlim = c(0, 3000))

# search for water neutral loss
ceramide_spectra$nl_18 <- containsNeutralLoss(ceramide_spectra,
                                              18.010565,
                                              tolerance = 0.005)

table(ceramide_spectra$nl_18)

# search for sphingoid base fragments from C17 iso sphingoid
ceramide_spectra$sph_c17_0 <- containsMz(ceramide_spectra,
                                         c(238.2530,
                                           250.2529,
                                           268.2635),
                                         tolerance = 0.005)

table(ceramide_spectra$sph_c17_0)
table(ceramide_spectra$sph_c17_0 & ceramide_spectra$nl_18)

# search for sphingoid base fragments from C17 iso sphingoid
ceramide_spectra$mz_184 <- containsMz(ceramide_spectra,
                                      c(184.0733),
                                      tolerance = 0.005)

# search for odd precursor mz
ceramide_spectra$odd_precursorMz <- as.integer(precursorMz(ceramide_spectra)) %% 2 == 1

# filter spectra potentially containing a c17 iso sphingoid
ceramide_spectra_filter1 <- ceramide_spectra[ceramide_spectra$sph_c17_0 & ceramide_spectra$nl_18]
ceramide_spectra_filter2 <- ceramide_spectra[ceramide_spectra$mz_184 & ceramide_spectra$odd_precursorMz]

# plot potential candidates
points(rtime(ceramide_spectra_filter1),
       precursorMz(ceramide_spectra_filter1),
       col = "red",
       cex = 2)

points(rtime(ceramide_spectra_filter2),
       precursorMz(ceramide_spectra_filter2),
       col = "green",
       cex = 2)

# create similarity plot
similarity <- compareSpectra(ceramide_spectra_filter1, tolerance = 0.005)
colnames(similarity) <-rownames(similarity) <- round(precursorMz(ceramide_spectra_filter1), 4)
pheatmap(similarity, color=myColor, breaks=myBreaks, main = "Dotproduct")

# search for specific spectra matching Cer(d17:1/22:0(OH))
ceramide_spectra_filter3 <- ceramide_spectra_filter1[abs(precursorMz(ceramide_spectra_filter1) - precursorMz(sphingo_lib[12])) < 0.005]

compareSpectra(ceramide_spectra_filter2, sphingo_lib[12])

plotSpectraMirror(ceramide_spectra_filter1[1], sphingo_lib[12])
plotSpectraMirror(ceramide_spectra_filter1[2], sphingo_lib[12])
