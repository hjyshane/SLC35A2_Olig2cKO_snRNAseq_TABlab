## Integration with Harmony
#

# Read previous object if needed.
# processed_list <- qs::qread(file.path(qsave_dir, "04_processed_rna_list.qs"))

# using Harmony library
integrated <- RunHarmony(integrated, "cov_har")

# Integration with Seruat
merged <- merge(processed_list[[1]], y = processed_list[-1])

DefaultAssay(merged) <- "RNA"

merged <- IntegrateLayers(
    object = merged, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = FALSE
)
