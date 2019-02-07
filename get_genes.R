
# Get gene choices for all samples

seurat_obj <- list.files("data/seurat", pattern = "*.seurat_small.Rda", full.names = TRUE)

for (i in seq_along(seurat_obj)) {
  
  load(seurat_obj[i])
  
}

samples <- list(s1_id, s2_id, s3_id, s4_id, s5_id)

names(samples) <- c("s1_id",
                    "s2_id",
                    "s3_id",
                    "s4_id",
                    "s5_id")

genes <- purrr::map(samples, ~ sort(rownames(.x@data)))
names(genes) <- names(samples)

save(genes, file = "data/seurat_genes.Rda")
