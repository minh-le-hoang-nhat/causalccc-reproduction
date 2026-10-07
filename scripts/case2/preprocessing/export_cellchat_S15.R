# Export de l'échantillon S15 (= SRR11396173) du jeu CellChat data_humanSkin (case study 2),
# lisible depuis Python (le .rda ne l'est pas).
# Usage : Rscript scripts/case2/preprocessing/export_cellchat_S15.R   (depuis la racine du projet)
suppressMessages(library(Matrix))

e <- new.env()
load("data/case2/cellchat/data_humanSkin_CellChat.rda", envir = e)
d <- e$data_humanSkin$data   # gènes x cellules, valeurs normalisées (log), pas des comptes
m <- e$data_humanSkin$meta

out <- "data/case2/preprocessing/cache/case_study2_S15"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

strip_bc <- function(x) sub("-1$", "", sub("^[^_]*_", "", x))

# méta de toutes les cellules (pour les tableaux croisés)
meta_all <- data.frame(cell = rownames(m), sample = sub("_.*", "", rownames(m)),
                       barcode = strip_bc(rownames(m)), m, row.names = NULL)
write.csv(meta_all, "data/case2/preprocessing/cache/case_study2_meta.csv", row.names = FALSE)

keep <- sub("_.*", "", colnames(d)) == "S15"
ds <- d[, keep]
writeMM(ds, file.path(out, "matrix.mtx"))
writeLines(rownames(ds), file.path(out, "genes.txt"))
writeLines(strip_bc(colnames(ds)), file.path(out, "barcodes.txt"))
write.csv(meta_all[meta_all$sample == "S15", ], file.path(out, "meta.csv"), row.names = FALSE)

cat("S15 :", ncol(ds), "cellules x", nrow(ds), "gènes ->", out, "\n")
print(table(meta_all$labels[meta_all$sample == "S15"]))
