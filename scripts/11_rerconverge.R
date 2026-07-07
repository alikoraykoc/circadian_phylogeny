#!/usr/bin/env Rscript
# 11_rerconverge.R
# Gene-level test: which of the 18 genes has relative evolutionary rates
# associated with the diel trait. Consumes the collected fixed-topology trees
# (step 05, uniform model Q.MAMMAL+F+R6) and the binary trait (step 07).

suppressMessages({
  library(RERconverge)
  library(ape)
})

proj <- Sys.getenv("PROJ", unset = getwd())
res  <- file.path(proj, "results", "rerconverge")
dir.create(res, showWarnings = FALSE, recursive = TRUE)

# rer_input.trees is "<gene>\t<newick>" per line, built in step 05.
trees <- readTrees(file.path(proj, "results", "branchlengths", "rer_input.trees"))

# relative evolutionary rates across all genes
RERmat <- getAllResiduals(trees, transform = "sqrt")

# ---- foreground: diurnal tips (the derived convergent state) ----
diel <- read.csv(file.path(proj, "data", "diel_activity.csv"))
fg_tips <- diel$species[diel$activity == "diurnal"]

# Use "terminal" for tip branches only (conservative);
# "all" would also include ancestral branches leading to diurnal clades.
fgTree <- foreground2Tree(fg_tips, trees, clade = "terminal")
paths  <- tree2Paths(fgTree, trees)

cor <- correlateWithBinaryPhenotype(RERmat, paths,
                                    min.sp = 10, min.pos = 2)

cor <- cor[order(cor$P), ]
write.csv(cor, file.path(res, "rer_gene_association.csv"))
cat("Top phenotype-associated genes:\n")
print(head(cor, 10))
