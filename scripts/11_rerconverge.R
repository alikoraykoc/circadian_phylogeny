#!/usr/bin/env Rscript
# 11_rerconverge.R
# Gene-level test: which of the 18 genes has relative evolutionary rates
# associated with the diel trait. Consumes the collected fixed-topology trees
# (04, one shared model) and the binary trait (07).

suppressMessages(library(RERconverge))   # install from GitHub: nclark-lab/RERconverge

proj <- Sys.getenv("PROJ")
res  <- file.path(proj, "results", "rerconverge")
dir.create(res, showWarnings = FALSE, recursive = TRUE)

# rer_input.trees is "<gene>\t<newick>" per line, built in step 04.
trees <- readTrees(file.path(proj, "results", "branchlengths", "rer_input.trees"))

# relative evolutionary rates across all genes
RERmat <- getAllResiduals(trees, transform = "sqrt", weighted = TRUE, scale = TRUE)

# ---- foreground: the diurnal (or reversal) tips/branches ----
# EDIT: list the foreground tip labels (independently derived state).
fg_tips <- c()  # e.g. c("Macaca_mulatta","Rattus_norvegicus_diurnal_relative", ...)
stopifnot(length(fg_tips) > 0)

fgTree <- foreground2Tree(fg_tips, trees, clade = "terminal")   # or "all"/"ancestral"
paths  <- tree2Paths(fgTree, trees)

cor <- correlateWithBinaryPhenotype(RERmat, paths,
                                    min.sp = 10, min.pos = 2)

cor <- cor[order(cor$P), ]
write.csv(cor, file.path(res, "rer_gene_association.csv"))
cat("Top phenotype-associated genes:\n")
print(head(cor, 10))
