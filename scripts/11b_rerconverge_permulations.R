#!/usr/bin/env Rscript
# 11b_rerconverge_permulations.R
# Permulation null for the RERconverge gene-level test.
#
# Step 11 used correlateWithBinaryPhenotype's parametric p-values. Those are
# known to be miscalibrated under phylogenetic structure: the authors of
# RERconverge published permulations (Saputra, Kowalczyk et al. MBE 2021)
# precisely because parametric p-values deviate from uniform on null phenotypes.
# This script recomputes the association with an empirical null built from
# permulated phenotypes, which preserve the number of foreground species and
# their phylogenetic relationships.
#
# sisters_list is NULL on purpose. With NULL, getForegroundInfoClades builds the
# foreground tree with clade = "terminal", which is exactly what step 11 used.
# Passing a sisters_list would switch it to clade = "all" and the permulation
# would then be the null for a different analysis than the one being tested.
#
# Env: NUMPERMS (default 1000), PERMMODE ("cc" or "ssm", default ssm)
# Out: results/rerconverge/rer_permulation_pvalues_<permmode>.csv

suppressMessages({
  library(RERconverge)
  library(ape)
})

proj <- Sys.getenv("PROJ", unset = getwd())
res  <- file.path(proj, "results", "rerconverge")
dir.create(res, showWarnings = FALSE, recursive = TRUE)

numperms <- as.integer(Sys.getenv("NUMPERMS", unset = "1000"))
permmode <- Sys.getenv("PERMMODE", unset = "ssm")
root_sp  <- Sys.getenv("ROOT_SP",  unset = "Sarcophilus_harrisii")
set.seed(20260902)

trees  <- readTrees(file.path(proj, "results", "branchlengths", "rer_input.trees"))
RERmat <- getAllResiduals(trees, transform = "sqrt")

diel    <- read.csv(file.path(proj, "data", "diel_activity.csv"))
fg_vec  <- diel$species[diel$activity == "diurnal"]

# Guard the interface assumption that has broken this pipeline before: every
# foreground label must exist as a tip in the master tree.
missing <- setdiff(fg_vec, trees$masterTree$tip.label)
if (length(missing) > 0) {
  stop("foreground species absent from master tree: ", paste(missing, collapse = ", "))
}
if (!(root_sp %in% trees$masterTree$tip.label)) {
  stop("root species not in master tree: ", root_sp)
}

cat(sprintf("genes=%d  foreground=%d/%d tips  permmode=%s  numperms=%d\n",
            nrow(RERmat), length(fg_vec), length(trees$masterTree$tip.label),
            permmode, numperms))

# ---- observed association, identical settings to step 11 ----
fgTree <- foreground2Tree(fg_vec, trees, clade = "terminal", plotTree = FALSE)
paths  <- tree2Paths(fgTree, trees)
realcor <- correlateWithBinaryPhenotype(RERmat, paths, min.sp = 10, min.pos = 2)

# ---- permulation null ----
t0 <- Sys.time()
perms <- getPermsBinary(numperms = numperms, fg_vec = fg_vec,
                        sisters_list = NULL, root_sp = root_sp,
                        RERmat = RERmat, trees = trees,
                        mastertree = trees$masterTree,
                        permmode = permmode, method = "k", min.pos = 2)
cat(sprintf("permulations done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# permpvalcor returns a data.frame with columns permpval and permstats, rows
# named by gene, not a bare vector. Index it explicitly.
permP <- permpvalcor(realcor, perms)
stopifnot("permpval" %in% colnames(permP))

out <- data.frame(
  gene          = rownames(realcor),
  Rho           = realcor$Rho,
  N             = realcor$N,
  parametric_P  = realcor$P,
  parametric_p_adj = p.adjust(realcor$P, method = "BH"),
  permulation_P = permP[rownames(realcor), "permpval"],
  stringsAsFactors = FALSE
)
out$permulation_p_adj <- p.adjust(out$permulation_P, method = "BH")
out <- out[order(out$permulation_P), ]

outfile <- file.path(res, sprintf("rer_permulation_pvalues_%s.csv", permmode))
write.csv(out, outfile, row.names = FALSE)
cat("\n"); print(out, row.names = FALSE)
cat(sprintf("\nwrote %s\n", outfile))
