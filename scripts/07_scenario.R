#!/usr/bin/env Rscript
# 07_scenario.R
# Reconstruct ancestral diel activity, identify transition branches, and emit:
#   (1) a PCOC scenario string,
#   (2) a HyPhy-labelled tree (foreground branches tagged {Foreground}),
#   (3) a plain list of transition branches for the discordance flagging (step 08),
#   (4) a phylogenetic-signal summary.
#
# Direction note: mammals are ancestrally nocturnal, so the biologically clean
# convergent state is DERIVED diurnality (plus any reversals to nocturnality).
# Read transitions off the ASR; do NOT hand-code from tip states alone.

suppressMessages({
  library(ape); library(phytools); library(corHMM)
})

proj <- Sys.getenv("PROJ", unset = normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")))
res  <- file.path(proj, "results", "scenario")
tree <- read.tree(file.path(proj, "data", "species_tree.nwk"))
diel <- read.csv(file.path(proj, "data", "diel_activity.csv"), stringsAsFactors = FALSE)
# expected columns: species, activity  (activity in {diurnal, nocturnal, cathemeral, crepuscular})

# ---- 1. code the binary trait (keep power up; middle states need a rule) ----
# EDIT this mapping to your decision. Example: diurnal = 1, everything else = 0.
diel$binary <- ifelse(diel$activity == "diurnal", 1L, 0L)
states <- setNames(diel$binary, diel$species)
states <- states[tree$tip.label]                 # align to tree order
stopifnot(!any(is.na(states)))                    # every tip must have a state

# ---- 2. phylogenetic signal (report; expect high for mammal activity) ----
D <- tryCatch(phylo.d.report <- NA, error = function(e) NA)  # placeholder
lambda <- phytools::phylosig(tree, states, method = "lambda", test = TRUE)
cat(sprintf("Pagel's lambda = %.3f (p = %.3g)\n", lambda$lambda, lambda$P))

# ---- 3. ancestral state reconstruction (corHMM, rate-class = 1 to start) ----
dat <- data.frame(sp = tree$tip.label, st = states[tree$tip.label] + 1L)  # corHMM wants 1/2
fit <- corHMM(phy = tree, data = dat, rate.cat = 1, model = "ARD",
              node.states = "marginal")
anc <- fit$states                                 # marginal probs at internal nodes

# ---- 4. identify transition branches (parent state != child state) ----
# For each edge, compare the MAP state of its parent and child nodes.
node_state <- function(probs) apply(probs, 1, which.max)  # 1 or 2
tip_state  <- states[tree$tip.label] + 1L
internal   <- node_state(anc)
all_states <- c(tip_state, internal)              # indexed 1..(Ntip+Nnode)

edges <- tree$edge                                # columns: parent, child
trans <- which(all_states[edges[,1]] != all_states[edges[,2]])
cat(sprintf("Transition branches detected: %d\n", length(trans)))
cat("If this is small (say < 4-5 independent events), PCOC power will be low;\n")
cat("revisit taxon sampling before committing (see pcoc_sim in step 09).\n")

# ---- 5. write outputs ----
dir.create(res, showWarnings = FALSE, recursive = TRUE)

# (3) transition branch list: parent_node child_node child_label(if tip)
tb <- data.frame(edge_id = trans,
                 parent  = edges[trans, 1],
                 child   = edges[trans, 2])
tb$child_label <- ifelse(tb$child <= Ntip(tree), tree$tip.label[tb$child], NA)
write.csv(tb, file.path(res, "transition_branches.csv"), row.names = FALSE)

# (2) HyPhy-labelled tree: tag foreground (transition) branches with {Foreground}
# EDIT: HyPhy branch labelling is by node; use the child node names. A robust way
# is to write the tree, then run `hyphy label-tree` or annotate with your own tags.
labelled <- tree
# Simple approach: name internal nodes, then tag transition children in the newick.
labelled$node.label <- paste0("N", seq_len(Nnode(labelled)))
write.tree(labelled, file.path(res, "species_labelled_base.nwk"))
cat("Foreground = transition branches in transition_branches.csv.\n")
cat("Convert to HyPhy {Foreground} tags with hyphy `label-tree` or a small script.\n")

# (1) PCOC scenario string: PCOC expects branch/node specifications for the
# convergent transitions. Build from the child-node ids of transition branches.
# EDIT to match your pcoc version's scenario syntax (node ids vs. tip clades).
scen <- paste(tb$child, collapse = "/")
writeLines(scen, file.path(res, "pcoc_scenario.txt"))
cat(sprintf("PCOC scenario written (%d transition nodes).\n", nrow(tb)))
