#!/usr/bin/env Rscript
# 07_scenario.R
# Reconstruct ancestral diel activity and emit, for each ASR model:
#   (1) node states (nocturnal/diurnal) for every node of the Upham tree,
#   (2) transition branches WITH DIRECTION (gain of diurnality vs reversal),
#   (3) convergent event membership: for each transition, the contiguous clade
#       of branches that remain in the derived state (this is what PCOC needs),
#   (4) a phylogenetic-signal summary.
#
# Direction note: mammals are ancestrally nocturnal, so the biologically clean
# convergent state is DERIVED diurnality; reversals to nocturnality are analysed
# as a SEPARATE convergent class, never merged into the same class. PCOC fits a
# single convergent amino-acid profile per run, so mixing branches that move
# toward opposite phenotypes into one class is incoherent and self-cancelling.
#
# Two models are fitted:
#   ER  (primary)     equal rates, root free. Best AIC; recovers a nocturnal
#                     placental ancestor without us imposing it.
#   ARD (sensitivity) all-rates-different, root free. Reported so the manuscript
#                     can show the convergence result does not hinge on the ASR
#                     model choice.
#
# States are coded 1 = nocturnal (ancestral), 2 = diurnal (derived).

suppressMessages({
  library(ape); library(phytools); library(corHMM)
})

proj <- Sys.getenv("PROJ", unset = getwd())
res  <- file.path(proj, "results", "scenario")
dir.create(res, showWarnings = FALSE, recursive = TRUE)

tree <- read.tree(file.path(proj, "data", "species_tree.nwk"))
diel <- read.csv(file.path(proj, "data", "diel_activity.csv"), stringsAsFactors = FALSE)

STATE_LAB <- c("nocturnal", "diurnal")
NOCT <- 1L
DIUR <- 2L

# ---- 1. code the binary trait ----
# Anything that is not scored diurnal is treated as nocturnal (the ancestral
# state). With the current data all 60 tips are strictly diurnal or nocturnal,
# so no intermediate states are being collapsed here.
diel$binary <- ifelse(diel$activity == "diurnal", 1L, 0L)
states <- setNames(diel$binary, diel$species)
states <- states[tree$tip.label]
stopifnot(!any(is.na(states)))

n_diur <- sum(states == 1L)
n_noct <- sum(states == 0L)
cat(sprintf("Tip states: %d diurnal, %d nocturnal (of %d)\n", n_diur, n_noct, Ntip(tree)))
if (min(n_diur, n_noct) / Ntip(tree) > 0.4) {
  cat("NOTE: the diel sample is close to balanced. Real mammals are predominantly\n")
  cat("nocturnal, so a balanced sample can pull the ASR toward a diurnal ancestor.\n")
  cat("This is why the ER/ARD comparison below is reported rather than assumed.\n")
}

# ---- 2. phylogenetic signal ----
lambda <- phytools::phylosig(tree, states, method = "lambda", test = TRUE)
cat(sprintf("Pagel's lambda = %.3f (p = %.3g)\n\n", lambda$lambda, lambda$P))

# ---- 3. helpers ----
tip_state <- states[tree$tip.label] + 1L
edges <- tree$edge

children_of <- function(node) edges[edges[, 1] == node, 2]

# Tips subtended by a node (the universal key for matching this node across
# trees with different internal numbering: PCOC, IQ-TREE, ape all differ).
tips_of <- function(node) {
  if (node <= Ntip(tree)) tree$tip.label[node] else extract.clade(tree, node)$tip.label
}

# The convergent event for a transition onto child node `c`: node c itself
# (the transition branch) plus every descendant branch that REMAINS in the
# derived state. Descent stops at any node that reverts, so a nested reversal
# is correctly excluded from the convergent group rather than being swept in.
convergent_clade <- function(c, all_states) {
  derived <- all_states[c]
  keep <- c(c)
  stack <- c(c)
  while (length(stack) > 0) {
    n <- stack[length(stack)]; stack <- stack[-length(stack)]
    for (ch in children_of(n)) {
      if (all_states[ch] == derived) {
        keep <- c(keep, ch)
        stack <- c(stack, ch)
      }
    }
  }
  keep  # first element is always the transition node, which PCOC requires
}

# ---- 4. fit both models and write outputs ----
fit_and_write <- function(model, tag) {
  set.seed(1)
  dat <- data.frame(sp = tree$tip.label, st = tip_state)
  fit <- corHMM(phy = tree, data = dat, rate.cat = 1, model = model,
                node.states = "marginal")

  internal <- apply(fit$states, 1, which.max)
  all_states <- c(tip_state, internal)  # indexed 1..(Ntip+Nnode)

  root <- Ntip(tree) + 1L
  cat(sprintf("== %s (%s) ==  AIC = %.2f\n", tag, model, fit$AIC))
  cat(sprintf("   root state: %s (P_noct = %.3f, P_diur = %.3f)\n",
              STATE_LAB[internal[1]], fit$states[1, 1], fit$states[1, 2]))

  trans <- which(all_states[edges[, 1]] != all_states[edges[, 2]])

  # node state table (every node, so downstream tools need not refit)
  ns <- data.frame(
    node  = seq_len(Ntip(tree) + Nnode(tree)),
    state = STATE_LAB[all_states],
    is_tip = seq_len(Ntip(tree) + Nnode(tree)) <= Ntip(tree),
    label = c(tree$tip.label, rep(NA, Nnode(tree)))
  )
  write.csv(ns, file.path(res, sprintf("node_states_%s.csv", tag)), row.names = FALSE)

  # Same states, keyed by DESCENDANT TIP SET rather than by node id. ape,
  # ete3/IQ-TREE and PCOC each number nodes differently, and a node id is
  # meaningless once a tree is pruned to a gene's taxa. The tip set is the one
  # key that survives both renumbering and pruning, so every downstream consumer
  # (HyPhy foreground labelling, TDG09 node labels, step 08) should join on this
  # rather than re-deriving states from the transition list.
  nts <- data.frame(
    node   = ns$node,
    state  = ns$state,
    is_tip = ns$is_tip,
    tips   = vapply(ns$node, function(n) paste(sort(tips_of(n)), collapse = ";"), ""),
    stringsAsFactors = FALSE
  )
  write.csv(nts, file.path(res, sprintf("node_tipsets_%s.csv", tag)), row.names = FALSE)

  # transition table with direction
  tb <- data.frame(
    edge_id = trans,
    parent  = edges[trans, 1],
    child   = edges[trans, 2],
    from    = STATE_LAB[all_states[edges[trans, 1]]],
    to      = STATE_LAB[all_states[edges[trans, 2]]]
  )
  tb$direction <- ifelse(tb$from == "nocturnal" & tb$to == "diurnal", "gain",
                  ifelse(tb$from == "diurnal" & tb$to == "nocturnal", "reversal", NA))
  tb$child_label <- ifelse(tb$child <= Ntip(tree), tree$tip.label[tb$child], NA)
  tb$n_desc_tips <- sapply(tb$child, function(c) length(tips_of(c)))

  write.csv(tb, file.path(res, sprintf("transition_branches_%s.csv", tag)), row.names = FALSE)

  # convergent event membership, one row per (event, node)
  ev_rows <- list()
  eid <- 0L
  for (i in seq_len(nrow(tb))) {
    eid <- eid + 1L
    nodes <- convergent_clade(tb$child[i], all_states)
    for (j in seq_along(nodes)) {
      n <- nodes[j]
      ev_rows[[length(ev_rows) + 1L]] <- data.frame(
        event_id      = eid,
        direction     = tb$direction[i],
        transition_edge = tb$edge_id[i],
        node          = n,
        node_role     = if (j == 1L) "transition" else "convergent",
        tips          = paste(sort(tips_of(n)), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
  ev <- do.call(rbind, ev_rows)
  write.csv(ev, file.path(res, sprintf("convergent_events_%s.csv", tag)), row.names = FALSE)

  n_gain <- sum(tb$direction == "gain")
  n_rev  <- sum(tb$direction == "reversal")
  gain_branches <- sum(ev$direction == "gain")
  rev_branches  <- sum(ev$direction == "reversal")
  cat(sprintf("   gains of diurnality:      %2d events, %3d convergent branches\n",
              n_gain, gain_branches))
  cat(sprintf("   reversals to nocturnality:%2d events, %3d convergent branches\n",
              n_rev, rev_branches))

  # PCOC power floor (see CLAUDE.md): fewer than ~4-5 independent events means
  # low power, which is a signal to revisit taxon sampling, not to push on.
  for (d in c("gain", "reversal")) {
    ne <- sum(tb$direction == d)
    if (ne < 4) {
      cat(sprintf("   WARNING: only %d independent %s events; PCOC power will be low.\n", ne, d))
    }
  }
  # A convergent class covering most of the BRANCHES leaves PCOC no ancestral
  # contrast to work against, which makes the run uninformative rather than
  # negative. Branches, not tips: an internal convergent node subtends tips that
  # may have reverted, so a tip-based count overstates the convergent fraction.
  nbranch <- nrow(edges)
  for (d in c("gain", "reversal")) {
    nb <- sum(ev$direction == d)
    frac <- nb / nbranch
    cat(sprintf("   %-9s class: %3d / %d branches convergent (%.0f%%)\n",
                d, nb, nbranch, 100 * frac))
    if (frac > 0.6) {
      cat(sprintf("   WARNING: %s class covers %.0f%% of branches; little ancestral\n",
                  d, 100 * frac))
      cat("            contrast remains, so PCOC would be uninformative, not negative.\n")
    }
  }
  cat("\n")
  invisible(tb)
}

er  <- fit_and_write("ER",  "ER")
ard <- fit_and_write("ARD", "ARD")

# Primary analysis is ER. Keep the historical filename pointing at the primary
# model so steps 08 and prep_tdg09_inputs.py keep working unchanged; they read
# edge_id/parent/child/child_label, which are still present.
file.copy(file.path(res, "transition_branches_ER.csv"),
          file.path(res, "transition_branches.csv"), overwrite = TRUE)

# HyPhy-labelled base tree (foreground tagging happens in step 12)
labelled <- tree
labelled$node.label <- paste0("N", seq_len(Nnode(labelled)))
write.tree(labelled, file.path(res, "species_labelled_base.nwk"))

cat("Primary model: ER. transition_branches.csv = transition_branches_ER.csv\n")
cat("Outputs per model: node_states_*.csv, transition_branches_*.csv, convergent_events_*.csv\n")
