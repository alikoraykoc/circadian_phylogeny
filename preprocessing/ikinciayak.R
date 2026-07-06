library(ape)
library(phangorn)
library(phytools)

base_dir <- "C:/Users/alper/Downloads/circadian_upham"

upham <- read.tree(file.path(base_dir, "circadian_60sp_upham2019_mcc.nwk"))

timetree <- read.tree(file.path(base_dir, "circadian_60sp_renamed.nwk"))

timetree$tip.label <- gsub(" ", "_", timetree$tip.label)
upham$tip.label    <- gsub(" ", "_", upham$tip.label)

cat("TimeTree'de olup Upham'da olmayan:\n")
print(setdiff(timetree$tip.label, upham$tip.label))

cat("Upham'da olup TimeTree'de olmayan:\n")
print(setdiff(upham$tip.label, timetree$tip.label))

rf <- RF.dist(unroot(timetree), unroot(upham), normalize = TRUE)

cat("\nNormalize RF mesafesi:", round(rf, 3), "\n")
cat("(0 = ayn?? topoloji, 1 = tamamen farkl??)\n")

assoc <- cbind(timetree$tip.label, timetree$tip.label)

pdf(file.path(base_dir, "tree_comparison.pdf"), width = 14, height = 14)

cophyloplot(
  timetree, upham,
  assoc = assoc,
  length.line = 4,
  space = 30,
  gap = 3,
  fsize = 0.6,
  lwd = 1
)

dev.off()

cat("Kar????la??t??rma g??rseli kaydedildi: tree_comparison.pdf\n")