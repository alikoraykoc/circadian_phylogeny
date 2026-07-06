# ============================================================
# Upham et al. (2019) Posterior Tree Pruning + MCC
# FINAL versiyon (6 sister species de??i??imi, do??ru isimlerle)
# ============================================================

library(ape)
library(phangorn)
library(phytools)

# ------------------------------------------------------------
# 1. Yollar
# ------------------------------------------------------------
base_dir     <- "C:/Users/alper/Downloads/circadian_upham"
tree_dir     <- file.path(base_dir, "DNAonly_4098sp_topoFree_NDexp")
output_file  <- file.path(base_dir, "circadian_60sp_upham2019_mcc.nwk")

# ------------------------------------------------------------
# 2. A??a?? dosyalar??n?? listele (subsample: 10000 -> 1000)
# ------------------------------------------------------------
tree_files <- list.files(tree_dir, pattern = "\\.tre$", full.names = TRUE)
cat("Toplam .tre dosyas??:", length(tree_files), "\n")

# T??m 10000'i kullanmak istersen bir sonraki sat??r?? yorum sat??r?? yap
tree_files <- tree_files[seq(1, length(tree_files), by = 10)]
cat("Kullan??lacak a??a?? say??s??:", length(tree_files), "\n")

# ------------------------------------------------------------
# 3. Hedef t??r listesi (6 sister species de??i??imi ile)
# ------------------------------------------------------------
focal_species <- c(
  "Pan_troglodytes", "Equus_caballus", "Bos_taurus", "Acinonyx_jubatus",
  "Gorilla_gorilla", "Homo_sapiens", "Ictidomys_tridecemlineatus",
  "Camelus_bactrianus", "Ochotona_princeps", "Arvicanthis_niloticus",
  "Sciurus_vulgaris", "Martes_flavigula", "Hylobates_pileatus",
  "Phacochoerus_africanus", "Rhynchocyon_petersi",
  "Cebus_capucinus",             # ger??ek: Cebus_imitator
  "Callicebus_cupreus",          # ger??ek: Plecturocebus_cupreus
  "Xerus_rutilus", "Nasua_narica", "Cynomys_ludovicianus",
  "Giraffa_camelopardalis",      # ger??ek: Giraffa_tippelskirchi
  "Suricata_suricatta", "Callithrix_jacchus", "Vicugna_pacos",
  "Nanger_dama", "Myrmecobius_fasciatus", "Eira_barbara",
  "Aepyceros_melampus", "Speothos_venaticus", "Octodon_degus",
  "Nycticebus_coucang",
  "Leopardus_guigna",            # ger??ek: Leopardus_geoffroyi
  "Aotus_nancymaae", "Acomys_cahirinus", "Microcebus_murinus",
  "Choloepus_hoffmanni", "Pipistrellus_pipistrellus", "Jaculus_jaculus",
  "Paguma_larvata", "Sylvilagus_floridanus", "Macrotis_lagotis",
  "Cynocephalus_volans", "Saccopteryx_bilineata", "Peromyscus_eremicus",
  "Lagorchestes_hirsutus", "Sarcophilus_harrisii", "Rattus_norvegicus",
  "Dasypus_novemcinctus", "Tapirus_terrestris", "Uromys_caudimaculatus",
  "Phalanger_gymnotis",
  "Pseudocheirus_peregrinus",    # ger??ek: Pseudocheirus_occidentalis
  "Puma_concolor", "Panthera_pardus", "Hyaena_hyaena", "Caracal_caracal",
  "Mirza_coquereli",             # ger??ek: Mirza_zaza
  "Tragelaphus_eurycerus", "Daubentonia_madagascariensis", "Orycteropus_afer"
)
cat("Hedef t??r say??s??:", length(focal_species), "\n")

# ------------------------------------------------------------
# 4. Tip e??le??tirme
# ------------------------------------------------------------
sample_tree <- read.tree(tree_files[1])
upham_tips <- sample_tree$tip.label
upham_short <- sapply(strsplit(upham_tips, "_"),
                      function(x) paste(x[1:2], collapse = "_"))

match_idx <- match(focal_species, upham_short)
present_mask <- !is.na(match_idx)
missing_species <- focal_species[!present_mask]

cat("\n--- T??r E??le??tirme ---\n")
cat("Bulunan:", sum(present_mask), "/", length(focal_species), "\n")

if (length(missing_species) > 0) {
  cat("\nEKS??K t??rler:\n")
  for (sp in missing_species) {
    cat("  -", sp, "\n")
    genus <- strsplit(sp, "_")[[1]][1]
    matches <- upham_short[grep(paste0("^", genus, "_"), upham_short)]
    if (length(matches) > 0) {
      cat("    Ayn?? cinsten mevcut:",
          paste(unique(matches), collapse = ", "), "\n")
    }
  }
  stop("Eksik t??rleri ????z.")
}

focal_tips_full <- upham_tips[match_idx]
cat("T??m t??rler bulundu.\n")

# ------------------------------------------------------------
# 5. Pruning
# ------------------------------------------------------------
cat("\nPruning ba??l??yor (", length(tree_files), "a??a??)...\n")
pruned_trees <- vector("list", length(tree_files))
for (i in seq_along(tree_files)) {
  tr <- read.tree(tree_files[i])
  pruned_trees[[i]] <- drop.tip(tr, setdiff(tr$tip.label, focal_tips_full))
  if (i %% 100 == 0) cat(".")
}
cat("\nPruning tamamland??.\n")
class(pruned_trees) <- "multiPhylo"

# ------------------------------------------------------------
# 6. MCC tree
# ------------------------------------------------------------
cat("\nMCC tree hesaplan??yor...\n")
mcc_tree <- maxCladeCred(pruned_trees, tree = TRUE, rooted = TRUE)
cat("MCC tree haz??r.\n")

# ------------------------------------------------------------
# 7. Tip isimlerini Genus_species format??na indir
# ------------------------------------------------------------
mcc_tree$tip.label <- sapply(strsplit(mcc_tree$tip.label, "_"),
                             function(x) paste(x[1:2], collapse = "_"))

# ------------------------------------------------------------
# 8. Sister species ??? ger??ek t??r ismi geri ??evirme
# ------------------------------------------------------------
rename_map <- c(
  "Cebus_capucinus"           = "Cebus_imitator",
  "Callicebus_cupreus"        = "Plecturocebus_cupreus",
  "Pseudocheirus_peregrinus"  = "Pseudocheirus_occidentalis",
  "Mirza_coquereli"           = "Mirza_zaza",
  "Giraffa_camelopardalis"    = "Giraffa_tippelskirchi",
  "Leopardus_guigna"          = "Leopardus_geoffroyi"
)
for (old_name in names(rename_map)) {
  idx <- which(mcc_tree$tip.label == old_name)
  if (length(idx) == 1) {
    mcc_tree$tip.label[idx] <- rename_map[old_name]
    cat("Yeniden adland??r??ld??:", old_name, "->", rename_map[old_name], "\n")
  }
}

# ------------------------------------------------------------
# 9. Kalite kontrol
# ------------------------------------------------------------
cat("\n--- Kalite Kontrol ---\n")
cat("Tip say??s??:", length(mcc_tree$tip.label), "(beklenen: 60)\n")
cat("Ultrametrik mi:", is.ultrametric(mcc_tree), "\n")
cat("K??k var m??:", is.rooted(mcc_tree), "\n")
cat("Binary mi:", is.binary(mcc_tree), "\n")

if (!is.ultrametric(mcc_tree)) {
  mcc_tree <- force.ultrametric(mcc_tree, method = "extend")
  cat("force.ultrametric uyguland??:", is.ultrametric(mcc_tree), "\n")
}

# ------------------------------------------------------------
# 10. Kaydet ve g??rsel
# ------------------------------------------------------------
write.tree(mcc_tree, file = output_file)
cat("\n>> MCC tree kaydedildi:", output_file, "\n")

pdf(file.path(base_dir, "upham2019_mcc_tree.pdf"), width = 10, height = 14)
plot(mcc_tree, cex = 0.7, label.offset = 0.5)
axisPhylo()
title("Upham et al. 2019 - 60 species MCC tree")
dev.off()
cat(">> G??rsel: upham2019_mcc_tree.pdf\n")

cat("\n=== TAMAMLANDI ===\n")
