# ============================================================
# SCRIPT 34: Five explanatory programmed cell death schematics
#
# Run after Script 32 and before Script 33.
# These conceptual pathway diagrams are literature-based; they
# are not figures estimated from TCGA, GEO, or cell-line data.
# Sources (primary research):
# Apoptosis: Li et al., Cell 1997, doi:10.1016/S0092-8674(00)80434-1;
#            Juo et al., Curr Biol 1998, doi:10.1016/S0960-9822(07)00420-4.
# Necroptosis: Sun et al., Cell 2012, doi:10.1016/j.cell.2011.11.031.
# Pyroptosis: Shi et al., Nature 2015, doi:10.1038/nature15514.
# Ferroptosis: Dixon et al., Cell 2012, doi:10.1016/j.cell.2012.03.042;
#              Yang et al., Cell 2014, doi:10.1016/j.cell.2013.12.010.
# Cuproptosis: Tsvetkov et al., Science 2022, doi:10.1126/science.abf0529.
# ============================================================

source("03_scripts/00_project_config.R")

out_dir_s34 <- "05_figures/manuscript/individual/pcd_pathways"
dir.create(out_dir_s34, recursive = TRUE, showWarnings = FALSE)

palette_s34 <- c(
  Trigger = "#E8F1FF",
  Signaling = "#E7F4F0",
  Effector = "#FFF0DF",
  Outcome = "#FBE8EA",
  Defense = "#EDEAF9"
)

node_s34 <- function(id, x, y, label, category, width = 2.05) {
  data.frame(
    id = id, x = x, y = y, label = label,
    category = category, width = width,
    stringsAsFactors = FALSE
  )
}

link_s34 <- function(from, to, label = "", dashed = FALSE) {
  data.frame(
    from = from, to = to, label = label, dashed = dashed,
    stringsAsFactors = FALSE
  )
}

draw_pathway_s34 <- function(title, subtitle, nodes, links, footer) {
  box_height <- 0.88
  graphics::par(mar = c(1.0, 0.45, 2.6, 0.45), xpd = NA)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 16), ylim = c(0, 5))
  graphics::title(main = title, cex.main = 1.35, line = 1.2)
  graphics::mtext(subtitle, side = 3, line = 0.05, cex = 0.91)

  for (i in seq_len(nrow(links))) {
    a <- nodes[nodes$id == links$from[i], , drop = FALSE]
    b <- nodes[nodes$id == links$to[i], , drop = FALSE]
    stopifnot(nrow(a) == 1L, nrow(b) == 1L)
    dx <- b$x - a$x
    dy <- b$y - a$y
    stopifnot(dx != 0 || dy != 0)
    clip <- function(w) {
      tx <- if (dx == 0) Inf else (w / 2) / abs(dx)
      ty <- if (dy == 0) Inf else (box_height / 2) / abs(dy)
      min(tx, ty)
    }
    ta <- clip(a$width)
    tb <- clip(b$width)
    graphics::arrows(
      a$x + ta * dx, a$y + ta * dy,
      b$x - tb * dx, b$y - tb * dy,
      length = 0.10, angle = 22, code = 2, lwd = 1.8,
      lty = if (links$dashed[i]) 2 else 1,
      col = if (links$dashed[i]) "#7558A6" else "#4B5969"
    )
    if (nzchar(links$label[i])) {
      graphics::text(
        (a$x + b$x) / 2, (a$y + b$y) / 2 + 0.25,
        links$label[i], cex = 0.74, col = "#684F94"
      )
    }
  }

  for (i in seq_len(nrow(nodes))) {
    n <- nodes[i, ]
    graphics::rect(
      n$x - n$width / 2, n$y - box_height / 2,
      n$x + n$width / 2, n$y + box_height / 2,
      col = palette_s34[[n$category]], border = "#536578", lwd = 1.25
    )
    graphics::text(n$x, n$y, n$label, cex = 0.82,
                   col = "#132537", font = 2)
  }
  graphics::text(0.35, 0.17, footer, adj = c(0, 0.5),
                 cex = 0.69, col = "#4D5B66")
  invisible(NULL)
}

save_pathway_s34 <- function(stem, title, subtitle, nodes, links, footer) {
  base <- file.path(out_dir_s34, stem)
  render <- function() draw_pathway_s34(title, subtitle, nodes, links, footer)

  grDevices::png(paste0(base, ".png"), width = 4200,
                 height = 1600, res = 300, bg = "white")
  tryCatch(render(), finally = grDevices::dev.off())

  grDevices::pdf(paste0(base, ".pdf"), width = 14,
                 height = 5.33, useDingbats = FALSE, bg = "white")
  tryCatch(render(), finally = grDevices::dev.off())

  grDevices::tiff(paste0(base, ".tiff"), width = 4200,
                  height = 1600, res = 300, compression = "lzw",
                  bg = "white")
  tryCatch(render(), finally = grDevices::dev.off())

  files <- paste0(base, c(".png", ".pdf", ".tiff"))
  stopifnot(all(file.exists(files)), all(file.info(files)$size > 0))
  files
}

all_pathway_files_s34 <- character()

# 1. Apoptosis: representative intrinsic and extrinsic routes.
nodes <- do.call(rbind, list(
  node_s34("damage", 1.25, 3.70, "DNA damage /\ncellular stress", "Trigger"),
  node_s34("bax", 3.70, 3.70, "BAX / BAK\nactivation", "Signaling"),
  node_s34("cytc", 6.15, 3.70, "Mitochondrial\ncytochrome c", "Signaling"),
  node_s34("apaf", 8.60, 3.70, "APAF1 apoptosome\n+ caspase-9", "Effector", 2.30),
  node_s34("ligand", 1.25, 1.62, "FASL / TRAIL\nligands", "Trigger"),
  node_s34("receptor", 3.70, 1.62, "Death receptor\nFAS / DR4 / DR5", "Signaling"),
  node_s34("disc", 6.15, 1.62, "FADD / DISC\nassembly", "Signaling"),
  node_s34("c8", 8.60, 1.62, "Caspase-8\nactivation", "Effector"),
  node_s34("exec", 11.35, 2.66, "Executioner\ncaspases-3 / -7", "Effector", 2.42),
  node_s34("death", 14.20, 2.66, "Orderly cellular\ndismantling", "Outcome", 2.40)
))
links <- do.call(rbind, list(
  link_s34("damage", "bax"), link_s34("bax", "cytc"),
  link_s34("cytc", "apaf"), link_s34("apaf", "exec"),
  link_s34("ligand", "receptor"), link_s34("receptor", "disc"),
  link_s34("disc", "c8"), link_s34("c8", "exec"),
  link_s34("exec", "death")
))
all_pathway_files_s34 <- c(all_pathway_files_s34, save_pathway_s34(
  "Figure_PCD_01_Apoptosis", "Apoptosis: intrinsic and extrinsic routes",
  "Two representative initiation routes converge on executioner caspases",
  nodes, links,
  "Conceptual schematic | Li et al., Cell 1997; Juo et al., Curr Biol 1998"
))

# 2. Necroptosis: representative TNF-triggered route, not a universal pathway.
nodes <- do.call(rbind, list(
  node_s34("tnf", 1.30, 2.55, "TNF / TNFR1\nsignaling", "Trigger"),
  node_s34("gate", 3.92, 2.55, "Caspase-8 activity\nrestrained", "Signaling", 2.30),
  node_s34("nec", 6.55, 2.55, "RIPK1 / RIPK3\nnecrosome", "Signaling", 2.25),
  node_s34("mlkl", 9.10, 2.55, "RIPK3 activates\nMLKL", "Effector"),
  node_s34("mem", 11.70, 2.55, "MLKL reaches\nplasma membrane", "Effector", 2.25),
  node_s34("lysis", 14.35, 2.55, "Membrane rupture\nand inflammatory release", "Outcome", 2.65)
))
links <- do.call(rbind, lapply(seq_len(nrow(nodes) - 1L),
  function(i) link_s34(nodes$id[i], nodes$id[i + 1L])))
all_pathway_files_s34 <- c(all_pathway_files_s34, save_pathway_s34(
  "Figure_PCD_02_Necroptosis", "Necroptosis: RIPK3 and MLKL execution",
  "Representative TNF-induced route when caspase-8-mediated apoptosis is restrained",
  nodes, links,
  "Conceptual schematic | Sun et al., Cell 2012"
))

# 3. Pyroptosis: canonical and human noncanonical routes.
nodes <- do.call(rbind, list(
  node_s34("sensor", 1.25, 3.68, "Damage / pathogen\nsensors", "Trigger"),
  node_s34("infl", 3.75, 3.68, "Inflammasome\n+ ASC", "Signaling"),
  node_s34("c1", 6.20, 3.68, "Caspase-1\nactivation", "Effector"),
  node_s34("lps", 1.25, 1.50, "Cytosolic\nLPS", "Trigger"),
  node_s34("c45", 3.75, 1.50, "Human caspase-4 / -5\n(noncanonical)", "Effector", 2.40),
  node_s34("gsd", 8.88, 2.62, "GSDMD cleavage\nN-terminal fragment", "Effector", 2.40),
  node_s34("pores", 11.55, 2.62, "GSDMD plasma\nmembrane pores", "Effector", 2.20),
  node_s34("lysis", 14.32, 2.62, "Cell lysis and\ncytokine release", "Outcome", 2.35)
))
links <- do.call(rbind, list(
  link_s34("sensor", "infl"), link_s34("infl", "c1"),
  link_s34("c1", "gsd"), link_s34("lps", "c45"),
  link_s34("c45", "gsd"), link_s34("gsd", "pores"),
  link_s34("pores", "lysis")
))
all_pathway_files_s34 <- c(all_pathway_files_s34, save_pathway_s34(
  "Figure_PCD_03_Pyroptosis", "Pyroptosis: gasdermin D pore formation",
  "Inflammasome and cytosolic-LPS routes converge on GSDMD cleavage",
  nodes, links,
  "Conceptual schematic | Shi et al., Nature 2015"
))

# 4. Ferroptosis: iron-driven oxidation opposed by GSH/GPX4 defense.
nodes <- do.call(rbind, list(
  node_s34("iron", 1.38, 4.08, "Redox-active\niron", "Trigger"),
  node_s34("pufa", 4.00, 2.95, "PUFA-containing\nmembrane lipids", "Signaling", 2.38),
  node_s34("perox", 7.15, 3.65, "Lipid peroxide\naccumulation", "Effector", 2.42),
  node_s34("failure", 10.45, 3.65, "Membrane\ndamage", "Effector"),
  node_s34("death", 13.72, 3.65, "Ferroptotic\ncell death", "Outcome"),
  node_s34("cystine", 1.38, 1.50, "SLC7A11\ncystine import", "Defense"),
  node_s34("gsh", 4.00, 1.50, "Glutathione\n(GSH)", "Defense"),
  node_s34("gpx", 7.15, 1.50, "GPX4 lipid\nperoxide repair", "Defense", 2.42)
))
links <- do.call(rbind, list(
  link_s34("iron", "perox"), link_s34("pufa", "perox"),
  link_s34("perox", "failure"), link_s34("failure", "death"),
  link_s34("cystine", "gsh"), link_s34("gsh", "gpx"),
  link_s34("gpx", "perox", label = "limits", dashed = TRUE)
))
all_pathway_files_s34 <- c(all_pathway_files_s34, save_pathway_s34(
  "Figure_PCD_04_Ferroptosis", "Ferroptosis: lipid peroxide accumulation",
  "Iron and oxidizable lipids promote damage; GSH/GPX4 supplies a protective defense",
  nodes, links,
  "Conceptual schematic | Dixon et al., Cell 2012; Yang et al., Cell 2014"
))

# 5. Cuproptosis: simplified mitochondrial lipoylated-protein model.
nodes <- do.call(rbind, list(
  node_s34("copper", 1.35, 2.55, "Copper\naccumulation", "Trigger"),
  node_s34("fdx", 3.95, 2.55, "FDX1-dependent\ncopper / lipoylation", "Signaling", 2.48),
  node_s34("tca", 6.65, 2.55, "Lipoylated TCA\nproteins (e.g., DLAT)", "Signaling", 2.55),
  node_s34("agg", 9.35, 2.55, "Protein aggregation\nand Fe-S protein loss", "Effector", 2.55),
  node_s34("stress", 12.02, 2.55, "Mitochondrial\nproteotoxic stress", "Effector", 2.35),
  node_s34("death", 14.45, 2.55, "Cuproptotic\ncell death", "Outcome", 2.12)
))
links <- do.call(rbind, lapply(seq_len(nrow(nodes) - 1L),
  function(i) link_s34(nodes$id[i], nodes$id[i + 1L])))
all_pathway_files_s34 <- c(all_pathway_files_s34, save_pathway_s34(
  "Figure_PCD_05_Cuproptosis", "Cuproptosis: copper-linked proteotoxicity",
  "A simplified mitochondrial pathway centered on lipoylated TCA-cycle proteins",
  nodes, links,
  "Conceptual schematic | Tsvetkov et al., Science 2022"
))

# Append the diagrams to the canonical figure manifest created by Script 32.
# Running Script 32 again overwrites that manifest; rerun Script 34 afterwards.
manifest_path_s34 <-
  "05_figures/manuscript/Publication_figure_manifest.csv"
stopifnot(file.exists(manifest_path_s34), length(all_pathway_files_s34) == 15L)
manifest_s34 <- read.csv(manifest_path_s34, stringsAsFactors = FALSE,
                         check.names = FALSE)
stopifnot(all(c("file", "category", "exists", "size_bytes") %in%
              names(manifest_s34)))
manifest_s34 <- manifest_s34[
  !manifest_s34$file %in% all_pathway_files_s34, , drop = FALSE
]
new_rows_s34 <- data.frame(
  file = all_pathway_files_s34,
  category = "Individual",
  exists = file.exists(all_pathway_files_s34),
  size_bytes = file.info(all_pathway_files_s34)$size,
  stringsAsFactors = FALSE
)
manifest_s34 <- rbind(manifest_s34, new_rows_s34)
write.csv(manifest_s34, manifest_path_s34, row.names = FALSE)
cat("Saved five pathway diagrams in PNG, PDF, TIFF at:", out_dir_s34,
    "\nCanonical individual figures:",
    sum(manifest_s34$category == "Individual") / 3L, "\n")
