rm(list=ls())
library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(tidytree)
library(ggtreeExtra)
library(ape)

tree_input_path <- "E:/2025.8-2026.7/iScience_revision/phylogenetic_tree/mOTUs.treefile"
tree <- read.tree(tree_input_path)
selected_outgroup <- "Candidatus_Nanosynbacter_sp._HMT-352"
cat("Rooting tree with outgroup: '", selected_outgroup, "'...\n", sep="")

if (selected_outgroup %in% tree$tip.label) {
    rooted_tree <- root(tree, outgroup = selected_outgroup, resolve.root = TRUE)
    cat("Rooting successful!\n")
} else {
    stop(paste("Error: Selected outgroup '", selected_outgroup, "' not found in the tree."))
}

output_nwk <- "E:/2025.8-2026.7/iScience_revision/phylogenetic_tree/rooted_tree.nwk"
write.tree(rooted_tree, file = output_nwk)
cat("Rooted tree saved to: ", output_nwk, "\n")

cat("\n=== Validation Results ===\n")
cat("Original tree rooted: ", is.rooted(tree), "\n")
cat("New tree rooted: ", is.rooted(rooted_tree), "\n")

root_node_id <- length(rooted_tree$tip.label) + 1
cat("Root Node ID: ", root_node_id, "\n")

info_path <- "E:/2025.8-2026.7/iScience_revision/phylogenetic_tree/species100.txt"
info_data <- read.table(info_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
info_data$species <- trimws(info_data$species)

custom_colors <- c(
    "Pseudomonadota" = "#87BEFA",
    "Actinomycetota" = "#FF9999",
    "Spirochaetota"  = "#7dc43a",
    "Bacteroidota"   = "#DE9AEA",
    "Chlamydiota"    = "#FF7F00",
    "Bacillota"      = "#FFE7A3",
    "Residuals"      = "#A6A6A6"
)

# Map phyla to groups
info_data <- info_data %>%
    mutate(phylum_group = ifelse(phylum %in% names(custom_colors), phylum, "Residuals"))

# Group the tree
phylum_group_list <- split(info_data$species, info_data$phylum_group)
tree_grouped <- groupOTU(rooted_tree, phylum_group_list)

max_x <- max(fortify(tree_grouped)$x, na.rm = TRUE)
alignment_x <- max_x * 1.15 

ggtree(tree_grouped, layout = "circular", aes(color = group), size = 0.5) +
    scale_color_manual(
        values = custom_colors,
        breaks = names(custom_colors),
        guide = guide_legend(title = "Phylum Classification")
    ) +
    geom_segment(
        aes(x = x, xend = alignment_x * 0.98, y = y, yend = y, color = group),
        data = td_filter(isTip),
        linetype = "dotted", size = 0.2, alpha = 0.5, show.legend = FALSE
    ) +
    geom_text(
        aes(
            x = alignment_x, 
            y = y, 
            label = label,
            color = group,
            angle = ifelse(angle > 90 & angle < 270, angle + 180, angle),
            hjust = ifelse(angle > 90 & angle < 270, 1, 0)
        ),
        data = td_filter(isTip),
        size = 2.5,
        show.legend = FALSE
    ) +
    xlim(0, max_x * 2.2) + 
    theme(
        legend.position = "right",
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
    )

species_labels <- data.frame(
    tip_label = rooted_tree$tip.label,
    species_name = gsub("_", " ", rooted_tree$tip.label),
    tip_number = 1:length(rooted_tree$tip.label)
)

