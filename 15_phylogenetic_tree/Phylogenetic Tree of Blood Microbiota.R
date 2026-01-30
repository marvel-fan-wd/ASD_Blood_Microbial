####------------------------------------------------------- (13).phylogenetic tree -------------------------------------------------------
rm(list=ls())
library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)
library(tidytree)
library(ggtreeExtra)
library(ape)

# ================= 1. Rooting the Tree =================

# Load the unrooted tree file
tree_input_path <- "E:/2025.8-2026.7/iScience_revision/15_phylogenetic_tree/full_analysis100.treefile"
tree <- read.tree(tree_input_path)

# Set the chosen outgroup
selected_outgroup <- "Candidatus_Nanosynbacter_sp._HMT-352"

cat("Rooting tree with outgroup: '", selected_outgroup, "'...\n", sep="")

if (selected_outgroup %in% tree$tip.label) {
    rooted_tree <- root(tree, outgroup = selected_outgroup, resolve.root = TRUE)
    cat("Rooting successful!\n")
} else {
    stop(paste("Error: Selected outgroup '", selected_outgroup, "' not found in the tree."))
}

# Save the rooted tree
output_nwk <- "E:/2025.8-2026.7/iScience_revision/15_phylogenetic_tree/rooted_tree_with_outgroup_Candidatus_Nanosynbacter_sp._HMT-352.nwk"
write.tree(rooted_tree, file = output_nwk)
cat("Rooted tree saved to: ", output_nwk, "\n")

# Basic Validation
cat("\n=== Validation Results ===\n")
cat("Original tree rooted: ", is.rooted(tree), "\n")
cat("New tree rooted: ", is.rooted(rooted_tree), "\n")

# FIXED: Standard way to get root node ID in 'ape'
root_node_id <- length(rooted_tree$tip.label) + 1
cat("Root Node ID: ", root_node_id, "\n")

# ================= 2. Data Preparation & Custom Coloring =================

info_path <- "E:/2025.8-2026.7/iScience_revision/15_phylogenetic_tree/phylum9_species100.txt"
info_data <- read.table(info_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
info_data$species <- trimws(info_data$species)

# Define custom color mapping
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

# ================= 3. Visualization =================

max_x <- max(fortify(tree_grouped)$x, na.rm = TRUE)
alignment_x <- max_x * 1.15 

p2 <- ggtree(tree_grouped, layout = "circular", aes(color = group), size = 0.5) +
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
    ) +
    labs(title = "Phylogenetic Tree of Blood Microbiota")

# Save the plot
pdf_output <- "E:/2025.8-2026.7/iScience_revision/15_phylogenetic_tree/Blood_Microbiota_Phylogenetic_Tree.pdf"
ggsave(pdf_output, p2, width = 12, height = 11)
print(p2)

# ================= 4. Generating Metadata Files =================

# 1. Basic Label Info
species_labels <- data.frame(
    tip_label = rooted_tree$tip.label,
    species_name = gsub("_", " ", rooted_tree$tip.label),
    tip_number = 1:length(rooted_tree$tip.label)
)
write.table(species_labels, "E:/2025.8-2026.7/iScience_revision/15_phylogenetic_tree/species_tip_labels_info.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 2. Complete Metadata
complete_species_info <- merge(species_labels, info_data, by.x = "tip_label", by.y = "species", all.x = TRUE)
write.table(complete_species_info, "E:/2025.8-2026.7/iScience_revision/15_phylogenetic_tree/complete_species_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)
