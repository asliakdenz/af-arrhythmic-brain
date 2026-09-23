# ==============================================================================
# 08b_circular_connectivity_figure.R
#
# Figure 1: circular connectogram of the NBS-identified hyperconnectivity
# subnetwork (AF > Controls, t = 4.5), nodes grouped by Schaefer 7-network
# assignment plus the subcortical block, chord colour scaled to the
# t-statistic.
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Network-based statistics -> Functional NBS and hub analysis
# Figure:  Main-text Figure 1
#
# Inputs:
#   - Edge list CSV from 04a_extract_edges.m for the primary NBS result
#   - Node label table from 04e_harvard_oxford_labels.py (RegionName,
#     Top_Overlap_HO_Label); if absent, Schaefer parcel names are used
#
# Output:
#   - SVG circular plot in OUTPUT_DIR
# ==============================================================================

suppressPackageStartupMessages({ library(dplyr); library(circlize) })

# ------------------------------------------------------------------------------
# Configuration  --- edit paths for your environment
# ------------------------------------------------------------------------------
EDGE_FILE <- paste0("/path/to/nbs/02_results/fc/design_1_total/",
                    "contrast_01_af_gt_control/T4.5/Extent/",
                    "nbs_result_T4.5_Extent_significant_edges.csv")
LABEL_FILE <- paste0("/path/to/nbs/02_results/fc/design_1_total/",
                     "contrast_01_af_gt_control/T4.5/Extent/",
                     "hubs_harvard_oxford.csv")          # from 04e
OUTPUT_DIR <- "/path/to/derivatives/figures/circular/"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
edges <- read.csv(EDGE_FILE)
if (!"Region1" %in% names(edges)) edges <- read.delim(EDGE_FILE)
hubs <- if (file.exists(LABEL_FILE)) read.csv(LABEL_FILE) else data.frame(RegionName = character(0))
if (!"Top_Overlap_HO_Label" %in% names(hubs)) {
  message("No Top_Overlap_HO_Label column available; using Schaefer parcel names.")
  hubs$Top_Overlap_HO_Label <- rep(NA_character_, nrow(hubs))
}

# ==============================================================================
# 2. LABEL PARSING
# ==============================================================================
get_network <- function(region) {
  sapply(region, function(x) {
    if(grepl("Vis", x)) return("Vis")
    if(grepl("SomMot", x)) return("SomMot")
    if(grepl("DorsAttn", x)) return("DorsAttn")
    if(grepl("SalVentAttn", x)) return("SalVentAttn")
    if(grepl("Limbic", x)) return("Limbic")
    if(grepl("Cont", x)) return("Cont")
    if(grepl("Default", x)) return("Default")
    return("Subcort") 
  })
}

get_hemi <- function(region) {
  sapply(region, function(x) {
    if(grepl("LH|lh", x)) return("L")
    if(grepl("RH|rh", x)) return("R")
    return("Unknown")
  })
}

get_anatomical_label <- function(region) {
  sapply(region, function(x) {
    
    # Subcortical
    if(grepl("-lh|-rh", x)) {
      raw_anat <- gsub("-lh|-rh", "", x)
      anat <- switch(raw_anat,
                     "HIP"  = "Hippocampus",
                     "AMY"  = "Amygdala",
                     "pTHA" = "Post.Thalamus",
                     "aTHA" = "Ant.Thalamus",
                     "NAc"  = "Accumbens",
                     "GP"   = "Pallidum",
                     "PUT"  = "Putamen",
                     "CAU"  = "Caudate",
                     raw_anat)
      hemi <- if(grepl("lh", x)) ".L" else ".R"
      return(paste0(anat, hemi))
    }
    
    # Cortical
    match_row <- hubs[hubs$RegionName == x, , drop = FALSE]
    if (nrow(match_row) > 0 && isTRUE(!is.na(match_row$Top_Overlap_HO_Label[1]))) {
      anat <- match_row$Top_Overlap_HO_Label[1]
    } else {
      anat <- gsub("7Networks_[LR]H_[A-Za-z]+_", "", x)
    }
    hemi <- if(grepl("LH", x)) ".L" else if(grepl("RH", x)) ".R" else ""
    
    # Aggressive abbreviations
    anat <- gsub("\\(formerly Supplementary Motor Cortex\\)", "", anat)
    anat <- gsub("\\(includes H1 and H2\\)", "", anat)
    anat <- gsub("Juxtapositional Lobule Cortex", "SMA", anat)
    anat <- gsub("Heschl's", "Heschls", anat)
    anat <- gsub("Frontal Orbital", "OFC", anat)
    anat <- gsub("Planum Temporale", "Planum.Temp", anat)
    anat <- gsub("Planum Polare", "Planum.Polare", anat)
    anat <- gsub(", anterior division", ".Ant", anat)
    anat <- gsub(", posterior division", ".Post", anat)
    anat <- gsub(", superior division", ".Sup", anat)
    anat <- gsub(", inferior division", ".Inf", anat)
    anat <- gsub(", pars opercularis", ".ParsOper", anat)
    anat <- gsub(", pars triangularis", ".ParsTri", anat)
    anat <- gsub("Superior ", "Sup.", anat)
    anat <- gsub("Middle ", "Mid.", anat)
    anat <- gsub("Inferior ", "Inf.", anat)
    anat <- gsub("Anterior ", "Ant.", anat)
    anat <- gsub("Posterior ", "Post.", anat)
    anat <- gsub("Lateral ", "Lat.", anat)
    anat <- gsub("Cingulate ", "Cingu.", anat)
    anat <- gsub("Paracingulate ", "Paracingu.", anat)
    anat <- gsub("Precentral ", "Precent.", anat)
    anat <- gsub("Postcentral ", "Postcent.", anat)
    anat <- gsub("Precuneous ", "Precun.", anat)
    anat <- gsub("Supramarginal ", "Supramarg.", anat)
    anat <- gsub("Parahippocampal ", "Parahipp.", anat)
    anat <- gsub("Temporal ", "Temp.", anat)
    anat <- gsub("Occipital ", "Occip.", anat)
    anat <- gsub("Parietal ", "Pariet.", anat)
    anat <- gsub("Opercular ", "Operc.", anat)
    anat <- gsub("Central ", "Cent.", anat)
    anat <- gsub("Frontal ", "Front.", anat)
    anat <- gsub(" Gyrus|Gyrus", ".Gyr", anat)
    anat <- gsub(" Cortex|Cortex", ".Cort", anat)
    anat <- gsub(" Lobule|Lobule", ".Lob", anat)
    anat <- gsub(" ", "", anat) 
    anat <- gsub("\\.\\.", "\\.", anat) 
    
    return(paste0(anat, hemi))
  })
}

all_nodes <- unique(c(edges$Region1, edges$Region2))

node_df <- data.frame(Raw_Region = all_nodes) %>%
  mutate(
    Network = get_network(Raw_Region),
    Hemi    = get_hemi(Raw_Region),
    Label   = get_anatomical_label(Raw_Region)
  ) %>%
  arrange(Network, desc(Hemi), Label) %>%
  mutate(Label = make.unique(Label, sep = "."))

links <- edges %>%
  select(Region1, Region2, t_value) %>%
  left_join(node_df, by=c("Region1"="Raw_Region")) %>% rename(from = Label) %>%
  left_join(node_df, by=c("Region2"="Raw_Region")) %>% rename(to = Label) %>%
  select(from, to, t_value)

# ==============================================================================
# 3. COLOUR & GAP CONFIGURATION
# ==============================================================================
net_colors <- c(
  "Cont"        = "#7A82E4",  
  "Default"     = "#9DB2F0",  
  "DorsAttn"    = "#B452CD",  
  "SalVentAttn" = "#E94CA1",  
  "SomMot"      = "#44C4D9",  
  "Vis"         = "#56008C",  
  "Limbic"      = "#F2ED96",  
  "Subcort"     = "#C5A1E8"   
)

grid_col <- net_colors[as.character(node_df$Network)]
names(grid_col) <- node_df$Label

net_runs <- rle(as.character(node_df$Network))
gaps <- rep(1, nrow(node_df))
gaps[cumsum(net_runs$lengths)] <- 6 

# Get min and max T-values to create a color gradient
min_t <- min(links$t_value)
max_t <- max(links$t_value)

# Function maps T-values to gradient (Light Pink -> Deep Dark Magenta)
edge_color_fun <- colorRamp2(c(min_t, max_t), c("#F8BBD0", "#880E4F"))
link_colors <- edge_color_fun(links$t_value)

# ==============================================================================
# 4. DRAW THE CHORD DIAGRAM
# ==============================================================================
circos.clear()
par(mar = c(2, 2, 2, 12)) # Leave room for TWO legends on the right

circos.par(gap.degree = gaps, start.degree = 90)

chordDiagram(
  x = links,
  order = node_df$Label,
  grid.col = grid_col,
  col = link_colors,         
  transparency = 0.3,
  annotationTrack = "grid", 
  preAllocateTracks = list(
    list(track.height = 0.25), 
    list(track.height = 0.05)  
  )
)

# Anatomical Labels (Outer)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.65)
}, bg.border = NA)

# Network Color Bands (Inner)
for(net in unique(node_df$Network)) {
  net_nodes <- node_df$Label[node_df$Network == net]
  highlight.sector(
    sector.index = net_nodes, 
    track.index = 2, 
    col = net_colors[net], 
    text = net,                  
    facing = "bending.inside",   
    text.vjust = 0.5,            
    text.col = "white",          
    niceFacing = TRUE,
    cex = 0.7,                   
    font = 2                     
  )
}

# ==============================================================================
# 5. LEGENDS
# ==============================================================================
# 1. Main Network Legend
legend(x = 1.1, y = 0.8, 
       legend = names(net_colors), 
       fill = net_colors, 
       title = "Network",
       box.lty = 0, cex = 0.9, xpd = TRUE)       

# 2. CONTINUOUS T-BAR LEGEND
# Define placement coordinates for the gradient box
xl <- 1.12
xr <- 1.18
yb <- -0.3
yt <-  0.1

# Generate gradient image
t_seq <- seq(min_t, max_t, length.out = 100)
grad_colors <- edge_color_fun(t_seq)
# Reverse to draw bottom (lowest t) to top (highest t)
legend_image <- as.raster(matrix(rev(grad_colors), ncol=1)) 

# Draw the Box and the Raster Image
rasterImage(legend_image, xl, yb, xr, yt, xpd = TRUE)
rect(xl, yb, xr, yt, xpd = TRUE)

# Add Text Labels (Min, Mid, Max)
text(x = xr + 0.02, y = yb, labels = round(min_t, 2), adj = c(0, 0.5), cex = 0.9, xpd = TRUE)
text(x = xr + 0.02, y = (yb+yt)/2, labels = round(mean(c(min_t, max_t)), 2), adj = c(0, 0.5), cex = 0.9, xpd = TRUE)
text(x = xr + 0.02, y = yt, labels = round(max_t, 2), adj = c(0, 0.5), cex = 0.9, xpd = TRUE)

# Add Title for Color Bar
text(x = xl - 0.02, y = yt + 0.05, labels = "t-statistic", adj = c(0, 0), font = 2, cex = 1.0, xpd = TRUE)

# ==============================================================================
# 6. EXPORT TO SVG
# ==============================================================================
out_svg <- file.path(OUTPUT_DIR, "Figure1_circular_subnetwork_T4.5.svg")
dev.copy(svg, file = out_svg, width = 12, height = 10)
dev.off()
cat(sprintf("Saved: %s\n", out_svg))