#!/usr/bin/env Rscript

trn_plot <- function(pc_a, pc_b) {
  ggplot2::ggplot(retx, ggplot2::aes_string(x = paste0("PC", pc_a), y = paste0("PC", pc_b))) +
    geom_jitter(aes(color = ancestry)) +
    scale_color_manual(values = brewer.pal(6, "Set2"), name = "Ancestry") +
    ggplot2::xlab(paste0("PC", pc_a, " (", round(trn_pve$pve[pc_a], 2) * 100, "% PVE)")) +
    ggplot2::ylab(paste0("PC3", pc_b," (", round(trn_pve$pve[pc_b], 2) * 100, "% PVE)")) +
    ggplot2::theme(text = element_text(size = 16), legend.text = element_text(size = 16),
                   legend.key = element_blank())
  }

pca_plot <- function(pc_a, pc_b, grp_by) {
  mplot <- ggplot2::ggplot(subset(pca_df, grp == grp_by), aes_string(x = pc_a, y = pc_b)) +
    geom_jitter(aes(color = val)) + theme(legend.position = "none") + xlim(-40, 40) +
    ylim(-40, 40)
  if (grp_by == "ancestry") {
    mplot <- mplot + scale_color_manual(values = brewer.pal(6, "Set2"), name = grp_by)
  } else if (grp_by == "grouping") {
    mplot <- mplot + scale_color_manual(values = brewer.pal(3, "Dark2")[1:2], name = "grouping")
  }
  return(mplot)
}
