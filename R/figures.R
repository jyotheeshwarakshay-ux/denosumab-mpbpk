# Generates two publication-quality figures for the project writeup, saved
# as PNG into figures/ (created if missing). Uses ggplot2 for Figure 1
# (already a project dependency's ecosystem-neighbor and available in this
# environment) and base R graphics for Figure 2 (diagram/DiagrammeR are not
# installed here, so a hand-laid-out base-R schematic is used instead, per
# the task's own fallback instruction).
#
# Both figures use white backgrounds and a single maroon accent color
# (#6E1A2E) reserved for emphasis, not applied indiscriminately.
#
# Run from the project root (denosumab-mpbpk/) so the relative output path
# resolves.

library(ggplot2)

dir.create("figures", showWarnings = FALSE)

MAROON <- "#6E1A2E"

# ===========================================================================
# FIGURE 1: Bone ABC validation - this model vs. published benchmark,
# with the 2-fold accuracy band shown as a shaded region.
# ===========================================================================

bone_abc_model     <- 5.70   # R/06_bone_abc_validation.R
benchmark          <- 7.27   # Shah & Betts 2013
band_low           <- benchmark / 2
band_high          <- benchmark * 2

fig1_df <- data.frame(
  label = factor(
    c("This model\n(distribution-only)", "Published benchmark\n(Shah & Betts 2013)"),
    levels = c("This model\n(distribution-only)", "Published benchmark\n(Shah & Betts 2013)")
  ),
  value = c(bone_abc_model, benchmark)
)

p1 <- ggplot(fig1_df, aes(x = label, y = value)) +
  geom_rect(
    aes(xmin = -Inf, xmax = Inf, ymin = band_low, ymax = band_high),
    inherit.aes = FALSE, fill = MAROON, alpha = 0.08
  ) +
  geom_col(aes(fill = label), width = 0.55, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.2f%%", value)),
    vjust = -0.7, size = 5.2, fontface = "bold", color = "grey20"
  ) +
  scale_fill_manual(values = c(MAROON, "grey45")) +
  scale_y_continuous(
    limits = c(0, band_high * 1.15),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Bone antibody biodistribution coefficient: model vs. benchmark",
    subtitle = sprintf(
      "Shaded band = 2-fold accuracy range around the benchmark (%.1f%%–%.1f%%)",
      band_low, band_high
    ),
    x = NULL, y = "Bone ABC (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    plot.title    = element_text(face = "bold", size = 15, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11, color = "grey35", margin = margin(b = 10)),
    axis.text.x   = element_text(size = 12, color = "grey20", lineheight = 1.1),
    axis.text.y   = element_text(size = 11, color = "grey30"),
    axis.title.y  = element_text(size = 12, color = "grey20"),
    plot.margin   = margin(12, 16, 10, 10)
  )

ggsave("figures/fig1_bone_abc.png", p1,
       width = 1200 / 150, height = 800 / 150, dpi = 150, bg = "white")

# ===========================================================================
# FIGURE 2: Model topology schematic - generic tissue distribution loop,
# plus the one structural exception (liver portal-vein routing).
# Hand-laid-out base R graphics (rect / arrows / text) on a 0-100 grid.
# ===========================================================================

png("figures/fig2_topology.png", width = 1200, height = 800, res = 150)
par(mar = c(3, 1, 3.5, 1))

plot(NA, xlim = c(-8, 106), ylim = c(0, 102),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", asp = NA)

draw_box <- function(x0, x1, y0, y1, label, border = MAROON, fill = "white",
                      lwd = 2, cex = 1, font = 2, textcol = "grey15") {
  rect(x0, y0, x1, y1, col = fill, border = border, lwd = lwd)
  text(mean(c(x0, x1)), mean(c(y0, y1)), label, cex = cex, font = font, col = textcol)
}

# --- Main distribution loop (generic tissue) --------------------------------

draw_box(6, 38, 86, 98, "Central Plasma")
draw_box(6, 38, 62, 74, "Tissue Vascular")
draw_box(6, 38, 38, 50, "Tissue Interstitial")
draw_box(6, 38, 14, 26, "Shared Lymph Node")

arrows(22, 86, 22, 74, length = 0.12, lwd = 2, col = "grey25")
arrows(22, 62, 22, 50, length = 0.12, lwd = 2, col = "grey25")
arrows(22, 38, 22, 26, length = 0.12, lwd = 2, col = "grey25")

text(42, 80, "PLQ\n(plasma flow)", adj = 0, cex = 0.78, col = "grey30")
text(42, 56, expression("Lymph" %*% "(1" - sigma[V] * ")"), adj = 0, cex = 0.78, col = "grey30")
text(42, 32, expression("Lymph" %*% "(1" - sigma[IS] * ")"), adj = 0, cex = 0.78, col = "grey30")

# Return path: Lymph Node -> Central Plasma, routed along the left margin
lines(c(6, 0), c(20, 20), lwd = 2, col = "grey25")
lines(c(0, 0), c(20, 92), lwd = 2, col = "grey25")
arrows(0, 92, 6, 92, length = 0.12, lwd = 2, col = "grey25")
text(-6, 56, "Lymph return\n(Σ LF)", srt = 90, cex = 0.78, col = "grey30")

# --- Liver portal routing (the one structural exception) -------------------
# Liver sits in the SAME row as Central Plasma (visual parity: it returns to
# plasma directly, same as the main loop's lymph node does) and is reached
# by two short horizontal connectors in the gap between the two boxes -
# deliberately NOT long diagonals, so nothing crosses the main loop's labels.

draw_box(64, 100, 86, 98, "Liver")

# Central Plasma -> Liver, arterial (grey dashed: NOT special - Liver
# receives ordinary arterial inflow exactly like any other tissue)
arrows(38, 95, 64, 95, length = 0.1, lwd = 1.5, col = "grey45", lty = 2)
text(51, 97.3, "PLQ_Liver (arterial, as usual)", cex = 0.62, col = "grey45", font = 3)

# Liver's combined outlet -> Central Plasma (maroon: part of the exception)
arrows(64, 90, 38, 90, length = 0.1, lwd = 2, col = MAROON)
text(51, 87.7, "Combined portal outlet", cex = 0.65, col = MAROON, font = 3)

# Gut organs sit BELOW Liver, aligned with the Tissue Vascular row for
# visual parity, with short vertical arrows feeding straight up into Liver.
gut_x <- list(Spleen = c(64, 73), Pancreas = c(74.5, 83.5), SI = c(85, 91.5), LI = c(93, 100))
for (nm in names(gut_x)) {
  xr <- gut_x[[nm]]
  draw_box(xr[1], xr[2], 64, 76, nm, cex = 0.66, lwd = 1.5)
}

for (nm in names(gut_x)) {
  xm <- mean(gut_x[[nm]])
  arrows(xm, 76, xm, 86, length = 0.08, lwd = 1.5, col = MAROON)
}
text(82, 81, "(PLQ − LF) portal inflow", cex = 0.68, col = MAROON, font = 3)

title(main = "Antibody distribution topology: generic tissue vs. liver portal routing",
      cex.main = 1.15, font.main = 2, line = 1.6)

mtext(
  paste0("Liver receives portal venous inflow from Spleen, Pancreas, SI and LI ",
         "before its combined outflow reaches plasma — the one\nstructural exception in the 15-tissue model. ",
         "All other tissues follow the generic loop on the left."),
  side = 1, line = 1, cex = 0.68, col = "grey40"
)

dev.off()

cat("Figures written:\n")
cat(" -", normalizePath("figures/fig1_bone_abc.png"), "\n")
cat(" -", normalizePath("figures/fig2_topology.png"), "\n")
