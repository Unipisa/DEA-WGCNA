library(VennDiagram)
library(grid)

# Function to create and save a Venn diagram
save_venn_diagram_old <- function(gene_lists, filename, format = "pdf", width = 6, height = 6, resolution = 300) {
  if (format == "pdf") {
    pdf(filename, width = width, height = height)
  } else if (format == "png") {
    png(filename, width = width * resolution, height = height * resolution, res = resolution)
  } else {
    stop("Unsupported format. Use 'pdf' or 'png'.")
  }
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  venn.plot <- venn.diagram(
    x = gene_lists,
    filename = NULL,
    fill = c("blue", "red"),
    alpha = 0.5,
    cat.col = c("blue", "red"),
    cat.cex = 1.2,
    margin = 0
  )
  grid.draw(venn.plot)
  dev.off()
}

save_venn_diagram <- function(gene_lists, filename, format = "pdf",
                              width = 6, height = 6, resolution = 300) {
  
  # --- input checks ---
  if (!is.list(gene_lists) || length(gene_lists) != 2) {
    stop("gene_lists must be a list of length 2 (two gene sets).")
  }
  
  # normalizza a vettori character unici
  x1 <- unique(as.character(gene_lists[[1]]))
  x2 <- unique(as.character(gene_lists[[2]]))
  
  a1 <- length(x1)
  a2 <- length(x2)
  ia <- length(intersect(x1, x2))
  
  # --- open device ---
  if (format == "pdf") {
    pdf(filename, width = width, height = height)
  } else if (format == "png") {
    png(filename, width = width * resolution, height = height * resolution, res = resolution)
  } else {
    stop("Unsupported format. Use 'pdf' or 'png'.")
  }
  
  on.exit({
    dev.off()
  }, add = TRUE)
  
  # silenzia log di VennDiagram
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  # --- handle degenerate cases (empty sets) ---
  if (a1 == 0 && a2 == 0) {
    plot.new()
    text(0.5, 0.60, "Venn diagram not drawn", cex = 1.2)
    text(0.5, 0.45, "Both gene sets are empty", cex = 1.1)
    text(0.5, 0.30, sprintf("|A|=%d  |B|=%d  |A∩B|=%d", a1, a2, ia), cex = 1.0)
    return(invisible(NULL))
  }
  
  if (a1 == 0 || a2 == 0) {
    plot.new()
    text(0.5, 0.60, "Venn diagram not drawn", cex = 1.2)
    text(0.5, 0.45, "One gene set is empty", cex = 1.1)
    text(0.5, 0.30, sprintf("|A|=%d  |B|=%d  |A∩B|=%d", a1, a2, ia), cex = 1.0)
    return(invisible(NULL))
  }
  
  # --- draw with protection ---
  res <- tryCatch({
    venn.plot <- VennDiagram::venn.diagram(
      x = list(x1, x2),
      filename = NULL,
      disable.logging = TRUE,
      fill = c("blue", "red"),
      alpha = 0.5,
      cat.col = c("blue", "red"),
      cat.cex = 1.2,
      margin = 0,
      scaled = FALSE  # aumenta robustezza (evita geometrie degeneri)
    )
    grid::grid.newpage()
    grid::grid.draw(venn.plot)
    TRUE
  }, error = function(e) {
    # fallback: scrive messaggio invece di far crashare tutto
    plot.new()
    text(0.5, 0.60, "VennDiagram error", cex = 1.2)
    text(0.5, 0.45, conditionMessage(e), cex = 0.9)
    text(0.5, 0.30, sprintf("|A|=%d  |B|=%d  |A∩B|=%d", a1, a2, ia), cex = 1.0)
    FALSE
  })
  
  invisible(res)
}


# Utility function to save PCA plots
save_PCA <- function(pca_scores, classes, labels, output_file, colors = c("red", "blue"), markers = c(19, 17), width = 6, height = 6) {
  # Define colors and markers based on classes
  color_map <- ifelse(classes == 1, colors[1], colors[2])
  marker_map <- ifelse(classes == 1, markers[1], markers[2])
  
  # Determine plot ranges with padding
  x_range <- range(pca_scores[, 1])
  y_range <- range(pca_scores[, 2])
  x_padding <- diff(x_range) * 0.05
  y_padding <- diff(y_range) * 0.05
  
  # Open the appropriate graphics device
  if (grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
    pdf(output_file, width = width, height = height)
  } else if (grepl("\\.png$", output_file, ignore.case = TRUE)) {
    png(output_file, width = width * 100, height = height * 100, res = 300)
  } else {
    stop("Unsupported file format. Use '.pdf' or '.png'.")
  }
  
  # Set plot margins
  par(mar = c(5, 5, 1, 1))
  
  # Create PCA plot
  plot(
    pca_scores[, 1], pca_scores[, 2],
    xlab = "Principal Component 1",
    ylab = "Principal Component 2",
    main = "",
    col = color_map,
    pch = marker_map,
    xlim = c(x_range[1] - x_padding, x_range[2] + x_padding),
    ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)
  )
  
  # Add point labels
  text(
    pca_scores[, 1], pca_scores[, 2],
    labels = labels,
    pos = 3,
    cex = 1
  )
  
  # Add legend
  legend(
    "bottomleft",
    legend = c("Cases", "Controls"),
    col = colors,
    pch = markers,
    title = "Classes"
  )
  
  # Close the device
  dev.off()
}
save_variance_plot <- function(eigenvalues, output_file, width = 6, height = 6) {
  variance_explained <- eigenvalues / sum(eigenvalues) * 100
  
  if (grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
    pdf(output_file, width = width, height = height)
  } else if (grepl("\\.png$", output_file, ignore.case = TRUE)) {
    png(output_file, width = width * 100, height = height * 100, res = 300)
  } else {
    stop("Unsupported file format. Use '.pdf' or '.png'.")
  }
  
  par(mar = c(5, 5, 2, 1))
  
  barplot(
    variance_explained,
    names.arg = seq_along(eigenvalues),
    xlab = "Principal Component",
    ylab = "Explained Variance (%)",
    main = "Explained Variance by Principal Components"
  )
  
  dev.off()
}
# Utility function to save hierarchical clustering plots
save_hierarchical_clustering <- function(sample_tree, output_file, width = 6, height = 4) {
  # Open the appropriate graphics device
  if (grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
    pdf(output_file, width = width, height = height)
  } else if (grepl("\\.png$", output_file, ignore.case = TRUE)) {
    png(output_file, width = width * 100, height = height * 100, res = 300)
  } else {
    stop("Unsupported file format. Use '.pdf' or '.png'.")
  }
  
  # Set graphical parameters
  par(cex = 0.6)
  par(mar = c(0, 4, 2, 0))
  
  # Plot the hierarchical clustering dendrogram
  plot(
    sample_tree,
    main = "Hierarchical Clustering Dendrogram",
    sub = "",
    xlab = "",
    ylab = "Euclidean distance (average linkage)",
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5
  )
  
  # Close the graphics device
  dev.off()
}


save_volcano_plot <- function(out_table, ALPHA_LIMMA, TH_FC, output_file, width, height) {
  if (grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
    pdf(output_file, width = width, height = height)
  } else if (grepl("\\.png$", output_file, ignore.case = TRUE)) {
    png(output_file, width = width * 100, height = height * 100, res = 300)
  } else {
    stop("Unsupported file format. Use '.pdf' or '.png'.")
  }
  par(mar = c(5, 5, 1, 1)) # Increase right margin for legend
  # Define significance criteria
  significant <- with(out_table, adj.P.Val < ALPHA_LIMMA & abs(logFC) > TH_FC)
  upregulated <- with(out_table, adj.P.Val < ALPHA_LIMMA & logFC > TH_FC)
  downregulated <- with(out_table, adj.P.Val < ALPHA_LIMMA & logFC < -TH_FC)
  # Set up colors for points
  colors <- ifelse(upregulated, "red", ifelse(downregulated, "blue", "grey"))
  # Create the volcano plot
  plot(
    out_table$logFC, 
    -log10(out_table$adj.P.Val), 
    col = colors, 
    pch = 20, 
    xlab = "Log2(FC)", 
    ylab = "-Log10(p_adj)", 
  )
  # Add horizontal and vertical dashed lines at thresholds
  abline(h = -log10(ALPHA_LIMMA), col = "black", lty = 2)  # Horizontal line for p-value threshold
  abline(v = TH_FC, col = "black", lty = 2)             # Vertical line for logFC threshold
  abline(v = -TH_FC, col = "black", lty = 2)            # Vertical line for logFC threshold
  dev.off()
}


# Utility function to save Scale-Free Topology (SFT) plots
save_sft_plots <- function(sft, powers, output_file, width = 9, height = 5) {
  # Generate the full output file name
  formatted_file <- output_file
  print(formatted_file)
  # Open the appropriate graphics device
  if (grepl("\\.pdf$", formatted_file, ignore.case = TRUE)) {
    pdf(formatted_file, width = width, height = height)
  } else if (grepl("\\.png$", formatted_file, ignore.case = TRUE)) {
    png(formatted_file, width = width, height = height, units = "in", res = 300)
  } else {
    stop("Unsupported file format. Use '.pdf' or '.png'.")
  }
  
  # Set plotting layout and parameters
  par(mfrow = c(1, 2))
  par(cex = 0.9)
  
  # Scale-free topology fit index plot
  plot(
    sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit, signed R^2",
    type = "n",
    main = "Scale independence"
  )
  text(
    sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = 0.9,
    col = "red"
  )
  abline(h = 0.90, col = "red")
  
  # Mean connectivity plot
  plot(
    sft$fitIndices[, 1], sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = "Mean connectivity"
  )
  text(
    sft$fitIndices[, 1], sft$fitIndices[, 5],
    labels = powers,
    cex = 0.9,
    col = "red"
  )
  
  # Close the graphics device
  dev.off()
}

# Utility function to save dendrogram plots
save_dendrogram_plot <- function(net, output_file, width = 12, height = 9) {
  # Generate the full output file name
  formatted_file <- output_file
  
  # Open the appropriate graphics device
  if (grepl("\\.pdf$", formatted_file, ignore.case = TRUE)) {
    pdf(formatted_file, width = width, height = height)
  } else if (grepl("\\.png$", formatted_file, ignore.case = TRUE)) {
    png(formatted_file, width = width, height = height, units = "in", res = 300)
  } else {
    stop("Unsupported file format. Use '.pdf' or '.png'.")
  }
  
  mergedColors <- labels2colors(net$colors)
  
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  dev.off()
}

save_ME_old <- function(MEs, output_file){
  nrows_plot = ncol(MEs)%/%4+1
  png(output_file, width = 12, height = 9, units = "in", res = 300) # Adjust resolution as needed
  par(mfrow = c(nrows_plot, 5), mar = c(3, 3, 2, 1)) # Adjust margins if needed
  for (i in 1:ncol(MEs)) {
    plot(MEs[, i], type = "l", main = colnames(MEs)[i], xlab = "Observation", ylab = "Value", xaxt = "n")
    axis(1, at = 1:nrow(MEs), labels = rownames(MEs), las = 2, cex.axis = 0.7)
  }
  dev.off()
}

save_ME <- function(MEs, output_file, ncol_plot = 4, nrow_plot = 3,
                    width_in = 12, height_in = 9, res = 300) {
  
  n_mod <- ncol(MEs)
  if (n_mod == 0) return(invisible(NULL))
  
  per_page <- ncol_plot * nrow_plot
  n_pages <- ceiling(n_mod / per_page)
  
  # base name: "xxx.png" -> "xxx_p01.png"
  base <- sub("\\.png$", "", output_file, ignore.case = TRUE)
  
  for (p in seq_len(n_pages)) {
    idx <- ((p - 1) * per_page + 1):min(p * per_page, n_mod)
    
    fn <- sprintf("%s_p%02d.png", base, p)
    png(fn, width = width_in, height = height_in, units = "in", res = res)
    on.exit(dev.off(), add = TRUE)
    
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    par(mfrow = c(nrow_plot, ncol_plot),
        mar = c(2.5, 2.5, 2, 0.5),
        oma = c(0, 0, 0, 0))
    
    for (i in idx) {
      y <- MEs[, i]
      if (length(y) == 0 || all(is.na(y))) {
        plot.new()
        title(main = colnames(MEs)[i])
        text(0.5, 0.5, "All NA / empty")
        next
      }
      
      plot(y, type = "l",
           main = colnames(MEs)[i],
           xlab = "Observation", ylab = "Value",
           xaxt = "n")
      
      n <- nrow(MEs)
      at <- pretty(seq_len(n), n = 6)
      at <- at[at >= 1 & at <= n]
      axis(1, at = at, labels = rownames(MEs)[at], las = 2, cex.axis = 0.6)
    }
    
    dev.off()
    on.exit(NULL, add = FALSE)  # evita di accumulare on.exit nel loop
  }
  
  invisible(NULL)
}


save_ME_turquoise <- function(MEs, output_file){
  pdf(output_file, width = 6, height = 4)
  par(mar = c(4, 4, 2, 1)) # Adjust margins if needed
  plot(MEs[,"MEturquoise"], type = "l", main = "MEturquoise", xlab = "Observation", ylab = "PC1", xaxt = "n")
  axis(1, at = 1:nrow(MEs), labels = rownames(MEs), las = 2, cex.axis = 0.7)
  dev.off()
}

# Utility function to save heatmap plots
save_heatmap_plot <- function(matrix, xLabels, yLabels, ySymbols, textMatrix, output_file,
                              height_adjustment = 0, width = 8, height_per_row = 0.35,
                              min_height = 4, max_height = 30,# clamp device height (inches)
                              mar_bottom = 6, mar_left = 6.5, mar_top = 2.5, mar_right = 2,
                              colors = blueWhiteRed(50), zlim = c(-1, 1), png_res = 300) {
  
  # Calculate dynamic height based on the number of modules
  height <- length(yLabels) + height_adjustment
  
  # Generate the full output file name
  formatted_file <- output_file
  
  # Open the appropriate graphics device
  if (grepl("\\.pdf$", formatted_file, ignore.case = TRUE)) {
    pdf(formatted_file, width = width, height = height)
  } else if (grepl("\\.png$", formatted_file, ignore.case = TRUE)) {
    png(formatted_file, width = width, height = height, units = "in", res = 300)
  } else {
    stop("Unsupported file format. Use '.pdf' or '.png'.")
  }
  
  # Plot the heatmap
  par(mfrow = c(1,1), mar = c(mar_bottom, mar_left, mar_top, mar_right), xpd = NA)
  
  labeledHeatmap(
    Matrix = matrix,
    xLabels = xLabels,
    yLabels = yLabels,
    ySymbols = ySymbols,
    colorLabels = FALSE,
    colors = colors,
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.75,
    zlim = zlim,
    main = paste("Module-trait relationships")
  )
  
  # Close the graphics device
  dev.off()
}


save_scatter_plot <- function(geneModuleMembership, geneTraitSignificance, moduleGenes, column, module, output_file, width = 12, height = 9) {
  
  # Generate the full output file name
  formatted_file <- output_file
  
  # Open the appropriate graphics device
  if (grepl("\\.pdf$", formatted_file, ignore.case = TRUE)) {
    pdf(formatted_file, width = width, height = height)
  } else if (grepl("\\.png$", formatted_file, ignore.case = TRUE)) {
    png(formatted_file, width = width, height = height, units = "in", res = 300)
  } else {
    stop("Unsupported file format. Use '.pdf' or '.png'.")
  }
  
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for trait GROUP",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
}

save_scatter_plot_with_selections <- function(
    geneModuleMembership, geneTraitSignificance, moduleGenes, column, module, output_file,
    MM_threshold, GS_threshold, r_squared, a, b,
    width = 12, height = 9
) {
  formatted_file <- output_file
  if (grepl("\\.pdf$", formatted_file, ignore.case = TRUE)) {
    pdf(formatted_file, width = width, height = height)
  } else if (grepl("\\.png$", formatted_file, ignore.case = TRUE)) {
    png(formatted_file, width = width, height = height, units = "in", res = 300)
  } else {
    stop("Unsupported file format. Use '.pdf' or '.png'.")
  }
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  
  # ---- Plot points (all) ----
  verboseScatterplot(x, y,
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for trait GROUP",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "cyan")
  
  # ---- Highlight selections ----
  # Ellipse
  ellipse_inside <- ((x - 1)^2)/(a^2) + ((y - 1)^2)/(b^2) <= 1
  points(x[ellipse_inside], y[ellipse_inside], col = "purple", pch = 19)
  # Circle
  radius <- sqrt(r_squared)
  circle_inside <- (x - 1)^2 + (y - 1)^2 <= r_squared
  points(x[circle_inside], y[circle_inside], col = "red", pch = 19)
  # Square
  square_inside <- (x >= MM_threshold) & (y >= GS_threshold)
  points(x[square_inside], y[square_inside], col = "blue", pch = 19)
  
  # ---- Draw threshold boundaries ----
  # Square (blue dashed)
  abline(h = GS_threshold, col = "blue", lty = 3, lwd = 2)
  abline(v = MM_threshold, col = "blue", lty = 3, lwd = 2)
  # Circle (red dashed)
  theta <- seq(0, 2 * pi, length.out = 200)
  lines(1 + radius * cos(theta), 1 + radius * sin(theta), col = "red", lwd = 2, lty = 2)
  # Ellipse (purple dashed)
  lines(1 + a * cos(theta), 1 + b * sin(theta), col = "purple", lwd = 2, lty = 2)
  
  # ---- Optional: Add a legend ----
  legend("bottomright",
         legend = c("All genes", "Ellipse", "Circle", "Square"),
         col = c("cyan", "purple", "red", "blue"),
         pch = 19, lty = c(NA, 2, 2, 3), lwd = 2, bty = "n")
  
  dev.off()
}
