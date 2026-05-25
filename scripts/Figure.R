
# ============================================================================
# Line chart of running profiles of LA-HM with its four variants (Figure 3)
# ============================================================================

# Load required libraries
library(ggplot2)
library(reshape2)
library(ggsci)

# Set the folder path containing the CSV files
folder_path <- "combined_csv"

# Retrieve paths of all CSV files in the specified folder
csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file to process data and generate charts
for (file_path in csv_files) {
  
  # Read the CSV file
  mydata2 <- read.csv(file_path)
  
  # Rename a specific column for consistency in plotting
  colnames(mydata2)[colnames(mydata2) == "LDDITSN"] <- "LA-HM"
  
  # Reshape data from wide format to long format for ggplot2 compatibility
  mydata2 <- melt(mydata2, id.vars = "Second", variable.name = "algorithm", value.name = "value")
  
  # Rename columns to standard names
  colnames(mydata2) <- c("Time", "algorithm", "value")
  
  # Determine the range of the dependent variable
  y_min <- min(mydata2$value)
  y_max <- max(mydata2$value)
  
  # Dynamically calculate the step size for the Y-axis to aim for 6 tick marks
  step1 <- (y_max - y_min) / 5  # 5 intervals yield 6 tick marks
  step1 <- round(step1, -2)     # Round step size to the nearest 100 (adjust -2 to -1 or others as needed)
  
  # Generate the ggplot visualization
  P6 <- ggplot(data = mydata2, aes(x = Time, y = value, group = algorithm, color = algorithm, shape = algorithm)) +
    geom_line(size = 1.1) +  # Set line width
    # geom_point(size = 1.1) +  # Optional: add points to lines
    xlab("Time") +  # Set X-axis label
    ylab("Average Objective Value") +  # Set Y-axis label
    theme_bw() +  # Use a black-and-white theme
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.background = element_rect(fill = "white", colour = NA),  # White panel background
      plot.background = element_rect(fill = "transparent", colour = NA),  # Transparent plot background
      axis.title = element_text(size = 14, face = "bold"),  # Format axis titles
      axis.text = element_text(size = 12),  # Format axis tick text
      legend.position = "bottom",  # Move legend below the plot
      legend.title = element_blank(),  # Remove the legend title
      legend.text = element_text(size = 11),  # Format legend text
      legend.key = element_rect(fill = "transparent"),  # Transparent legend key background
      legend.box = "horizontal",  # Align legend items horizontally
      legend.box.margin = margin(-5, 0, -10, 0),  # Adjust spacing around the legend box
      plot.margin = margin(10, 10, 5, 10)  # Reduce whitespace at the bottom of the plot
    ) +
    scale_x_continuous(limits = c(0, 200), breaks = seq(0, 200, 20)) +  # Set static X-axis limits and ticks
    scale_y_continuous(breaks = seq(0, max(mydata2$value), by = step1)) +  # Apply dynamic Y-axis ticks
    # scale_color_rickandmorty()  # Alternative color scheme (commented out)
    scale_color_jco()  # Apply JCO color palette from ggsci package
  # scale_color_brewer(palette = "RdBu")  # Alternative color scheme (commented out)
  # scale_color_viridis_d(option = "plasma")  # Alternative color scheme (commented out)
  
  # Extract the base file name (without path and extension)
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Construct the output PDF path
  pdf_file_path <- file.path(folder_path, paste0(file_name, "_JCO.pdf"))
  
  # Save the generated chart to a PDF file with specified dimensions
  ggsave(pdf_file_path, plot = P6, width = 6, height = 5, units = "in", device = "pdf")
  
  # Print a confirmation message to the console
  cat("PDF chart saved to:", pdf_file_path, "\n")
}



# ============================================================================
# Boxplot f percent gaps among LA-HM, LA-ITS and MLNS for best and average
# objective values (Figure OS1)
# ============================================================================

library(ggplot2)
library(ggsci)
library(dplyr)
library(cowplot)
library(grid) 

# =======================================================
# Data Preparation
# =======================================================
SetBgap <- read.csv("SetBgapLNSVerify.csv")
SetBgap$Algorithm <- factor(SetBgap$Algorithm, levels = c("LA-HM", "LA-ITS", "MLNS"))
SetBgap$Type <- factor(SetBgap$Type, levels = c("Best", "Average"))

pos <- position_dodge(width = 0.8) 
my_fill <- scale_fill_jco()
bg_color <- "#E5E5E5" 

# =======================================================
# Create the Inset Plot
# =======================================================
data_zoom <- SetBgap %>% filter(Algorithm %in% c("LA-HM", "LA-ITS"))

p_inset <- ggplot(data_zoom, aes(x = Algorithm, y = Value, fill = Type)) +
  geom_boxplot(alpha = 0.8, width = 0.7, position = pos, outlier.size = 1) +
  stat_summary(aes(group = Type), fun = mean, geom = "point",
               shape = 18, size = 2.5, color = "red", position = pos) +
  my_fill + 
  scale_y_continuous(name = NULL) + 
  scale_x_discrete(name = NULL) +   
  theme_bw() +
  theme(
    legend.position = "none",
    # Set plot background to white with a dashed black border
    plot.background = element_rect(fill = "white", color = "black", linetype = "dashed", size = 0.5), 
    
    # Set inner panel background to white
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 10, color = "black"), 
    
    # Use light gray grid lines for visibility
    panel.grid.major = element_line(size = 0.2, color = "gray90"),
    plot.margin = margin(5, 5, 5, 5) 
  )

# =======================================================
# Create the Main Plot
# =======================================================
p_main <- ggplot(SetBgap, aes(x = Algorithm, y = Value, fill = Type)) +
  
  # Initialize axes to prevent potential annotation rendering issues
  geom_blank() + 
  
  # 1. Background rectangle (aligned with data coordinates)
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -0.15, ymax = 0.8, 
           fill = "white", alpha = 0.5) + 
  
  # 2. Dashed border for the zoom-in target area
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -0.15, ymax = 0.8, 
           alpha = 0, color = "black", linetype = "dashed", size = 0.3) +
  
  # 3. Standard plot layers
  geom_boxplot(alpha = 0.7, width = 0.7, position = pos, outlier.size = 1.5) +
  stat_summary(aes(group = Type), fun = mean, geom = "point",
               shape = 18, size = 3, color = "red", position = pos) +
  my_fill +
  scale_y_continuous(name = "Percent gap (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    text = element_text(size = 14),
    axis.title = element_text(face="bold"),
    axis.title.x = element_blank(),
    legend.position = c(0.1, 0.9),
    legend.title = element_blank(),
    legend.background = element_rect(fill="white", color="white")
  )

# =======================================================
# Define Coordinate Configurations (Precisely Adjusted)
# =======================================================
# Notes:
# box_top_y: Controls the vertical base position of the connector polygon. Decrease this value if it does not touch the target.
# box_top_left_x / right_x: Controls the width of the connection base. Adjust to align seamlessly.

zoom_conf <- list(
  # Inset dimensions and positions
  # 1. Scale up the inset plot size
  inset_w = 0.55,   
  inset_h = 0.54,   
  
  # 2. Adjust position to accommodate the new size
  inset_x = 0.112, 
  inset_y = 0.27,   
  
  # Coordinates for the base of the connection area
  box_top_left_x  = 0.103,   
  box_top_right_x = 0.673,  
  box_top_y       = 0.177    
)

# =======================================================
# Assemble Plots (Low-level Grid Layout)
# =======================================================

# 1. Prepare the connection polygon
poly_grob <- polygonGrob(
  x = unit(c(zoom_conf$box_top_left_x, zoom_conf$box_top_right_x, 
             zoom_conf$inset_x + zoom_conf$inset_w, zoom_conf$inset_x), "npc"),
  y = unit(c(zoom_conf$box_top_y, zoom_conf$box_top_y, 
             zoom_conf$inset_y, zoom_conf$inset_y), "npc"),
  gp = gpar(fill = "white", col = NA, alpha = 0.5) 
)

# 2 Prepare the connecting dashed lines
line_left <- segmentsGrob(
  x0 = unit(zoom_conf$box_top_left_x, "npc"), y0 = unit(zoom_conf$box_top_y, "npc"),
  x1 = unit(zoom_conf$inset_x, "npc"),        y1 = unit(zoom_conf$inset_y, "npc"),
  gp = gpar(col = "black", lty = "dashed", lwd = 1, size = 0.5)
)

line_right <- segmentsGrob(
  x0 = unit(zoom_conf$box_top_right_x, "npc"),               y0 = unit(zoom_conf$box_top_y, "npc"),
  x1 = unit(zoom_conf$inset_x + zoom_conf$inset_w, "npc"),   y1 = unit(zoom_conf$inset_y, "npc"),
  gp = gpar(col = "black", lty = "dashed", lwd = 1, size = 0.5)
)

# 3 Combine the elements
final_plot <- ggdraw() +
  draw_plot(p_main) +             # Main plot
  draw_grob(poly_grob) +          # Connection polygon
  draw_grob(line_left) +          # Left dashed line
  draw_grob(line_right) +         # Right dashed line
  draw_plot(p_inset,              # Inset plot
            x = zoom_conf$inset_x, 
            y = zoom_conf$inset_y, 
            width = zoom_conf$inset_w, 
            height = zoom_conf$inset_h)


print(final_plot)

# Save the final visualization as a PDF
ggsave("SetBGapLNSVerify.pdf", plot = final_plot, width = 6, height = 4.5)

