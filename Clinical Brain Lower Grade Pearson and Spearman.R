#The dataset was manually imported from the cBioPortal platform by selecting the TCGA Brain Lower Grade Glioma study, with source data provided by GDAC Firehose.

# Load required libraries
install.packages("tidyverse")
install.packages("plotly")
install.packages("MASS")
installed.packages("fields")
install.packages("viridis")

library(tidyverse)       # Data wrangling and ggplot2
library(plotly)          # For interactive plots
library(MASS)            # For 2D kernel density estimation
library(fields)          # For interpolating density values at point locations
library(viridis)         # For perceptually uniform color palettes

# Load your pre-loaded TCGA clinical data
data <- lgg_tcga_clinical_data

# Select relevant columns and remove rows with missing values
data <- data %>%
  dplyr::select(Patient.ID, Fraction.Genome.Altered, Mutation.Count) %>%
  drop_na()

# Calculate 2D density for each point
dens <- kde2d(data$Fraction.Genome.Altered, data$Mutation.Count, n = 200)
data$density <- interp.surface(dens, data.frame(
  x = data$Fraction.Genome.Altered,
  y = data$Mutation.Count
))

# Compute Pearson and Spearman correlations
pearson <- cor.test(data$Fraction.Genome.Altered, data$Mutation.Count, method = "pearson")
spearman <- cor.test(data$Fraction.Genome.Altered, data$Mutation.Count, method = "spearman")

# Create a density-colored scatter plot
p <- ggplot(data, aes(x = Fraction.Genome.Altered, y = Mutation.Count)) +
  geom_point(aes(color = density,
                 text = paste0("Patient ID: ", Patient.ID,
                               "<br>Fraction Genome Altered: ", Fraction.Genome.Altered,
                               "<br>Mutation Count: ", Mutation.Count)),
             size = 1.3, alpha = 0.85) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom"
  ) +
  labs(
    title = "Mutation Count vs Fraction Genome Altered",
    x = "Fraction Genome Altered",
    y = "Mutation Count",
    color = "Density"
  )


# Convert to interactive Plotly object with a custom annotation for correlations
interactive_plot <- ggplotly(p, tooltip = "text") %>%
  layout(
    hoverlabel = list(bgcolor = "white"),
    annotations = list(
      list(
        x = 0.92, y = 0.05, xref = 'paper', yref = 'paper',
        text = paste0(
          "<span style='font-size:10px'><b>",
          "n = ", nrow(data),
          "<br>Pearson: ", round(pearson$estimate, 3),
          " (p < 0.001)",
          "<br>Spearman: ", round(spearman$estimate, 3),
          " (p < 0.001)</b></span>"
        ),
        showarrow = FALSE,
        xanchor = "left",
        yanchor = "bottom",
        align = "left",
        bgcolor = "white",
        bordercolor = "white",
        borderwidth = 0.5,
        borderpad = 1
      )
    ),
    plot_bgcolor = "white",
    paper_bgcolor = "white"
  )

# Display the final interactive plot
interactive_plot

