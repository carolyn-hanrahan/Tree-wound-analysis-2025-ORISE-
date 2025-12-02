# Tree wound analysis R code 
# Carolyn Hanrahan 
# December 2nd, 2025 

library(vegan)
library(dplyr)
install.packages("readxl")   # only needed once
library(readxl)
library(tidyr)
install.packages(c("ggplot2", "vegan", "ggrepel"))
library(ggplot2)
library(vegan)
library(ggrepel)


df <- read_excel("Analysis - Summer 2025 Chestnut Experiment (1).xlsx")

# Select numeric variables + treatment label
wound_data <- df %>%
  select(Treatment,
         `Wound length (mm)`,
         `Total discoloration length (mm)`,
         `discolor - above wound`,
         `discolor - below wound`,
         `stem diam - at wound`,
         `stem diam - 1 cm above wound`,
         `stem diam - 1 cm below wound`) %>%
  drop_na()  # remove rows with missing measurements



str(wound_data)


# Convert all numeric variables to numeric explicitly
numeric_cols <- c(
  "Wound length (mm)",
  "Total discoloration length (mm)",
  "discolor - above wound",
  "discolor - below wound",
  "stem diam - at wound",
  "stem diam - 1 cm above wound",
  "stem diam - 1 cm below wound"
)

# Convert columns to numeric (coercing non-numeric entries to NA)
wound_data[numeric_cols] <- lapply(wound_data[numeric_cols], \(x) as.numeric(as.character(x)))

wound_data <- wound_data %>% drop_na(all_of(numeric_cols))


# Create matrix 
X <- wound_data %>% select(all_of(numeric_cols))
X_scaled <- scale(X)

# Run NMDS 
set.seed(123)
nmds_result <- metaMDS(X_scaled, distance = "bray", k = 2, trymax = 200)

# Test treatment group differences (PERMANOVA)
adonis2(X_scaled ~ Treatment, data = wound_data, method = "euclidean")

# Fit environmental vectors onto NMDS 
env <- envfit(nmds_result, X_scaled)
env
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot NMDS with treatment groups 

plot(nmds_result, type = "n")
points(nmds_result, display = "sites",
       col = as.factor(wound_data$Treatment),
       pch = 19, cex = 1.2)

ordiellipse(nmds_result, wound_data$Treatment, draw = "polygon",
            alpha = 0.2, col = 1:length(unique(wound_data$Treatment)))

plot(env, col = "black")   # add variable vectors like stem diameter

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Better version of plot 

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$Treatment <- wound_data$Treatment

env_vectors <- scores(env, display = "vectors")
env_vectors <- as.data.frame(env_vectors)

# rescale arrows for visibility
env_vectors$NMDS1 <- env_vectors$NMDS1 * 0.7
env_vectors$NMDS2 <- env_vectors$NMDS2 * 0.7

# extract environmental vectors 
env_vectors$label <- rownames(env_vectors)

# covex hulls for treatment groups 
library(dplyr)

hulls <- nmds_scores %>%
  group_by(Treatment) %>%
  slice(chull(NMDS1, NMDS2))

# plot in ggplot 

ggplot() +
  # hull polygons
  geom_polygon(data = hulls,
               aes(x = NMDS1, y = NMDS2, fill = Treatment,
                   group = Treatment),
               alpha = 0.25, color = NA) +
  
  # points
  geom_point(data = nmds_scores,
             aes(x = NMDS1, y = NMDS2, color = Treatment),
             size = 3, alpha = 0.8) +
  
  # envfit arrows
  geom_segment(data = env_vectors,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 0.7, color = "black") +
  
  # labels with repelling
  geom_text_repel(data = env_vectors,
                  aes(x = NMDS1, y = NMDS2, label = label),
                  size = 4, color = "black") +
  
  theme_minimal(base_size = 14) +
  labs(title = "NMDS of Chestnut Wound Characteristics",
       x = "NMDS1", y = "NMDS2") +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())


