#############
# Proteomic analysis: Adult vs Aged, treatment: ethanol vs. water
# Mina Peyton
# 24DEC2021
#############
# Heatmaps

# Adult 15 common proteins
# Aged 15 common proteins

library(pheatmap)

library(readr)
Adult_common <- read_csv("Adult_common15.csv")
Aged_common <- read_csv("Aged_common15.csv")

colnames(Adult_common)

Etoh <-as.matrix(Adult_common[,14:19])
H20 <- as.matrix(Adult_common[,20:26])

rownames(H20) <-Adult_common$`Gene Symbol`
rownames(Etoh) <-Adult_common$`Gene Symbol`

Adult_hp <- cbind(H20,Etoh)

annotation_col = data.frame(
  Group = factor(rep(c("H20", "Etoh"), c(7,6))), Time = 1:13)

annotation_groups = list(Group = c(H20 = "blue", Etoh="red"))

p1 <- pheatmap(Adult_hp, scale = "row", cluster_cols = FALSE,
               clustering_distance_rows = "correlation",cluster_rows = TRUE,
               clustering_method = "average",
               show_rownames = TRUE,
               annotation_col = annotation_col[1],
               annotation_colors = annotation_groups[1],
               cellwidth = 15,
               cellheight = 15,
               show_colnames = FALSE)

# Open a pdf file
pdf("Adult_heatmap.pdf", width = , height = ) 
# 2. Create a plot
p1
# Close the pdf file
dev.off() 


#### Aged

colnames(Aged_common)

Etoh <-as.matrix(Aged_common[,13:18])
H20 <- as.matrix(Aged_common[,19:24])

rownames(H20) <-Aged_common$Gene
rownames(Etoh) <-Aged_common$Gene

Aged_hp <- cbind(H20,Etoh)

annotation_col = data.frame(
  Group = factor(rep(c("H20", "Etoh"), c(6,6))), Time = 1:12)

annotation_groups = list(Group = c(H20 = "blue", Etoh="red"))

p1 <- pheatmap(Aged_hp, scale = "row", cluster_cols = FALSE,
         clustering_distance_rows = "correlation",cluster_rows = TRUE,
         clustering_method = "average",
         show_rownames = TRUE,
         annotation_col = annotation_col[1],
         annotation_colors = annotation_groups[1],
         cellwidth = 15,
         cellheight = 15,
         show_colnames = FALSE)

# Open a pdf file
pdf("Aged_heatmap.pdf", width = , height = ) 
# 2. Create a plot
p1
# Close the pdf file
dev.off() 
