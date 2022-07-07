#############
# Proteomic analysis: Adult vs Aged, treatment: ethanol vs. water
# Mina Peyton
# 24DEC2021
#############
# Volcano Plots
#
# Adult etoh/h20 & Aged etoh/h20

library(ggplot2)
library(ggrepel)
library(ggpubr)

library(readxl)
Adult_4190 <- read_excel("Adult_4190.xlsx")
Aged_4192 <- read_excel("Aged_4192.xlsx")

# rename columns needed
colnames(Adult_4190)
names(Adult_4190)[1]<- "Accession"
names(Adult_4190)[5]<- "log2ratio"
names(Adult_4190)[6]<-"pvalue"

colnames(Aged_4192)
names(Aged_4192)[1] <- "Accession"
names(Aged_4192)[5]<- "log2ratio"
names(Aged_4192)[6]<-"pvalue"

# Add a column to the data frame to specify if they are UP- or DOWN- regulated 
# (log2FoldChange respectively positive or negative)

# add a column of NAs
Adult_4190$diffexpressed <- "NO"
# if log2Foldchange > 0.485 and pvalue < 0.05, set as "UP" 
Adult_4190$diffexpressed[Adult_4190$log2ratio > 0 & Adult_4190$pvalue< 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Adult_4190$diffexpressed[Adult_4190$log2ratio < -0 & Adult_4190$pvalue < 0.05] <- "DOWN"

# make basic scatter plot and conver to -log10pvalue
p <- ggplot(data=Adult_4190, aes(x=log2ratio, 
                              y= -log10(pvalue),
                              col = diffexpressed))+
  geom_point(aes(fill = diffexpressed),pch=21, colour="Black", size=2) +
  theme_classic()+
  scale_fill_manual(values= mycolors)

# Add vertical lines for log2FoldChange thresholds and one horizontal line for 
# the p-value threshold + color differentially expressed phosphopeptides
p2 <- p + geom_vline(xintercept=c(-0.485, 0.485), col="gray") +
  geom_hline(yintercept=-log10(0.05), col="gray") + 
  scale_x_continuous(limits=c(-1, 2), breaks = seq(-1,2,1)) +
  scale_y_continuous(limits=c(0,4), breaks = seq(0,4,1))

# scale_color_manual(values=c("blue", "red", "gray"))

# Open a pdf file
pdf("20220307_AdultProteins_volcano.pdf", width = , height = ) 
# 2. Create a plot
p2
# Close the pdf file
dev.off() 


############ Aged volcano plot

# add a column of NAs
Aged_4192$diffexpressed <- "NO"
# if log2Foldchange > 0.485 and pvalue < 0.05, set as "UP" 
Aged_4192$diffexpressed[Aged_4192$log2ratio > 0 & Aged_4192$pvalue< 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Aged_4192$diffexpressed[Aged_4192$log2ratio < -0 & Aged_4192$pvalue < 0.05] <- "DOWN"

# make basic scatter plot and conver to -log10pvalue
p <- ggplot(data=Aged_4192, aes(x=log2ratio, 
                                 y= -log10(pvalue),
                                 col = diffexpressed))+
  geom_point(aes(fill = diffexpressed),pch=21, colour="Black", size=2)+
  theme_classic()+
  scale_fill_manual(values= mycolors)

# automate a bit: create a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "gray")
names(mycolors) <- c("DOWN", "UP", "NO")

# Add vertical lines for log2FoldChange thresholds and one horizontal line for 
# the p-value threshold + color differentially expressed phosphopeptides
p2 <- p + geom_vline(xintercept=c(-0.485, 0.485), col="gray") +
  geom_hline(yintercept=-log10(0.05), col="gray") + 
  scale_x_continuous(limits=c(-.50, 1), breaks = seq(-.50,1,.5)) +
  scale_y_continuous(limits=c(0,3), breaks = seq(0,3,1))

scale_color_manual(values=c("green", "red", "gray"))

# Open a pdf file
pdf("20220307_AgedProteins_volcano.pdf", width = , height = ) 
# 2. Create a plot
p2
# Close the pdf file
dev.off() 

# automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "gray")
names(mycolors) <- c("DOWN", "UP", "NO")
# p3 <- p2 + scale_colour_manual(values = mycolors)

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
Ovx_PGs$delabel <- NA
Ovx_PGs$delabel[Ovx_PGs$diffexpressed != "NO"] <- Ovx_PGs$Accession[Ovx_PGs$diffexpressed != "NO"]

ggplot(data=Old_PGs, aes(x=log2ratio, 
                         y=-log10(pvalue), 
                         col=diffexpressed, label=delabel))+
  geom_point() + 
  theme_minimal() +
  geom_text()

# organize the labels 
library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=Old_PGs, aes(x=log2ratio, 
                         y=-log10(pvalue), 
                         col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel()+
  scale_color_manual(values=c("green", "black", "red")) +
  geom_vline(xintercept=c(-0.485, 0.485), col="gray") +
  geom_hline(yintercept=-log10(0.05), col="gray") +
  scale_x_continuous(limits=c(-4, 4), breaks = seq(-4,4,2))

##### make sure x and y axis scales can fit data points == always set
# or else will use default and eliminate points automatically

####### make volcano plot function with labels
volcanoPlot <- function(data){
  ggplot(data, aes(x=log2ratio,
                   y= -log10(pvalue),
                   fill = diffexpressed,
                   label=delabel)) +
    geom_point(aes(fill = diffexpressed),pch=21, colour="Black", size=2) +
    scale_fill_manual(values= mycolors) +
    theme_classic() +
    geom_text_repel()+
    geom_vline(xintercept=c(-0.485, 0.485), col="gray") +
    geom_hline(yintercept=-log10(0.05), col="gray") +
    scale_x_continuous(limits=c(-8, 8), breaks = seq(-8,8,2)) +
    scale_y_continuous(limits=c(0,20), breaks = seq(0,20,5))
}

# make volcano pot without labels
volcanoPlot2 <- function(data){
  ggplot(data, aes(x=log2ratio,
                   y= -log10(pvalue),
                   col = diffexpressed))+
    geom_point(aes(fill = diffexpressed),pch=21, colour="Black", size=2) +
    theme_classic()+
    scale_fill_manual(values= mycolors)+
    geom_vline(xintercept=c(-0.485, 0.485), col="gray") +
    geom_hline(yintercept=-log10(0.05), col="gray") +
    scale_x_continuous(limits=c(-8, 8), breaks = seq(-8,8,2)) +
    scale_y_continuous(limits=c(0,20), breaks = seq(0,20,5))
}
