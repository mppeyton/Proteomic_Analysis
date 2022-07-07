#############
# Proteomic analysis: Adult vs Aged, treatment: ethanol vs. water
# Mina Peyton
# 24DEC2021
#############
# Overrepresentation Analysis

# Aged etoh/h20 == modsites threshold: pvalue < 0.05 & 1FC
# Adult etoh/h20 == modsites threshold: pvalue < 0.05 & 1FC 

library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(ReactomePA)
library(writexl)
library(ggplot2)
library(forcats)

BiocManager::install("org.Rn.eg.db", force = TRUE)

# convert AccessionIDs into EntrezID
Adultsig <- Adult_sig$Accession
Adultdown <- Adult_down$Accession
Adultup <- Adult_up$Accession

Adultsig <- bitr(Adultsig, fromType="UNIPROT", toType="ENTREZID", 
                  OrgDb="org.Rn.eg.db") #288
Adultsig <- Adultsig$ENTREZID

Adultdown <- bitr(Adultdown, fromType="UNIPROT", toType="ENTREZID", 
               OrgDb="org.Rn.eg.db") #142
Adultdown <- Adultdown$ENTREZID

Adultup <- bitr(Adultup, fromType="UNIPROT", toType="ENTREZID", 
                 OrgDb="org.Rn.eg.db") #146
Adultup <- Adultup$ENTREZID

##### Aged
Agedsig <-Aged_sig$Accession
Ageddown <- Aged_down$Accession
Agedup <- Aged_up$Accession

Agedsig <- bitr(Agedsig, fromType="UNIPROT", toType="ENTREZID", 
                 OrgDb="org.Rn.eg.db") #135
Agedsig <- Agedsig$ENTREZID 

Ageddown <- bitr(Ageddown, fromType="UNIPROT", toType="ENTREZID", 
                  OrgDb="org.Rn.eg.db") #76
Ageddown <- Ageddown$ENTREZID

Agedup <- bitr(Agedup, fromType="UNIPROT", toType="ENTREZID", 
                OrgDb="org.Rn.eg.db") #59
Agedup <- Agedup$ENTREZID

# GO analysis of Adult proteome
# compare up and down clusters
Adultproteome <- list("Down" = Adultdown, 
                   "Up" = Adultup)

AdultBP <- compareCluster(geneCluster = Adultproteome, fun = enrichGO,
                       ont = "BP",
                       OrgDb = org.Rn.eg.db,
                       pvalueCutoff = 2, readable = TRUE)
# ck <- setReadable(ck, OrgDb = org.Rn.eg.db, keyType="ENTREZID")

write.csv(AdultMF@compareClusterResult, "AdultMFresults.csv") 
write.csv(AdultCC@compareClusterResult, "AdultCCresults.csv") 
write.csv(AdultBP@compareClusterResult, "AdultBPresults.csv") 

AdultBPp <- dotplot(AdultBP, showCategory = 5)
# Open a pdf file
pdf("MF.pdf", width = 9, height = ) 
# 2. Create a plot
MF
# Close the pdf file
dev.off()

# GO analysis of Aged proteome
# compare up and down clusters

Agedproteome <- list("Down" = Ageddown, 
                      "Up" = Agedup)

AgedMF <- compareCluster(geneCluster = Agedproteome, fun = enrichGO,
                          ont = "MF",
                          OrgDb = "org.Rn.eg.db",
                          pvalueCutoff = 2, readable = TRUE)
# ck <- setReadable(ck, OrgDb = org.Rn.eg.db, keyType="ENTREZID")

write.csv(AgedMF@compareClusterResult, "AgedMFresults.csv") 
write.csv(AgedCC@compareClusterResult, "AgedCCresults.csv") 
write.csv(AgedBP@compareClusterResult, "AgedBPresults.csv") 

AgedMFp <- dotplot(AgedMF, showCategory = 5)

enrichGO(gene=Ageddown, OrgDb = org.Rn.eg.db, keyType = "ENTREZID",
         pvalueCutoff = 2, ont = "MF")

# compare Adult vs Aged proteomes GO terms
# for BP, CC, MF

Adult_Aged <- list("Adult" = Adultsig, 
                  "Aged" = Agedsig)

ckBP <- compareCluster(geneCluster = Adult_Aged, fun = enrichGO,
                     ont = "BP",
                     OrgDb = org.Rn.eg.db,
                     pvalueCutoff = 2, readable = TRUE)
# ck <- setReadable(ck, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
head(ckMF, 10) 

BP <- dotplot(ckBP, showCategory = 5)

ckMF <- pairwise_termsim(ckMF)  
emapplot(ckMF)
cp <- cnetplot(ckMF, shadowtext = 'none')

# Open a pdf file
pdf("compareAdultAged_BP.pdf", width = 9, height = ) 
# 2. Create a plot
BP
# Close the pdf file
dev.off() 

###### compare Adult up/down proteomes for KEGG
Adultdown_KEGG <- enrichKEGG(gene = Adultdown, organism = "rno",
                             pvalueCutoff = 2)
Adultdown_KEGG <- setReadable(Adultdown_KEGG, OrgDb = org.Rn.eg.db, keyType="ENTREZID")

AD_KEGGp <-dotplot(Adultdown_KEGG, showCategory = 12)
write.csv(Adultdown_KEGG@result,"AD_KEGGresults.csv")

Adultdowncp <- cnetplot(Adultdown_KEGG, shadowtext = 'none',
                      foldChange = fold_changesAdultdown)

fold_changesAdultdown <- Adult_down$log2ratio
names(fold_changesAdultdown) <-Adult_down$`Gene Symbol`

## upregulated Adult proteins

Adultup_KEGG <- enrichKEGG(gene = Adultup, organism = "rno",
                             pvalueCutoff = 2)
Adultup_KEGG <- setReadable(Adultup_KEGG, OrgDb = org.Rn.eg.db, keyType="ENTREZID")

AU_KEGGp <-dotplot(Adultup_KEGG, showCategory = 10, x = "count")

write.csv(Adultup_KEGG@result,"AU_KEGGresults.csv")
View(Adultup_KEGG@result)

Adultupcp <- cnetplot(Adultup_KEGG, shadowtext = 'none',
                    foldChange = fold_changesAdultup)

fold_changesAdultup <- Adult_up$log2ratio
names(fold_changesAdultup) <-Adult_up$`Gene Symbol`

###### compare Aged up/down proteomes for KEGG
Ageddown_KEGG <- enrichKEGG(gene = Ageddown, organism = "rno",
                             pvalueCutoff = 2)

AgedD_KEGGp <-dotplot(Ageddown_KEGG, showCategory = 10)
write.csv(Ageddown_KEGG@result,"AgedD_KEGGresults.csv")
View(Ageddown_KEGG@result)

Agedup_KEGG <- enrichKEGG(gene = Agedup, organism = "rno",
                           pvalueCutoff = 2)
AgedU_KEGGp <-dotplot(Agedup_KEGG, showCategory = 10)
write.csv(Agedup_KEGG@result,"AgedU_KEGGresults.csv")
View(Agedup_KEGG@result)

####### KEGG for Adult and Aged proteomes
Adult_KEGG <- enrichKEGG(gene = Adultsig, organism = "rno",
                            pvalueCutoff = 2)
Adult_KEGG <- setReadable(Adult_KEGG, OrgDb = org.Rn.eg.db, keyType="ENTREZID")

Adult_KEGGp <-dotplot(Adult_KEGG, showCategory = 10)
write.csv(Adult_KEGG@result,"Adult_KEGG.csv")
View(Adult_KEGG@result)

Adultcp <- cnetplot(Adult_KEGG, shadowtext = 'none',
                   foldChange = fold_changesAdult)

ggsave(Adultcp, file="AdultcpKEGG.pdf", width=10, height=7, dpi=300)


Adult_KEGGep <- pairwise_termsim(Adult_KEGG)  
Adult_KEGGeplot <- emapplot(Adult_KEGGep, shadowtext = FALSE,
                            showCategory =10)

ggsave(Adult_KEGGeplot, file="Adult_KEGGeplot.pdf", width=10, height=7, dpi=300)

cnetplot(Adult_KEGG, shadowtext = "none", showCategory = selected_pathways,
         foldChange = fold_changesAdult)

selected_pathways <- (Adult_KEGG$Description[3])

##### edit cnetplot colors and select specific nodes 1-3
Adult_KEGGcp1 <- cnetplot(Adult_KEGG, shadowtext = 'none', 
                          showCategory = Adult_KEGG$Description[1], 
                          foldChange = fold_changesAdult)
Adult_KEGGcp2 <- cnetplot(Adult_KEGG, shadowtext = 'none', 
                          showCategory = Adult_KEGG$Description[2], 
                          foldChange = fold_changesAdult)
Adult_KEGGcp3 <- cnetplot(Adult_KEGG, shadowtext = 'none', 
                          showCategory = Adult_KEGG$Description[3], 
                          foldChange = fold_changesAdult)

ggsave(Adult_KEGGcp1, file="Adult_KEGGcp1.pdf", width=5, height=5, dpi=300)
ggsave(Adult_KEGGcp2, file="Adult_KEGGcp2.pdf", width=5, height=5, dpi=300)
ggsave(Aged_KEGGcp3, file="Aged_KEGGcp3.pdf", width=5, height=5, dpi=300)

################### Aged
Aged_KEGG <- enrichKEGG(gene = Agedsig, organism = "rno",
                         pvalueCutoff = 2)
Aged_KEGG <- setReadable(Aged_KEGG, OrgDb = org.Rn.eg.db, keyType="ENTREZID")

Aged_KEGGp <-dotplot(Aged_KEGG, showCategory = 10)
write.csv(Aged_KEGG@result,"Aged_KEGG.csv")
View(Aged_KEGG@result)

Aged_KEGGcp <- cnetplot(Aged_KEGG, shadowtext = 'none',
                    foldChange = fold_changesAged)

ggsave(Aged_KEGGcp, file="Aged_KEGGcp.pdf", width=10, height=7, dpi=300)

##### edit cnetplot colors and select specific nodes 1-3
Aged_KEGGcp1 <- cnetplot(Aged_KEGG, shadowtext = 'none', 
                          showCategory = Aged_KEGG$Description[1], 
                          foldChange = fold_changesAged)
Aged_KEGGcp2 <- cnetplot(Aged_KEGG, shadowtext = 'none', 
                          showCategory = Aged_KEGG$Description[5], 
                          foldChange = fold_changesAged)
Aged_KEGGcp3 <- cnetplot(Aged_KEGG, shadowtext = 'none', 
                          showCategory = Aged_KEGG$Description[3], 
                          foldChange = fold_changesAged)

ggsave(Aged_KEGGcp1, file="Aged_KEGGcp1.pdf", width=5, height=5, dpi=300)
ggsave(Aged_KEGGcp2, file="Aged_KEGGcp2.pdf", width=5, height=5, dpi=300)
ggsave(Aged_KEGGcp3, file="Aged_KEGGcp3.pdf", width=5, height=5, dpi=300)

##### KEGG eplot
Aged_KEGGep <- pairwise_termsim(Aged_KEGG)  
Aged_KEGGeplot <- emapplot(Aged_KEGGep, shadowtext = FALSE,
                           showCategory = 10)

ggsave(Aged_KEGGeplot, file="Aged_KEGGeplot.pdf", width=10, height=7, dpi=300)


###### compare Adult vs. Aged KEGG enrichment
KEGG <- compareCluster(geneCluster = Adult_Aged, fun = enrichKEGG,
                       organism = "rat",
                       pvalueCutoff = 2)
KEGG <- setReadable(KEGG, OrgDb = "org.Rn.eg.db", keyType="ENTREZID")
head(KEGG, 10) 

write.csv(KEGG@compareClusterResult,"compareKEGGresults.csv")

KGp <- dotplot(KEGG, showCategory = 5)

KEGG <- pairwise_termsim(KEGG)  
emapplot(KEGG)

cp <- cnetplot(KEGG, shadowtext = 'none')

fold_changesAll <- c(fold_changesAdult,fold_changesAged)
head(fold_changesAll)

# Open a pdf file
pdf("compareAdult_AgedKEGGcp.pdf", width = 9, height = ) 
# 2. Create a plot
cp
# Close the pdf file
dev.off()

###### compare Adult up/down proteomes for Reactome
AdultdownR <- enrichPathway(gene=Adultdown, organism = "rat",
                                pvalueCutoff = 0.05, pAdjustMethod = "BH",
                                readable=TRUE)

AdultdownRp <-dotplot(AdultdownR, showCategory = 10)
write.csv(AdultdownR@result,"AdultdownR.csv")

AdultupR <- enrichPathway(gene=Adultup, organism = "rat",
                              pvalueCutoff = 0.05, pAdjustMethod = "BH",
                              readable=TRUE)
AdultupRp <-dotplot(AdultupR, showCategory = 10)
write.csv(AdultupR@result,"AdultupR.csv")

###### compare Aged up/down proteomes for Reactome
AgeddownR <- enrichPathway(gene=Ageddown, organism = "rat",
                            pvalueCutoff = 0.05, pAdjustMethod = "BH",
                            readable=TRUE)

AgeddownRp <-dotplot(AgeddownR, showCategory = 10)
write.csv(AgeddownR@result,"AgeddownR.csv")

AgedupR <- enrichPathway(gene=Agedup, organism = "rat",
                          pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          readable=TRUE)
AgedupRp <-dotplot(AgedupR, showCategory = 10)
write.csv(AgedupR@result,"AgedupR.csv")
View(AgedupR@result)

###### Reactome for Adult
AdultR<- enrichPathway(gene=Adultsig, organism = "rat",
                      pvalueCutoff = 0.05, pAdjustMethod = "BH",
                      readable=TRUE)
AdultRp <-dotplot(AdultR, showCategory = 10)
write.csv(AdultR@result,"AdultR.csv")
View(AdultR@result)

AdultRcp <- cnetplot(AdultR, shadowtext = 'none',
                   foldChange = fold_changesAdult)

fold_changesAdult <- Adult_sig$log2ratio
names(fold_changesAdult) <- Adult_sig$`Gene Symbol`
head(fold_changesAdult)

ggsave(AdultRcp, file="AdultcpR.pdf", width=10, height=7, dpi=300)


###### Reactome for Aged
AgedR<- enrichPathway(gene=Agedsig, organism = "rat",
                         pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         readable=TRUE)
AgedRp <-dotplot(AgedR, showCategory = 10)
write.csv(AgedR@result,"AgedR.csv")
View(AgedR@result)

Agedcp <- cnetplot(AgedR, shadowtext = 'none',
                   foldChange = fold_changesAged)

fold_changesAged <- Aged_sig$log2ratio
names(fold_changesAged) <- Aged_sig$Gene
head(fold_changesAged)

ggsave(Agedcp, file="AgedcpR.pdf", width=10, height=7, dpi=300)


# compare Reactome between Adult and Aged
cReactome <- compareCluster(geneCluster = Adult_Aged, fun = enrichPathway,
                            organism = "mouse", pvalueCutoff = 1, readable = TRUE)
# no enrichment found in any gene clusters

Reactp <- dotplot(cReactome, showCategory = 5)








