library(maftools)
library(grDevices)
library(VennDiagram)
#download the maf file from TCGA into your working directory
laml = read.maf(maf = skcm.maf)
#Exploratory analysis of the maf file
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#Comparing mutation load aross TCGA cohort
laml.mutload = tcgaCompare(maf = laml, cohortName = 'TCGA-SKCM')
#Finding cancer driver genes using OncodriveCLUST
laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
#Known oncogenic pathways
OncogenicPathways(maf = laml)
#intersection between the output of the three tools skcm_candra, skcm_mut,skcm_onco and their overlap with the CGC
overlap_3=intersect(skcm_candra$HGNC,intersect(skcm_mut$gene,skcm_onco$Hugo_Symbol))
overlap_CGC=intersect(overlap_3,census_GRCH37_v90$Gene.name
#Venn diagram
temp <- venn.diagram(list(MutSigCV=skcm_mut$gene, CanDrA=skcm_candra$HGNC,OncoCLUST=skcm_onco$Hugo_Symbol),fill = c("red", "green","blue"),cex = 2,cat.fontface = 4,filename=NULL)
pdf(file="skcm.pdf")
grid.draw(temp)
dev.off()

