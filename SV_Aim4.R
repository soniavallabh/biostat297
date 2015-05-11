options(stringsAsFactors = FALSE)

setwd("~/Desktop/Biostat297_final_project")

## READ IN DATA -- full gene list from Samocha 2014, icelandic LOF list, and MacArthur Lab gene lists (autosomal dominant disease
## genes Berg 2013, autosomal recessive disease genes Berg 2013, drug targets Nelson 2012, and FDA-approved drug targets)
Samocha = read.table("Samocha.txt", sep = "\t", header = TRUE)
fisherTable = Samocha
FDAdrugList = read.table("fda_approved_drug_targets.tsv", sep = "\t", header = FALSE)
drugList = read.table("drug_targets_nelson.tsv", sep = "\t", header = FALSE)
domList = read.table("berg_dominant.tsv", sep = "\t", header = FALSE)
recList = read.table("berg_recessive.tsv", sep = "\t", header = FALSE)
iceland = read.table("IcelandKO.txt", skip=3, sep="\t", header = TRUE)
colnames(iceland) = tolower(colnames(iceland))
colnames(iceland)[15] = "imputed"

# FDA DRUG LIST GENES ARE NOT SIGNIFICANTLY ENRICHED IN ICELAND HOMOZYGOUS KNOCKOUT GENES
fisherTable$FDAdrug = fisherTable$gene %in% FDAdrugList$V1
fisherTable$homLOF = fisherTable$gene %in% iceland$gene[iceland$imputed > 0]
table(fisherTable$FDAdrug)  
table(fisherTable[,c("FDAdrug", "homLOF")])

fisher.test(table(fisherTable[,c("FDAdrug", "homLOF")]))
# p = 0.4314, odds ratio = 1.168

# # DRUG LIST GENES ARE NOT SIGNIFICANTLY ENRICHED IN ICELAND HOMOZYGOUS KNOCKOUT GENES
fisherTable$drugList = fisherTable$gene %in% drugList$V1
fisherTable$homLOF = fisherTable$gene %in% iceland$gene[iceland$imputed > 0]
table(fisherTable$drugList)  
table(fisherTable[,c("drugList", "homLOF")])

fisher.test(table(fisherTable[,c("drugList", "homLOF")]))
# p = 0.345, odds ratio = 0.706

# # DOMINANT DISEASE GENES ARE SIGNIFICANTLY DEPLETED IN ICELAND HOMOZYGOUS KNOCKOUT GENES
fisherTable$dominant = fisherTable$gene %in% domList$V1
fisherTable$homLOF = fisherTable$gene %in% iceland$gene[iceland$imputed > 0]
table(fisherTable$dominant)  
table(fisherTable[,c("dominant", "homLOF")])

fisher.test(table(fisherTable[,c("dominant", "homLOF")]))
# p = 0.000172, odds ratio = 0.49

# # RECESSIVE DISEASE GENES ARE NOT SIGNIFICANTLY DEPLETED IN ICELAND HOMOZYGOUS KNOCKOUT GENES
fisherTable$recessive = fisherTable$gene %in% recList$V1
fisherTable$homLOF = fisherTable$gene %in% iceland$gene[iceland$imputed > 0]
table(fisherTable$recessive)  
table(fisherTable[,c("recessive", "homLOF")])

fisher.test(table(fisherTable[,c("recessive", "homLOF")]))
# p = 0.86, odds ratio = 1.02
