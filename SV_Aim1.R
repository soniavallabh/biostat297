# command enter to execute a line of code.  Code will run in the console.
# use ? before function name to get help
# r is 1-based
# In r data frames, rows are 1st dimension, columns are 2nd
# use dollar signs to specify a column in a data structure e.g. masterTable$gene
# use square brackets to index into a vector or data frame
# ! means not.  for example !is.na means "is not NA"

# set stringsAsFactors = FASLE to keep strings being from being read in as values
options(stringsAsFactors = FALSE)

# set our working directory to where the data is
setwd("~/Desktop/Biostat297_final_project")

# ASSEMBLE DATA FRAME
# Reading in data -- requires iteraction. Default is space-delineated.
iceland = read.table("IcelandKO.txt", skip=3, sep="\t", header = TRUE)
# change column names to all lower-case using tolower
colnames(iceland) = tolower(colnames(iceland))
# table is a built-in function that counts the number of times each
# element occurs in a table or vector.
lofAlleleCount = table(iceland$gene)
# R uses vectors instead of lists. Initialize with c(). for example:
x = c(42, 547, 65, 1, 1)
nIcelanders = 2636
head(iceland$sequencing.maf....*nIcelanders*2)
?aggregate
# ~ means "as a function of"
# create a "combined allele frequency LOF" score for each gene.  We'll call this lofCAF.
lofCAF = aggregate(sequencing.maf.... ~ gene, data = iceland, FUN = sum)
plot(lofCAF$sequencing.maf...., as.numeric(lofAlleleCount))
#read in second dataset, Samocha.txt, that contains number of coding bp per gene for 18,000 genes
samocha = read.table("Samocha.txt", sep="\t", header = TRUE)
#count number of genes in dataset and number of unique genes
length(samocha$gene)
length(unique(samocha$gene))
# "merge" combines two data frames. "all.x" keeps the size of the first one "x")
masterTable = merge(samocha, lofCAF, by = "gene", all.x = TRUE)
colnames(samocha)
colnames(masterTable)
dim(masterTable)
# index into data frames with brackets. Rows are 1st dimension, columns are 2nd.
# Overwrite masterTable data frame with a reduced version of itself containing
# only the columns we want. To do this, create a vector of desired columns.
# If we wanted to, we could do this by column name as well.
masterTable = masterTable[,c(1,3:9,13)]
#make the column name "sequencing.maf...." more friendly using colnames and indexing
colnames(masterTable)[9] = "lofcaf"
# We need to replace all of the "NA" values in the lofcaf column with "0"
# First we use is.na to create a vector of TRUE/FALSE telling us locations of NA
# Then we index into masterTable$lofcaf using this vector, and set the trues = 0.
masterTable$lofcaf[is.na(masterTable$lofcaf)] = 0 
# Now we're going to merge more data into master Table -- the lofAlleleCount table
# First we'll convert it from a table to a data frame
lofAlleleCount = as.data.frame(lofAlleleCount)
colnames(lofAlleleCount)[1] = "gene"
masterTable = merge(masterTable, lofAlleleCount, by = "gene", all.x = TRUE)
colnames(masterTable)[10] = "number_of_alleles"
masterTable$number_of_alleles[is.na(masterTable$number_of_alleles)] = 0

# plot LOFCAF vs. BP

# RELATIONSHIP BETWEEN LOFCAF, NUMBER_OF_ALLELES, AND BP
# Compare combined allele frequency of LOFs versus number of LOF alleles. Not as correlated as one might expect.
plot(masterTable$lofcaf, masterTable$number_of_alleles)
# Compare number of coding bp with number of LOF alleles
plot(masterTable$bp, masterTable$number_of_alleles, log = "x", xlab = "length of coding sequence (bp)", pch = 20, ylab = "number of LOF alleles", main = "Number of LOF alleles vs. length of coding sequence", xaxt = "n")
axis(side = 1, at = 10^(2:5), labels = formatC(10^(2:5), format = "fg", big.mark = ","))
# Check for a correlation -- default is Pearson's correlation.
# Put in parens the same thing that you plotted. Comes out to ~0.25
cor.test(masterTable$bp, masterTable$number_of_alleles)
# spearman uses rank instead of numerical values; less sensitive to outliers. Comes out to ~0.18
cor.test(masterTable$bp, masterTable$number_of_alleles, method = "spearman")
# looking at the y-axis outlier 
#which.max will give you the row number in column number_of_alleles
# use that row number to index into masterTable, returning the whole row
# The gene is SSPO -- it has lots of LOFs depsite not being very long
masterTable[which.max(masterTable$number_of_alleles),]
# x-axis outlier must be Titin because it's super huge.  Let's test this
masterTable[which.max(masterTable$bp),]
# Plots bp vs. lofcaf
plot(masterTable$bp, masterTable$lofcaf)
# Test correlations b/t bp and lofcaf
cor.test(masterTable$bp, masterTable$lofcaf)
cor.test(masterTable$bp, masterTable$lofcaf, method = "spearman")


# RELATIONSHIP BETWEEN LOFCAF, NUMBER_OF_ALLELES AND CONSTRAINT
# Import Kaitlin's list of 1003 most constrained genes.
constraint = read.table("Samocha_constraint.txt", sep="\t", header = TRUE)
constraint = constraint[,c(2,13)]
masterTable = merge(masterTable, constraint, by = "gene", all.x = TRUE)
# use wilcoxon rank sum to compare lofcaf for constrained versus uncontrained genes. 
# wilcoxon rank sum is for unpaired data (signed rank is for paired data)
wilcox.test(masterTable$lofcaf[!is.na(masterTable$mis_z_sign)], masterTable$lofcaf[is.na(masterTable$mis_z_sign)])
# how much does mean lofcaf differ between constrained and non-constrained genes? Comes out to a lot (0.03 vs. 0.45)
mean(masterTable$lofcaf[!is.na(masterTable$mis_z_sign)])
mean(masterTable$lofcaf[is.na(masterTable$mis_z_sign)])
sd(masterTable$lofcaf[!is.na(masterTable$mis_z_sign)])
sd(masterTable$lofcaf[is.na(masterTable$mis_z_sign)])
mean(masterTable$number_of_alleles[!is.na(masterTable$mis_z_sign)])
mean(masterTable$number_of_alleles[is.na(masterTable$mis_z_sign)])
median(masterTable$number_of_alleles[!is.na(masterTable$mis_z_sign)])
median(masterTable$number_of_alleles[is.na(masterTable$mis_z_sign)])
IQR(masterTable$number_of_alleles[!is.na(masterTable$mis_z_sign)])
IQR(masterTable$number_of_alleles[is.na(masterTable$mis_z_sign)])


# plot lofcaf versus contraint, for genes with a constraint score. LOF and missense constraint are not perfectly correlated.
# Look at outliers with high lofcaf despite being constrained.
plot(masterTable$lofcaf[!is.na(masterTable$mis_z_sign)], masterTable$mis_z_sign[!is.na(masterTable$mis_z_sign)], xlab = "LOF combined allele frequency", ylab = "constraint", main = "contraint vs. LOF combined allele frequency", pch = 20)
plot(jitter(masterTable$number_of_alleles[!is.na(masterTable$mis_z_sign)], 0.5), masterTable$mis_z_sign[!is.na(masterTable$mis_z_sign)], xlab = "number of LOF alleles", ylab = "constraint", main = "Constraint vs. number of LOF alleles", pch = 20)
hist(masterTable$lofcaf[!is.na(masterTable$mis_z_sign)], xlab = "LOF allele count", col = "lightblue", ylab = "Number of constained genes", main = "LOF allele counts for genes reported as constrained")

# CALCULATE LOFs/PERSON
sum(iceland$sequencing.maf....)/100
