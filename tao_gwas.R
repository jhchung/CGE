
setwd("~/GitHub/cge_final/cge_final")
# Read ped file with genotype data
ped <- read.csv("./questiondat_tao.csv", header=TRUE)

## check missing along columns. Use inline function to calculate the missingness
## for each SNP
missing <- apply(ped,2,function(x) sum(is.na(x))/length(x))

## remove SNPs with missing rate > 10%
missing <- subset(missing, (missing > 0.1))

ped <- ped[,(!names(ped) %in% names(missing))]

### Create a function to performe Hardy-Weinberg test
HWtest <- function(x){
  # Remove NA values from input
  xnna=x[!is.na(x)]
  # Calculate some ratio... what is it?
  pbar<-sum(xnna)/(length(xnna)*2)                    
  E=c((1-pbar)^2,2*pbar*(1-pbar),pbar^2)*length(xnna)
  # get a count of the different genotypes
  O=table(xnna)
  t=sum((O-E)^2/E)
  p=1-pchisq(t,df=1)
  return(c(pbar,p))
}

# Calculate HWE pbar and pvalue for each SNP
alleles <- ped[,(2:dim(ped)[2])]
hwe_test <- apply(alleles,2,function(x)HWtest(x))    ### for cases and controls

## HWE RESULTS: snp1 and snp2 fail HWE test. It seems that control samples are
## severely out of HWE for both while case only snp1 is out of HWE. This could
## be due to genotype errors or rare alleles or that snp1 and snp2 are the
## significant SNPs

###test association
# Add additive model genotypes to the table
testdat <- as.data.frame(ped)
# Convert phenotype to a factor
testdat$PHENOTY <- factor(testdat$PHENOTY)
# Use a general linearized model with binomial distribution to test
total <- dim(testdat)[1]
pb <- txtProgressBar(min=0,max=total, style=3)
glm_tests <- list()
for (i in 2:length(names(testdat))){
  setTxtProgressBar(pb, i)
  snp_glm <- summary(glm(as.formula(paste('PHENOTY~',
                                          names(testdat)[i])), 
                 family = 'binomial', data = testdat))$coef
  glm_tests[[snp]] <- snp_glm
  i = i+1
}
close(pb)


glm_pval_dat <- as.data.frame(matrix(nrow = length(glm_tests), ncol = 2))
names(glm_pval_dat) <- c('snp', 'pvalue')
for (i in 1:length(glm_tests)){
  glm_pval_dat[i,'snp'] <- names(glm_tests)[i]
  if (dim(glm_tests[[i]])[1] == 1 ){
    glm_pval_dat[i,'pvalue'] <- NA
  } else {
    glm_pval_dat[i,'pvalue'] <- glm_tests[[i]][2,4]
  }
}

# Remove NA values
glm_pval_dat <- glm_pval_dat[complete.cases(glm_pval_dat),]

top_snp <- subset(glm_pval_dat, pvalue == min(glm_pval_dat$pvalue))
# Top SNP is rs661891 in gene PRKCQ

# annotate SNP locations
annot <- read.table('illumina550_chr10.txt', header = TRUE)
snp_pval <- merge(annot, glm_pval_dat, by.x = 'name', by.y = 'snp')

# Gene where top SNP is located
gene <- data.frame(gene = 'PKRCQ',
          chr = 'chr10',
          txStart = 6469105,
          txEnd = 6622254)

## extract SNPs within gene
gene_snps <- subset(snp_pval, 
                    (chromStart >= (gene$txStart-10000) &
                      chromEnd <= (gene$txEnd + 10000)))

# Combine pvalues using Fisher's test which assumes the combined pvalues
# follow a chi-square distribution with 2*n degrees of freedom
t_fisher <- -2*sum(log(gene_snps$pvalue))
gene_chisq <- pchisq(t_fisher, 2*dim(gene_snps)[1],lower.tail = FALSE)
