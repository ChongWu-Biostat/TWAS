library(data.table)
library(matlib)
library(Rcpp)
library(RcppArmadillo)
library(bigmemory)
library(mvtnorm)
library(MASS)

suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))
suppressMessages(library(aSPU2))
source("dist_support.R")
# Some codes are directly copied from TWAS source codes
option_list = list(
make_option("--sumstats", action="store", default=NA, type='character',
help="summary statistics (rds file and must have SNP and Z column headers) [required]"),
make_option("--out", action="store", default=NA, type='character',
help="Path to output files [required]"),
make_option("--weights", action="store", default=NA, type='character',
help="File listing molecular weight (rds files and must have columns WGT,ID,CHR,P0,P1) [required]"),
make_option("--weights_dir", action="store", default=NA, type='character',
help="Path to directory where weight files (WGT column) are stored [required]"),
make_option("--ref_ld", action="store", default=NA, type='character',
help="Reference LD files in binary PLINK format [required]"),
make_option("--gene_list", action="store", default=NA, type='character',
help="Gene list we want to analyze, currently only single chromosome analyses are performed [required]"),
make_option("--chr", action="store", default=NA, type='character',
help="Chromosome to analyze, currently only single chromosome analyses are performed [required]"),
make_option("--max_nperm", action="store", default= 1000000, type='integer',
help="maximum number of permutation for aSPU or daSPU. [default: 1000000]"),
make_option("--force_model", action="store", default=NA, type='character',
help="Force specific predictive model to be used, no flag (default) means select most significant cross-val. Options: blup,lasso,top1,enet")
)

opt = parse_args(OptionParser(option_list=option_list))


# This function (allele.qc) is downloaded from TWAS (http://gusevlab.org/projects/fusion/#typical-analysis-and-output)
allele.qc = function(a1,a2,ref1,ref2) {
    ref = ref1
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip1 = flip
    
    ref = ref2
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip2 = flip;
    
    snp = list()
    snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    
    return(snp)
}

sumstat = readRDS(opt$sumstats)

wgtlist = read.table(opt$weights,head=T,as.is=T)
wgtlist = wgtlist[as.character(wgtlist$CHR) == as.character(opt$chr) , ]
chr = unique(wgtlist$CHR)


gene.list = read.table(opt$gene_list)
gene.list[,1] = as.character(gene.list[,1])
gene.list = gene.list[,1]

wgtlist = wgtlist[wgtlist[,"ID"] %in% gene.list,]

genos = read_plink(paste(opt$ref_ld,chr,sep=''),impute="avg")

# Match summary data to input, record NA where summary data is missing
m = match( genos$bim[,2] , sumstat$SNP )
sum.missing = is.na(m)
sumstat = sumstat[m,]
sumstat$SNP = genos$bim[,2]
sumstat$A1[ sum.missing ] = genos$bim[sum.missing,5]
sumstat$A2[ sum.missing ] = genos$bim[sum.missing,6]

# QC / allele-flip the input and output
qc = allele.qc( sumstat$A1 , sumstat$A2 , genos$bim[,5] , genos$bim[,6] )

# Flip Z-scores for mismatching alleles
sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
sumstat$A1[ qc$flip ] = genos$bim[qc$flip,5]
sumstat$A2[ qc$flip ] = genos$bim[qc$flip,6]

# Remove strand ambiguous SNPs (if any)
if ( sum(!qc$keep) > 0 ) {
    genos$bim = genos$bim[qc$keep,]
    genos$bed = genos$bed[,qc$keep]
    sumstat = sumstat[qc$keep,]
}

out.res = as.data.frame(matrix(NA,nrow(wgtlist),15))
colnames(out.res) = c("gene","CHR","P0","P1","#nonzero_SNPs","TWAS_asy","SSU_asy","SPU(1)","SPU(2)","SPU(3)","SPU(4)","SPU(5)","SPU(6)","SPU(Inf)","aSPU")

## For each wgt file:
for ( w in 1:nrow(wgtlist) ) {
    tryCatch({
        #cat( unlist(wgtlist[w,]) , '\n' )
        # Load weights
        wgt.file = paste(opt$weights_dir,"/",wgtlist$WGT[w],sep='')
        load(wgt.file)
        # Remove NAs (these should not be here)
        wgt.matrix[is.na(wgt.matrix)] = 0
        
        # Match up the SNPs and weights
        m = match( snps[,2] , genos$bim[,2] )
        m.keep = !is.na(m)
        snps = snps[m.keep,]
        wgt.matrix = wgt.matrix[m.keep,]
        cur.genos = scale(genos$bed[,m[m.keep]])
        cur.bim = genos$bim[m[m.keep],]
        # Flip WEIGHTS for mismatching alleles
        qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
        wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
        rm(snps)
        
        cur.FAIL = FALSE
        
        # Match up the SNPs and the summary stats
        m = match(cur.bim[,2] , sumstat$SNP)
        cur.Z = sumstat$Z[m]
        
        # Identify the best model
        if ( !is.na(opt$force_model) ) {
            mod.best = which( colnames(wgt.matrix) == opt$force_model )
            if ( length(mod.best) == 0 ) {
                cat( "WARNING : --force_model" , mod.best ,"does not exist for", unlist(wgtlist[w,]) , "\n")
                cur.FAIL = TRUE
            }
        } else {
            mod.best = (which.max(cv.performance[1,]))
            if ( names(mod.best) == "top1" ) {
                # cat( "WARNING: top eQTL is the best predictor for this gene, continuing with 2nd-best model\n" )
                mod.best = names( which.max(cv.performance[1,colnames(cv.performance)!="top1"]) )
                mod.best = which( colnames(cv.performance) == mod.best )
            }
        }
        
        if ( sum(wgt.matrix) == 0 ) {
            cat( "WARNING : " , unlist(wgtlist[w,]) , "had", length(cur.Z) , "overlapping SNPs, but non with non-zero expression weights, try more SNPS or a different model\n")
            cur.FAIL = TRUE
        }
        
        # Compute LD matrix
        if ( length(cur.Z) == 0 ) {
            cat( "WARNING : " , unlist(wgtlist[w,]) , " had no overlapping SNPs\n")
            cur.FAIL = TRUE
        } else {
            cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
            cur.miss = is.na(cur.Z)
            # Impute missing Z-scores
            if ( sum(cur.miss) != 0 ) {
                if ( sum(!cur.miss) == 0 ) {
                    cat( "WARNING : " , unlist(wgtlist[w,]) , " had no overlapping GWAS Z-scores\n")
                    cur.FAIL = TRUE
                } else {
                    cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
                    cur.impz = cur.wgt %*% cur.Z[!cur.miss]
                    cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
                    cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
                }
            }
            
            if ( !cur.FAIL ) {
                non_zero_ind <- (wgt.matrix[,mod.best] !=0)
                diag_element <- as.vector(wgt.matrix[non_zero_ind,mod.best])
                diag_sd <- diag_element/sum(abs(diag_element))
                weight_diag <- diag(diag_sd,nrow = length(diag_sd))
                
                Zstat.w <- weight_diag %*% cur.Z[non_zero_ind]
                corSNP.w <- weight_diag %*% cur.LD[non_zero_ind, non_zero_ind] %*% t(weight_diag)
                
                U = cur.Z
                U = as.matrix(U)
                
                V = cur.LD
                weight = wgt.matrix[,mod.best]
                
                name = rownames(V)
                rownames(U) = name
                ### Sum, SSU and UminP tests
                pSum = Sum(U=Zstat.w, CovS=corSNP.w)
                pSSU = SumSqU(U=Zstat.w, CovS=corSNP.w)
                
                res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e3)
                
                if(min(res$pvs) < 5e-3 & opt$max_nperm >=1e4) {
                    res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e4)
                }
                
                if(min(res$pvs) < 5e-4  & opt$max_nperm >=1e5) {
                    res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e5)
                }
                
                if(min(res$pvs) < 5e-5  & opt$max_nperm >=1e6) {
                    res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e6)
                }
                
                if(min(res$pvs) < 5e-6  & opt$max_nperm >=1e7) {
                    res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e7)
                }
                
                if(min(res$pvs) < 5e-7  & opt$max_nperm >=1e8) {
                    res = aSPU(U,V,weight,pow=c(1:6,Inf),n.perm = 1e8)
                }
                
                out.res[w,] = c(as.character(wgtlist[w,"ID"]),wgtlist[w,"CHR"],wgtlist[w,"P0"],wgtlist[w,"P1"],sum(weight!=0),pSum,pSSU,res$pvs)
            }
        }
        print(w)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

saveRDS(out.res,opt$out)

write.table(out.res, "output.txt")

