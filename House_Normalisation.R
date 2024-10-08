##Date: 8/10/2024
##Author; Sarah Grant
##Description: Calculating the mean of 100 genes with lowest variance across all samples 

exprs_var<-apply(assay(dds_TA),1,var)
exprs_mean<-apply(assay(dds_TA),1,mean)

hist(exprs_var)
hist(exprs_mean, breaks=50)

dds_TA

colMeans(raw_counts[rownames(raw_counts) %in% names(exprs_var[order(exprs_var)][1:100]),])
