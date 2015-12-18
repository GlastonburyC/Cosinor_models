library(data.table)
library(lme4)
library(gtools)

file<-as.numeric(Sys.getenv("SGE_TASK_ID"))
f_path=paste("EB_F_scaled_expr_",file,".txt",sep="")
rna_seq=read.table(f_path,head=T,sep="\t",stringsAsFactors=F)
eigen_vector<-read.table("eigen.eigenU.txt",head=F)
eigen_vector=eigen_vector[,716:720]

not_diabetic=read.table("~/not_diabetic.IDs",head=T,sep="\t",stringsAsFactors=F)
covs=read.table("qc_F_freezev1.txt",head=T,sep="\t",stringsAsFactors=F)
row.names(eigen_vector)=covs$SampleID
eigen_vector=eigen_vector[row.names(eigen_vector) %in% not_diabetic[,2],]
not_diabetic$BiopsyDate=as.POSIXct(not_diabetic$BiopsyDate,format="%d/%m/%Y")
rna_seq_exon=rna_seq[,1]


rna_seq=rna_seq[,colnames(rna_seq) %in% covs[,1]]
not_diabetic=not_diabetic[not_diabetic[,2] %in% colnames(rna_seq),]
rna_seq=rna_seq[,colnames(rna_seq) %in% not_diabetic[,2]]
covs=covs[covs[,1] %in% colnames(rna_seq),]

index=NULL
rna_seq=rna_seq[,mixedsort(colnames(rna_seq))]
index=mixedorder(not_diabetic$EBTwinID)
not_diabetic=not_diabetic[index,]
index=NULL
index=mixedorder(covs$SampleID)
covs=covs[index,]


  family <- as.factor(as.matrix(covs['Family']))
  zygosity <- as.factor(as.matrix(covs['Zygosity']))
  INSERT_SIZE_MODE <- as.numeric(as.matrix(covs['INSERT_SIZE_MODE']))
  GC_mean <- as.numeric(as.matrix(covs['GC_mean']))
  PrimerIndex <- as.factor(as.matrix(covs['PrimerIndex']))
  age <- as.numeric(as.matrix(covs['AGE']))
# batch <- as.factor(as.matrix(covs['Set']))
  BMI <- as.numeric(as.matrix(covs['BMI']))

t=NULL
day_of_year=NULL
results=NULL
gene=NULL
CosTime=NULL
SinTime=NULL
# Calculate the day of the year i.e. Jan 1st =  1.
not_diabetic$BiopsyDate=as.POSIXct(not_diabetic$BiopsyDate,format="%d/%m/%Y")
day_of_year=as.numeric(strftime(not_diabetic$BiopsyDate,format="%j"))


# Divide by number of days in the year to get t.
t=day_of_year/365

#Transform t into a cosinor relationship

CosTime=cos(2*pi*t)
SinTime=sin(2*pi*t)

# This allows you to then fit the linear model: expression ~ CosTime + SinTime + fixed_effect + ( 1 | random_effects). Compare to Null.
# Intercept = MESOR - the mean around which the expression of gene[i] oscillates
# Ampitude can be calculated from The coefficients of SinTime (B1) + CosTime (B2) i.e. sqrt(B1^2+B2^2)
# Acrophase (peak expression level in cycle) can be calulated by computing arctan(-B2/B1)

cosinor_lm = function(x) {
for(i in 1:dim(rna_seq)[1]){

	gene=as.matrix(as.numeric(rna_seq[i,1:length(rna_seq)]))

full=lmer(gene ~ scale(CosTime) + scale(SinTime) + scale(age)+ scale(INSERT_SIZE_MODE)+ scale(GC_mean) + (1 | family) + (1 | zygosity)  + scale(as.matrix(eigen_vector[,1]))+ scale(as.matrix(eigen_vector[,2])) + scale(as.matrix(eigen_vector[,3]))+ scale(as.matrix(eigen_vector[,4]))+ scale(as.matrix(eigen_vector[,5])) + (1|PrimerIndex),REML=FALSE)
null=lmer(gene ~ scale(age) + scale(INSERT_SIZE_MODE) + scale(GC_mean) + (1| PrimerIndex) + (1 | family) + (1 | zygosity) + scale(as.matrix(eigen_vector[,1]))+ scale(as.matrix(eigen_vector[,2])) + scale(as.matrix(eigen_vector[,3]))+ scale(as.matrix(eigen_vector[,4]))+ scale(as.matrix(eigen_vector[,5])),REML=FALSE)


full_bmi=lmer(gene ~ scale(CosTime) + scale(SinTime) + scale(age) + scale(BMI)+ scale(INSERT_SIZE_MODE) + scale(GC_mean) + (1|PrimerIndex)+ (1 | family) + (1 | zygosity) + scale(as.matrix(eigen_vector[,1]))+ scale(as.matrix(eigen_vector[,2])) + scale(as.matrix(eigen_vector[,3]))+ scale(as.matrix(eigen_vector[,4]))+ scale(as.matrix(eigen_vector[,5])),REML=FALSE)
null_bmi=lmer(gene ~ scale(age) + scale(BMI) + scale(INSERT_SIZE_MODE) + scale(GC_mean) + (1|PrimerIndex)+ (1 | family) + (1 | zygosity)+ scale(as.matrix(eigen_vector[,1]))+ scale(as.matrix(eigen_vector[,2])) + scale(as.matrix(eigen_vector[,3]))+ scale(as.matrix(eigen_vector[,4]))+ scale(as.matrix(eigen_vector[,5])),REML=FALSE)


pval_BMI=anova(full_bmi,null_bmi)$Pr[2]

pval=anova(full,null)$Pr[2]

#results=rbind(results,data.frame(gene=rna_seq_exon[i],mesor=fixef(full)[1],CosTime=fixef(full)[2],SinTime=fixef(full)[3],age=fixef(full)[4],p_value=pval,p_value_BMIadj=pval_BMI))
cat(rna_seq_exon[i],'\t',fixef(full)[1],'\t',fixef(full)[2],'\t',fixef(full)[3],'\t',fixef(full)[4],fixef(full_bmi)[3],'\t',pval,'\t',pval_BMI,'\n')
full=null=NULL

}
}

sink(file=paste("/home/glastonc/cosinor/Season_F_expression_",file,".tab",sep=""))
if(file==1){cat('Exon','\t','MESOR','\t','CosTime','\t','SinTime','\t','age','\t','BMI','\t','p_value','\t','BMI_adj_pvalue','\n')}
cosinor_lm()
sink()
