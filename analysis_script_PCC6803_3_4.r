



ascomb<-function(out,outFC){

genes<-outFC[[16]]$genes
dfeat<-outFC[[14]]


as<-grep("as",genes$SystematicName)
int<-grep("int",genes$SystematicName)
UTR<-grep("5'UTR",genes$SystematicName)
us<-grep("us",genes$SystematicName)
nc<-grep("NC",genes$SystematicName,ignore.case=TRUE)
crispr<-grep("Spacer",genes$SystematicName)
mRNA<-seq(1:length(genes$SystematicName))
mRNA<-mRNA[-c(as,int,us,nc,UTR,crispr)]

r<-list()

for(i in 1:length(as)){

n<-grep(as.character(genes$Name[as[i]]),as.character(genes$Name[mRNA]))
if (length(n) >1){

}
if(length(n)>0){

r[[i]]<-n
}
else{

r[[i]]<-NA
}
}

asname<-c()
mname<-c()
anno<-c()

for (i in 1:length(r)){
 for(j in 1: length(r[[i]])){
  if(is.na(r[[i]][j]) == FALSE){
 asname<-c(asname,as.character(genes$SystematicName[as[i]]))
 mname<-c(mname,as.character(genes$SystematicName[mRNA[r[[i]][1]]]))
  anno<-c(anno,as.character(genes$Genenames[as[i]]))
     
}
else{
asname<-c(asname,as.character(genes$SystematicName[as[i]]))
 mname<-c(mname,NA)
  anno<-c(anno,as.character(genes$Genenames[as[i]]))
}
}
}

anno<-cbind(asname,mname,anno)

out<-list(anno, as, mRNA, r)

as_sig<-matrix(,nrow(anno),ncol(dfeat))
colnames(as_sig)<-colnames(dfeat)
for(i in 1:ncol(dfeat)){
as_sig[,i]<-assig2(outFC, out, i)
}
asoverview<-matrix(colSums(as_sig,na.rm = TRUE),1,ncol(as_sig))
colnames(asoverview)<-colnames(dfeat)

as_sig<-cbind(as_sig,anno)
asout<-list(as_sig,asoverview)
asout
}

assig2<-function(outFC, outanno, count){

anno<-outanno[[1]]
as<-outanno[[2]]
mRNA<-outanno[[3]]
r<-outanno[[4]]
genes<-outFC[[16]]$genes
dfeat<-outFC[[14]]

assig<-c()


for (i in 1:length(r)){
 for(j in 1: length(r[[i]])){
  if(is.na(r[[i]][j]) == FALSE){
  
  if(dfeat[as[i],count] !=0 & dfeat[mRNA[r[[i]][j]],count] !=0){
  
  assig<-c(assig,1)
    }
else{
assig<-c(assig,0)}
	
}


else{
assig<-c(assig,NA)
}
}
}


assig
}






tableout<-function(x){

genes<-x[[16]]$genes
dfeat<-x[[9]]
dfeatcl<-x[[14]]

pars<-matrix(,8,1)
colnames(pars)<-"Parameter"
cloess<-matrix(,8,ncol(dfeat))
quant<-matrix(,8,ncol(dfeat))
colnames(cloess)<-paste(colnames(dfeat),"cyclicloess")
colnames(quant)<-paste(colnames(dfeat),"quantile")
rownames(cloess)<-c("total","mRNA","5'UTR","asRNA","ncRNA","intern","unsure","CRISPR_Spacer")
rownames(quant)<-c("total","mRNA","5'UTR","asRNA","ncRNA","intern","unsure","CRISPR_Spacer")

as<-grep("as",genes$SystematicName)
int<-grep("int",genes$SystematicName)
UTR<-grep("5'UTR",genes$SystematicName)
us<-grep("us",genes$SystematicName)
nc<-grep("NC",genes$SystematicName,ignore.case=TRUE)
crispr<-grep("Spacer",genes$SystematicName)
mRNA<-seq(1:length(genes$SystematicName))
mRNA<-mRNA[-c(as,int,us,nc,UTR,crispr)]

for(i in 1:ncol(dfeat)){
cloess[1,i]<-length(which(dfeatcl[,i]!=0))
cloess[2,i]<-length(which(dfeatcl[mRNA,i]!=0))
cloess[3,i]<-length(which(dfeatcl[UTR,i]!=0))
cloess[4,i]<-length(which(dfeatcl[as,i]!=0))
cloess[5,i]<-length(which(dfeatcl[nc,i]!=0))
cloess[6,i]<-length(which(dfeatcl[int,i]!=0))
cloess[7,i]<-length(which(dfeatcl[us,i]!=0))
cloess[8,i]<-length(which(dfeatcl[crispr,i]!=0))

quant[1,i]<-length(which(dfeat[,i]!=0))
quant[2,i]<-length(which(dfeat[mRNA,i]!=0))
quant[3,i]<-length(which(dfeat[UTR,i]!=0))
quant[4,i]<-length(which(dfeat[as,i]!=0))
quant[5,i]<-length(which(dfeat[nc,i]!=0))
quant[6,i]<-length(which(dfeat[int,i]!=0))
quant[7,i]<-length(which(dfeat[us,i]!=0))
quant[8,i]<-length(which(dfeat[crispr,i]!=0))
}

pars[1,1]<-paste("LogFC",x[[17]], sep="=")
pars[2,1]<-paste("adj_pValue",x[[18]], sep="=")
pars[3,1]<-paste("Method",x[[19]], sep="=")
#pars<-list(LogFC=x[[17]],adj_pValue=x[[18]],method=x[[19]])
final<-cbind(quant,cloess,pars)
final
}
#out<-list(RGf,RG,fitprob,fitfeat,targets,G3,kprob,kfeat,G2prob,G2feat)



createFC<-function(x, design,lev){

kprob<-x[[8]]
kfeat<-x[[9]]
G2prob<-x[[14]]
G2feat<-x[[15]]
G2probcl<-x[[16]]
G2featcl<-x[[17]]

fitprob <- lmFit(G2prob, design)
fitfeat <- lmFit(G2feat, design)

fitprobcl <- lmFit(G2probcl, design)
fitfeatcl <- lmFit(G2featcl, design)

fitprob$genes<-kprob$genes
fitfeat$genes<-kfeat$genes

fitprobcl$genes<-kprob$genes
fitfeatcl$genes<-kfeat$genes

#print(design)
con<-matrix(,length(lev),length(lev))
rownames(con)<-lev
colnames(con)<-lev
for(i in 1 :length(lev)){
   for(j in 1 :length(lev)){
con[i,j]<-paste(colnames(con)[i],rownames(con)[j], sep="-")}} 

con1<-select.list(sort(con), multiple=TRUE)

con2<-matrix(,length(con1),length(con1))
rownames(con2)<-paste("(",con1,")", sep="")
colnames(con2)<-paste("(",con1,")", sep="")
for(i in 1 :length(con1)){
   for(j in 1 :length(con1)){
con2[i,j]<-paste(colnames(con2)[i],rownames(con2)[j], sep="-")}} 

con3<-select.list(sort(con2), multiple=TRUE)

con4<-c(con1,con3)

cont.wt <- makeContrasts(contrasts=con4,levels=design)


fitprob2 <- contrasts.fit(fitprob, cont.wt)
fitprob2 <- eBayes(fitprob2)

fitfeat2 <- contrasts.fit(fitfeat, cont.wt)
fitfeat2 <- eBayes(fitfeat2)

prob<-topTable(fitprob2, sort="none",n=Inf)
feat<-topTable(fitfeat2, sort="none",n=Inf)

p.value<-select.list(as.character(seq(0.01,0.05, by=0.01)), title="significant_p.Value")
p.value<-as.numeric(p.value)

lfc<-select.list(as.character(seq(0.5,3, by=0.1)), title="significant_log2_foldchange")
lfc<-as.numeric(lfc)  

met<-select.list(c("separate","global"), title="decideTests_method (see manual)")

fitprob2cl <- contrasts.fit(fitprobcl, cont.wt)
fitprob2cl <- eBayes(fitprob2cl)

fitfeat2cl <- contrasts.fit(fitfeatcl, cont.wt)
fitfeat2cl <- eBayes(fitfeat2cl)

probcl<-topTable(fitprob2cl, sort="none",n=Inf)
featcl<-topTable(fitfeat2cl, sort="none",n=Inf)

dprob<-decideTests(fitprob2,method=met,adjust.method="BH",p.value=p.value,lfc=lfc)
dfeat<-decideTests(fitfeat2,method=met,adjust.method="BH",p.value=p.value,lfc=lfc)

dprob2<-cbind(dprob, fitprob2$genes)
dprob2<-avesums(dprob,  ID=dprob2$TU.ID)
dprob2<-dprob2[,1:length(con4)]


dprobcl<-decideTests(fitprob2cl,method=met,adjust.method="BH",p.value=p.value,lfc=lfc)
dfeatcl<-decideTests(fitfeat2cl,method=met,adjust.method="BH",p.value=p.value,lfc=lfc)

dprob2cl<-cbind(dprobcl, fitprob2cl$genes)
dprob2cl<-avesums(dprobcl,  ID=dprob2cl$TU.ID)
dprob2cl<-dprob2cl[,1:length(con4)]

prob2<-cbind(prob,dprob[,1:length(con4)])
feat2<-cbind(feat,dfeat[,1:length(con4)])

prob2cl<-cbind(probcl,dprobcl[,1:length(con4)])
feat2cl<-cbind(featcl,dfeatcl[,1:length(con4)])

outFC<-list(prob2, feat2, prob2cl, feat2cl, prob, feat, dprob, dprob2, dfeat, probcl, featcl, dprobcl, dprob2cl, dfeatcl, con4, fitfeat2, lfc, p.value, met)

outFC
}



selectreps<-function(targets){



l<-paste(targets$Cy3, targets$rep)


dismiss<-select.list(l, multiple=TRUE)

dis2<-c()

for(i in 1:length(dismiss)){
dis2<-c(dis2,grep(dismiss[i], l))
}


targets$Cy3[dis2]<-"excluded_replicate"
#targets<-targets[-dis2,]


f <- factor(targets$Cy3)
lev<-levels(f)
design <- model.matrix(~0+f)
colnames(design) <- lev

x<-list(design, dismiss,lev)
x
}


avesums<-function (x, ID = rownames(x)) 
{
    if (is.null(x)) 
        return(NULL)
    x <- as.matrix(x)
    nspots <- nrow(x)
    narrays <- ncol(x)
    ID <- as.character(ID)
    iu <- !duplicated(ID)
    if (mode(x) == "character") 
        return(x[iu, , drop = FALSE])
    t <- ID[iu]
    nprobes <- length(t)
    y <- x[iu, , drop = FALSE]
    for (i in 1:length(t)) y[i, ] <- colSums(x[ID == t[i], , 
        drop = FALSE], na.rm = TRUE)
    y
}


analyze<-function(name="array", s=0, e=10000, offset=6.7, whole=TRUE, makeline=TRUE, linie=TRUE){

n<-paste("results_",name,sep="")
d<-dir()
dd<-grep(n,d)
di<-paste(getwd(),paste("results_",name,sep=""),sep="/")
if(length(dd) < 1){

dir.create(di)
}
count=1
pack<-loadedNamespaces()
q="TRUE"
#if( length(grep("amap", pack))==0){
#library(amap)}

#if( length(grep("R.utils", pack))==0){
#library(R.utils)}
#if( length(grep("R.graphics", pack))==0){
#library(R.graphics)}
#if( length(grep("gplots", pack))==0){
#library(gplots)}
require(amap)
require(scales)
require(limma)
require(Biobase)
require(arrayQualityMetrics)
direc<-list.files()
d<-grep("out.Rdata",direc)
if(length(d)==0){


targets<-readTargets()
targets$Cy3<-make.names(as.character(targets$Cy3))

b<-sort.list(targets$Cy3)
tar<-targets
for(i in 1:length(b)){
tar[i,]<-targets[b[i],]}
targets<-tar


#anno<-read.csv("annofile.txt", sep="\t")
load("anno.Rdata")
nas<-which(is.na(anno$SystematicName))
nas2<-seq(1,length(nas))
levels(anno$SystematicName) <- c(levels(anno$SystematicName),nas2)
anno$SystematicName[nas]<-nas2

a<-preprocess(targets,anno)
quality(a,di=di)
}



else{
load("out.Rdata")
a<-out

}


#load("out.Rdata")
#a<-out
#print(is(a))


#print(targets)



#if(length(a) == 1){
#a<-preprocess(targets,anno)

 #quality(a)}

 q<-select.list(c("Genomeplot","Arrayauswertung"), multiple=F, title="Anwendung")
 if(q == "Genomeplot"){
 print(offset)
 
  drawgenomeplot(a,name=name, s=s, e=e, whole=whole, offset=offset, linie=linie, di=di)}
 else{
 por=TRUE
analyze2(a, name=name, di=di)
while(por==TRUE){
por<-select.list(c("TRUE","FALSE"), multiple=F, title="another comparison?")
if(q==TRUE){
count<-count+1
name2<-paste(name, count, sep="-")
analyze2(a, name=name2, di=di)
}}}
#a
}


#--------------------------------------------------------------
preprocess<-function(x=targets, y=anno){


targets=x
b<-sort.list(targets$Cy3)
tar<-targets
for(i in 1:length(b)){
tar[i,]<-targets[b[i],]}
targets<-tar

anno=y


RG<- read.maimages(targets, source="agilent", green.only=TRUE)

new_names<-paste(targets$Cy3, targets$rep, sep=".-.")
colnames(RG$E)<-new_names
colnames(RG$Eb)<-new_names

RG$genes<-anno

RG_normexp<-backgroundCorrect(RG, method="normexp",normexp.method="mle", offset=50)

sub1<-list(RG, RG_normexp)
sub2<-list()
count<-0

weights_vect<-rep(1, nrow(RG$E))

wei<-select.list(c("FALSE", "TRUE"), multiple=F, title="use_weights_for_cyclicloess?")

if(wei=="TRUE"){

load("weights_vect.Rdata")
}


for(i in 1:2){
temp1<- normalizeBetweenArrays(sub1[[i]], method="quantile")
sub2[[i+count]]<-temp1
count<-count+1
temp2<- normalizeBetweenArrays(sub1[[i]], method="cyclicloess", weights=weights_vect, cyclic.method="affy")
sub2[[i+count]]<-temp2

}
sub2[[5]]<-normalizeBetweenArrays(RG, method="none")
names(sub2)<-c("quantile","cyclicloess","normexp_quantile","normexp_cyclicloess","raw")

sub3<-sub2
selec<-select.list(c("TRUE","FALSE"),title = "exclude non-expressed probes?")
for(i in 1:5){
neg95 <- apply(sub2[[i]]$E[sub2[[i]]$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95))

if(selec==TRUE){
cutoff <- matrix(1.05*neg95,nrow(sub2[[i]]),ncol(sub2[[i]]),byrow=TRUE)
isexpr <- rowSums(sub2[[i]]$E > cutoff) >= 2#ncol(sub2[[i]])/length(unique(targets$Cy3))
sub3[[i]] <- sub2[[i]][sub2[[i]]$genes$ControlType==0 & isexpr,]
}
sub3[[i]] <- sub2[[i]][sub2[[i]]$genes$ControlType==0,]
}


sub_probes<-sub3
sub_features<-sub3
sub_TU<-sub3



for(i in 1:5){
sub_probes[[i]]<-avereps(sub_probes[[i]], ID=sub_probes[[i]]$genes$ProbeName)
sub_features[[i]]<-avereps(sub_probes[[i]], ID=sub_probes[[i]]$genes$SystematicName)
sub_TU[[i]]<-avereps(sub_probes[[i]], ID=sub_probes[[i]]$genes$TU.ID)
}




###


out<-list(targets=targets, raw=RG, all_normalized=sub2, all_probes_excluded=sub3, probes=sub_probes, features=sub_features, TUs=sub_TU)



geo<-cbind(RG$genes$FeatureNum,(RG$E))

write.table(geo, file="file_for_geo_submission_raw.txt", sep="\t")
geo<-cbind(RG$genes$FeatureNum,out[[3]][[1]]$E)
write.table(geo, file="file_for_geo_submission_quantile.txt", sep="\t")
geo<-cbind(RG$genes$FeatureNum,out[[3]][[2]]$E)
write.table(geo, file="file_for_geo_submission_cyclicloess.txt", sep="\t")
geo<-cbind(RG$genes$FeatureNum,out[[3]][[3]]$E)
write.table(geo, file="file_for_geo_submission_normexp_quantile.txt", sep="\t")
geo<-cbind(RG$genes$FeatureNum,out[[3]][[4]]$E)
write.table(geo, file="file_for_geo_submission_normexo_cyclicloess.txt", sep="\t")


save(out, file="out.Rdata")
print("preprocess")
out
}

#---------------------------------------


quality<-function(x, name="array", di){



targets<-x[[1]]
raw1<-x[[2]]
sub2<-x[[3]]
nab<-paste(di,"Densities_boxplots_PCA.pdf", sep="/")
pdf(file=nab, paper="a4r", width=0, height=0)

for(i in 1:5){
	plotDensities(sub2[[i]][[1]], main=names(sub2)[i], legend="topright")
}

for(i in 1:5){
	boxplot(sub2[[i]][[1]], main=names(sub2)[i])
}
for(i in 1:5){
	plotMDS(sub2[[i]][[1]], main=names(sub2)[i])
}

dev.off()

nab<-paste(di,"raw_data_quality", sep="/")
p<-matrix(targets$Cy3,length(colnames(raw1$E)),1, dimnames=list(colnames(raw1$E), "condition"))
p<-as.data.frame(p) 
m<-data.frame(labelDescription="condition", row.names="condition")
phenoData <- new("AnnotatedDataFrame", data=p, varMetadata=m)
a<-ExpressionSet(assayData=raw1$E, phenoData=phenoData)
arrayQualityMetrics(a, outdir = nab, spatial=TRUE, intgroup="condition",do.logtransform=TRUE)

for(i in 1:4){
nab<-paste(di,names(sub2)[i], sep="/")
p<-matrix(targets$Cy3,length(colnames(sub2[[i]]$E)),1, dimnames=list(colnames(sub2[[i]]$E), "condition"))
p<-as.data.frame(p) 
m<-data.frame(labelDescription="condition", row.names="condition")
phenoData <- new("AnnotatedDataFrame", data=p, varMetadata=m)
a<-ExpressionSet(assayData=sub2[[i]]$E, phenoData=phenoData)
arrayQualityMetrics(a, outdir = nab, spatial=TRUE, intgroup="condition",do.logtransform=FALSE)
}





for(i in 1:4){
dat<-sub2[[i]]$E
nab<-paste(di,names(dat)[i], sep="/")
jpeg(file=paste(di,nab,sep="/"), width = 1920, height = 1080)
par(mfrow = c(1, 1))
plot(as.data.frame(dat), main="cyclicloess")
dev.off()
}









}


#--------------------------------------------------


analyze2<-function(x, name="array", di){

targets<-x[[1]]

f <- factor(targets$Cy3)
lev<-levels(f)
design <- model.matrix(~0+f)
colnames(design) <- lev

con<-matrix(,length(lev),length(lev))
rownames(con)<-lev
colnames(con)<-lev
for(i in 1 :length(lev)){
   for(j in 1 :length(lev)){
con[i,j]<-paste(colnames(con)[i],rownames(con)[j], sep="-")}} 

con1<-select.list(sort(con), multiple=TRUE)

con2<-matrix(,length(con1),length(con1))
rownames(con2)<-paste("(",con1,")", sep="")
colnames(con2)<-paste("(",con1,")", sep="")
for(i in 1 :length(con1)){
   for(j in 1 :length(con1)){
con2[i,j]<-paste(colnames(con2)[i],rownames(con2)[j], sep="-")}} 

con3<-select.list(sort(con2), multiple=TRUE)

con4<-c(con1,con3)

cont.wt <- makeContrasts(contrasts=con4,levels=design)

sub_probes_fit<-x[[5]]
sub_features_fit<-x[[6]]
sub_TU_fit<-x[[7]]
Nas<-which(is.na(x[[7]][[1]]$genes$TU.ID))

for(i in 1:5){
sub_TU_fit[[i]]$E<-sub_TU_fit[[i]]$E[-Nas,]
sub_TU_fit[[i]]$Eb<-sub_TU_fit[[i]]$Eb[-Nas,]
sub_TU_fit[[i]]$genes<-sub_TU_fit[[i]]$genes[-Nas,]
}


re<-select.list(choices=c("Yes","No"), multiple=F, title="Average identities required?")

if(re=="Yes"){
	variants<-c("quantile","cyclicloess","normexp_qunatile","normexp_cyclicloess","raw")
	ave_probes_fit<-x[[5]]
	ave_features_fit<-x[[6]]
	for(i in 1:5){
		ave_probes_fit[[i]]<- lmFit(ave_probes_fit[[i]], design)
		ave_features_fit[[i]]<- lmFit(ave_features_fit[[i]], design)
	}
	what<-select.list(c("Probes","Features"), multiple=F, title="Individual average identities required?")
	if(what=="Probes"){
		for(i in 1:5){
			naf<-paste("average_intensity", variants[i],"probes.txt", sep="_")
			naf<-paste(di,naf, sep="/")
			write.table(cbind(ave_probes_fit[[i]]$coefficients,ave_probes_fit[[i]]$genes), file=naf, sep="\t", row.names=F, quote=F)
		}
	}
	if(what=="Features"){
		for(i in 1:5){
			naf<-paste("average_intensity", variants[i],"features.txt", sep="_")
			naf<-paste(di,naf, sep="/")
			write.table(cbind(ave_features_fit[[i]]$coefficients,ave_features_fit[[i]]$genes), file=naf, sep="\t", row.names=F, quote=F)
		}
	}
}



for(i in 1:5){
sub_probes_fit[[i]]<- eBayes(contrasts.fit(lmFit(sub_probes_fit[[i]], design), cont.wt))
sub_features_fit[[i]]<- eBayes(contrasts.fit(lmFit(sub_features_fit[[i]], design), cont.wt))
sub_TU_fit[[i]]<- eBayes(contrasts.fit(lmFit(sub_TU_fit[[i]], design), cont.wt))
}

pdf(file="Sigma vs A plot.pdf", paper="a4r", width=0, height=0)
par(mfrow=c(3,2))
for( i in 1:5){
	plotSA(sub_probes_fit[[i]], main=names(sub_probes_fit)[i])
}
par(mfrow=c(1,1))


par(mfrow=c(3,2))
for( i in 1:5){
	plotSA(sub_features_fit[[i]], main=names(sub_features_fit)[i])
}
par(mfrow=c(1,1))

for( i in 1:5){
	plotSA(sub_TU_fit[[i]], main=names(sub_TU_fit)[i])
}
par(mfrow=c(1,1))
dev.off()

results_probes<-list()
results_features<-list()
results_TU<-list()

p.value<-select.list(as.character(seq(0.01,0.05, by=0.01)), title="significant_p.Value")
p.value<-as.numeric(p.value)

lfc<-select.list(as.character(seq(0.5,3, by=0.1)), title="significant_log2_foldchange")
lfc<-as.numeric(lfc)  

met<-select.list(c("separate","global"), title="decideTests_method (see manual)")

for(i in 1:5){
results_probes[[i]]<-topTable(sub_probes_fit[[i]], sort="none",n=Inf)
results_features[[i]]<-topTable(sub_features_fit[[i]], sort="none",n=Inf)
results_TU[[i]]<-topTable(sub_TU_fit[[i]], sort="none",n=Inf)
}

decide_probes<-list()
decide_features<-list()
decide_TU<-list()

for(i in 1:5){
decide_probes[[i]]<-decideTests(sub_probes_fit[[i]],method=met,adjust.method="BH",p.value=p.value,lfc=lfc)
decide_features[[i]]<-decideTests(sub_features_fit[[i]],method=met,adjust.method="BH",p.value=p.value,lfc=lfc)
decide_TU[[i]]<-decideTests(sub_TU_fit[[i]],method=met,adjust.method="BH",p.value=p.value,lfc=lfc)
}


summary_probes<-c()
summary_features<-c()
summary_TU<-c()

for(i in 1:5){

summary_probes<-rbind(summary_probes, colSums(abs(decide_probes[[i]])))
summary_features<-rbind(summary_features, colSums(abs(decide_features[[i]])))
summary_TU<-rbind(summary_TU, colSums(abs(decide_TU[[i]])))

write.table(summary_probes, file="summary_probes.txt", sep="\t")
write.table(summary_features, file="summary_probes.txt", sep="\t")
write.table(summary_TU, file="summary_TU.txt", sep="\t")
}


rownames(summary_probes)<-c("quantile","cyclicloess","normexp_qunatile","normexp_cyclicloess","raw")
rownames(summary_features)<-c("quantile","cyclicloess","normexp_quantile","normexp_cyclicloess", "raw")
rownames(summary_TU)<-c("quantile","cyclicloess","normexp_quantile","normexp_cyclicloess", "raw")

summary_probes<-cbind(summary_probes, rowSums(summary_probes))
summary_features<-cbind(summary_features, rowSums(summary_features))
summary_TU<-cbind(summary_TU, rowSums(summary_TU))

naf<-paste(di,"preprocessing_impact_summary.jpg", sep="/" )
jpeg(file=naf)
barplot(t(summary_features[,1:(ncol(summary_features)-1)]),legend.text=TRUE,args.legend = list(x = "topleft", bty="n", cex=0.5),cex.names=0.7,main=paste(paste("adj. P.value =",p.value),paste("Foldchange =",lfc),paste("Method =",met),sep=", "))
dev.off()

for(i in 1:5){
results_probes[[i]]<-cbind(results_probes[[i]],decide_probes[[i]])
an<-rep(NA,nrow(results_probes[[i]]))
an[1:3]<-c(paste("adj. P.value =",p.value),paste("Foldchange =",lfc),paste("Method =",met))
results_probes[[i]]<-cbind(results_probes[[i]],an)

results_features[[i]]<-cbind(results_features[[i]],decide_features[[i]])
an<-rep(NA,nrow(results_features[[i]]))
an[1:3]<-c(paste("adj. P.value =",p.value),paste("Foldchange =",lfc),paste("Method =",met))
results_features[[i]]<-cbind(results_features[[i]],an)

results_TU[[i]]<-cbind(results_TU[[i]],decide_TU[[i]])
an<-rep(NA,nrow(results_TU[[i]]))
an[1:3]<-c(paste("adj. P.value =",p.value),paste("Foldchange =",lfc),paste("Method =",met))
results_TU[[i]]<-cbind(results_TU[[i]],an)
}




variants<-c("quantile","cyclicloess","normexp_qunatile","normexp_cyclicloess","raw")
for(i in 1:5){
naf<-paste("results", variants[i],"features.txt", sep="_")
naf<-paste(di,naf, sep="/")
write.table(results_features[[i]], file=naf, sep="\t")

naf<-paste("results", variants[i],"probes.txt", sep="_")
naf<-paste(di,naf, sep="/")
write.table(results_probes[[i]], file=naf, sep="\t")

naf<-paste("results", variants[i],"TU.txt", sep="_")
naf<-paste(di,naf, sep="/")
write.table(results_TU[[i]], file=naf, sep="\t")
}

}
#################
ave_probes_fit[[i]]<- lmFit(ave_probes_fit[[i]], design)
		ave_features_fit[[i]]<- lmFit(ave_features_fit[[i]], design)
#-------------------------------------------
genomeplotdata<-function(x,y, makeline=TRUE){

#out<-list(targets=targets, raw=RG, all_normalized=sub2, all_probes_excluded=sub3, probes=sub_probes, features=sub_features)
sub2<-x[[3]]
#sel<-select.list(names(sub2))
#sel<-grep(paste("i",sel,"i"), paste("i",names(sub2),"i"))

temp_out<-list()

for(ii in 1:length(sub2)){
	targets<-x[[1]]
	conditions<-unique(targets$Cy3)
	
	f <- factor(targets$Cy3)
	lev<-levels(f)
	design <- model.matrix(~0+f)
	colnames(design) <- lev
	
	dat<-lmFit(sub2[[ii]], design)
	
	data1<-cbind(dat$coefficients,dat$genes)
  # data1<-sub2[[ii]]
  # data1<-avereps(data1, ID=data1$genes$ProbeName)
  # data1<-cbind(data1$E,data1$genes)
  
  replicons<-unique(data1$Replicon)
  replicons<-(na.omit(replicons))

  replicon<-list()
  for(i in 1:length(replicons)){
    temp<-which(data1$Replicon==replicons[i])
	temp_data<-data1[temp,]
	replicon[[i]]<-temp_data
	names(replicon[i])<-replicons[i] 
   }

 
 
 # load("NGS_data.Rdata")
  
  output<-list()
  
  for(jj in 1:length(replicons)){
  
  
  data1<-replicon[[jj]]
  
  # clist<-numeric(0)
  # data2<-matrix(, nrow(data1), length(conditions))
  # colnames(data2)<-conditions
  
  # for(i in 1:length(conditions)){
  
    # clist<-grep(paste(conditions[i],".-.", sep=""), colnames(data1))
  # dat<-matrix(,nrow(data1),length(clist))
   
  # for( j in 1:length(clist)){
   
   # dat[,j]<-data1[,clist[j]]
   
   
   # }
   # data2[,i]<-rowMeans(dat)
  
  
  # }
  
  data2<-cbind(data1[,1:length(conditions)], data1)
  
  b<-sort.list(data2$start, na.last=T)
  
  data2<-as.matrix(data2)
  data3<-matrix(,length(b),ncol(data2))
  for(i in 1:length(b)){
  data3[i,]<-data2[b[i],]}
  colnames(data3)<-colnames(data2)
  
  b<-sort.list(data3[,ncol(data3)], na.last=NA)
 # data4<-matrix(,length(b),ncol(data2))
  #for(i in 1:length(b)){
  #data4[i,]<-data3[b[i],]}
  #colnames(data4)<-colnames(data3)
  data5<-as.data.frame(data3)
  orientation<-data5$orientation
  or<-paste(orientation,1,sep="")
  or<-as.numeric(or)
  dd<-data5[,1:length(conditions)]
  ddd<-as.matrix(dd)
  ddd<-as.numeric(ddd)
  ddd<-matrix(ddd,nrow(dd),ncol(dd))
 colnames(ddd)<-colnames(data5[,1:length(conditions)])
 klj<-grep("orientation", colnames(data5))
 
 tab<-cbind(ddd, data5[,(klj-2):klj])
 tab[,ncol(tab)]<-paste(tab[,ncol(tab)],1,sep="")
 # naf<-which(is.na(data5$SystematicName))
 # naf2<-seq(1,length(naf))
 # naf3<-data5$SystematicName
  #naf3<-as.character(naf3)
  #naf3[naf]<-naf2

  features<-unique(data5$SystematicName)
 
  
line<-rep( list(list()), length(features) )
lin<-rep( list(list()), length(features) )
g<-as.matrix(tab)
g<-as.numeric(g)
g<-matrix(g,nrow(tab),ncol(tab))
colnames(g)<-colnames(tab)
posx<-(g[,(ncol(g))-2] + g[,(ncol(g))] /2)

len<-as.numeric(as.character(tab[,(ncol(tab)-1)]))-as.numeric(as.character(tab[,(ncol(tab)-2)]))
len[which(len<=10)]<-50 # workaround für anno fehler
tab[,ncol(tab)]<-len

ma<-max(ddd)
mi<-min(ddd)
offset<-round(mi-3.55, digits=0)
if (offset >= (mi-3.55)) {
offset=offset-0.5
}




if (makeline==TRUE){
for(i in 1:length(features)){


b<-which(paste("i",features[i],"i") == paste("i",data5$SystematicName,"i") )
#print(b)
line[[i]]<-rep(list(list()), length(conditions)+1)

  line[[i]][[1]]<-posx[b]
  
    for(j in 1: length(conditions)){
	  line[[i]][[j+1]]<-(ddd[b,j]-(offset))*(or[b[1]])*(1)
	  
	  }}

}





jff<-paste(replicons[jj],".gff", sep="")
jffTU<-paste(replicons[jj],"_TU",".gff", sep="")

rep_length<-read.delim(jff, header=FALSE)
rep_length<-rep_length[4,]
rep_length<-strsplit(as.character(rep_length)," ")

lll<-length(rep_length[[1]])
rep_length<-as.numeric(rep_length[[1]][(lll-1):lll])


#gff2<-read.csv(jff, sep="\t")
gff2<-gfftojff(jff)

gffTU<-gfftojffTU(jffTU) 


  
graph_file<- paste(replicons[jj],"_rev.grp", sep="")
t_file<-paste(replicons[jj],"_fwd.grp", sep="")

di<-dir()

hh<-grep(graph_file, di)
hhh<-grep(t_file, di)


if(length(hh) > 0 & length(hhh) > 0){
	graph<-read.table(graph_file)
	t<-read.table(t_file)
	
	l<-seq(1,nrow(t))
	t<-cbind(l,t)
	graph<-cbind(l,graph)
	graph<-as.matrix(graph)
	v<-as.matrix(t)
}

temp<-grep("sequencing_conditions.txt", di)
if(length(temp) > 0){
seqset<-read.table("sequencing_conditions.txt", sep="\t", header=FALSE)
seqset<-as.matrix(seqset)
}

if(length(hh) == 0 | length(hhh) == 0){
graph<-matrix(,rep_length[2],(ncol(seqset)-1))
t<-matrix(,rep_length[2],(ncol(seqset)-1))
t[,]<-0
graph[,]<-0

l<-seq(1,nrow(t))
t<-cbind(l,t)
graph<-cbind(l,graph)
graph<-as.matrix(graph)
v<-as.matrix(t)
  
}

seqset<-c()
for(jjj in 1:(ncol(v)-1)){
seqset<-c(seqset,paste("cond_", jjj, sep=""))
}
seqset<-matrix(seqset,1,length(seqset))
temp<-grep("sequencing_conditions.txt", di)
if(length(temp) > 0){
seqset<-read.table("sequencing_conditions.txt", sep="\t", header=FALSE)
seqset<-as.matrix(seqset)
}
colnames(v)[2:ncol(v)]<-as.character(seqset[1,])
colnames(graph)[2:ncol(graph)]<-as.character(seqset[1,])

  temp_output<-list(conditions,posx,tab,orientation,line,gff2,graph,v,lin,gffTU,offset)
  output[[jj]]<-temp_output
  }
  print(replicons)
  names(output)<-replicons
  print(names(sub2)[ii])
  temp_out[[ii]]<-output
 } 
  
  output<-temp_out
  names(output)<-names(sub2)
  save(output, file="genomeplotd.Rdata")
  output
  }
  
  
  
  
  
  
  
  
  
drawgenomeplot<-function(x,name="array", s=0, e=10000,offset=offset, whole=TRUE, makeline=TRUE, linie=linie, di){
targets<-x[[1]]
conditions<-unique(targets$Cy3)

print(targets$Cy3)
ppp<-select.list(c("yes","no"), multiple=F, title="resort the conditions?")
if(ppp=="yes"){
new<-"la"
con<-conditions
count<-1
for(i in 1:length(con)){
ö<-grep(paste("i",new[length(new)],"i"),paste("i",con,"i"))
if(length(ö)>0){
con<-con[-ö]}
ne<-select.list(con, multiple=F, title=paste("choose condition",count, sep=" "))
count<-count+1
new<-c(new,ne)
}


new<-new[-1]
conditions<-new
}

q<-select.list(conditions, multiple=T, title="conditions_included_in genomeplot")
conditionsused<-q





color<-colorselect(conditionsused)


direc<-list.files()
d<-grep("genomeplotd.Rdata",direc)
if(length(d)==0){
out_data<-genomeplotdata(x,conditions, makeline=makeline)

}
else{
load("genomeplotd.Rdata")
out_data<-output
}



sel<-select.list(names(out_data))
sel<-which(names(out_data)==sel)
norm_name<-names(out_data)[sel]

out_data<-out_data[[sel]]

seqset1<-colnames(out_data[[1]][[8]])
seqset1<-seqset1[2:length(seqset1)]
seqset<-select.list(seqset1, multiple=T, title="conditions_included_in genomeplot")
see<-c()
for(i in 1:length(seqset)){
temp<-grep(paste("i",seqset[i], "i"), paste("i", seqset1, "i"))
see<-c(see,temp)


}
seqset<-see+1
seqcol<-colorselect(seqset1[seqset-1], ind=2, trans=TRUE)
if(whole==TRUE){
for(jj in 1:length(out_data)){

output<-out_data[[jj]]

ma<-max(output[[3]][,1:(length(output[[3]])-3)])
#mi<-min(output[[3]][,1:(length(output[[3]])-3)])
#offset<-round(mi-3.55, digits=0)
#if (offset >= (mi-3.55)) {
#offset=offset-0.5
#}

offset<-output[[11]]
yli=ma-offset+0.5
ylim=c(-yli,yli)

print(offset)
print(ylim)

"pdf_ausgabe"
nab<-paste(names(out_data)[jj],"_",norm_name,"_", name, ".pdf", sep="")
st<-1
en2<-nrow(output[[7]])
en<-nrow(output[[7]])%/%10000
en<-en
pdf(file=paste(di,nab,sep="/"), paper="a4r", width=0, height=0)

for (i in 1:en){
if(en > 0){
genomeplot((i-1)*10000,(i*10000),output,conditionsused,color,offset=offset,ylim=ylim,  linie=linie, repliconsize=en2, seqcol=seqcol, seqset=seqset)}}
genomeplot(max(0,(en2-10000)),en2,output,conditionsused,color,offset=offset,ylim=ylim,  linie=linie, repliconsize=en2,seqcol=seqcol, seqset=seqset)
dev.off()


}
}
else{


replicons1<-select.list(names(out_data))
replicons<-which(names(out_data)==replicons1)
output<-out_data[[replicons]]

ma<-max(output[[3]][,1:(length(output[[3]])-3)])
#mi<-min(output[[3]][,1:(length(output[[3]])-3)])
#offset<-round(mi-3.55, digits=0)
#if (offset >= (mi-3.55)) {
#offset=offset-0.5
#}

offset<-output[[11]]
yli=ma-offset+0.5
ylim=c(-yli,yli)

print(offset)
print(ylim)

"pdf_ausgabe"
nab<-paste(names(out_data)[replicons1],"_",norm_name,"_", name,"_",s,"_",e, ".pdf", sep="")
st<-1
en2<-nrow(output[[7]])
en<-nrow(output[[7]])%/%10000
en<-en
pdf(file=paste(di,nab,sep="/"), paper="a4r", width=0, height=0)


genomeplot(s,e,output,conditionsused,color,offset=offset,ylim=ylim,  linie=linie,repliconsize=en2, seqcol=seqcol, seqset=seqset)
dev.off()
}
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

#new<-"la"
#con<-conditionsused
#count<-1
#for(i in 1:length(con)){
#ö<-grep(new[length(new)],con)
#if(length(ö)>0){
#con<-con[-ö]}
#ne<-select.list(con, multiple=F, title=paste("choose condition",count, sep=" "))
#count<-count+1
#new<-c(new,ne)
#}
#conditionsused<-new[2:length(new)]





colorselect<-function(x, ind=1, trans=FALSE){
conditions<-x




k<-paste(length(conditions),"own_colors_in_file", sep="_")
qq<-select.list(c(k, "select_from_pulldown_menu"), multiple=F, title="colors")
if(qq==k) {
colo<-paste("color",ind,".txt", sep="")
colo<-read.table(colo, sep="\t",header=TRUE)
q<-colnames(colo)

}
else{

new<-numeric(0)
for(i in 1:length(conditions)){
t<-paste(conditions[i],"color", sep=" ") 
q<-select.list(colors(), multiple=F, title=t)
new<-c(new,q)
}
q<-new
}
#if(length(q)!=length(conditions)){
#print("!!wrong number of colors!!")
#q<-colorselect(conditions)}
print(q)
if(trans==TRUE){
alph<-select.list(as.character(seq(0,1,by=0.1)))
q<-alpha(q,alph)
par(mfrow=c(1,2))
}
print(q)
barplot(1:length(q)/1:length(q),col=q, axes=F, legend.text=conditions, ylim=c(0,5))
if(trans==TRUE){
barplot(1:length(q),col=q, axes=F, legend.text=conditions, ylim=c(0,5), beside=FALSE, space=-0.8)
}
like<-select.list(c("I like the colors","I don't like the colors"), multiple=F)
if(like!="I like the colors"){
par(mfrow=c(1,1))
q<-colorselect(conditions, ind=ind ,trans=trans)}
par(mfrow=c(1,1))
q

}
#----------------------------






  
  #points
  #---------------------------------------------------------
  genomeplot<-function(x,y, z,zz, xy,offset=6.7, wid=0.2, ylim=c(-13,13),rest=FALSE,  seq=TRUE, names=TRUE, linie=TRUE, versatz=FALSE,punkte=FALSE, log=TRUE, nobox=FALSE, lt=1,tres=11.21, array=TRUE, gene=FALSE,  Achse=TRUE, anno=TRUE, cs=0.2, pch=18, repliconsize, seqcol, seqset){
  
 conditions<-z[[1]]
 posx<-z[[2]]
 tab<-z[[3]]
 orientation<-z[[4]]
 line<-z[[5]]
 gff2<-z[[6]]
 graph<-z[[7]][,c(1,seqset)]
 v<-z[[8]][,c(1,seqset)]
  gffTU<-z[[10]]
 color<-xy
 conditionsused<-zz
 
 #print=c(offset, ylim)
 
 b<-numeric(0)
 for(i in 1:length(conditionsused)){
 b<-c(b,which(conditionsused[i] == conditions))
 }
 
 ta<-matrix(,nrow(tab),length(conditionsused)+3)
 l<-ncol(ta)
 ll<-ncol(tab)
 r<-nrow(tab)
 cl<-colnames(tab)
 tab<-as.matrix(tab)
 tab<-as.numeric(tab)
 tab<-matrix(tab,r,ll)
 colnames(tab)<-cl
 ta[,(l-2):(l)]<-tab[,(ll-2):(ll)]
 for(i in 1:length(conditionsused)){
 ta[,i]<-tab[,b[i]]
 }
  
  colnames(ta)[(l-2):(l)]<-c("start","end","len")
  colnames(ta)[1:(l-3)]<-conditionsused
  tab<-ta
if(log==TRUE){
#ylim=c(-13.35,13.35)
}

 "Plotfenster zeichnen"
plot(0,0 , xlim=c(x,y),ylim=ylim, xlab="Position (nt)", ylab="log2 expression",pch=20, type="n",frame.plot=FALSE, axes=FALSE, new=TRUE)



ll<-numeric(0)
for (i in 1:nrow(gff2)){
if	  (as.numeric(gff2[i,4])>x-20000 & as.numeric(gff2[i,4])<y+20000)
ll<-c(ll, i)
}
if(length(ll)>0){
if (names == TRUE){
for (j in 1:length(ll)){
s=(gff2[ll[j],5]-gff2[ll[j],4])


text(((s/2)+gff2[ll[j],4]),gff2[ll[j],10]*2.4*2, gff2[ll[j],2], cex=0.2, col="white", font=4)
text(((s/2)+gff2[ll[j],4]),gff2[ll[j],8]*3.6*2, gff2[ll[j],9], cex=0.2, col="white", font=4)
}
}






}
"454 daten"
if(seq==TRUE){

#grap<-as.matrix(graph)
#ki<-(v[x:y,3]+v[x:y,2])

#k<-max(v[x:y,2],ki)
#ai<-(grap[x:y,3]+grap[x:y,2])

#a<-min(grap[x:y,2],ai)
#f<-abs(a)/max(ylim)
#q<-k/max(ylim)
#gr<-graph
#vv<-v
#grap[,2]<-grap[,2]/f
#grap[,3]<-grap[,3]/f
#v[,2]<-v[,2]/q
#v[,3]<-v[,3]/q
#v[,4]<-2.3
#grap[,4]<-(-2.3)
vv<-v
gr<-graph
if (log==TRUE){



#points(vv[x:y,1],(log2(vv[x:y,2]+vv[x:y,3]))+3.3, type="h",  col="gray86")
#points(vv[x:y,1],(log2(vv[x:y,2]))+3.3, type="h",  col="gray60")

xx<-x
yy<-y
if(xx-500 > 0){
xx<-xx-550
}
if(yy+500 < nrow(vv)){
yy<-yy+500
}



maxfor<-max(sqrt((vv[x:y,2:ncol(vv)])))
maxrv<-max(sqrt((gr[x:y,2:ncol(gr)])*(-1)))
minrv<--maxrv

normfor<-maxfor/(ylim[2]-3.55)
normrv<-maxrv/(ylim[2]-3.55)
 
normfor<-max(1,normfor)
normrv<-max(1,normrv)

for(i in 1:(ncol(vv)-1)){



#seqfor<-(log2(vv[(xx):(yy),i+1]))+3.3
#seqrv<-(-1*(log2(-1*(gr[xx:yy,i+1]))))-3.3









seqfor<-(((vv[(xx):(yy),i+1])^(1/2))/normfor)  +3.55
seqrv<-((-1*((-1*(gr[xx:yy,i+1])))^(1/2))/normrv) -3.55


seqfor[c(1,length(seqfor))]<-3.55
seqrv[c(1,length(seqrv))]<--3.55

 infsfor<-which(seqfor==Inf | seqfor==-Inf)
 infsrv<-which(seqrv==Inf | seqrv==-Inf)
 
 seqfor[infsfor]<-3.55
 seqrv[infsrv]<--3.55
 
 
 

 
 
polygon(vv[(xx):(yy),1],seqfor, col=seqcol[i], border=NA)
polygon(vv[(xx):(yy),1],seqrv, col=seqcol[i], border=NA)

 
#polygon(vv[(xx):(yy),1],seqfor, col=alpha(seqcol[i],0.4), border=NA)
#polygon(vv[(xx):(yy),1],seqrv, col=alpha(seqcol[i],0.4), border=NA)

#polygon(vv[(xx):(yy),1],seqfor,border=alpha(seqcol[i],0.3), col=NA)
#polygon(vv[(xx):(yy),1],seqrv,border=alpha(seqcol[i],0.3), col=NA)

#points(vv[xx:yy,1],seqfor, type="l",  col=seqcol[i])
#points(vv[xx:yy,1],seqrv, type="l",  col=seqcol[i])

#points(gr[x:y,1],(-1*(log2(-1*(gr[x:y,2]+gr[x:y,3]))))-3.3, type="l",  col="gray86")
#points(gr[x:y,1],(-1*(log2(-1*(gr[x:y,2]))))-3.3, type="l",  col="gray60")


#symbols((xx+(yy-xx)/2),-1.15, rectangles=matrix(c(yy-xx,2.3),1,2), bg="white", add=TRUE, inches=FALSE, lty=0)
#symbols((xx+(yy-xx)/2),1.15, rectangles=matrix(c(yy-xx,2.3),1,2), bg="white", add=TRUE, inches=FALSE, lty=0)



}

legend("topright", legend=colnames(vv)[2:ncol(vv)],bty="n", cex=0.8, fill=alpha(seqcol,0.2))

}

else
{
points(v[x:y,1],(v[x:y,2]+v[x:y,3])+3.3, type="h",  col="gray86")
points(v[x:y,1],v[x:y,2]+3.3, type="h",  col="gray60")
symbols((x+(y-x)/2),1.15, rectangles=matrix(c(y-x,3.3),1,2), bg="white", add=TRUE, inches=FALSE, lty=0)

points(grap[x:y,1],(grap[x:y,2]+grap[x:y,3])-3.3, type="h",  col="gray86")
points(grap[x:y,1],grap[x:y,2]-3.3, type="h",  col="gray60")
symbols((x+(y-x)/2),-1.15, rectangles=matrix(c(y-x,3.3),1,2), bg="white", add=TRUE, inches=FALSE, lty=0)
legend("topright", legend=c("(-) total RNA","(+) with primary 5'ends" ),bty="n", cex=0.8, fill=c("gray86","gray60"), title="log2 read number")
}



}







#restrictionssites einzeichnen für TMDH
if(rest==TRUE){
ll<-numeric(0)
for (i in 1:length(hpaII)){
if	  (hpaII[i]>x-10000 & hpaII[i]<y+10000)
ll<-c(ll, i)
}
for(i in 1:length(ll)){

lines(c(hpaII[ll[i]],hpaII[ll[i]]),c(3.55,7), col="grey60")
lines(c(hpaII[ll[i]],hpaII[ll[i]]),c(-3.55,-7), col="grey60")
}
}


"threshold"
#abline(h=tres-offset, lty=1,col="grey70")
#abline(h=-tres+offset, lty=1,col="grey70")






"Lines"
if(linie==TRUE){

ll<-numeric(0)
for (i in 1:length(line)){
if	  (line[[i]][[1]][1]>x-20000 & line[[i]][[1]][1]<y+20000)
ll<-c(ll, i)
}
b


for(i in 1:length(ll)){
  for(jj in 1:length(b)){
  
  
      lines(line[[ll[i]]][[1]],(line[[ll[i]]][[b[jj]+1]])*-1,col=color[jj])
	
	  }}


}
"Arraydaten"

if(array==TRUE){


#ll<-numeric(0)
#for (i in 1:nrow(tab)){
#if	  (tab[i,(ncol(tab)-2)]>x-10000 & tab[i,(ncol(tab)-2)]<y+10000){
#ll<-c(ll, i)}
#}
ll<-which(as.numeric(as.character(tab[,(ncol(tab)-2)]))>x-1000 & as.numeric(as.character(tab[,(ncol(tab)-2)]))<y+1000)

l<-length(ll)
if(length(ll)>0){
dimen<- matrix(,l,2)
dimen[,1]<-tab[ll,length(conditionsused)+3]
dimen[,2]<-wid
#for(j in 1:length(ll)){
if(punkte==TRUE){
for(i in 1:length(conditionsused)){
    assign(colnames(tab)[i],((((tab[ll,i]))-(offset))*as.numeric(paste(orientation[ll],1, sep=""))*(-1)))
	
	points(c(posx[ll]), c(get(colnames(tab)[i])), col=color[i] ,pch=pch,cex=cs)
	
  }

}

else{

for(i in 1:length(conditionsused)){
    assign(colnames(tab)[i],((((tab[ll,i]))-(offset))*as.numeric(paste(orientation[ll],1, sep=""))*(-1)))
	symbols( c(posx[ll]), c(get(colnames(tab)[i])), rectangles=dimen, inches=FALSE, add=TRUE, fg=FALSE, bg=color[i])
  }
#}
}
}
}
"Achsen"
#Axis(, side=2, at=c(2.3,3.3,4.3,5.3,6.3,7.3,8.3,9.3,10.3,11.3,12.3,13.3), labels=c(9,10,11,12,13,14,15,16,17,18,19,20), cex=0.5, col="black") 
#Axis(, side=2, at=c(-2.3,-3.3,-4.3,-5.3,-6.3,-7.3,-8.3,-9.3,-10.3,-11.3,-12.3,-13.3), labels=c(9,10,11,12,13,14,15,16,17,18,19,20), cex=0.5, col="black") 

#att<-seq(2.3,ylim[2],by=1)
#lab<-seq(2.3+offset,2.3+offset+length(att),by=1)
#lab<-lab[1:length(att)]
#Axis(, side=2, at=att, labels=lab, cex=0.5, col="black") 
#att<-seq(ylim[1],-2.3,by=1)
#lab<-seq(-2.3-offset,-2.3-offset-length(att),by=-1)
#lab<-lab[1:length(att)]
#att<-sort(att, decreasing=TRUE)
#Axis(, side=2, at=att, labels=lab, cex=0.5, col="black") 

att<-seq(1,30)
att<-att-offset
s<-which(att >= 3.55)
att<-att[s]
lab<-att+offset
Axis(, side=2, at=att, labels=lab, cex=0.5, col="black") 
Axis(, side=2, at=att*(-1), labels=lab, cex=0.5, col="black") 



if(seq==TRUE){

if(log==TRUE){

ticksfor<-(sqrt(c(1,10,25,50,100,250,500,1000,2500,5000,10000,25000,50000,100000,500000,10^6))/normfor)+3.55
ticksrv<-(sqrt(c(1,10,25,50,100,250,500,1000,2500,5000,10000,25000,50000,100000,500000,10^6))/normrv)*(-1)-3.55

labfor<-c(1,10,25,50,100,250,500,1000,2500,5000,10000,25000,50000,100000,500000,10^6)
labrv<-c(1,10,25,50,100,250,500,1000,2500,5000,10000,25000,50000,100000,500000,10^6)

#Axis(, side=4, at=c(3.3,4.3,5.3,6.3,7.3,8.3,9.3,10.3,11.3,12.3,13.3,14.3,15.3,16.3,17.3), labels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14), cex=0.5, col="lightgrey") 
#Axis(, side=4, at=c(-3.3,-4.3,-5.3,-6.3,-7.3,-8.3,-9.3,-10.3,-11.3,-12.3,-13.3,-14.3,-15.3,-16.3,-17.3), labels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14), cex=0.5, col="lightgrey")} 
Axis(, side=4, at=ticksfor, labels=labfor, cex=0.3, col="lightgrey") 
Axis(, side=4, at=ticksrv, labels=labrv, cex=0.3, col="lightgrey") 
}

else
{
Axis(, side=4, at=c(2.3,3.3,4.3,5.3,6.3,7.3,8.3,9.3,10.3,11.3,12.3,13.3,14.3), labels=c(0,round(k/12),round(k*2/12),round(k*3/12),round(k*4/12),round(k*5/12),round(k*6/12),round(k*7/12),round(k*8/12),round(k*9/12),round(k*10/12),round(k*11/12),round(k)), cex=0.5, col="lightgrey") 
Axis(, side=4, at=c(-2.3,-3.3,-4.3,-5.3,-6.3,-7.3,-8.3,-9.3,-10.3,-11.3,-12.3,-13.3,-14.3), labels=c(0,round(abs(a)/12),round(abs(a)*2/12),round(abs(a)*3/12),round(abs(a)*4/12),round(abs(a)*5/12),round(abs(a)*6/12),round(abs(a)*7/12),round(abs(a)*8/12),round(abs(a)*9/12),round(abs(a)*10/12),round(abs(a)*11/12),round(abs(a))), cex=0.5, col="lightgrey")} }




abline(h=0.5*3, lty=1, col="gray86")
abline(h=(-0.5)*3, lty=1, col="gray86")


if(nobox==TRUE){
lt<-0}
if(Achse==TRUE){
"x-Achse"
num<- seq(0,repliconsize, by=2500)
n1m<- (seq(0,repliconsize, by=2500)+1250)
num1<-cbind( num, -0.35)
num2<-cbind( num, 0.15)
num3<-cbind( num, -0.15)
n1m1<-cbind( n1m, 0.35)
n1m2<-cbind( n1m, 0.15)
n1m3<-cbind( n1m, -0.15)
text(num1, label=(num), cex=0.75)
arrows(num3[,1],num3[,2], num2[,1], num2[,2], length=0)
arrows(n1m3[,1],n1m3[,2], n1m2[,1], n1m2[,2], length=0)
text(n1m1, label=(n1m), cex=0.75)
abline(h=0, col="black")}

colo<-as.matrix(gff2[,7])
#color<-as.matrix(gff_genes[,7])





"Annotation"

# genes 
ll<-numeric(0)
for (i in 1:nrow(gff2)){
if	  (as.numeric(gff2[i,4])>x-20000 & as.numeric(gff2[i,4])<y+20000)
ll<-c(ll, i)
}
if(length(ll)>0){

for (j in 1:length(ll)){
s=(gff2[ll[j],5]-gff2[ll[j],4])

symbols(((s/2)+as.numeric(gff2[ll[j],4])),as.numeric(gff2[ll[j],8])*3,rectangles=matrix(c(s,1.3),1,2), bg=colo[ll[j],], add=TRUE, inches=FALSE, lty=lt)
}



if (names == TRUE)
for (j in 1:length(ll)){
s=(gff2[ll[j],5]-gff2[ll[j],4])


text(((s/2)+gff2[ll[j],4]),gff2[ll[j],10]*2.4, gff2[ll[j],2], cex=0.7, col="black", font=4)
text(((s/2)+gff2[ll[j],4]),gff2[ll[j],8]*3.6, gff2[ll[j],9], cex=0.7, col="black", font=4)

text(((s/2)+gff2[ll[j],4]),gff2[ll[j],10]*2*(ylim[2]-1), gff2[ll[j],2], cex=0.01, col="white", font=4)
text(((s/2)+gff2[ll[j],4]),gff2[ll[j],8]*2*(ylim[2]-1), gff2[ll[j],9], cex=0.01, col="white", font=4)
}
}
# TU annotation

if(is.na(gffTU)==FALSE){
if(length(ll)>0){
ll<-numeric(0)
for (i in 1:nrow(gffTU)){
if	  (as.numeric(gffTU[i,4])>x-20000 & as.numeric(gffTU[i,4])<y+20000)
ll<-c(ll, i)
}

for (j in 1:length(ll)){
s=(gffTU[ll[j],5]-gffTU[ll[j],4])

symbols(((s/2)+as.numeric(gffTU[ll[j],4])),as.numeric(gffTU[ll[j],8])*5.5,rectangles=matrix(c(s,0.65),1,2), bg="#ffedcc", add=TRUE, inches=FALSE, lty=lt)
}

for (j in 1:length(ll)){
s=(gffTU[ll[j],5]-gffTU[ll[j],4])


text(((s/2)+gffTU[ll[j],4]),gffTU[ll[j],10]*5.5, gffTU[ll[j],2], cex=0.7, col="black", font=4)
text(((s/2)+gffTU[ll[j],4]),gffTU[ll[j],10]*2*(ylim[2]-1), gffTU[ll[j],2], cex=0.01, col="white", font=4)
}
}

}


legend("topleft", 19, conditionsused, text.col = color,bty="n",  cex=0.7)



"x-Achse"
if(Achse==FALSE){
text(0,0.35, label="0", cex=0.75)
text(repliconsize,0.35, label=repliconsize, cex=0.75)
arrows(0,0.15, 0, -0.15, length=0)
arrows(repliconsize,0.15, repliconsize, -0.15, length=0)
num<- (seq(0,repliconsize, by=1000000)+250000)
n1m<- (seq(0,repliconsize, by=1000000)-250000)
num1<-cbind( num, -0.35)
num2<-cbind( num, 0.15)
num3<-cbind( num, -0.15)
n1m1<-cbind( n1m, 0.35)
n1m2<-cbind( n1m, 0.15)
n1m3<-cbind( n1m, -0.15)
text(num1, label=as.character(num), cex=0.75)
arrows(num3[,1],num3[,2], num2[,1], num2[,2], length=0)
arrows(n1m3[,1],n1m3[,2], n1m2[,1], n1m2[,2], length=0)
text(n1m1, label=as.character(n1m), cex=0.75)
abline(h=0, col="black")}




}  


#genomeplot(60000,70000,output[[2]],conditionsused,color,offset=offset,ylim=ylim,  linie=linie, repliconsize=en2, seqcol=seqcol, seqset=seqset)



#genomeplot(0,10000,output,output[[1]],c("black","dodgerblue2"),offset=offset,ylim=ylim,  linie=FALSE, seq=FALSE,tres=6.7,rest=TRUE)
 #genomeplot(0,3573470,output,output[[1]][2],c("dodgerblue2"),offset=offset,ylim=ylim,  linie=FALSE,seq=FALSE, names=FALSE, nobox=TRUE, Achse=FALSE,punkte=TRUE,tres=7.3)

#genomeplot((0,3573470,output,output[[1]][2],c("dodgerblue2"),offset=offset,ylim=ylim,  linie=FALSE)
#analyze(whole=FALSE, s=40000, e=50000, linie=FALSE)


#pdf(file="test.pdf", paper="a4r", width=0, height=0)

#for (i in 1:10){
#genomeplot((i-1)*10000,(i*10000),output,output[[1]],c("black","dodgerblue2"),offset=offset,ylim=ylim,  linie=FALSE, seq=FALSE,tres=6.7,rest=TRUE)}
#genomeplot(3563470,3573470,output,conditionsused,color,offset=offset,ylim=ylim,  linie=TRUE)
#dev.off()

#pdf(file="TMDH_hpaII_5000bp.pdf", paper="a4r", width=0, height=0)

#for (i in 1:714){
#genomeplot((i-1)*5000,(i*5000),output,output[[1]],c("black","dodgerblue2"),offset=offset,ylim=ylim,  linie=FALSE, seq=FALSE,tres=6.7,rest=TRUE)}
#genomeplot(3573470-5000,3573470,output,output[[1]],c("black","dodgerblue2"),offset=offset,ylim=ylim,  linie=FALSE, seq=FALSE,tres=6.7,rest=TRUE)
#dev.off()

gfftojff<-function(x){
genome<-read.csv(x,sep="\t", header=FALSE,comment.char="#")

genes<-grep("gene", genome[,3])
asRNA<-grep("asRNA", genome[,3])
UTR<-grep("5'UTR", genome[,3])
ncRNA<-grep("ncRNA", genome[,3])
rRNA<-grep("rRNA", genome[,3])
tRNA<-grep("tRNA", genome[,3])
intern<-grep("intern", genome[,3])

al<-c(genes,asRNA,UTR,ncRNA,tRNA,rRNA)

genome<-genome[al,]


genes<-grep("gene", genome[,3])
asRNA<-grep("asRNA", genome[,3])
UTR<-grep("5'UTR", genome[,3])
ncRNA<-grep("ncRNA", genome[,3])
rRNA<-grep("rRNA", genome[,3])
tRNA<-grep("tRNA", genome[,3])
intern<-grep("intern", genome[,3])
noname<-c(UTR,intern)
anno<-strsplit(as.character(genome[,9]),";")

anno_gene<-matrix(,length(anno),1)
anno_tag<-matrix(,length(anno),1)
for(i in 1:length(anno)){
temp<-grep("gene=", anno[[i]])
temp2<-grep("locus_tag=", anno[[i]])

if(length(temp>0)){
anno_gene[i,]<-gsub("gene=", "",anno[[i]][temp])
}

if(length(temp2>0)){
anno_tag[i,]<-gsub("locus_tag=", "",anno[[i]][temp2])
}
}

temp<-which(anno_gene=="NA")
anno_gene[temp,]<-""

jff<-matrix(,nrow(genome),10)

jff[,1]<-as.character(genome[,3])

if(length(noname)>0){
jff[-noname,2]<-anno_tag[-noname]
}

if(length(noname)==0){
jff[,2]<-anno_tag
}


jff[,3]<-as.character(genome[,3])
jff[asRNA,2]<-""
jff[,4]<-genome[,4]
jff[,5]<-genome[,5]
jff[,6]<-paste(genome[,7],1, sep="")
jff[,7]<-"#4682B4"
jff[genes,7]<-"#4682B4"
jff[asRNA,7]<-"#EE0000"
jff[UTR,7]<-"#F0F8FF"
jff[ncRNA,7]<-"#FFB90F"
jff[intern,7]<-"#87CEFA"
jff[rRNA,7]<-"#8C8C8C"
jff[tRNA,7]<-"#90EE90"

jff[,8]<-paste(genome[,7],0.5, sep="")
jff[,9]<-anno_gene
jff[,10]<-jff[,8]

jff<-jff[order(as.numeric(jff[,4])),]


write.table(jff, file="temp.txt", row.names = FALSE, sep="\t")
jff<-read.csv("temp.txt", sep="\t")
unlink("temp.txt")

jff


}

gfftojffTU<-function(x){



try(genome<-read.csv(x,sep="\t", header=FALSE,comment.char="#"))

if(exists("genome")==TRUE){

anno<-strsplit(as.character(genome[,9]),";")

anno_type<-matrix(,length(anno),1)
anno_index<-matrix(,length(anno),1)

for(i in 1:length(anno)){
temp<-grep("type=", anno[[i]])
temp2<-grep("index=", anno[[i]])

if(length(temp>0)){
anno_type[i,]<-gsub("TU_type=", "",anno[[i]][temp])
}

if(length(temp2>0)){
anno_index[i,]<-gsub("index=", "",anno[[i]][temp2])
}

}

jff<-matrix(,nrow(genome),10)

jff[,1]<-as.character(genome[,3])

jff[,2]<-anno_index
jff[,3]<-as.character(genome[,3])
jff[,4]<-genome[,4]
jff[,5]<-genome[,5]
jff[,6]<-paste(genome[,7],1, sep="")
jff[,7]<-"#ffe4b2"
jff[,8]<-paste(genome[,7],0.5, sep="")
jff[,9]<-anno_type
jff[,10]<-jff[,8]

jff<-jff[order(as.numeric(jff[,4])),]


write.table(jff, file="temp.txt", row.names = FALSE, sep="\t")
jff<-read.csv("temp.txt", sep="\t")
unlink("temp.txt")
}

else{
jff<-NA
}
jff


}
