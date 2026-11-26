# Descriptive NMDS complement using immutable primary Bray-Curtis matrices.
# Run: Rscript --vanilla run_nmds.R INPUT_DIRECTORY NEW_OUTPUT_DIRECTORY
suppressPackageStartupMessages(library(vegan))
stopifnot(as.character(packageVersion('vegan')) == '2.6.8')
args <- commandArgs(TRUE);stopifnot(length(args)==2)
inp<-normalizePath(args[1]);out<-args[2]
if(dir.exists(out)) stop('Output directory exists; choose a new execution directory.')
dir.create(out,recursive=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_SESSION_INFO.txt'))
write_tsv<-function(x,name)write.table(x,file.path(out,name),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
RNGkind(kind='Mersenne-Twister',normal.kind='Inversion',sample.kind='Rejection')
allpoints<-list();allshepard<-list();alldiag<-list()
for(m in c('16S','ITS')) {
  seed<-if(m=='16S')2026091501L else 2026091502L
  a<-as.matrix(read.delim(file.path(inp,paste0(m,'_bray.tsv')),row.names=1,check.names=FALSE))
  stopifnot(nrow(a)==9,ncol(a)==9,identical(rownames(a),colnames(a)),all(is.finite(a)),all(a>=0),all(a<=1),max(abs(a-t(a)))<1e-12,all(diag(a)==0))
  stopifnot(setequal(rownames(a),c('B1','B2','B3','M1','M2','M3','A1','A2','A3')))
  set.seed(seed)
  fit<-metaMDS(as.dist(a),k=2,engine='monoMDS',try=100,trymax=200,
     autotransform=FALSE,noshare=FALSE,wascores=FALSE,expand=FALSE,
     parallel=1,trace=1,plot=FALSE,halfchange=FALSE,maxit=1000,
     model='global',weakties=TRUE,stress=1,smin=1e-4,sfgrmin=1e-7,sratmax=0.999999)
  saveRDS(fit,file.path(out,paste0(m,'_NMDS_NATIVE.rds')))
  z<-scores(fit,display='sites',choices=1:2)
  allpoints[[m]]<-data.frame(marker=m,sample=rownames(z),NMDS1=z[,1],NMDS2=z[,2],point_stress_contribution=goodness(fit))
  allshepard[[m]]<-data.frame(marker=m,sample_1=rownames(a)[fit$iidx],sample_2=rownames(a)[fit$jidx],observed_dissimilarity=fit$diss,ordination_distance=fit$dist,monotonic_fitted_distance=fit$dhat)
  stopifnot(nrow(allshepard[[m]])==36)
  alldiag[[m]]<-data.frame(marker=m,n=9,k=2,seed=seed,stress_1=fit$stress,
     random_starts=fit$tries,repeated_best=fit$converged,best_try=fit$bestry,
     engine_iterations=fit$iters,engine_stop_code=fit$icause,
     nonmetric_fit_R2=1-fit$stress^2,linear_fit_R2=cor(fit$dist,fit$dhat)^2,
     vegan_version=as.character(packageVersion('vegan')),
     input=if(m=='16S')'Bray-Curtis; rarefied retained prokaryotic counts' else 'Bray-Curtis; non-rarefied fungal proportions')
  pdf(file.path(out,paste0(m,'_NATIVE_SHEPARD.pdf')),width=6,height=5)
  stressplot(fit,main=paste(m,'NMDS Shepard plot'));dev.off()
  writeLines(capture.output(str(fit,max.level=1)),file.path(out,paste0(m,'FIT_STRUCTURE.txt')))
}
write_tsv(do.call(rbind,allpoints),'NMDS_COORDINATES.tsv')
write_tsv(do.call(rbind,allshepard),'NMDS_SHEPARD.tsv')
write_tsv(do.call(rbind,alldiag),'NMDS_DIAGNOSTICS.tsv')
cat('\nNMDS execution completed.\n');print(do.call(rbind,alldiag))
