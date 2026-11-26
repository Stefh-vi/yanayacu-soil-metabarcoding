# Refit current beta endpoints from the exact retained distance matrices.
# Rscript --vanilla replay_beta.R PACKAGE_ROOT NEW_OUTPUT_DIRECTORY
# Native fits use vegan 2.6.8. The two alpha P values per marker are imported
# from the current seven-row category table, independently refitted by analysis/.
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(permute))
suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(packageVersion('vegan')) == '2.6.8')
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
root <- normalizePath(args[1],mustWork=TRUE)
out <- args[2]
if(dir.exists(out)) stop('Choose a new output directory; existing results are not overwritten.')
inp <- file.path(root,'sequence_workflow','native_analysis','ecology','inputs')
meta <- read.delim(file.path(root,'data','metadata','soil_metadata_source_with_preparation.tsv'),check.names=FALSE)
ids <- c('B1','B2','B3','M1','M2','M3','A1','A2','A3')
meta <- meta[match(ids,meta$SampleID),,drop=FALSE]
stopifnot(identical(as.character(meta$SampleID),ids))
meta$zone <- factor(meta$zone,levels=c('Low degradation','Moderate degradation','High degradation'))
meta$locality <- factor(meta$locality,levels=c('J','SI'))
P <- as.matrix(read.delim(file.path(inp,'permutations_within_locality_4319.tsv'),check.names=FALSE))
storage.mode(P) <- 'integer'
stopifnot(identical(dim(P),c(4319L,9L)),identical(colnames(P),ids),
          all(apply(P,1,function(x) identical(sort(unname(x)),1:9))),
          !any(apply(P,1,function(x) identical(unname(x),1:9))),
          all(apply(P,1,function(x) all(meta$locality[x]==meta$locality))))
current <- read.delim(file.path(inp,'current_category_tests.tsv'),check.names=FALSE)
disp_expected <- read.delim(file.path(inp,'current_dispersion_tests.tsv'),check.names=FALSE)
stopifnot(nrow(current)==7L,nrow(disp_expected)==3L)
dir.create(out,recursive=TRUE)
write_tsv <- function(x,name) write.table(x,file.path(out,name),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
specs <- data.frame(marker=c('16S','ITS','ITS'),endpoint=c('Bray-Curtis','Bray-Curtis proportions','Hellinger-Euclidean'),
                    file=c('16S_bray_lingoes.tsv','ITS_bray.tsv','ITS_hellinger.tsv'))
betas <- list(); dispersions <- list()
for(i in seq_len(nrow(specs))) {
  s <- specs[i,]
  x <- as.matrix(read.delim(file.path(inp,s$file),row.names=1,check.names=FALSE))
  x <- x[ids,ids,drop=FALSE]; storage.mode(x) <- 'double'
  stopifnot(all(is.finite(x)),all(x>=0),max(abs(x-t(x)))<1e-12,all(diag(x)==0))
  d <- as.dist(x)
  fit <- dbrda(d ~ zone + Condition(locality),data=meta)
  an <- anova(fit,permutations=P,model='reduced',parallel=1)
  stopifnot(an$Df[1]==2L,an$Df[2]==5L)
  b <- data.frame(Marker=s$marker,Endpoint=s$endpoint,F=an$F[1],P_raw=an$`Pr(>F)`[1],
          R2_total=fit$CCA$tot.chi/fit$tot.chi,
          partial_R2=fit$CCA$tot.chi/(fit$CCA$tot.chi+fit$CA$tot.chi))
  e <- current[current$Marker==s$marker & current$Endpoint==s$endpoint,]
  stopifnot(nrow(e)==1L,abs(b$F-e$`F statistic`)<1e-9,abs(b$P_raw-e$`P raw`)<1e-12,
            abs(b$partial_R2-e$`Partial R2 after locality`)<1e-9)
  bd <- betadisper(d,meta$zone,type='median',bias.adjust=TRUE,sqrt.dist=FALSE,add=FALSE)
  dp <- permutest(bd,permutations=P,pairwise=FALSE,parallel=1)
  ds <- data.frame(Marker=s$marker,Endpoint=s$endpoint,F=dp$tab$F[1],P_raw=dp$tab$`Pr(>F)`[1])
  de <- disp_expected[i,]
  stopifnot(abs(ds$F-de$`F statistic`)<1e-8,abs(ds$P_raw-de$`P raw`)<1e-12)
  stem <- paste(s$marker,tools::file_path_sans_ext(s$file),sep='_')
  saveRDS(list(beta=fit,anova=an,dispersion=bd,dispersion_test=dp),file.path(out,paste0(stem,'_MODELS.rds')))
  write_tsv(data.frame(Sample=ids,Distance=bd$distances),paste0(stem,'_PERMDISP_DISTANCES.tsv'))
  betas[[i]] <- b; dispersions[[i]] <- ds
}
b <- do.call(rbind,betas); ds <- do.call(rbind,dispersions)
for(marker in c('16S','ITS')) {
  ix <- which(current$Marker==marker)
  for(j in which(b$Marker==marker)) {
    current$`P raw`[current$Marker==marker & current$Endpoint==b$Endpoint[j]] <- b$P_raw[j]
  }
  q <- p.adjust(current$`P raw`[ix],method='holm',n=length(ix))
  stopifnot(max(abs(q-current$`P Holm`[ix]))<1e-12)
  current$`P Holm`[ix] <- q
}
ds$P_Holm <- NA_real_
ds$P_Holm[ds$Marker=='ITS'] <- p.adjust(ds$P_raw[ds$Marker=='ITS'],method='holm',n=2L)
write_tsv(b,'REFITTED_BETA.tsv')
write_tsv(ds,'REFITTED_PERMDISP.tsv')
write_tsv(current,'CURRENT_CATEGORY_FAMILIES.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'SESSION_INFO.txt'))
write_json(list(status='PASS_NATIVE_BETA_AND_DISPERSION_REPLAY',beta_fits=3L,dispersion_fits=3L,
  alpha_refitted_here=FALSE,alpha_source='current_category_tests.tsv; see independent numerical alpha replay',
  category_Holm_families=list('16S'=3L,'ITS'=4L),dispersion_family='ITS Holm n2; 16S single diagnostic',
  limitation='Permutation restriction within locality does not make PERMDISP a locality-adjusted dispersion model.'),
  file.path(out,'REPLAY_RECEIPT.json'),auto_unbox=TRUE,pretty=TRUE)
