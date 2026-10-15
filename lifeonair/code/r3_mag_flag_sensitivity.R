#!/usr/bin/env Rscript
# Focused R3 sensitivity: baseline and exclusion-union fits, without repeating
# the main 2,000-resample climate analysis. Run from repository root in eq-lifeonair.
# Implements the published genus, genomic-trait and candidate CBB models on
# identical frozen inputs, removing flagged MAGs before constructing predictors.
suppressPackageStartupMessages(library(lme4))
IN <- 'lifeonair/results/inputs'; OUT <- 'lifeonair/results/review_r3'
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
rd <- function(f,...) read.delim(f,quote='',check.names=FALSE,...)
wt <- function(x,n) write.table(x,file.path(OUT,n),sep='\t',quote=FALSE,row.names=FALSE)
meta <- rd(file.path(IN,'samples.tsv')); meta <- meta[meta$trip!='' & meta$compartment!='' & meta$site %in% 1:60,]
cls0 <- rd(file.path(IN,'atmo_classification_v2.tsv'))
gr0 <- rd(file.path(IN,'rep_growth_v2.tsv')); gr0 <- gr0[is.finite(gr0$d_hours)&gr0$d_hours>0,]
gs0 <- rd(file.path(IN,'rep_genome_size.tsv'))
ko0 <- rd(file.path(IN,'rep_cbb_ko.tsv'))
ab0 <- rd(gzfile(file.path(IN,'coverm_rep_abundance_long.tsv.gz')),header=FALSE,col.names=c('sample','rep','abund'))
ab0 <- ab0[ab0$rep!='unmapped' & is.finite(ab0$abund)&ab0$abund>0,]
tax <- rbind(rd(file.path(IN,'r232_bac120.summary.tsv'))[,1:2],rd(file.path(IN,'r232_ar53.summary.tsv'))[,1:2])
rk <- function(cl,p) ifelse(grepl(paste0(p,'__[^;]'),cl),sub(paste0('.*',p,'__([^;]*).*'),'\\1',cl),'unassigned')
cl <- tax$classification; tax$phylum <- rk(cl,'p'); tax$order <- paste(tax$phylum,rk(cl,'o'),sep=';');tax$genus <- paste(tax$order,rk(cl,'g'),sep=';')
orph <- readLines('lifeonair/results/reads/orphan_sub_reps.txt')
flags <- readLines(file.path(OUT,'mag_exclusion_union.txt'))
models <- counts <- NULL
record <- function(m,variant,family,model,terms) {
 s <- summary(m)$coefficients
 for (term in intersect(terms,rownames(s))) models <<- rbind(models,data.frame(variant,family,model,term,n=nobs(m),estimate=s[term,1],se=s[term,2],singular=if(inherits(m,'merMod')) isSingular(m) else NA))
}
for (variant in c('baseline','excluded')) {
 drop <- if(variant=='baseline') orph else union(orph,flags)
 cls <- cls0[!cls0$rep_id %in% drop,]; gr <- gr0[!gr0$rep %in% drop,]; gs <- gs0[!gs0$genome_id %in% drop,]
 carrier <- function(g) unique(cls$rep_id[cls$marker=='hydB' & grepl(paste0('Group_',g,'$'),cls$assignment)])
 sets <- list(c1l=carrier('1l'),c1h=carrier('1h'),c2a=carrier('2a'),cfI=setdiff(readLines(file.path(IN,'guild_coxL_reps_v2.txt')),drop))
 for (k in names(sets)) {gr[[k]] <- gr$rep %in% sets[[k]];gs[[k]] <- gs$genome_id %in% sets[[k]]}
 gr$logd <- log(gr$d_hours)
 for(k in c('phylum','order','genus')) gr[[k]] <- tax[[k]][match(gr$rep,tax$user_genome)]
 gs$order <- tax$order[match(gs$genome_id,tax$user_genome)];gs$logd <- gr$logd[match(gs$genome_id,gr$rep)]
 gs$log_est_mb <- log(gs$length_bp/(gs$completeness/100)/1e6);gs$gc100 <- 100*gs$gc
 for (resp in c('log_est_mb','gc100','logd')) {
  m <- lmer(as.formula(paste(resp,'~ c1l+c1h+c2a+cfI+(1|order)')),gs)
  record(m,variant,'trait',resp,paste0(names(sets),'TRUE'))
 }
 ab <- ab0[!ab0$rep %in% drop,];ab$logd <- gr$logd[match(ab$rep,gr$rep)];ab <- ab[!is.na(ab$logd),]
 ab$genus <- gr$genus[match(ab$rep,gr$rep)]
 gm <- xtabs(abund~sample+genus,ab);gm <- sweep(as.matrix(unclass(gm)),1,rowSums(gm),'/')*100
 gm <- as.data.frame(gm);gm$sample_id <- rownames(gm);gm <- merge(meta,gm,by='sample_id')
 gen <- setdiff(names(gm),names(meta));dm <- colMeans(gm[gm$compartment=='Deep',gen]);gen <- gen[dm>=.05]
 k <- paste(gm$site,gm$trip);keys <- intersect(k[gm$compartment=='PlantRoot'],k[gm$compartment=='Deep'])
 L <- t(sapply(keys,function(z) {
  a <- apply(as.matrix(gm[k==z & gm$compartment=='PlantRoot',gen,drop=FALSE]),2,median)
  b <- apply(as.matrix(gm[k==z & gm$compartment=='Deep',gen,drop=FALSE]),2,median)
  log2((a+.01)/(b+.01))
 }))
 gt <- data.frame(genus=gen,deep_mean_pct=dm[gen],log2fc_root=apply(L,2,median),median_logd=tapply(gr$logd,gr$genus,median)[gen])
 for (pair in list(c('c1l','frac_1l'),c('c1h','frac_1h'),c('c2a','frac_2a'),c('cfI','frac_formI'))) gt[[pair[2]]] <- tapply(gr[[pair[1]]],gr$genus,mean)[gen]
 gt$phylum <- sub(';.*','',gt$genus);gt$order <- sub('^([^;]*;[^;]*);.*','\\1',gt$genus)
 for (rank in c('phylum','order')) {
  m <- lm(as.formula(paste('log2fc_root~median_logd+frac_1l+frac_1h+frac_2a+frac_formI+',rank)),gt,weights=sqrt(deep_mean_pct))
  record(m,variant,'genus',rank,c('median_logd','frac_1l','frac_1h','frac_2a','frac_formI'))
 }
 kt <- ko0[!ko0$rep %in% drop,];kt$cbb <- as.integer(kt$K01601==1 & kt$K01602==1 & kt$K00855==1)
 kt$is1l <- kt$rep %in% sets$c1l;kt$is1h <- kt$rep %in% sets$c1h
 kt <- kt[xor(kt$is1l,kt$is1h),];kt$is1l <- as.integer(kt$is1l)
 # h123_v2 uses named order alone and excludes genomes lacking an order.
 clkt <- tax$classification[match(kt$rep,tax$user_genome)]
 kt$order <- ifelse(grepl('o__[^;]',clkt),sub('.*o__([^;]*).*','\\1',clkt),NA)
 kt <- kt[!is.na(kt$order),]
 kt$comp <- gs$completeness[match(kt$rep,gs$genome_id)]
 for (model in c('within_order','completeness_adjusted','hc90')) {
  z <- if(model=='hc90') kt[kt$comp>=90,] else kt
  f <- if(model=='completeness_adjusted') cbb~is1l+scale(comp)+(1|order) else cbb~is1l+(1|order)
  m <- glmer(f,z,family=binomial);record(m,variant,'CBB',model,'is1l')
 }
 counts <- rbind(counts,data.frame(variant,excluded=length(drop),genomes=nrow(gs),growth_genomes=nrow(gr),genera=nrow(gt),carriers_1l=length(sets$c1l),carriers_1h=length(sets$c1h)))
}
wt(models,'mag_sensitivity_models.tsv');wt(counts,'mag_sensitivity_counts.tsv')
# Baseline reproduction guards against accidental differences in the sensitivity implementation.
original <- rd('lifeonair/results/h1_mechanism_v2/m2d_genus_models.tsv')
expected <- original$Estimate[original$model=='growth + carriage + phylum' & original$term=='frac_1l']
actual <- models$estimate[models$variant=='baseline' & models$family=='genus' & models$model=='phylum' & models$term=='frac_1l']
stopifnot(abs(actual-expected)<1e-8)
orig <- rd('lifeonair/results/a5_genome_traits/trait_models_within_order.tsv')
z <- models[models$variant=='baseline' & models$family=='trait',]
ix <- match(paste(z$model,z$term),paste(orig$response,orig$term));stopifnot(max(abs(z$estimate-orig$estimate[ix]))<1e-7)
orig_cbb <- rd('lifeonair/results/h123_v2/A4_cbb_cocarriage_by_group.tsv')[1,]
cols <- c(within_order='within_order_log_odds_1l_vs_1h',completeness_adjusted='completeness_adjusted_log_odds',hc90='hc90_log_odds')
for (model in names(cols)) {
 actual <- models$estimate[models$variant=='baseline' & models$family=='CBB' & models$model==model]
 stopifnot(length(actual)==1,abs(actual-orig_cbb[[cols[[model]]]])<1e-7)
}
writeLines(capture.output(sessionInfo()),file.path(OUT,'R-session.txt'))
print(counts);print(models[models$term %in% c('frac_1l','is1l'),])
