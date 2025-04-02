#dd <- "/research/bsi/projects/staff_analysis/m216453/2023_06_22_BMDD/Data/microbiomeHD_cleaned_data/"
dd <- "/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Code/mPower_paper/Data/EffectSize/"

setwd(dd)
files <- list.files()
otu_filter <- function(feature.dat, prev = 0.1, dep = 100){ # otu * sam
  idx <- apply(feature.dat, 1, function(x) sum(x>0)>(ncol(feature.dat)*prev))
  idx2 <- colSums(feature.dat) > dep
  return(feature.dat[idx,idx2])
} 
setwd(dd)
files <- list.files(pattern = 'Rdata$')
res <- NULL
xx <- files[!(files %in% unique(res$data))]
for(file in files){
  tmp <- otu_table <- meta.dat <- tax_table <- NULL
  cat('---------')
  cat(file,':')
  load(file)
  case.tmp <- unique(metadata$DiseaseState)[-grep('H',unique(metadata$DiseaseState))]
  print(table(metadata$DiseaseState,metadata$grp))
  xx <- cbind(table(metadata$DiseaseState,metadata$grp))
  xx2 <- xx[,colnames(xx)[!(colnames(xx) %in% "H")]]
  case.tmp <- paste0(names(xx2[xx2!=0]),collapse=',')
  control.tmp  <- paste0(names(xx2[xx2==0]),collapse=',')
  if(file =="autism_kb.Rdata"){formu <- '~grp'}
  if(file =="crc_baxter.Rdata"){metadata$bmi <- as.numeric(metadata$BMI_s); formu <- '~grp+bmi'}
  if(file =="crc_zackular.Rdata"){
    metadata$age <- iconv(metadata$age, from = "latin1", to = "UTF-8")
    metadata$age <- as.numeric(gsub("\u00A0", "", metadata$age))
    # metadata$age <- as.numeric(gsub('\xa0\xa0','',metadata$age))
    formu <- '~grp+sex+age'
  }
  if(file =="crc_zhao.Rdata"){formu <- '~grp'}
  if(file =="hiv_lozupone.Rdata"){formu <- '~grp+age'}
  if(file =="hiv_noguerajulian.Rdata"){metadata$sex <- metadata$host_sex_s; metadata$age <- metadata$host_age_s; formu <- '~grp+sex+age'}
  if(file =="mhe_zhang.Rdata"){formu <- '~grp'}
  if(file =="nash_chan.Rdata"){formu <- '~grp'}
  if(file =="ob_goodrich.Rdata"){formu <- '~grp+sex+age'}
  if(file =="par_scheperjans.Rdata"){formu <- '~grp'}
  if(file =="ra_littman.Rdata"){formu <- '~grp'}
  if(file =="crc_zeller.Rdata"){metadata$age <- metadata$`Age (years)`;metadata$bmi <- metadata$`BMI (kg/m)`; formu <- '~grp+age+bmi'}
  # if(file =="asd_son.Rdata"){metadata$sex <- metadata$host_sex_s; formu <- '~grp+sex'}
  # if(file =="cdi_schubert.Rdata"){metadata$sex <- metadata$gender; formu <- '~grp+sex+age'}
  # if(file =="cdi_vincent_v3v5.Rdata"){formu <- '~grp'}
  # if(file =="cdi_youngster.Rdata"){formu <- '~grp'}
  # if(file =="crc_xiang.Rdata"){formu <- '~grp'}
  # if(file =="crc_zeller.Rdata"){metadata$age <- metadata$`Age (years)`;metadata$bmi <- metadata$`BMI (kg/m)`; formu <- '~grp+age+bmi'}
  # if(file =="edd_singh.Rdata"){formu <- '~grp'}
  # if(file =="hiv_dinh.Rdata"){metadata$sex <- metadata$sex_s;metadata$bmi <- metadata$BMI_s;metadata$age <- metadata$age_s; formu <- '~grp+sex+bmi+age'}
  # if(file =="ibd_alm.Rdata"){metadata$sex <- metadata$gender; formu <- '~grp+sex+age'}
  # if(file =="ibd_engstrand_maxee.Rdata"){formu <- '~grp'}
  # if(file =="ibd_gevers_2014.Rdata"){formu <- '~grp'}
  # if(file =="ibd_huttenhower.Rdata"){metadata$sex <- metadata$gender; formu <- '~grp+age+gender'}
  # if(file =="nash_ob_baker.Rdata"){formu <- '~grp+sex+age'}
  # if(file =="ob_gordon_2008_v2.Rdata"){formu <- '~grp'}
  # if(file =="ob_ross.Rdata"){metadata$age <- metadata$age_at_visit; formu <- '~grp+sex+age+BMI'}
  # if(file =="ob_zupancic.Rdata"){metadata$sex <- metadata$sex_s; formu <- '~grp+sex'}
  # if(file =="t1d_alkanani.Rdata"){formu <- '~grp'}
  # if(file =="t1d_mejialeon.Rdata"){metadata$age <- metadata$age_in_years; formu <- '~grp+age+BMI'}

  otu.tab <- inner_join(otu_table %>% rownames_to_column('otu'), tax_table[,c('otu','genus')]) %>% column_to_rownames('otu')
  otu.tab <- aggregate(.~genus, data = otu.tab, function(x) sum(x)) %>% column_to_rownames('genus')
  # otu.tab <- otu_table #[If use raw OTU table]
  meta.dat <- metadata
  meta.dat$grp <- as.factor(meta.dat$grp)
  feature.dat <- otu_filter(otu.tab)
  meta.dat <- meta.dat[colnames(feature.dat),]
  idx <- names(which(colSums(feature.dat) > 0))
  if(length(idx)>0){
    meta.dat <- meta.dat[idx,, drop =F]
    print(table(meta.dat$grp))
    cat(';')
    print(min(colSums(feature.dat)))
    cat('\n')
    colnames(meta.dat)
    table(meta.dat$DiseaseState)
    if('H' %in% names(table(meta.dat$grp))){
      n.H.tmp <- table(meta.dat$grp)['H']
      n.D.tmp <- table(meta.dat$grp)[!(names(table(meta.dat$grp)) %in% 'H')]
    }else{
      n.H.tmp <- table(meta.dat$grp)['nonIBD']
      n.D.tmp <- table(meta.dat$grp)[!(names(table(meta.dat$grp)) %in% 'nonIBD')]
    }
    try({linda.obj  <- MicrobiomeStat::linda(feature.dat = feature.dat, meta.dat = meta.dat,
                                             formula = formu, feature.dat.type = 'count', prev.filter = 0, max.abund.filter = 0)})
    try({tmp <- linda.obj$output[[1]] %>% dplyr::select(log2FoldChange, pvalue) %>% 
      mutate(n.H =n.H.tmp, 
             n.D =n.D.tmp, 
             Cases = case.tmp,
             Controls = control.tmp,
             m = nrow(feature.dat),
             formula = gsub('\\~grp\\+','',formu),
             data = file, 
             pi0 = pi0est(linda.obj$output[[1]][,'pvalue'],lambda=seq(0.1, max(linda.obj$output[[1]][,'pvalue']) - 0.05, 0.05))$pi0)})
    res <- rbind(res, tmp)
  }else{
    cat('[====',file,' does not have sufficient samples ====]\n')
  }
}

# res$data <- gsub("asd_son.Rdata","Son 2015, ASD",res$data)
# res$data <- gsub("cdi_schubert.Rdata","Schubert 2014, CDI",res$data)
# res$data <- gsub("cdi_vincent_v3v5.Rdata","Vincent 2013, CDI",res$data)
# res$data <- gsub("cdi_youngster.Rdata","Youngster 2014, CDI",res$data)
# res$data <- gsub("edd_singh.Rdata","Singh 2015, EDD",res$data)
# res$data <- gsub("hiv_dinh.Rdata","Dinh 2015, HIV",res$data)
# res$data <- gsub("ibd_gevers_2014.Rdata","Gevers 2014, IBD",res$data)
# res$data <- gsub("ob_ross.Rdata","Ross 2015, OB",res$data)
# res$data <- gsub("ob_zupancic.Rdata","Zupancic 2012, OB",res$data)
# res$data <- gsub("t1d_alkanani.Rdata","Alkanani 2015, T1D",res$data)
# res$data <- gsub("t1d_mejialeon.Rdata" ,"Mejia-Leon 2014, T1D",res$data)
# res$data <- gsub("ibd_alm.Rdata","Papa 2012, IBD",res$data)
# res$data <- gsub("ibd_huttenhower.Rdata","Morgan 2012, IBD",res$data)
# res$data <- gsub("crc_xiang.Rdata","Chen 2012, CRC",res$data)
# res$data <- gsub("ibd_engstrand_maxee.Rdata","Willing 2010, IBD",res$data)
# res$data <- gsub("nash_ob_baker.Rdata","Zhu 2013, NASH",res$data)
# res$data <- gsub("ob_gordon_2008_v2.Rdata","Turnbaugh 2009, OB",res$data)

res$data <- gsub("autism_kb.Rdata","Kang 2013, ASD",res$data)
res$data <- gsub("crc_baxter.Rdata","Baxter 2016, CRC",res$data)
res$data <- gsub("crc_zackular.Rdata","Zackular 2014, CRC",res$data)
res$data <- gsub("crc_zhao.Rdata","Wang 2012, CRC",res$data)
res$data <- gsub("hiv_lozupone.Rdata","Lozupone 2013, HIV",res$data)
res$data <- gsub("hiv_noguerajulian.Rdata","Noguera-Julian 2016, HIV",res$data)
res$data <- gsub("mhe_zhang.Rdata","Zhang 2013, LIV",res$data)
res$data <- gsub("nash_chan.Rdata","Wong 2013, NASH",res$data)
res$data <- gsub("ob_goodrich.Rdata","Goodrich 2014, OB",res$data)
res$data <- gsub("par_scheperjans.Rdata" ,"Scheperjans 2015, PAR",res$data)
res$data <- gsub("ra_littman.Rdata","Scher 2013, ART",res$data)
res$data <- gsub("crc_zeller.Rdata","Zeller 2014, CRC",res$data)

# save(res, file = "/research/bsi/projects/staff_analysis/m216453/2023_05_06_PowerShiny/res.Rdata")




res2 <- res %>% 
  group_by(data) %>% 
  mutate(max = max(abs(log2FoldChange))) %>% 
  mutate(dataset = data, 
         data = paste0(data,'[differential(%)=',round(1-pi0,3),';max(LFC)=',round(max,3),']')) 



dt <- as.data.frame(unique(res2[,-c(1:2)])) %>% dplyr::select(-c(data)) %>% 
  mutate(adjusted = gsub('\\~grp','',formula), `Differential taxa(proportion)`=round(1-pi0, 3)) %>% 
  dplyr::rename(`max(log2FoldChange)`=max, `Control(n)` = n.H, `Case(n)` = n.D, `#tested feature`= m) %>% dplyr::select(-c('formula','pi0')) %>% 
  dplyr::select(dataset, Controls,`Control(n)`, Cases, `Case(n)`, `#tested feature`, adjusted, `max(log2FoldChange)`,`Differential taxa(proportion)`) %>% 
  mutate(dataset, `max(log2FoldChange)` = round(`max(log2FoldChange)`, 3)) 
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
colnames(dt) <- firstup(colnames(dt))
# save(dt, file = "/research/bsi/projects/staff_analysis/m216453/2023_05_06_PowerShiny/dt.Rdata")
dt2 <- dt %>% 
  dplyr::rename(`Percentage of differential taxa`=`Differential taxa(proportion)`, `Max log2 fold change` =`Max(log2FoldChange)`)

head(res)
p0 <- ggplot(res[res$data %in% dt2$Dataset,], aes(x=log2FoldChange))+
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity", 
                 binwidth = 0.1, color = brewer.pal(8,'Set1')[4], fill = 'white')+
  geom_density(fill = brewer.pal(8,'Set1')[3], alpha=0.7) + theme_classic() + 
  facet_wrap(~data, scale='free') + 
theme(legend.position = 'none',
      axis.text = element_text(color='black', size=14),
      axis.title = element_text(color='black', size=14),
      strip.text = element_text(color='black', size=10)) +
  labs(y = 'Density')
p0

for(file in files){
  tmp <- otu_table <- meta.dat <- tax_table <- NULL
  cat('---------')
  cat(file,':')
  load(file)
  case.tmp <- unique(metadata$DiseaseState)[-grep('H',unique(metadata$DiseaseState))]
  print(table(metadata$DiseaseState,metadata$grp))
  xx <- cbind(table(metadata$DiseaseState,metadata$grp))
  xx2 <- xx[,colnames(xx)[!(colnames(xx) %in% "H")]]
  case.tmp <- names(xx2[xx2!=0])
}



setwd(dd)
file = 'crc_zackular.Rdata'
load(file)
case.tmp <- paste0(unique(metadata$DiseaseState)[-grep('H',unique(metadata$DiseaseState))], collapse = ';')
print(table(metadata$DiseaseState,metadata$grp))
xx <- cbind(table(metadata$DiseaseState,metadata$grp))
if(file =="crc_zackular.Rdata"){
  metadata$age <- iconv(metadata$age, from = "latin1", to = "UTF-8")
  metadata$age <- as.numeric(gsub("\u00A0", "", metadata$age))
  # metadata$age <- as.numeric(gsub('\xa0\xa0','',metadata$age))
  formu <- '~grp+sex+age'
}
otu.tab <- inner_join(otu_table %>% rownames_to_column('otu'), tax_table[,c('otu','genus')]) %>% column_to_rownames('otu')
otu.tab <- aggregate(.~genus, data = otu.tab, function(x) sum(x)) %>% column_to_rownames('genus')
meta.dat <- metadata
meta.dat$grp <- as.factor(meta.dat$grp)
feature.dat <- otu_filter(otu.tab)
meta.dat <- meta.dat[colnames(feature.dat),]
idx <- names(which(colSums(feature.dat) > 0))
meta.dat <- meta.dat[idx,, drop =F]
print(table(meta.dat$grp))
cat(';')
print(min(colSums(feature.dat)))
cat('\n')
colnames(meta.dat)
table(meta.dat$DiseaseState)
if('H' %in% names(table(meta.dat$grp))){
  n.H.tmp <- table(meta.dat$grp)['H']
  n.D.tmp <- table(meta.dat$grp)[!(names(table(meta.dat$grp)) %in% 'H')]
}else{
  n.H.tmp <- table(meta.dat$grp)['nonIBD']
  n.D.tmp <- table(meta.dat$grp)[!(names(table(meta.dat$grp)) %in% 'nonIBD')]
}
try({linda.obj  <- MicrobiomeStat::linda(feature.dat = feature.dat, meta.dat = meta.dat,
                                         formula = formu, feature.dat.type = 'count', prev.filter = 0, max.abund.filter = 0)})
try({tmp <- linda.obj$output[[1]] %>% dplyr::select(log2FoldChange, pvalue) %>% 
  mutate(n.H =n.H.tmp, 
         n.D =n.D.tmp, 
         Cases = case.tmp,
         Controls = control.tmp,
         m = nrow(feature.dat),
         formula = gsub('\\~grp\\+','',formu),
         data = file, 
         pi0 = pi0est(linda.obj$output[[1]][,'pvalue'],lambda=seq(0.1, max(linda.obj$output[[1]][,'pvalue']) - 0.05, 0.05))$pi0
         )})
density(tmp$log2FoldChange)



idx <- rownames(tmp)[which(tmp$pvalue<0.05)]
Y <- feature.dat + 0.5
logY <- log2(Y)
W <- t(t(logY) - colMeans(logY))
prop <- t(t(feature.dat)/colSums(feature.dat))
p1 <- merge(W[idx,] %>% t(), meta.dat[,'grp', drop =F], by = 0) %>% column_to_rownames('Row.names')%>% melt() %>% 
  ggplot(aes(x=grp,y = value, color=grp)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size=0.8, alpha=0.5) + 
  facet_wrap(~variable, scale='free', nrow= 2) + 
  scale_color_brewer(palette = 'Dark2') + 
  labs(x='', y ='clr(abundance)') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text(color='black', size=14),
        axis.title = element_text(color='black', size=14),
        strip.text = element_text(face='italic',color='black', size=10)) 
p1
p2 <- merge(prop[idx,] %>% t(), meta.dat[,'grp', drop =F], by = 0) %>% column_to_rownames('Row.names')%>% melt() %>% 
  ggplot(aes(x=grp,y = sqrt(value), color=grp)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size=0.8, alpha=0.5) + 
  facet_wrap(~variable, scale='free', nrow= 2) + 
  scale_color_brewer(palette = 'Dark2') + 
  labs(x='', y ='sqrt(proportion)') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text(color='black', size=14),
        axis.title = element_text(color='black', size=14),
        strip.text = element_text(face='italic',color='black', size=10)) 
p12 <- ggarrange(p1, p2, labels = c('b','c'), nrow = 1)

ggarrange(p0, p12, labels = c('a',''), nrow = 2)
ggarrange(p0, p1,p2, labels = c('a','b','c'), nrow = 3,heights = c(1.5, 1, 1))
ggsave(file=paste0(rd,'Figure3.pdf'), width = 10, height = 12)
