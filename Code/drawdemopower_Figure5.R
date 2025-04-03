source('/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Code/mPower/R/utils.R')
source('/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Code/mPower/R/mPower.R')
tm <- load('/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Data/mforge_output/Daniel2018_AGP_UK(0.05prev).RData')

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(modeest)
library(tibble)
library(MASS)

## for shinyR demo
demo1 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
                test = 'Community', design = 'CaseControl',
                nSams = c(50, 100, 150, 200), grp.ratio = 0.5,
                iters = 50, alpha = 0.05, distance = 'BC',
                diff.otu.pct = 0.082, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
                covariate.eff.min = 0, covariate.eff.maxs = 0.978,
                confounder = 'no', depth.mu = 100000, depth.sd = 40000)
demo1
ggsave("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/demo1.png",
       width = 12, height = 4)

demo2 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
                test = 'Community', design = 'CaseControl',
                nSams = c(50, 100, 150, 200), grp.ratio = 0.5,
                iters = 50, alpha = 0.05, distance = 'BC',
                diff.otu.pct = 0.082, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
                covariate.eff.min = 0, covariate.eff.maxs = 0.978, conf.cov.cor = 0.6,
                confounder = 'yes', depth.mu = 100000, depth.sd = 40000)

demo2
ggsave("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/demo2.png",
       width = 12, height = 4)

## For case studies
## case1: community-level power
res1.noconf <- mPower(feature.dat = feature.dat, model.paras = model.paras,
                      test = 'Community', design = 'CaseControl',
                      nSams = c(20,40,60,80,100), grp.ratio = 0.5,
                      iters = 1000, alpha = 0.05, distance = 'BC',
                      diff.otu.pct = 0.082, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
                      covariate.eff.min = 0, covariate.eff.maxs = 0.978,
                      confounder = 'no', depth.mu = 100000, depth.sd = 40000)

res1.conf <- mPower(feature.dat = feature.dat, model.paras = model.paras,
               test = 'Community', design = 'CaseControl',
               nSams = c(20,40, 60, 80,100), grp.ratio = 0.5,
               iters = 1000, alpha = 0.05, distance = 'BC',
               diff.otu.pct = 0.082, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
               covariate.eff.min = 0, covariate.eff.maxs = 0.978, conf.cov.cor = 0.6,
               confounder = 'yes', depth.mu = 100000, depth.sd = 40000)
ggarrange(res1.noconf$plot, res1.conf$plot, nrow=2)
ggsave("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case1.pdf", width = 14, height = 12)
save(res1.noconf, res1.conf, file = '/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case1.RData')
P1 <- generate_plot3(data = res1.noconf$TPR.beta, effect.size = 'covariate.eff.max', ylab= 'Power', matched.pair = F)
P2 <- generate_plot3(data = res1.conf$TPR.beta, effect.size = 'covariate.eff.max', ylab= 'Power', matched.pair = F)
ggarrange(P1,P2, nrow =1)

summary((res1.noconf$power$power-res1.conf$power$power)/res1.noconf$power$power)

## case2: no confounder vs. strong confounder
res2.noconf <- mPower(feature.dat = feature.dat, model.paras = model.paras,
               test = 'Taxa', design = 'CaseControl',
               nSams = 50, grp.ratio = 0.5,# 20 vs. 30
               iters = 100, alpha = 0.1,
               diff.otu.pct = 0.27, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
               covariate.eff.min = 0, covariate.eff.maxs = c(1,2,3,4,5),
               confounder = 'no', depth.mu = 100000, depth.sd = 40000)
res2.conf <- mPower(feature.dat = feature.dat, model.paras = model.paras,
               test = 'Taxa', design = 'CaseControl',
               nSams = 50, grp.ratio = 0.5,
               iters = 100, alpha = 0.1,
               diff.otu.pct = 0.27, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
               covariate.eff.min = 0, covariate.eff.maxs = c(1,2,3,4,5),
               confounder.eff.max = 1, conf.cov.cor = 0.6, conf.diff.otu.pct = 0.1, # confounder effect
               confounder = 'yes', depth.mu = 100000, depth.sd = 40000)
ggarrange(res2.noconf$plot, res2.conf$plot, nrow=2)
ggsave("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case2.pdf", width = 12, height = 8)
save(res2.noconf, res2.conf, file = '/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case2.RData')

summary((res2.noconf$aTPR$aTPR-res2.conf$aTPR$aTPR)/res2.noconf$aTPR$aTPR)

res2.noconf.genus <- mPower(feature.dat = feature.dat.genus, model.paras = model.paras.genus,
                      test = 'Taxa', design = 'CaseControl',
                      nSams = 50, grp.ratio = 0.5,# 20 vs. 30
                      iters = 100, alpha = 0.1,
                      prev.filter = 0.2, max.abund.filter = 0.002,
                      diff.otu.pct = 0.27, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
                      covariate.eff.min = 0, covariate.eff.maxs = c(0.5,1,1.5,2,2.5),
                      confounder = 'no', depth.mu = 100000, depth.sd = 40000)
res2.conf.genus <- mPower(feature.dat = feature.dat, model.paras = model.paras,
                    test = 'Taxa', design = 'CaseControl',
                    nSams = 50, grp.ratio = 0.5,
                    iters = 100, alpha = 0.1,
                    prev.filter = 0.2, max.abund.filter = 0.002,
                    diff.otu.pct = 0.27, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
                    covariate.eff.min = 0, covariate.eff.maxs = c(0.5,1,1.5,2,2.5),
                    confounder.eff.max = 1, conf.cov.cor = 0.6, conf.diff.otu.pct = 0.1, # confounder effect
                    confounder = 'yes', depth.mu = 100000, depth.sd = 40000)
ggarrange(res2.noconf.genus$plot, res2.conf.genus$plot, nrow=2)
ggsave("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case2.genus.pdf", width = 12, height = 8)
save(res2.noconf.genus, res2.conf.genus, file = '/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case2.genus.RData')


## case3: paired design
res3.match <- mPower(feature.dat = feature.dat, model.paras = model.paras,
                     test = 'Taxa', design = 'MatchedPair',
                     nSams = 25, iters = 100, alpha = 0.1,
                     diff.otu.pct = 0.27, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
                     covariate.eff.min = 0, covariate.eff.maxs = c(1,2,3,4,5),
                     confounder = 'no', depth.mu = 100000, depth.sd = 40000)
res3.casecontrol <- mPower(feature.dat = feature.dat, model.paras = model.paras,
                           test = 'Taxa', design = 'CaseControl',
                           nSams = 50, grp.ratio = 0.5,
                           iters = 100, alpha = 0.1,
                           diff.otu.pct = 0.27, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
                           covariate.eff.min = 0, covariate.eff.maxs = c(1,2,3,4,5),
                           confounder.eff.max = 1, conf.cov.cor = 0.6, conf.diff.otu.pct = 0.1, # confounder effect
                           confounder = 'yes', depth.mu = 100000, depth.sd = 40000)
gc()
res3.match$aTPR;res3.casecontrol$aTPR
mean((res3.match$aTPR$aTPR-res3.casecontrol$aTPR$aTPR)/res3.match$aTPR$aTPR)
ggarrange(res3.match$plot, res3.casecontrol$plot, nrow=2)
ggsave("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case3.pdf", width = 12, height = 8)
save(res3.match,res3.casecontrol, file = '/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case3.RData')

res3.match.genus <- mPower(feature.dat = feature.dat.genus, model.paras = model.paras.genus,
                     test = 'Taxa', design = 'MatchedPair',
                     nSams = 25, iters = 100, alpha = 0.1,
                     prev.filter = 0.2, max.abund.filter = 0.002,
                     diff.otu.pct = 0.27, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
                     covariate.eff.min = 0, covariate.eff.maxs = c(0.5,1,1.5,2,2.5),
                     confounder = 'no', depth.mu = 100000, depth.sd = 40000)
res3.casecontrol.genus <- mPower(feature.dat = feature.dat.genus, model.paras = model.paras.genus,
                           test = 'Taxa', design = 'CaseControl',
                           nSams = 50, grp.ratio = 0.5,
                           iters = 100, alpha = 0.1,
                           prev.filter = 0.2, max.abund.filter = 0.002,
                           diff.otu.pct = 0.27, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
                           covariate.eff.min = 0, covariate.eff.maxs = c(0.5,1,1.5,2,2.5),
                           confounder.eff.max = 1, conf.cov.cor = 0.6, conf.diff.otu.pct = 0.1, # confounder effect
                           confounder = 'yes', depth.mu = 100000, depth.sd = 40000)
res3.casecontrol.noconf.genus <- mPower(feature.dat = feature.dat.genus, model.paras = model.paras.genus,
                                 test = 'Taxa', design = 'CaseControl',
                                 nSams = 50, grp.ratio = 0.5,
                                 iters = 100, alpha = 0.1,
                                 prev.filter = 0.2, max.abund.filter = 0.002,
                                 diff.otu.pct = 0.27, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
                                 covariate.eff.min = 0, covariate.eff.maxs = c(0.5,1,1.5,2,2.5),
                                 confounder.eff.max = 1, conf.cov.cor = 0.6, conf.diff.otu.pct = 0.1, # confounder effect
                                 confounder = 'no', depth.mu = 100000, depth.sd = 40000)
gc()
res3.match.genus$aTPR;res3.casecontrol.genus$aTPR
mean((res3.match.genus$aTPR$aTPR-res3.casecontrol.genus$aTPR$aTPR)/res3.match.genus$aTPR$aTPR)
mean((res3.match.genus$aTPR$aTPR-res3.casecontrol.noconf.genus$aTPR$aTPR)/res3.match.genus$aTPR$aTPR)
ggarrange(res3.match.genus$plot, res3.casecontrol.genus$plot,res3.casecontrol.noconf.genus$plot, nrow=3)
ggsave("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case3.genus.pdf", width = 12, height = 8)
save(res3.match.genus,res3.casecontrol.genus, res3.casecontrol.noconf.genus, file = '/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case3.genus.RData')


## case4: depth
res4.low <- mPower(feature.dat = feature.dat, model.paras = model.paras,
                   test = 'Taxa', design = 'CaseControl',
                   nSams = 50, grp.ratio = 0.5,
                   iters = 100, alpha = 0.1,
                   prev.filter = 0, max.abund.filter = 0,
                   diff.otu.pct = 0.2, diff.otu.direct = 'balanced',diff.otu.mode = 'rare',
                   covariate.eff.min = 0, covariate.eff.maxs = c(1,2,3,4,5),
                   confounder.eff.max = 1, conf.cov.cor = 0.8, conf.diff.otu.pct = 0.05,
                confounder = 'no', depth.mu = 10000, depth.sd = 4000)
res4.high <- mPower(feature.dat = feature.dat, model.paras = model.paras,
                    test = 'Taxa', design = 'CaseControl',
                    nSams = 50, grp.ratio = 0.5,
                    iters = 100, alpha = 0.1,
                    prev.filter = 0, max.abund.filter = 0,
                    diff.otu.pct = 0.2, diff.otu.direct = 'balanced',diff.otu.mode = 'rare',
                    covariate.eff.min = 0, covariate.eff.maxs = c(1,2,3,4,5),
                    confounder.eff.max = 1, conf.cov.cor = 0.8, conf.diff.otu.pct = 0.05,
                    confounder = 'no', depth.mu = 100000, depth.sd = 40000)
res4.highest <- mPower(feature.dat = feature.dat, model.paras = model.paras,
                    test = 'Taxa', design = 'CaseControl',
                    nSams = 50, grp.ratio = 0.5,
                    iters = 100, alpha = 0.1,
                    prev.filter = 0, max.abund.filter = 0,
                    diff.otu.pct = 0.2, diff.otu.direct = 'balanced',diff.otu.mode = 'rare',
                    covariate.eff.min = 0, covariate.eff.maxs = c(1,2,3,4,5),
                    confounder.eff.max = 1, conf.cov.cor = 0.8, conf.diff.otu.pct = 0.05,
                    confounder = 'no', depth.mu = 1000000, depth.sd = 400000)

# (depth.mu^2) / (depth.sd^2 - depth.mu)
ggarrange(res4.low$plot, res4.high$plot,res4.highest$plot,nrow=3)
ggsave("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case4.pdf", width = 12, height = 10)
save(res4.low, res4.high, res4.highest, file = '/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case4.RData')


res4.low.genus <- mPower(feature.dat = feature.dat.genus, model.paras = model.paras.genus,
                   test = 'Taxa', design = 'CaseControl',
                   nSams = 50, grp.ratio = 0.5,
                   iters = 100, alpha = 0.1,
                   prev.filter = 0, max.abund.filter = 0,
                   diff.otu.pct = 0.2, diff.otu.direct = 'balanced',diff.otu.mode = 'rare',
                   covariate.eff.min = 0, covariate.eff.maxs = c(0.5,1,1.5,2,2.5),
                   confounder.eff.max = 1, conf.cov.cor = 0.8, conf.diff.otu.pct = 0.05,
                   confounder = 'no', depth.mu = 10000, depth.sd = 4000)
res4.high.genus <- mPower(feature.dat = feature.dat.genus, model.paras = model.paras.genus,
                    test = 'Taxa', design = 'CaseControl',
                    nSams = 50, grp.ratio = 0.5,
                    iters = 100, alpha = 0.1,
                    prev.filter = 0, max.abund.filter = 0,
                    diff.otu.pct = 0.2, diff.otu.direct = 'balanced',diff.otu.mode = 'rare',
                    covariate.eff.min = 0, covariate.eff.maxs = c(0.5,1,1.5,2,2.5),
                    confounder.eff.max = 1, conf.cov.cor = 0.8, conf.diff.otu.pct = 0.05,
                    confounder = 'no', depth.mu = 100000, depth.sd = 40000)
res4.highest.genus <- mPower(feature.dat = feature.dat.genus, model.paras = model.paras.genus,
                       test = 'Taxa', design = 'CaseControl',
                       nSams = 50, grp.ratio = 0.5,
                       iters = 100, alpha = 0.1,
                       prev.filter = 0, max.abund.filter = 0,
                       diff.otu.pct = 0.2, diff.otu.direct = 'balanced',diff.otu.mode = 'rare',
                       covariate.eff.min = 0, covariate.eff.maxs = c(0.5,1,1.5,2,2.5),
                       confounder.eff.max = 1, conf.cov.cor = 0.8, conf.diff.otu.pct = 0.05,
                       confounder = 'no', depth.mu = 1000000, depth.sd = 400000)

# (depth.mu^2) / (depth.sd^2 - depth.mu)
ggarrange(res4.low.genus$plot, res4.high.genus$plot,res4.highest.genus$plot,nrow=3)
ggsave("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case4.genus.pdf", width = 12, height = 10)
save(res4.low.genus, res4.high.genus, res4.highest.genus, file = '/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case4.genus.RData')




library(RColorBrewer)
load('/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case1.RData')
df <- rbind(res1.noconf$TPR.beta %>% mutate(Group='No confounder'),
            res1.conf$TPR.beta %>% mutate(Group='With confounder'))
p1 <- ggplot(df, aes(x = nSam, y = value, color = Group)) +
  geom_line(aes(group = Group), size=1) +  # Line plot
  geom_point(size = 3, aes(alpha = 0.5), show.legend = F) +  # Points
  scale_color_manual(values = brewer.pal(8,'Dark2')[1:2]) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +  # Error bars
  labs(title = "",
       x = "Sample Size",
       y = "Community-level power", color = '') +
  theme_classic() +
  theme(axis.title = element_text(color = 'black', size=16),
        legend.text = element_text(color = 'black', size=16),
        axis.text = element_text(color = 'black', size=16))
p1


load('/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case2.genus.RData')
res2.noconf <- res2.noconf.genus
res2.conf <- res2.conf.genus
df <- rbind(res2.noconf$aTPR %>% mutate(Group='No confounder'),
            res2.conf$aTPR %>% mutate(Group='With confounder'))
p2 <- ggplot(df, aes(x = `max log2 fold change`, y = aTPR, color = Group)) +
  geom_line(aes(group=Group), size=1) +
  geom_point(size = 3, aes(alpha = 0.5), show.legend = F) +
  scale_color_manual(values = brewer.pal(8,'Dark2')[3:4]) +
  labs(title = "", x = "max log2 fold change", y = "aTPR", color = '') +
  theme_classic() +
  theme(axis.title = element_text(color = 'black', size=16),
        legend.text = element_text(color = 'black', size=16),
        axis.text = element_text(color = 'black', size=16)
  )

p2


load('/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case3.genus.RData')
res3.match <- res3.match.genus
res3.casecontrol <- res3.casecontrol.genus
df <- rbind(res3.match$aTPR %>% mutate(Group='Matched Pair'),
            res3.casecontrol$aTPR %>% mutate(Group='Case-Control'))
p3 <- ggplot(df, aes(x = `max log2 fold change`, y = aTPR, color = Group)) +
  geom_line(aes(group=Group), size=1) +
  geom_point(size = 3, aes(alpha = 0.5), show.legend = F) +
  scale_color_manual(values = brewer.pal(8,'Dark2')[5:6]) +
  labs(title = "", x = "max log2 fold change", y = "aTPR", color = '') +
  theme_classic() +
  theme(axis.title = element_text(color = 'black', size=16),
        legend.text = element_text(color = 'black', size=16),
        axis.text = element_text(color = 'black', size=16)
        )
p3

load('/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/case4.genus.RData')
res4.low <- res4.low.genus
res4.high <- res4.high.genus
res4.highest <- res4.highest.genus
df <- rbind(res4.low$aTPR %>% mutate(Group='Mean depth=10k'),
            res4.high$aTPR %>% mutate(Group='Mean depth=100k'),
            res4.highest$aTPR %>% mutate(Group='Mean depth=1000k'))
p4 <- ggplot(df, aes(x = `max log2 fold change`, y = aTPR, color = Group)) +
  geom_line(aes(group=Group), size=1) +
  geom_point(size = 3, aes(alpha = 0.5), show.legend = F) +
  scale_color_manual(values = brewer.pal(8,'Dark2')[6:8]) +
  labs(title = "", x = "max log2 fold change", y = "aTPR", color = '') +
  theme_classic() +
  theme(axis.title = element_text(color = 'black', size=16),
        legend.text = element_text(color = 'black', size=16),
        axis.text = element_text(color = 'black', size=16)
  )
p4

library(ggpubr)
ggarrange(p1, p2, p3, p4, labels = c("a", "b", "c", "d"))
ggsave("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/Figure_cases_genus.pdf", width = 12, height = 7)





setwd('/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Data/mforge_output/')
files <- list.files()
files <- files[-grep('^para4',files)]
for(file in files){
  tm <- load(file)
  model.paras2 <- list(feature.dat=feature.dat,  ref.otu.tab=model.paras$ref.otu.tab)
  model.paras.genus2 <- list(feature.dat=feature.dat.genus,  ref.otu.tab=model.paras.genus$ref.otu.tab)
  model.paras <- model.paras2
  model.paras.genus <- model.paras.genus2
  save(model.paras, model.paras.genus, file=paste0('../../mPower/data/',file))
}
