require(MASS)
require(dirmult)

fit.NB<-function(x){
  x0=as.data.frame(x)
  x0=as.matrix(x0)
  notu=ncol(x0)
  sfs=rowSums(x0)
  params=matrix(NA,notu,ncol=2)
  colnames(params)=c('mu','phi')
  fits=list()
  for (i in 1:notu) {
    if(sum(x0[,i])>2){
      fit=try(glm.nb(x0[,i]~1+offset(log(sfs))),silent=T)
      fits[[i]]=fit
      if(methods::is(fit,"try-error")) {
        fit=try(glm(x0[,i]~1+offset(log(sfs)),family = poisson),silent=T)
        mu=exp(coef(fit))
        params[i,]=c(mu,0)
      }else{
        phi=1/fit$theta
        mu=exp(coef(fit))
        params[i,]=c(mu, phi)
      }
    }else{
      fit=try(glm(x0[,i]~1+offset(log(sfs)),family = poisson))
      mu=exp(coef(fit))
      params[i,]=c(mu,0)
    }
  }
  list(mu=params[,'mu'],phi=params[,'phi'],fits=fits)
}

fit.ZINB<-function(x){
  x0=as.data.frame(x)
  x0=as.matrix(x0)
  notu=ncol(x0)
  sfs=rowSums(x0)
  params=matrix(NA,notu,ncol=3)
  colnames(params)=c('mu','phi','p0')
  fits=list()
  for (i in 1:notu) {
    if(sum(x0[,i])>5){
      ctrl=pscl::zeroinfl.control(method = "L-BFGS-B")
      ctrl$reltol=NULL
      ctrl$factr=1e-3/.Machine$double.eps
      fit=try(pscl::zeroinfl(x0[,i]~1+offset(log(sfs))|1, dist = "negbin", link = "logit" , control=ctrl), silent=T)
      fits[[i]]=fit
      if(methods::is(fit,"try-error")) {
        fit=try(glm.nb(x0[,i]~1+offset(log(sfs))))
        phi=1/fit$theta
        mu=exp(coef(fit))
        params[i,]=c(mu,phi,0)
      }else{
        mu=exp(fit$coefficients$count) 
        theta=fit$theta
        p0=plogis(fit$coefficients$zero)
        params[i,]=c(mu,1/theta,p0)
      }
    }else{
      fit=try(glm.nb(x0[,i]~1+offset(sfs)))
      phi=1/fit$theta
      mu=exp(coef(fit))
      params[i,]=c(mu,phi,0)
    }
  }
  list(mu=params[,'mu'],phi=params[,'phi'],p0=params[,'p0'],fits=fits)
}

fit.ZIP<-function(x){
  x0=as.data.frame(x)
  x0=as.matrix(x0)
  notu=ncol(x0)
  sfs=rowSums(x0)
  params=matrix(NA,notu,ncol=2)
  colnames(params)=c('mu','p0')
  logliks=aics=rep(NA,notu)
  fits=list()
  
  for (i in 1:notu) {
    if(sum(x0[,i])>5){
      ctrl =pscl::zeroinfl.control(method = "L-BFGS-B")
      ctrl$reltol= NULL
      ctrl$factr=1e-3/.Machine$double.eps
      fit=try(pscl::zeroinfl(x0[,i]~1+offset(log(sfs))|1, dist = "poisson", link = "logit" , control=ctrl), silent=T)
      fits[[i]]=fit
      if(methods::is(fit,"try-error")) {
        fit=try(glm(x0[,i]~1+offset(log(sfs)),family = quasipoisson))
        mu=exp(coef(fit))
        params[i,]=c(mu,0)
      }else{
        mu=exp(fit$coefficients$count) 
        p0=plogis(fit$coefficients$zero)
        params[i,]=c(mu,p0)
      }
    }else{
      fit=try(glm(x0[,i]~1+offset(log(sfs)),family = quasipoisson))
      mu=exp(coef(fit))
      params[i,]=c(mu,0)
    }
  }
  list(mu=params[,'mu'],p0=params[,'p0'],fits=fits)
}

EstParaP <- function (otu.tab,  distrib=c('NB','ZINB','DM','ZIP')) { 
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste0('OTU', 1 : nrow(otu.tab))
  }
  otu.tab <- t(otu.tab)
  if(distrib=='DM'){
    dirmult.obj <- dirmult(otu.tab)
    phi <- dirmult.obj$theta
    pi <- dirmult.obj$pi
    ord <- order(pi, decreasing = TRUE) # order pi
    otu.tab <- otu.tab[, ord] # order the otu.tab based on OTU's pi
    pi <- pi[ord] # order pi 
  }
  
  if(distrib=='NB'){
    params <- fit.NB(otu.tab)
    phi <- params$phi
    pi <- params$mu
    names(pi) <- names(phi) <- colnames(otu.tab)
    ord <- order(pi, decreasing = TRUE) # order pi
    otu.tab <- otu.tab[, ord] # order the otu.tab based on OTU's pi
    pi <- pi[ord] # order pi 
    phi <- phi[ord]
  }
  
  if(distrib=='ZINB'){
    params <- fit.ZINB(otu.tab)
    phi <- params$phi
    pi <- params$mu
    p0 <- params$p0
    names(pi) <- names(phi) <- names(p0) <- colnames(otu.tab)
    ord <- order(pi, decreasing = TRUE) 
    otu.tab <- otu.tab[, ord] 
    pi <- pi[ord] 
    phi <- phi[ord]
    p0 <- p0[ord]
  }
  
  if(distrib=='ZIP'){
    params <- fit.ZIP(otu.tab)
    pi <- params$mu
    p0 <- params$p0
    names(pi) <- names(p0) <- colnames(otu.tab)
    ord <- order(pi, decreasing = TRUE) 
    otu.tab <- otu.tab[, ord] 
    pi <- pi[ord] 
    p0 <- p0[ord]
  }
  
  if(distrib=='ZINB'){
    return(list(lmu = log(pi), lphi = log(phi), lp0 = log(p0), otu.tab = otu.tab))
  }else if(distrib=='ZIP'){
    return(list(lmu = log(pi), lp0 = log(p0), otu.tab = otu.tab))
  }else{
    return(list(lmu = log(pi), lphi = log(phi), otu.tab = otu.tab))
  }
}

rzinbinom <- function(n, mu, theta, size, pi) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta=size
  rval=rnbinom(n, mu = mu, size = theta)
  rval[runif(n) < pi]=0
  rval
}

rzipois <- function(n, lambda, pi) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  rval=rpois(n, lambda = lambda)
  rval[runif(n) < pi]=0
  rval
}

simulate.data <- function (
  # Input data
  otu.tab, 
  nSam, nOTU, 
  
  # Count generation model
  model.paras = NULL, data.gen.model = c('DM','NB','ZINB','ZIP'), 
  
  # Signal structure
  diff.otu.pct = 0.1,  

  # Covariate effect
  grp.ratio = 1,  covariate.eff.mean = 1, covariate.eff.sd = 0, 

  # Sequence depth
  depth.mu = 5000, depth.theta = 5,  depth.conf.factor = 0
  
) {
  
  data.gen.model <- match.arg(data.gen.model)

  if (is.null(model.paras)) {
    otu.tab0 <- otu.tab
    otu.tab <- t(otu.tab)
    otu.tab.p <- otu.tab / rowSums(otu.tab)
    if(nOTU > ncol(otu.tab.p)){
      idx <- sample(ncol(otu.tab.p),nOTU, replace = T)
      otu.tab <- otu.tab[,idx]
    }else{
      otu.tab <- t(otu.tab[, order(colMeans(otu.tab.p), decreasing = TRUE)[1 : nOTU]])
    }
    
    if (data.gen.model == 'DM') {
      model.paras <- EstParaP(otu.tab, distrib = 'DM')
      }
    if (data.gen.model == 'ZINB') { 
      model.paras <- EstParaP(otu.tab, distrib = 'ZINB')
      }
    if (data.gen.model == 'NB') {
      model.paras <- EstParaP(otu.tab, distrib = 'NB')
    }
    if (data.gen.model == 'ZIP') {
      model.paras <- EstParaP(otu.tab, distrib = 'ZIP')
    }
  } else {
    if (nOTU <  length(model.paras$lmu)) {
      warning('Some rare OTUs will be dropped to achieve the requested OTU number!\n')
      model.paras$lmu <- model.paras$lmu[1 : (nOTU)]
      if(data.gen.model == 'ZINB'){
        model.paras$lphi <- model.paras$lphi[1 : (nOTU)]
        model.paras$lp0 <- model.paras$lp0[1 : (nOTU)]
      }
      if(data.gen.model == 'ZIP'){
        model.paras$lp0 <- model.paras$lp0[1 : (nOTU)]
      }
      if(data.gen.model == 'NB'){
        model.paras$lphi <- model.paras$lphi[1 : (nOTU)]
      }
    }else{
      model.paras$lmu <- sample(model.paras$lmu, nOTU, replace = T)
      if(data.gen.model == 'ZINB'){
        model.paras$lphi <- sample(model.paras$lphi, nOTU, replace = T)
        model.paras$lp0 <- sample(model.paras$lp0, nOTU, replace = T)
      }
      if(data.gen.model == 'ZIP'){
        model.paras$lp0 <- sample(model.paras$lp0, nOTU, replace = T)
      }
      if(data.gen.model == 'NB'){
        model.paras$lphi <- sample(model.paras$lphi, nOTU, replace = T)
      }
      
    }
  }
  
  X <- rnorm(nSam)
  X <- cbind(ifelse(X <= quantile(X, grp.ratio / (1 + grp.ratio)), 0, 1))

  # Generate log fold change for covariate effect- assume balanced change only
  eta.diff <- sample(c(rnorm(floor(nOTU / 2), mean = -covariate.eff.mean, sd = covariate.eff.sd), 
                       rnorm(nOTU - floor(nOTU / 2), mean = covariate.eff.mean, sd = covariate.eff.sd)))  %*% t(scale(X))

  # Get the differential OTU id
  otu.ord <- 1 : (nOTU)
  diff.otu.ind <- NULL
  diff.otu.num <- round(diff.otu.pct * nOTU)
  
  # assume random change of diff.otu
  diff.otu.ind <- c(diff.otu.ind, sample(otu.ord, diff.otu.num)) 

  eta.diff[setdiff(1 : nOTU, diff.otu.ind), ] <- 0

  if(data.gen.model == 'DM'){
    eta.exp <- exp(t(model.paras$lmu + eta.diff))
    otu.tab.prop <- eta.exp / rowSums(eta.exp)
    otu.tab.gamma <- otu.tab.prop * (1 - exp(model.paras$lphi)) / exp(model.paras$lphi)
    otu.tab.prop <- t(rdirichlet.m(t(otu.tab.gamma)))
    
    # Renormalize
    otu.tab.prop <- otu.tab.prop / rowSums(otu.tab.prop)
    otu.tab.prop <- t(otu.tab.prop)
    
    # Generate the sequence depth
    nSeq <- rnegbin(nSam, mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta) 
    otu.tab.sim <- sapply(1:ncol(otu.tab.prop), function (i) rmultinom(1, nSeq[i], otu.tab.prop[, i]))
  }

  if(data.gen.model =='NB'){
    # Generate the sequence depth
    lnSeq <- log(rnegbin(nSam, mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta))
    nSeq.m <- matrix(rep(lnSeq, nOTU), nrow = nOTU, byrow = TRUE)
    lmu0 <- matrix(rep(model.paras$lmu, nSam), ncol = nSam) 
    
    otu.exp <- exp(lmu0 + lnSeq + eta.diff) 
    otu.tab.sim <- rnbinom(n=length(otu.exp), mu=otu.exp, size=1/exp(model.paras$lphi))
    otu.tab.sim <- matrix(otu.tab.sim, ncol = nSam)
    }
  
  
  if(data.gen.model =='ZINB'){
    # Generate the sequence depth
    lnSeq <- log(rnegbin(nSam, mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta))
    nSeq.m <- matrix(rep(lnSeq, nOTU), nrow = nOTU, byrow = TRUE)
    lmu0 <- t(matrix(rep(model.paras$lmu, nSam), nrow = nSam, byrow = TRUE))
    
    otu.exp <- exp(lmu0 + lnSeq + eta.diff) 
    otu.tab.sim <- rzinbinom(n=length(otu.exp), mu=otu.exp,theta=1/exp(model.paras$lphi), pi=exp(model.paras$lp0))
    otu.tab.sim <- matrix(otu.tab.sim,ncol=nSam) 
  }
  
  if(data.gen.model =='ZIP'){
    # Generate the sequence depth
    lnSeq <- log(rnegbin(nSam, mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta))
    nSeq.m <- matrix(rep(lnSeq, nOTU), nrow = nOTU, byrow = TRUE)
    lmu0 <- t(matrix(rep(model.paras$lmu, nSam), nrow = nSam, byrow = TRUE))
    
    otu.exp <- exp(lmu0 + lnSeq + eta.diff) 
    otu.tab.sim <- rzipois(n=length(otu.exp),lambda=otu.exp,pi=exp(model.paras$lp0))
    otu.tab.sim <- matrix(otu.tab.sim,ncol=nSam) 
  }
  
  colnames(otu.tab.sim) <- paste0('S', 1 : ncol(otu.tab.sim))
  rownames(otu.tab.sim) <- paste0('O', 1 : nrow(otu.tab.sim))#names(model.paras$lmu)
  
  meta.sim <- cbind.data.frame(grp = X)
  rownames(meta.sim) <- colnames(otu.tab.sim)
  
  diff.otu.ind <- (1 : nOTU) %in% diff.otu.ind
  names(diff.otu.ind) <- rownames(otu.tab.sim)

  return(list(call = match.call(), model.paras = model.paras, otu.tab.sim = otu.tab.sim, otu.tab.used=otu.tab,
              diff.otu.ind = diff.otu.ind, meta.sim = meta.sim))
}

## copied from GUniFrac package, add model.paras
SimulateMSeq <- function(ref.otu.tab, model.paras, nSam = 100, nOTU = 500, diff.otu.pct = 0.1, 
          diff.otu.direct = c("balanced", "unbalanced"), diff.otu.mode = c("abundant", "rare", "mix"), 
          covariate.type = c("binary", "continuous"), 
          grp.ratio = 1, covariate.eff.mean = 1, covariate.eff.sd = 0, 
          confounder.type = c("none", "binary", "continuous", "both"), 
          conf.cov.cor = 0.6, conf.diff.otu.pct = 0, conf.nondiff.otu.pct = 0.1, 
          confounder.eff.mean = 0, confounder.eff.sd = 0, error.sd = 0, 
          depth.mu = 10000, depth.theta = 5, depth.conf.factor = 0) 
{
  diff.otu.direct <- match.arg(diff.otu.direct)
  diff.otu.mode <- match.arg(diff.otu.mode)
  covariate.type <- match.arg(covariate.type)
  confounder.type <- match.arg(confounder.type)
  ref.otu.tab0 <- ref.otu.tab
  #model.paras <- EstPara(ref.otu.tab = ref.otu.tab)
  sample.names <- colnames(model.paras$ref.otu.tab)
  ref.otu.tab <- model.paras$ref.otu.tab[(1:(nOTU)), ]
  idx.otu <- rownames(ref.otu.tab)
  idx.sample <- sample(sample.names, nSam)
  idx.nonsample <- colnames(ref.otu.tab)[!(colnames(ref.otu.tab) %in% idx.sample)]
  ref.otu.tab = ref.otu.tab[, idx.sample]
  ref.otu.tab.unselect = ref.otu.tab0[c(1:(nOTU)), ][, idx.nonsample]
  if (confounder.type == "none") {
    confounder.type <- "continuous"
    confounder.eff.mean <- 0
    confounder.eff.sd <- 0
    Z <- NULL
  }
  if (confounder.type == "continuous") 
    Z <- cbind(rnorm(nSam))
  if (confounder.type == "binary") 
    Z <- cbind(c(rep(0, nSam%/%2), rep(1, nSam - nSam%/%2)))
  if (confounder.type == "both") 
    Z <- cbind(rnorm(nSam), c(rep(0, nSam%/%2), rep(1, nSam - 
                                                      nSam%/%2)))
  rho <- sqrt(conf.cov.cor^2/(1 - conf.cov.cor^2))
  if (covariate.type == "continuous") {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + rnorm(nSam)
  }
  if (covariate.type == "binary") {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + rnorm(nSam)
    X <- cbind(ifelse(X <= quantile(X, grp.ratio/(1 + grp.ratio)), 
                      0, 1))
  }
  rownames(X) <- colnames(ref.otu.tab)
  covariate.eff.mean1 = covariate.eff.mean
  covariate.eff.mean2 = covariate.eff.mean
  if (diff.otu.direct == "balanced") {
    if (diff.otu.mode == "abundant") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd),
                           rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*% t(scale(X))
    }
    else if (diff.otu.mode == "rare") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, 
                                 sd = covariate.eff.sd), rnorm(nOTU - floor(nOTU/2), 
                                                               mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*% 
        t(scale(X))
    }
    else {
      eta.diff <- c(sample(c(rnorm(floor(nOTU/4), mean = -covariate.eff.mean1, 
                                   sd = covariate.eff.sd), rnorm(floor(nOTU/2) - 
                                                                   floor(nOTU/4), mean = covariate.eff.mean1, sd = covariate.eff.sd))), 
                    sample(c(rnorm(floor((nOTU - floor(nOTU/2))/2), 
                                   mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                             rnorm(nOTU - floor(nOTU/2) - floor((nOTU - 
                                                                   floor(nOTU/2))/2), mean = covariate.eff.mean2, 
                                   sd = covariate.eff.sd)))) %*% t(scale(X))
    }
  }
  if (diff.otu.direct == "unbalanced") {
    if (diff.otu.mode == "abundant") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, 
                        sd = covariate.eff.sd) %*% t(scale(X))
    }
    else if (diff.otu.mode == "rare") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, 
                        sd = covariate.eff.sd) %*% t(scale(X))
    }
    else {
      eta.diff <- c(sample(c(rnorm(floor(nOTU/2), mean = covariate.eff.mean1, 
                                   sd = covariate.eff.sd))), sample(c(rnorm(nOTU - 
                                                                              floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*% 
        t(scale(X))
    }
  }
  eta.conf <- sample(c(rnorm(floor(nOTU/2), mean = -confounder.eff.mean, 
                             sd = confounder.eff.sd), rnorm(nOTU - floor(nOTU/2), 
                                                            mean = confounder.eff.mean, sd = confounder.eff.sd))) %*% 
    t(scale(scale(Z) %*% rep(1, ncol(Z))))
  otu.ord <- 1:(nOTU)
  diff.otu.ind <- NULL
  diff.otu.num <- round(diff.otu.pct * nOTU)
  if (diff.otu.mode == "mix") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord, diff.otu.num))
  if (diff.otu.mode == "abundant") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[1:round(length(otu.ord)/4)], 
                                           diff.otu.num))
  if (diff.otu.mode == "rare") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[round(3 * 
                                                           length(otu.ord)/4):length(otu.ord)], diff.otu.num))
  if (length(diff.otu.ind) >= round(nOTU * conf.diff.otu.pct)) {
    conf.otu.ind1 <- sample(diff.otu.ind, round(nOTU * conf.diff.otu.pct))
  }
  else {
    conf.otu.ind1 <- diff.otu.ind
  }
  conf.otu.ind <- c(conf.otu.ind1, sample(setdiff(1:(nOTU), 
                                                  diff.otu.ind), round(conf.nondiff.otu.pct * nOTU)))
  eta.diff[setdiff(1:(nOTU), diff.otu.ind), ] <- 0
  eta.conf[setdiff(1:(nOTU), conf.otu.ind), ] <- 0
  eta.error <- matrix(rnorm(nOTU * nSam, 0, error.sd), nOTU, 
                      nSam)
  eta.exp <- exp(t(eta.diff + eta.conf + eta.error))
  eta.exp <- eta.exp * t(ref.otu.tab)
  ref.otu.tab.prop <- eta.exp/rowSums(eta.exp)
  ref.otu.tab.prop <- t(ref.otu.tab.prop)
  nSeq <- rnegbin(nSam, mu = depth.mu * exp(scale(X) * depth.conf.factor), 
                  theta = depth.theta)
  otu.tab.sim <- sapply(1:ncol(ref.otu.tab.prop), function(i) rmultinom(1, 
                                                                        nSeq[i], ref.otu.tab.prop[, i]))
  colnames(otu.tab.sim) <- rownames(eta.exp)
  rownames(otu.tab.sim) <- rownames(ref.otu.tab)
  diff.otu.ind = (1:nOTU) %in% diff.otu.ind
  conf.otu.ind = (1:nOTU) %in% conf.otu.ind
  return(list(otu.tab.sim = otu.tab.sim, covariate = X, confounder = Z, 
              diff.otu.ind = diff.otu.ind, otu.names = idx.otu, conf.otu.ind = conf.otu.ind))
}


## ==== Compare simulation models ====
pkg <- c('dplyr','ggplot2','RColorBrewer','tibble','reshape2','GUniFrac')
sapply(pkg, function(x) require(x, character.only = T))
dd <- '/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Data/mforge_output/'
rd <- '/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/'

setwd(dd)
file.names <- gsub('para4_','',list.files(pattern = '^para4'))
for(file.name in file.names){
  setwd(dd)
  para4 <- NULL
  load(paste0('para4_',file.name))
  
  load(file.name)
  data.gen.models <- c('DM','NB','ZINB','ZIP','Semi')
  iters <- 100
  nSams <- c(20, 40, 60)
  TPR <- FDR <- array(NA, dim = c(length(nSams), iters, length(data.gen.models)),
                      dimnames = list(nSams = paste0('S',nSams), iters = paste0('iter',1:iters), methods = data.gen.models))
  for(data.gen.model in data.gen.models){
    for(i in 1:iters){
      for(nSam in nSams){
        nOTU <- 200 #nrow(feature.dat)
        if(data.gen.model %in% c('DM','NB','ZINB','ZIP')){
          sim.obj <- simulate.data(otu.tab = t(para4[[data.gen.model]]$otu.tab), model.paras = para4[[data.gen.model]],
                                   nSam = nSam, nOTU = nOTU, data.gen.model = data.gen.model,
                                   diff.otu.pct = 0.1, depth.mu = 10000,
                                   covariate.eff.mean = 1)
          otu.tab.sim <- otu_filter(feature.dat = sim.obj$otu.tab.sim)
          truth <- sim.obj$diff.otu.ind
          fit <- NULL
          try({fit <- MicrobiomeStat::linda(feature.dat = otu.tab.sim, meta.dat = sim.obj$meta.sim, formula = '~ grp',prev.filter = 0.2, max.abund.filter = 0.002)})
          try({temp <- cal.tpr.fdr(p.adj.fdr = fit$output$grp$padj, truth = truth[rownames(fit$output$grp)], cutoff = 0.05)})
        }else{
          sim.obj <- SimulateMSeq(ref.otu.tab = feature.dat, model.paras = model.paras,
                                   nSam = nSam, nOTU = nOTU, diff.otu.pct = 0.1,
                                    diff.otu.direct = "balanced", diff.otu.mode ="mix",
                                    covariate.type = "binary", grp.ratio = 1, covariate.eff.mean = 1, covariate.eff.sd = 0,
                                    confounder.type = "none",
                                    conf.cov.cor = 0.6, conf.diff.otu.pct = 0, conf.nondiff.otu.pct = 0.1,
                                    confounder.eff.mean = 0, confounder.eff.sd = 0, error.sd = 0,
                                    depth.mu = 10000, depth.theta = 5, depth.conf.factor = 0)
          try({truth <- sim.obj$diff.otu.ind; names(truth) <- sim.obj$otu.names})
          fit <- NULL
          try({fit <- MicrobiomeStat::linda(feature.dat = sim.obj$otu.tab.sim, 
                                            meta.dat = as.data.frame(sim.obj$covariate) %>% 
                                              dplyr::rename(grp =1) %>% 
                                              mutate(grp= as.factor(grp)), 
                                            prev.filter = 0.2, max.abund.filter = 0.002,
                                            formula = '~grp')})
          try({temp <- cal.tpr.fdr(p.adj.fdr = fit$output$grp$padj, truth = truth[rownames(fit$output$grp)], cutoff = 0.05)})
        }

        TPR[paste0('S',nSam),paste0('iter',i),data.gen.model] <- temp['TPR']
        FDR[paste0('S',nSam),paste0('iter',i),data.gen.model] <- temp['FDR']
      }
    }
  }

  setwd(rd)
  save(TPR, FDR, file = paste0('model_comparison_',gsub('\\(0.05prev\\).RData','',file.name),'(0.05prev)_prevmaxabund_nOTU200.Rdata'))
}




setwd(rd)
file.names <- list.files(pattern = '_prevmaxabund_nOTU200.Rdata$')
DF <- NULL
col <- c(brewer.pal(8,'Set1')[1],brewer.pal(8,'Paired')[1:4])
names(col) <- c('Semi','DM','ZINB','NB','ZIP')
for(file.name in file.names){
  setwd(rd)
  load(file.name)
  if(length(grep("Callahan2017_Vaginal", file.name))==1) name = 'Human vaginal'
  if(length(grep("Daniel2018_AGP_UK", file.name))==1) name = 'Human stool'
  if(length(grep("GermanLakeWater", file.name))==1) name = 'Lake water'
  if(length(grep("Huttenhower_HMPv13OralCavity_ID1927", file.name))==1) name = 'Human oral cavity'
  
  if(length(grep("Huttenhower_HMPv35Skin_ID1928", file.name))==1) name = 'Human skin'
  if(length(grep("NoahFierer_Soil_ID2140", file.name))==1) name = 'Soil'
  if(length(grep("PlantSourface_ID2229", file.name))==1) name = 'Plant surface'
  if(length(grep("Wu2019_WWTP", file.name))==1) name = 'Wastewater'
  cat(name)
  tmp <- as.data.frame(t(as.data.frame(TPR))) %>% 
    mutate(method = gsub('.*\\.','',rownames(.))) %>% 
    reshape2::melt() %>% mutate(variable = as.numeric(gsub('S','', variable))) %>% mutate(source = name)
  DF <- rbind(DF, tmp)
}

DF$variable <- paste0('n=',DF$variable)
p1 <- ggplot(DF[DF$source %in% unique(DF$source)[grep('Human',unique(DF$source))],], aes(x = reorder(method, -value), y = value, fill = method)) + 
  stat_boxplot(geom = "errorbar", width = 0.2) + 
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values=col) + 
  facet_grid(source~variable, scales = 'free', space = 'free') + 
  theme_bw() + 
  labs(x = '', y = 'Power', fill = '') + 
  theme(axis.text = element_text(size = 14, colour = 'black'), 
        axis.title = element_text(size = 20, colour = 'black'), 
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, colour = 'black', angle = 90, hjust = 1), 
        legend.position = 'none',
        axis.ticks.x = element_blank(), 
        strip.text = element_text(size = 14, colour = 'black'))


p2 <- ggplot(DF[!(DF$source %in% unique(DF$source)[grep('Human',unique(DF$source))]),], aes(x = reorder(method, -value), y = value, fill = method)) + 
  stat_boxplot(geom = "errorbar", width = 0.2) + 
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values=col) + 
  facet_grid(source~variable, scales = 'free', space = 'free') + 
  theme_bw() + 
  labs(x = '', y = 'Power', fill = '') + 
  theme(axis.text = element_text(size = 14, colour = 'black'), 
        axis.title = element_text(size = 20, colour = 'black'), 
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, colour = 'black', angle = 90, hjust = 1), 
        legend.position = 'none',
        axis.ticks.x = element_blank(), 
        strip.text = element_text(size = 14, colour = 'black'))
ggpubr::ggarrange(p1, p2, labels = c('A','B'), font.label = list(size = 24, color = "black"))
setwd(rd)
ggsave(file = 'ModelsComp_200(0.05prev)_linda.pdf', width = 16, height = 8)



### Table for reference dataset info
load("/Users/luyang1/myicloud/Documents/Mayo_Research/2023_05_06_PowerShiny/Result/rstudio-export/Callahan2017_Vaginal(0.05prev)2.RData")
dim(feature.dat.genus)
dim(feature.dat)

setwd(dd)
file.names <- list.files(pattern = 'RData$')
load('Callahan2017_Vaginal(0.05prev).RData')
dim(feature.dat.genus)
dim(feature.dat.genus2)
dim(feature.dat)
file.names <- file.names[-grep('^para4_',file.names)]
sam <- matrix(NA, nrow=length(file.names), ncol=6)
rownames(sam) <- gsub('\\(.*','',file.names)
colnames(sam) <- c('Number of features(Species)', 'Number of samples(Species)','Reads depth(Species)',
                   'Number of features(Genus)', 'Number of samples(Genus)','Reads depth(Genus)')
for(file.name in file.names){
  cat('----',file.name,'---\n')
  rm(feature.dat);rm(feature.dat.genus)
  tm <- load(file.name)
  print(min(colSums(feature.dat)))
  if(min(rowSums(feature.dat>0))!=0){
    sam[gsub('\\(.*','',file.name),'Number of features(Genus)'] <- nrow(feature.dat.genus)
    sam[gsub('\\(.*','',file.name),'Number of samples(Genus)'] <- ncol(feature.dat.genus)
    sam[gsub('\\(.*','',file.name),'Reads depth(Genus)'] <- paste0(min(colSums(feature.dat.genus)),'-',max(colSums(feature.dat.genus)),'(',round(median(colSums(feature.dat.genus))),')')
    
    sam[gsub('\\(.*','',file.name),'Number of features(Species)'] <- nrow(feature.dat)
    sam[gsub('\\(.*','',file.name),'Number of samples(Species)'] <- ncol(feature.dat)
    sam[gsub('\\(.*','',file.name),'Reads depth(Species)'] <- paste0(min(colSums(feature.dat)),'-',max(colSums(feature.dat)),'(',round(median(colSums(feature.dat))),')')
  }
}

rownames(sam) <- gsub("Callahan2017_Vaginal",'Human vaginal (Callahan, Benjamin J., et al.,2017)',rownames(sam))
rownames(sam) <- gsub("Daniel2018_AGP_UK",'Human stool (McDonald, Daniel, et al.,2018)',rownames(sam))
rownames(sam) <- gsub("GermanLakeWater",'Lake water (Gonzalez, A., et al., 2018)',rownames(sam))
rownames(sam) <- gsub("Huttenhower_HMPv13OralCavity_ID1927",'Human oral cavity (The Human Microbiome Project Consortium.,2012)',rownames(sam))
rownames(sam) <- gsub("Huttenhower_HMPv35Skin_ID1928",'Human skin (The Human Microbiome Project Consortium.,2012)',rownames(sam))
rownames(sam) <- gsub("NoahFierer_Soil_ID2140",'Soil (Ramirez, K.S., et al.,2014)',rownames(sam))
rownames(sam) <- gsub("PlantSourface_ID2229",'Plant surface (Marzinelli, E.M., et al., 2015)',rownames(sam))
rownames(sam) <- gsub("Wu2019_WWTP",'Wastewater (Wu, L.W., et al., 2019)',rownames(sam))
sam <- as.data.frame(sam) %>% rownames_to_column('Study')
DT::datatable(sam)
# write.csv(sam,file='../Table2.csv', row.names = F)
