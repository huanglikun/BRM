#!/usr/bin/env RScript

# version 0.1
##############################################
# command
##############################################
opt <- commandArgs(TRUE);
# help
if (length(opt) < 1) {cat("Usage: Rscript cal_ua_rk.R ua_rk_conf.txt \n");q();}
#
file_conf <- opt[1]; # arguments configure file

# fixed global options
tol   <- 1e-13; # tolerance threshold
rstt  <- 0.001;
rend  <- 0.083;
cmstt <- 1;
cmend <- 10;
rk <- 0.083;

#########################################
# functions
#########################################
#
read_conf <- function(file_conf){
  conf <- list();
  #
  FH <- file(file_conf,"r");
  while(TRUE){
    rline <- readLines(FH,n=1);
    if (length(rline) == 0) break;
    if (length (grep ("^#", rline, perl=TRUE) ) > 0 ) next;
    if (length (grep ("^\\s", rline, perl=TRUE) ) > 0) next;
    if (nchar (rline) == 0) next;
    rline <- sub ("#.*", "", rline, perl=TRUE);
    rline <- sub ("\\s+$", "", rline, perl=TRUE);
    rline <- sub ("\\s*?=\\s*", "=", rline, perl=TRUE);
    arr <- strsplit (rline, "=")[[1]];
    conf[[ arr[1] ]] <- arr[2];
  }
  close(FH);
  return(conf);
}
# recursion
f <- function(k,r1){
  if(k==1){
    rnext <- r1;
  }else{
    rnext <- f(k-1,r1) + 0.5*(1-f(k-1,r1))^2*r1;
  }
}
# 
f_bisect <- function(x,parms){
  return(f(parms[1],x)-parms[2]);
}

# Bisection function is from GuangchuangYu
bisect <- function(f, u, v, eps, parms){    
  if ( f(u,parms) * f(v,parms) > 0)
    stop (" error, select another value for u or v...");
  if ( f(u,parms) < 0) {
    u1 <- u;
    v1 <- v;
  } else {
    u1 <- v;
    v1 <- u;
  }
  cnt <- 1;
  step <- log2(abs(u1-v1)/eps);
  while ( cnt < step ){
    m <- (u1+v1)/2;
    if ( f(m,parms) == 0)
      break;
    if ( f(m,parms) < 0 ) {
      u1 <- m;
    } else {
      v1 <- m;
    }
    cnt <- cnt + 1;
  }     
  return (m);
}

cal_cm <- function(x,parms){
  return(0.5 * (exp(x / 25) - 1) / (exp(x / 25) + 1) - parms[1])
}

##############################################
# main
##############################################
cat("Read global parameters from", file_conf, ".\n");
conf <- read_conf(file_conf);
for(name in names(conf) ){
  cat("Parameter ", name, "=", conf[[name]],  "\n", sep="");
}
cat("\n");

# parameters
if (!is.null(conf$precision)) tol <- as.numeric(conf$precision);
if (!is.null(conf$k)) k <- as.numeric(conf$k);
if (!is.null(conf$chrNum)){
  chrNum <- as.numeric(conf$chrNum);
}else{
  stop("Undefined chrNum.");
}
if (!is.null(conf$geneticLength)){
  geneticLength <- as.numeric(conf$geneticLength);
}else{
  if (!is.null(conf$genomeSize) && !is.null(conf$ratio)){
    genomeSize    <- as.numeric(conf$genomeSize);
    ratio         <- as.numeric(conf$ratio);
    geneticLength <- genomeSize/ratio;
  }else{
    stop("Undefined genomeSize or ratio.");
  }
}

r_root <- bisect(f_bisect,rstt,rend,tol,parms=c(k,rk));
r_root <- round(r_root,digits = 5);

cat("recombination rate is ",r_root,"\n");

cm_root <- bisect(cal_cm,cmstt,cmend,tol,parms = c(r_root));
cm_root <- round(cm_root,digits = 5);
cat("Genetic distance is ",cm_root," cM\n");

#
a5 <- 0.05/(chrNum+geneticLength/cm_root);
#a1 <- 0.01/(chrNum+geneticLength/cm_root);
ua <- round(abs(qnorm(a5/2)),digits = 6);
cat("ua is ",ua,"\n");
