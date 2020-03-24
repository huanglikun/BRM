#!/usr/bin/env RScript

# version 0.1
##############################################
# command
##############################################
opt <- commandArgs(TRUE);
# help
if (length(opt) < 1) {cat("Usage: Rscript cal_ua_fk.R ua_fk_conf.txt \n");q();}
#
file_conf <- opt[1]; # arguments configure file

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

r_root  <- 0.084;
cm_root <- 8.4;
if (k == 0 | k > 4){
  r_root  <- 0.0459;
  cm_root <- 4.5;
}else if (k == 1 | k == 2){
  r_root  <- 0.084;
  cm_root <- 8.4;
}else if (k == 3){
  r_root  <- 0.0582;
  cm_root <- 5.8;
}else if (k == 4){
  r_root  <- 0.0508;
  cm_root <- 5.0;
}else{
  cat("Not proper k?\n");
  q();
}

#
a5 <- 0.05/(chrNum+geneticLength/cm_root);
#a1 <- 0.01/(chrNum+geneticLength/cm_root);
ua <- round(abs(qnorm(a5/2)),digits = 6);
cat("ua is ",ua,"\n");
