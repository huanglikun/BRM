#!/usr/bin/env RScript

# version 0.1
##############################################
# command
##############################################
opt <- commandArgs(TRUE)
# help
if (length(opt) < 3) {cat("Usage: Rscript BRMstep3.peak_and_CI.R AFD.xls threshold.xls allpeaks.xls\n");q();}
#
file_afd <- opt[1];
file_var <- opt[2];
file_out <- opt[3];

# fixed global options
tol <- 0.01 # tolerance threshold

#########################################
# functions
#########################################
#
opt.span <- function(model, criterion = c("aicc", "gcv"), span.range = c(0.05, 0.95)) {
  as.crit <- function(x) {
    span <- x$pars$span;
    traceL <- x$trace.hat;
    sigma2 <- sum(x$residuals^2)/(x$n - 1);
    aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL - 2);
    gcv <- x$n * sigma2/(x$n - traceL)^2;
    result <- list(span = span, aicc = aicc, gcv = gcv);
    return(result);
  }
  criterion <- match.arg(criterion);
  fn <- function(span) {
    mod <- update(model, span = span);
    as.crit(mod)[[criterion]];
  }
  result <- optimize(fn, span.range);
  return(list(span = result$minimum, criterion = result$objective));
}

#
locfit_by_loess <- function(x, y){
  # loess+AICc
  fit0  <- loess(y~x, degree = 2)
  span1 <- opt.span(fit0, criterion="aicc")$span
  fit1  <- loess(y~x, degree=2, span=span1)
}

#
return_cross_pos <- function(all_plots, x1, x2, y, tol){
	subregion  <- all_plots[which(all_plots[,1]>x1 & all_plots[,1]<x2),];
	y_distance <- abs(subregion[,2] - y);
	pos        <- subregion[which(y_distance==min(y_distance) & y_distance<tol),1];
}

#########################################
# main
#########################################

# read into memory, as data frame
dat_afd <- read.table(file_afd, as.is=TRUE);
dat_var <- read.table(file_var, as.is=TRUE);

# one-by-one for every chromosome
chrom   <- c("#Chromosome");
peak_x  <- c("Peak Position");
peak_y  <- c("Peak Value");
type    <- c("Peak Direction");
left_x  <- c("Interval Start");
right_x <- c("Interval End");

# title write to file
out <- data.frame(chrom, peak_x, peak_y, type, left_x, right_x);

write.table(out, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=FALSE);

#
dat_chr <- dat_afd[,1];
dat_pos <- dat_afd[,2];
dat_avg <- dat_afd[,3];

dat_std <- sqrt(dat_var[,6]);
qtl_std <- sqrt(dat_var[,8]);

chr     <- data.frame(unique(dat_chr));
#
for(i in 1:length(chr[,1])){
  thischr <- chr[i,1];
  idx <- which(dat_chr == chr[i,1] & !is.na(dat_avg));
  total <- length(idx);
  cat("data size of", chr[i,1], "is", total, "\n");
  if(total==0) next;

  x <- dat_pos[idx];
  y <- dat_avg[idx];
  std  <- dat_std[idx];
  qstd <- qtl_std[idx];
  fit1 <- locfit_by_loess(x, y);

  newx <- seq(x[1],x[length(x)],100);
  newx_len <- length(newx);
  plo  <- data.frame(x=newx, y=predict(fit1, newx));
  peakxs   <- newx[which(diff(diff(plo$y)>0)!=0L)+1L];
  len  <- length(peakxs);
  if (len==0){
    k  <- ifelse(diff(c(plo$y[1],plo$y[length(plo)]))>0L,1L,-1L);
    ifelse(k>0L,k <- c(1L,-1L),k <- c(-1L,1L));
    pks <- data.frame(peak_x=c(plo$x[1],plo$x[length(plo$x)]),peak_y=c(plo$y[1],plo$y[length(plo$y)]),type=k);
  }else{
    k    <- diff(diff(plo$y)>0L);
    k    <- k[which(k!=0L)];
    stat     <- k[1]>0;
    stat0    <- stat;
    
    peakxs   <- (peakxs[1:len-1]+peakxs[2:len])/2;
    
    pks  <- sapply(as.data.frame(rbind(c(newx[1],peakxs), c(peakxs,newx[length(newx)]))),
                   function(k){ stat <<- !stat;
                   optimize(f=function(x){
                     predict(fit1, newdata=data.frame(x))
                   }, maximum=stat, interval=k)});
    pks  <- as.data.frame(t(pks));
    names(pks) <- c("peak_x","peak_y");
    pks$peak_x <- as.numeric(pks$peak_x);
    pks$peak_y <- as.numeric(pks$peak_y);
    pks$type   <- k;
 
  # add the first and the last position.
  firstPos   <- c(plo$x[1],plo$y[1],-1*pks$type[1]);
  lastPos    <- c(plo$x[length(plo$x)],plo$y[length(plo$y)],-1*pks$type[length(pks$type)]);
  pks        <- rbind(firstPos, pks, lastPos);
  k          <- pks$type;
  }
  
  all_peak_plots <- pks$peak_x;
  all_ident_y    <- pks$peak_y + k*1.65*pks$nearby_std;
  all_x1         <- all_peak_plots[1:(length(all_peak_plots)-1)];
  all_x3         <- floor(all_peak_plots[1:(length(all_peak_plots)-1)]);
  all_x2         <- all_peak_plots[2:(length(all_peak_plots))];
  all_x4         <- all_peak_plots[2:length(all_peak_plots)];
  pos1_list <- c(x[1])
  pos2_list <- c()
  all_x1_len<- length(all_x1);
  for (i in 1:all_x1_len){
	pos1         <- return_cross_pos(plo, all_x1[i], all_x2[i], all_ident_y[i+1], tol);
	j <- i ;
	while(j > 1 & length(pos1) == 0){
		j <- j - 1;
		pos1 <- return_cross_pos(plo, all_x1[j], all_x2[i], all_ident_y[i+1], tol);
	}
	pos2         <- return_cross_pos(plo, all_x3[i], all_x4[i], all_ident_y[i], tol);
	j <- i;
	while(j < all_x1_len & length(pos2) == 0){
		j <- j + 1;
		pos2 <- return_cross_pos(plo, all_x3[i], all_x4[j], all_ident_y[i], tol);
	}
    if (length(pos1) == 0) pos1 <- c(all_x1[1]);
    if (length(pos2) == 0) pos2 <- c(all_x4[length(all_x4)]);
    pos1_list  <- c(pos1_list,pos1);
    pos2_list  <- c(pos2_list,pos2);
  }
  pos2_list    <- c(pos2_list,x[length(x)]);
  pks$type   <- ifelse(k==1,"-","+");
  pks <-data.frame(thischr, pks,left_x=c(pos1_list), right_x=c(pos2_list));

  write.table(pks, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE);
}
