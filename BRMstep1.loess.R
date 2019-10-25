#!/usr/bin/env RScript

# version 0.1
##############################################
# command
##############################################
opt <- commandArgs(TRUE);
# help
if (length(opt) < 5) {cat("Usage: Rscript BRMstep1.loess.R loess_conf.txt chr_length.txt markers.bsa AFD.xls|AF1.xls|AF2.xls|AAF.xls AFD|AF1|AF2|AAF\n");q();}
#
file_conf <- opt[1]; # arguments configure file
file_chr  <- opt[2]; # chr name and length list
file_tab  <- opt[3]; # input(bsa file)
file_out  <- opt[4]; # output
# STATISTIC support AFD, absAFD, AF1, AF2, MAF, AAF(AFM)
STATISTIC <- as.character(opt[5]);

##############################################
# functions
##############################################
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

# code from fANCOVO package
# best smoothing parameter 
opt.span <- function(model, criterion = c("aicc", "gcv"), span.range = c(0.05, 0.95)) {
    as.crit <- function(x) {
        span   <- x$pars$span;
        traceL <- x$trace.hat;
        sigma2 <- sum(x$residuals^2)/(x$n - 1);
        aicc   <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL - 2);
        gcv    <- x$n * sigma2/(x$n - traceL)^2;
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

# block's middle position 
def_block_pos <- function(chr, size){
	# block number
	n <- as.integer(2*chr/size)+2;
	if( n%%2 != 0 ) n <- n+1;
	n <- as.integer(n/2);
	# block index and the middle position of each block
	i <- c(1:n);
	pos <- (i-1)*size;
	if(pos[n]>chr) pos[n]<-chr;
	return(pos);
}

# meta value of block
cal_block_meta <- function(loc, val, chr, size, depth){
	# input: location vector, value vector
	# input: chr length, block size, location depth vector
	pos <- def_block_pos(chr, size);
	idx <- as.integer(0.5+loc/size)+1;
	#
	num <- c();
	avg <- c();
	blockDepth <- c();
	for(i in 1:length(pos)){
		k  <- which(idx==i);
		no <- length(k);
		a <- NA;
		n <- 0;
		if (no > 0) {n <- sum(depth[k])};
		if( n >= MIN ) a <- sum(val[k])/n;
		num <- c(num, no);
		avg <- c(avg, a);
		blockDepth <- c(blockDepth, n);
	}
	return( list(pos=pos, avg=avg, num=num, blockDepth=blockDepth) );
}

check_variable <- function(A,B,C,D){
	idx <- c();
	# not missing
	idx <- c(idx, which(is.na(A)) );
	idx <- c(idx, which(is.na(B)) );
	idx <- c(idx, which(is.na(C)) );
	idx <- c(idx, which(is.na(D)) );
	#
	idx <- c(idx, which(!is.integer(A)) );
	idx <- c(idx, which(!is.integer(B)) );
	idx <- c(idx, which(!is.integer(C)) );
	idx <- c(idx, which(!is.integer(D)) );
	#
	idx <- c(idx, which(A+B==0) );
	idx <- c(idx, which(A+C==0) );
	idx <- c(idx, which(B+D==0) );
	idx <- c(idx, which(C+D==0) );
	#
	idx <- unique(idx);
	return(idx);
}

##############################################
# main
##############################################
ptm <- proc.time(); # begin to record runtime 

# read configure file
cat("Read global parameters from", file_conf, ".\n");
conf <- read_conf(file_conf);
for(name in names(conf) ){
	cat("Parameter ", name, "=", conf[[name]],  "\n", sep="");
}
cat("\n");

# global options
BLK <- as.numeric(conf$BLK) * as.numeric(conf$UNIT); # block size/step
MIN <- as.numeric(conf$MIN); # least total depth in one block 
DEG <- as.numeric(conf$DEG); # degree of polynomial
MINVALID <- as.numeric(conf$MINVALID); # least valid blocks in one chromosome

if(MIN==0) MIN <- 1;
if(MINVALID < 10) MINVALID <- 10;

# read into memory
tab <- read.table(file_tab, as.is=TRUE);
chr <- read.table(file_chr, as.is=TRUE);

# one-by-one for every chromosome
dat_chr <- c("#Chr.");
dat_pos <- c("Position");
dat_avg <- c("Block Average");
fit_avg <- c("Fitted Average");
fit_ser <- c("Std. Error");
dat_num <- c("No. of Markers");
blk_dph <- c("Depth in block");
# title write to file
cat("Export data to", file_out, ".\n");
if(STATISTIC=="AF1" || STATISTIC=="AF2"){
	out <- data.frame(dat_chr, dat_pos, dat_avg, fit_avg, fit_ser, dat_num, blk_dph);
}else{
	out <- data.frame(dat_chr, dat_pos, dat_avg, fit_avg, fit_ser, dat_num);
}
	
write.table(out, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=FALSE);
#
dat_chr <- tab[,1];
dat_pos <- tab[,2];
A 		<- tab[,3];
B 		<- tab[,4];
C 		<- tab[,5];
D 		<- tab[,6];
##
idx <- check_variable(A, B, C, D);
cat("Number of missing data:", length(idx), "\n");
dat_all_chr <- dat_chr;
dat_all_pos <- dat_pos;
if(length(idx)>0){
	A <- A[-idx];
	B <- B[-idx];
	C <- C[-idx];
	D <- D[-idx];
	dat_all_chr <- dat_chr[-idx];
	dat_all_pos <- dat_pos[-idx];
}
## 
pool1Depth <- A + B;
pool2Depth <- C + D;
#
cat("Number of imported chromosomes is", length(chr[,1]), ".\n");
for(i in 1:length(chr[,1])){
	dat_chr <- dat_all_chr;
	dat_pos <- dat_all_pos;
	idx <- which(dat_chr == chr[i,1]);
	total <- length(idx);
	cat("Data size of", chr[i,1], "is", total, ".\n");
	if(total==0) next;
	#
	cat("Calculate average value for each block.\n");
	block_af1 <- cal_block_meta(dat_pos[idx], A[idx], chr[i, 2], BLK, pool1Depth[idx]);
	block_af2 <- cal_block_meta(dat_pos[idx], C[idx], chr[i, 2], BLK, pool2Depth[idx]);

	# initialize
	block <- list();
	#
	if(STATISTIC=="AF1") {
		block <- block_af1;
	}else if (STATISTIC=="AF2") {
		block <- block_af2;
	}else if (STATISTIC=="AFD") {
		block$avg <- block_af1$avg - block_af2$avg;
		block$pos <- block_af1$pos;
		block$num <- block_af1$num;
	}else if (STATISTIC=="absAFD") {
		block$avg <- abs(block_af1$avg - block_af2$avg);
		block$pos <- block_af1$pos;
		block$num <- block_af1$num;
	}else if (STATISTIC=="MAF"){
		block$avg <- 0.5 - abs(0.5 - block_af1$avg);
		block$pos <- block_af1$pos;
		block$num <- block_af1$num;
	}else if (STATISTIC=="AAF" || STATISTIC=="AFM"){
		block$avg <- (block_af1$avg + block_af2$avg) / 2;
		block$pos <- block_af1$pos;
		block$num <- block_af1$num;
	}else{
		cat("No STATISTIC definded.\n");
		q();
	}
	#
	x <- as.numeric(block$pos);
	y <- as.numeric(block$avg);
	jdx <- which(!is.na(y));
	cat("Total blocks in ", chr[i,1], ": ", length(x), "\n", sep="");
	cat("Number of valid block:", length(jdx), "\n");
	# only consider those chromosomes that have at least MINVALID valid blocks
	if(length(jdx)<MINVALID) next;
	#
	fit0  <- loess(y[jdx]~x[jdx], degree=DEG);
	span1 <- opt.span(fit0, criterion="aicc")$span;
	fit1  <- loess(y[jdx]~x[jdx], degree=DEG, span=span1);
	cat("Span/Smoothing parameter (alpha):", sprintf("%.3f", span1), ".\n");
	cat("Number of spanned blocks:", sprintf("%.0f", span1*length(jdx)), ".\n");
	#
	plo <- predict(fit1, x, se=TRUE);
	value <- plo$fit;
	###################################
	# manually correct fitted value
	if(STATISTIC=="AFD")    value <- ifelse(value < -1,    -1, value);
	if(STATISTIC=="AFD")    value <- ifelse(value >  1,     1, value);
	if(STATISTIC=="absAFD") value <- ifelse(value <  0,     0, value);
	if(STATISTIC=="absAFD") value <- ifelse(value >  1,     1, value);
	if(STATISTIC=="AAF" || STATISTIC=="AFM") value <- ifelse(value <  0,     0, value);
	if(STATISTIC=="AAF" || STATISTIC=="AFM") value <- ifelse(value >  1,     1, value);
	if(STATISTIC=="AF1")     value <- ifelse(value <  0,     0, value);
	if(STATISTIC=="AF1")     value <- ifelse(value >  1,     1, value);
	if(STATISTIC=="AF2")     value <- ifelse(value <  0,     0, value);
	if(STATISTIC=="AF2")     value <- ifelse(value >  1,     1, value);
	if(STATISTIC=="MAF")    value <- ifelse(value <  0,     0, value);
	if(STATISTIC=="MAF")    value <- ifelse(value >  0.5, 0.5, value);
	####################################
	dat_chr <- rep(chr[i,1], length(block$num));
	dat_pos <- block$pos;
	dat_avg <- block$avg;
	fit_avg <- value;
	fit_ser <- plo$s;
	dat_num <- block$num;
	blk_dph <- block$blockDepth;
	# 
	cat("Writing\n");
	if(STATISTIC=="AF1" || STATISTIC=="AF2"){
		out <- data.frame(dat_chr, dat_pos, dat_avg, fit_avg, fit_ser, dat_num, blk_dph);
	}else{
		out <- data.frame(dat_chr, dat_pos, dat_avg, fit_avg, fit_ser, dat_num);
	}
	write.table(out, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE);
	cat("\n");
}
cat("\n");

# runtime
proc.time() - ptm; 
