#!/usr/bin/env RScript

# version 0.1
##############################################
# command
##############################################
opt <- commandArgs(TRUE);
# help
if (length(opt) < 5) {cat("Usage: Rscript BRMstep2.threshold.R threshold_conf.txt AAF.xls AF1.xls AF2.xls threshold.xls\n");q();}
#
file_conf <- opt[1]; # arguments configure file
file_tab  <- opt[2]; # input AAF or total pool AF
file_af1  <- opt[3]; # input AF1
file_af2  <- opt[4]; # input AF2
file_out  <- opt[5]; # output

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
N1 <- as.numeric(conf$n1); # number of pool 1
N2 <- as.numeric(conf$n2); # number of pool 2
T  <- as.numeric(conf$t);  # level of population. For DH or RI etc., T=0; F2 or F3 etc., T=1
ua <- as.numeric(conf$ua); # μα

# read into memory
tab <- read.table(file_tab, as.is=TRUE);
af1 <- read.table(file_af1, as.is=TRUE);
af2 <- read.table(file_af2, as.is=TRUE);

# one-by-one for every chromosome
dat_chr   <- c("#Chr.");
dat_pos   <- c("Position");
dat_val   <- c("Allele Frequency");
fit_val   <- c("Fitted Average");
dat_varua <- c("Sample threshold");
dat_var   <- c("Variance of sample");
thr_fix   <- c("Theoretical threshold");
thr_qtl   <- c("Variance if is QTL");

# title write to file
cat("Export data to", file_out, "\n");
out <- data.frame(dat_chr, dat_pos, dat_val, fit_val, dat_varua, dat_var, thr_fix, thr_qtl);
write.table(out, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=FALSE);
#
dat_chr <- tab[,1];
dat_pos <- tab[,2];
dat_val <- tab[,3];
fit_val <- tab[,4];
fit_af1 <- af1[,4];
fit_af2 <- af2[,4];
# sample threshold
dat_var <- (N1 + N2) / (2^T * N1 * N2) * fit_val * (1 - fit_val);
dat_stdua <- ua * sqrt(dat_var);
# theoretical threshold
var_fix <- (N1 + N2) / (2^T * N1 * N2) * 0.25; # Fix AF = 0.5
stdua_fix <- ua * sqrt(var_fix);
# variance if this location is QTL
var_qtl <- (fit_af1 * (1 - fit_af1) / (2^T * N1)) + (fit_af2 * (1 - fit_af2) / (2^T * N2));

# 
cat("Writing\n");
out <- data.frame(dat_chr, dat_pos, dat_val, fit_val, dat_stdua, dat_var, stdua_fix, var_qtl);
write.table(out, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE);
cat("\n");