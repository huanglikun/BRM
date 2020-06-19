#!/usr/bin/env RScript

# version 0.2
##############################################
# command
##############################################

opt <- commandArgs(TRUE);
# help
if (length(opt) < 3) {cat("Usage: Rscript BRM.R BRM_conf.txt chr_length.txt markers.bsa\n");q();}
#
file_conf   <- opt[1]; # arguments configure file
file_chrlen <- opt[2]; # chromosome length file
file_bsa    <- opt[3]; # bsa file

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

#
create_dir <- function(dir){
  if (file.exists(dir)){
    # cat(dir," already exists.\n");
  }else{
    cat("Creating ",dir," ...\n");
    tryCatch({dir.create(dir,recursive=T)},warning=function(w){s <- as.character(w);stop(s)},error=function(e){s <- as.character(e);stop(s)});
  }
}

#
chk_num_para <- function(x,name){
	if (is.null(x)){
		stop("Undefined ",name,".");
	}
	return(as.numeric(x));
}

#
chk_dir <- function(path){
	dir <- dirname(path);
	if (dir != '.') create_dir(dir);
}

#
chk_conf <- function(conf){
	if (is.null(conf$Design) || (conf$Design != "A" && conf$Design != "BH" && conf$Design != "BL")) stop("Incorrect Design.");
	if (is.null(conf$t) || (conf$t != 0 && conf$t != 1)) stop("Incorrect t.");
	conf$n1       <- chk_num_para(conf$n1,"n1");
	conf$n2       <- chk_num_para(conf$n2,"n2");
	conf$ua       <- chk_num_para(conf$ua,"ua");
	conf$UNIT     <- chk_num_para(conf$UNIT,"UNIT");
	conf$DEG      <- chk_num_para(conf$DEG,"DEG");
	conf$BLK      <- chk_num_para(conf$BLK,"BLK");
	conf$MIN      <- chk_num_para(conf$MIN,"MIN");
	conf$MINVALID <- chk_num_para(conf$MINVALID,"MINVALID");
	if (is.null(conf$Result1_File)) conf$Result1_File <- "result/result1.xls";
	if (is.null(conf$Result2_File)) conf$Result2_File <- "result/result2.xls";
	chk_dir(conf$Result1_File);
	chk_dir(conf$Result2_File);
	return(conf);
}

#####################
# step 1
#####################

# block's middle position 
def_block_pos <- function(chr, size){
	# block number
	n <- as.double(2*chr/size)+2;
	if( n%%2 != 0 ) n <- n+1;
	n <- as.double(n/2);
	# block index and the middle position of each block
	i <- c(1:n);
	pos <- (i-1)*size+floor(size/2); # middle position
	if(pos[n]>chr) pos[n]<-chr;
	return(pos);
}

# meta value of block
cal_block_meta <- function(loc, val, chr, size, depth, MIN){
	# input: location vector, value vector
	# input: chr length, block size, location depth vector
	pos <- def_block_pos(chr, size);
	idx <- as.double(0.5+loc/size)+1;
	#
	avg <- c();
	blockDepth <- c();
	for(i in 1:length(pos)){
		k  <- which(idx==i);
		no <- length(k);
		a <- NA;
		n <- 0;
		if (no > 0) {n <- sum(depth[k])};
		if( n >= MIN ) a <- sum(val[k])/n;
		avg <- c(avg, a);
	}
	return( list(pos=pos, avg=avg) );
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

#
run_step1 <- function(conf,file_chrlen,file_bsa,file_out,statistic){
	#
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
	# return
	fit_model <- list();
	# global options
	BLK <- as.numeric(conf$BLK) * as.numeric(conf$UNIT); # block size/step
	MIN <- as.numeric(conf$MIN); # least total depth in one block 
	DEG <- as.numeric(conf$DEG); # degree of polynomial
	MINVALID <- as.numeric(conf$MINVALID); # least valid blocks in one chromosome
	#
	if(MIN==0) MIN <- 1;
	if(MINVALID < 10) MINVALID <- 10;
	# read into memory
	tab <- read.table(file_bsa, as.is=TRUE);
	chro <- read.table(file_chrlen, as.is=TRUE);
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
	cat("Number of imported chromosomes is", length(chro[,1]), ".\n");
	for(i in 1:length(chro[,1])){
		dat_chr <- dat_all_chr;
		dat_pos <- dat_all_pos;
		idx <- which(dat_chr == chro[i,1]);
		total <- length(idx);
		cat("Data size of", chro[i,1], "is", total, ".\n");
		if(total==0) next;
		# initialize
		block <- list();
		#
		if (statistic=="AF1") {
			block <- cal_block_meta(dat_pos[idx], A[idx], chro[i, 2], BLK, pool1Depth[idx], MIN);
		}else if (statistic=="AF2") {
			block <- cal_block_meta(dat_pos[idx], C[idx], chro[i, 2], BLK, pool2Depth[idx], MIN);
		}else{
			block_af1 <- cal_block_meta(dat_pos[idx], A[idx], chro[i, 2], BLK, pool1Depth[idx], MIN);
			block_af2 <- cal_block_meta(dat_pos[idx], C[idx], chro[i, 2], BLK, pool2Depth[idx], MIN);
		}
		#
		if (statistic=="AFD") {
		    if (conf$Design == "BL"){
		      block$avg <- block_af2$avg - block_af1$avg;
		    }else{
		      block$avg <- block_af1$avg - block_af2$avg;
		    }
			block$pos <- block_af1$pos;
			block$num <- block_af1$num;
		}else if (statistic=="AAF"){
			block$avg <- (block_af1$avg + block_af2$avg) / 2;
			block$pos <- block_af1$pos;
			block$num <- block_af1$num;
		}else if (statistic!="AF1" & statistic!="AF2"){
			stop("No statistic definded.");
		}
		#
		x <- as.numeric(block$pos);
		y <- as.numeric(block$avg);
		jdx <- which(!is.na(y));
		cat("Total blocks in ", chro[i,1], ": ", length(x), "\n", sep="");
		cat("Number of valid block:", length(jdx), "\n");
		# only consider those chromosomes that have at least MINVALID valid blocks
		if(length(jdx)<MINVALID) next;
		#
		fits <- list();
		fit0  <- loess(y[jdx]~x[jdx], degree=DEG);
		span1 <- opt.span(fit0, criterion="aicc")$span;
		fit1  <- loess(y[jdx]~x[jdx], degree=DEG, span=span1);
		fit_model[[ i ]] <- fit1;
		#
		plo <- predict(fit1, x, se=TRUE);
		value <- plo$fit;
		###################################
		# manually correct fitted value
		if(statistic=="AFD")    value <- ifelse(value < -1,    -1, value);
		if(statistic=="AFD")    value <- ifelse(value >  1,     1, value);
		if(statistic=="AAF") value <- ifelse(value <  0,     0, value);
		if(statistic=="AAF") value <- ifelse(value >  1,     1, value);
		if(statistic=="AF1")     value <- ifelse(value <  0,     0, value);
		if(statistic=="AF1")     value <- ifelse(value >  1,     1, value);
		if(statistic=="AF2")     value <- ifelse(value <  0,     0, value);
		if(statistic=="AF2")     value <- ifelse(value >  1,     1, value);
		####################################
		dat_chr <- rep(chro[i,1], length(block$pos));
		dat_pos <- block$pos;
		dat_avg <- round(block$avg,4);
		fit_avg <- round(value,4);
		out <- data.frame(dat_chr, dat_pos, dat_avg, fit_avg);
		if (i == 1){
			allout <- out;
		}else{
			allout <- rbind(allout,out);
		}
	}
	cat("\n");
	return(list(fits=fit_model,data=allout));
}

#
run_step2 <- function(conf,data_aaf,data_af1,data_af2,file_out){
	#
	N1 <- as.numeric(conf$n1); # number of pool 1
	N2 <- as.numeric(conf$n2); # number of pool 2
	T  <- as.numeric(conf$t);  # level of population. For DH or RI etc., T=0; F2 or F3 etc., T=1
	ua <- as.numeric(conf$ua); # uα
	#
	dat_chr <- data_aaf[,1];
	dat_pos <- data_aaf[,2];
	dat_val <- data_aaf[,3];
	fit_val <- data_aaf[,4];
	fit_af1 <- data_af1[,4];
	fit_af2 <- data_af2[,4];
	# sample threshold
	dat_var <- (N1 + N2) / (2^T * N1 * N2) * fit_val * (1 - fit_val);
	dat_stdua <- ua * sqrt(dat_var);
	# theoretical threshold
	var_fix <- (N1 + N2) / (2^T * N1 * N2) * 0.25; # Fix AF = 0.5
	stdua_fix <- ua * sqrt(var_fix);
	# variance if this location is QTL
	var_qtl <- (fit_af1 * (1 - fit_af1) / (2^T * N1)) + (fit_af2 * (1 - fit_af2) / (2^T * N2));
	# 
	out <- data.frame(dat_chr, dat_pos, dat_val, fit_val, dat_stdua, dat_var, stdua_fix, var_qtl);
	return(list(theo=stdua_fix,data=out));
}

#
run_step3 <- function(conf,dat_afd,dat_var,file_out,fit_model_afd,fit_model_af1,fit_model_af2,threshold){
	#
	threshold <- as.numeric(threshold);
  	#
  	N1 <- as.numeric(conf$n1); # number of pool 1
  	N2 <- as.numeric(conf$n2); # number of pool 2
  	T  <- as.numeric(conf$t);  # level of population. For DH or RI etc., T=0; F2 or F3 etc., T=1
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
	# write title to file
	chrom   <- c("#Chr.");
	peak_x  <- c("Pos.");
	peak_y  <- c("Val.");
	type    <- c("Peak Dir.");
	left_x  <- c("Start");
	right_x <- c("End");
	out <- data.frame(chrom, peak_x, peak_y, type, left_x, right_x);
	write.table(out, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=FALSE);
	#
	dat_chr <- dat_afd[,1];
	dat_pos <- dat_afd[,2];
	dat_avg <- dat_afd[,3];
	fit_avg <- dat_afd[,4];
	#
	dat_std <- sqrt(dat_var[,6]);
	qtl_std <- sqrt(dat_var[,8]);
	#
	chro     <- data.frame(unique(dat_chr));
	# one-by-one for every chromosome
	#
	for(i in 1:length(chro[,1])){
	    thischr <- chro[i,1];
	    idx <- which(dat_chr == chro[i,1] & !is.na(dat_avg));
	    total <- length(idx);
	    cat("data size of", chro[i,1], "is", total, "\n");
	    if(total==0) next;
	    #
	    x <- dat_pos[idx];
	    y <- dat_avg[idx];
		fity <- fit_avg[idx];
	    std  <- dat_std[idx];
	    qstd <- qtl_std[idx];
	    fit2 <- locfit_by_loess(x, y);
		#
		peakxs_allidx <- which(diff(diff(fity)>0)!=0L)+1L;
		peakxs_flt0<- which(abs(fity[peakxs_allidx])>threshold);
		if(length(peakxs_flt0)==0) next;
		peakxs_idx <- peakxs_allidx[peakxs_flt0];
		peakallxs<- x[peakxs_allidx];
	    peakxs   <- x[peakxs_idx];
	    len  <- length(peakallxs);
		approx_peaks <- peakxs;
		appr_peaks_y <- fity[peakxs_idx];
		x_pool <- approx_peaks;
		y_pool <- appr_peaks_y;
		#
		ylength <- length(fity);
		xlength <- length(x);
		# left
		lefts <- (c(x[1],peakallxs[1:len-1]) + peakallxs) / 2;
		# right
		rights <- (peakallxs + c(peakallxs[2:len],x[xlength])) / 2;
		#
		idx_pass <- c();
		# 
	    if (len==0){
	      k         <- ifelse(diff(c(fity[1],fity[ylength]))>0L,1L,-1L);
	      k         <- ifelse(k>0L,c(1L,-1L),c(-1L,1L));
	      pks       <- data.frame(peak_x=c(x[1],x[xlength]),peak_y=c(fity[1],fity[ylength]),type=k);
	      peak_std  <- c(qstd[1],qstd[xlength]);
	    }else{
	      k         <- diff(fity)>0L;
		  x1k       <- ifelse(k[1]==TRUE,1L,-1L);
		  xlastk    <- ifelse(k[length(k)],-1L,1L);
		  k         <- diff(k);
		  #
	      k         <- k[which(k!=0L)];
		  k         <- k[peakxs_flt0];
		  k_flt_idx <- which(appr_peaks_y*k<0);
		  if(length(k_flt_idx)==0) next;
	      #
		  peakxs_idx    <- peakxs_idx[k_flt_idx];
		  appr_peaks_y  <- appr_peaks_y[k_flt_idx];
		  y_pool        <- appr_peaks_y;
		  x_pool        <- x_pool[k_flt_idx];
		  k             <- k[k_flt_idx];
		  m             <- k<0;
		  lefts         <- lefts[peakxs_flt0[k_flt_idx]];
		  rights        <- rights[peakxs_flt0[k_flt_idx]];
	      #
	      pks  <- sapply(as.data.frame(rbind(lefts, rights,m)),
	                     function(r){ 
	                        optimize(f=function(x){
	                        	predict(fit2, newdata=data.frame(x))},  # if let fit_afd = fit(obs_af1-obs_af2)
	                        	# predict(fit_model_af1[[ i ]],x) - predict(fit_model_af2[[ i ]],x)},  # if let fit_afd = fit_af1 - fit_af2
	                     	maximum=r[3], interval=r[c(1,2)])});
	      pks  <- as.data.frame(t(pks));
	      names(pks) <- c("peak_x","peak_y");
	      pks$peak_x <- trunc(as.numeric(pks$peak_x));
	      pks$peak_y <- round(as.numeric(pks$peak_y),4);
		  pks$peak_y <- ifelse(pks$peak_y > 1, 1, pks$peak_y);
		  pks$peak_y <- ifelse(pks$peak_y < -1, -1, pks$peak_y);
	      pks$type   <- k;
	      peak_af1   <- predict(fit_model_af1[[ i ]],pks$peak_x);
	      peak_af2   <- predict(fit_model_af2[[ i ]],pks$peak_x);
		  #
		  peak_af1   <- ifelse(peak_af1 < 0, 0, peak_af1);
		  peak_af2   <- ifelse(peak_af2 < 0, 0, peak_af2);
		  #
	      peak_std   <- as.numeric(sqrt(peak_af1 * (1 - peak_af1) / (2^T * N1) + peak_af2 * (1 - peak_af2) / (2^T * N2)));
		  # check if the first pos is a peak
		  if (abs(fity[1])>abs(threshold) & fity[1] * x1k < 0){
			  p1       <- c(x[1],fity[1],x1k);
			  pks      <- rbind(p1,pks);
			  peak_std <- c(qstd[1],peak_std);
			  y_pool   <- c(fity[1],y_pool);
			  x_pool   <- c(x[1],x_pool);
			  peakxs_idx <- c(1,peakxs_idx);
		  }
		  # check if the last pos is a peak
		  if (abs(fity[ylength])>abs(threshold) & fity[ylength] * xlastk < 0){
			  plast    <- c(x[ylength],fity[ylength],xlastk);
			  pks      <- rbind(pks,plast);
			  peak_std <- c(peak_std,qstd[ylength]);
			  y_pool   <- c(y_pool,fity[ylength]);
			  x_pool   <- c(x_pool,x[xlength]);
			  peakxs_idx <- c(peakxs_idx,ylength);
		  }
		  #
		  all_ident_y    <- pks$peak_y + pks$type*1.65*peak_std;
		}
		#
		peak_flt1   <- c();
		pos1_list   <- c();
		pos2_list   <- c();
		fltlen      <- length(pks$peak_x);
		while(length(idx_pass) < fltlen){
			# find max
			max_idx    <- which(abs(y_pool)==max(abs(y_pool)));
			peak_flt1  <- c(peak_flt1, max_idx);
			interval_i <- peakxs_idx[max_idx];
			# pos1
			if(interval_i == 1){
				pos1_x     <- x[1];
			}else{
				while(interval_i > 2 & (fity[interval_i]-all_ident_y[max_idx])*(fity[interval_i-1]-all_ident_y[max_idx]) > 0){
					interval_i <- interval_i - 1;
					if((fity[interval_i]-all_ident_y[max_idx])*(fity[interval_i-1]-all_ident_y[max_idx]) < 0) break;
				}
				distance_x1x2  <- x[interval_i] - x[interval_i-1];
				dis1_y1        <- abs(fity[interval_i-1]-all_ident_y[max_idx]);
				dis1_y2        <- abs(fity[interval_i]-all_ident_y[max_idx]);
				pos1_x         <- x[interval_i-1]+(distance_x1x2)*(dis1_y1/(dis1_y1+dis1_y2));
			}
			pos1_list      <- c(pos1_list,pos1_x);
			interval_i     <- peakxs_idx[max_idx];
			# pos2
			if(interval_i == xlength){
				pos2_x     <- x[xlength];
			}else{
				while(interval_i < xlength - 2 & (fity[interval_i]-all_ident_y[max_idx])*(fity[interval_i+1]-all_ident_y[max_idx]) > 0){
					interval_i <- interval_i + 1;
					if((fity[interval_i]-all_ident_y[max_idx])*(fity[interval_i+1]-all_ident_y[max_idx]) < 0) break;
				}
				distance_x1x2  <- x[interval_i+1] - x[interval_i];
				dis2_y1        <- abs(fity[interval_i]-all_ident_y[max_idx]);
				dis2_y2        <- abs(fity[interval_i+1]-all_ident_y[max_idx]);
				pos2_x         <- x[interval_i]+(distance_x1x2)*(dis2_y1/(dis2_y1+dis2_y2));
			}
			pos2_list      <- c(pos2_list,pos2_x);
			# loop until pass all candidate peaks
			flt_region_idx <- which(x_pool >= pos1_x & x_pool <= pos2_x);
			idx_pass       <- unique(c(idx_pass,flt_region_idx));
			y_pool[idx_pass] <- 0;
		}
		pks <- pks[sort(peak_flt1),];
		pks <- data.frame(thischr, pks,left_x=c(trunc(pos1_list)), right_x=c(trunc(pos2_list)));
		peak_flt2 <- which((pks$peak_y >= threshold | pks$peak_y <= -threshold) & pks$peak_y * pks$type < 0);
		pks <- pks[peak_flt2,];
		pks$type <- ifelse(pks$type==1,"-","+");
		write.table(pks, file_out, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE);
	}
}

#########################################
# main
#########################################
ptm <- proc.time(); # begin to record runtime 
#
# read configure file
cat("Read global parameters from", file_conf, ".\n");
conf <- read_conf(file_conf);
# check parameter
conf <- chk_conf(conf);
# print parameter
for(name in names(conf) ){
  cat("Parameter ", name, "=", conf[[name]],  "\n", sep="");
}

cat("\n");
# step 1
af1s <- run_step1(conf,file_chrlen,file_bsa,conf$AF1_File,"AF1");
af2s <- run_step1(conf,file_chrlen,file_bsa,conf$AF2_File,"AF2");
afds <- run_step1(conf,file_chrlen,file_bsa,conf$AFD_File,"AFD");
if (conf$Design == "A"){
  aafs <- run_step1(conf,file_chrlen,file_bsa,conf$AAF_File,"AAF");
}else if (conf$Design == "BH" || conf$Design == "BL"){
  aafs <- af2s;
  file.copy(conf$AF2_File,conf$AAF_File);
}else{
  stop("Is design not A or B?");
}

af1_model <- af1s$fits;
af1_data  <- af1s$data;
af2_model <- af2s$fits;
af2_data  <- af2s$data;
aaf_model <- aafs$fits;
aaf_data  <- aafs$data;
afd_model <- afds$fits;
afd_data  <- afds$data;

# step 2
step2s <- run_step2(conf,aaf_data,af1_data,af2_data,conf$STEP2OUT);
theoretical_threshold <- round(step2s$theo,4);
var_data              <- step2s$data;
cat("Theoretical threshold is ±",theoretical_threshold,".\n");

# threshold and title
cat("Export data to", conf$Result1_File, "\n");
thr_title <- "#Theoretical threshold is ±";
out <- data.frame(thr_title, theoretical_threshold);
write.table(out, conf$Result1_File, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=FALSE);

dat_chr  <- c("#Chr.");
dat_pos  <- c("Pos.");
dat_avg1 <- c("AF1-Observed");
fit_avg1 <- c("AF1-Expected");
dat_avg2 <- c("AF2-Observed");
fit_avg2 <- c("AF2-Expected");
dat_avg3 <- c("AFD-Observed");
fit_avg3 <- c("AFD-Expected");
dat_avg4 <- c("AFP-Observed");
fit_avg4 <- c("AFP-Expected");
dat_varua <- c("Sample threshold");
out <- data.frame(dat_chr, dat_pos, dat_avg1, fit_avg1, dat_avg2, fit_avg2, dat_avg3, fit_avg3, dat_avg4, fit_avg4, dat_varua);
write.table(out, conf$Result1_File, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE);

cat("Writing\n");
out <- data.frame(af1_data[,1], af1_data[,2], af1_data[,3], af1_data[,4], af2_data[,3], af2_data[,4], afd_data[,3], afd_data[,4], aaf_data[,3], aaf_data[,4], round(var_data[,5],4));
write.table(out, conf$Result1_File, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE);
cat("\n");

# step 3
run_step3(conf,afd_data,var_data,conf$Result2_File,afd_model,af1_model,af2_model,theoretical_threshold);

# runtime
proc.time() - ptm; 
#
cat("finished.\n")