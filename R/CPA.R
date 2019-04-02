# Seed a time series
# ------------------
df <- data.frame(time=seq(1:49),
                 event=c(rnorm(12, 35, 3), #13 ideal changepoint
                         rnorm(4, 45, 5), #17
                         rnorm(8, 53, 8), #25
                         rnorm(15, 30, 3), #40
                         rnorm(5, 40, 5), #45
                         rnorm(5, 80, 2))) #46

# Parameter list
# --------------
alpha=0.05 # Alpha level for changepoint determination
cpmax=100 # max number of changepoints allowed. This supersedes the theoretical max based on epochs.
min_seglen=6 # Minimum required length of consecutive measurements where a changepoint search is allowed
epochs=NULL # Number of epochs allowed where 2^epochs is the max theoretical changepoints findable.
# When epochs=NULL, epochs is estimated from number of observations and min_seglen
bootstrap_iter=1000 # Lowest recommended is 1000. Increasing also gives more p-value precision
replace=T # For bootstrapping, sampling with replacement or not

# Mean-shift changepoint algorithm
# --------------------------------
# Single changepoint based on mean shift
# Original code by Xu, Kass-Hout 2015
# Updated by Chung 2019
cp_mean_single <- function(
  df,
  bootstrap_iter,
  alpha,
  replace
){
  # Changepoint estimate based on max absolute CUSUM difference
  x_ctr <- cumsum(df$event - mean(df$event))
  cp <- which.max(abs(x_ctr)) + 1
  diff <- max(x_ctr) - min(x_ctr)
  # Bootstrap null distribution
  n <- ifelse(length(df$event) > 7, bootstrap_iter, factorial(length(df$event)))
  cpwo <- mat.or.vec(n, 1)
  cpdiff <- mat.or.vec(n, 1)
  for (n1 in 1:n){
    x_random <- sample(df$event, replace=replace)
    x_random_ctr <- cumsum(x_random - mean(x_random))
    cpwo[n1] <- which.max(abs(x_random_ctr)) + 1
    cpdiff[n1] <- max(x_random_ctr) - min(x_random_ctr)
  }
  conflevel <- sum(diff > cpdiff) / n
  # Assemble output
  out <- data.frame(
    cp=df$time[cp],
    cpindex=cp,
    cpfound=(1 - conflevel) <= alpha,

    alpha=alpha,
    pvalue=(1 - conflevel))
  return(out)
}

# Code start
# ----------
# Inspired by code written by Xu, Kass-Hout 2015
# Rewritten search start/stop strategy

# Initial changepoint estimate
a <- cp_mean_single(
  df=df,
  bootstrap_iter=bootstrap_iter,
  alpha=alpha,
  replace=replace)
if (a$cpfound){
  cpct <- epochct <- 1
  # Max epochs based on number of measurements
  if (is.null(epochs)) epochs <- floor(log2(nrow(df) / (min_seglen / 2)))
  # Constrain maximum number of findable changepoints
  while (cpct < cpmax & epochct <= epochs){
    cpct_now <- cpct
    a <- a[order(a$cpindex), ]
    # Identify current segments
    cps_now <- c(1, a$cpindex, length(df$event) + 1)
    begin <- cps_now[1]
    bhold <- data.frame()
    # Iterate over current segments
    for (j in c(2:(length(cps_now)))){
      end <- cps_now[j] - 1
      if ((end - begin + 1) >= min_seglen){
        subx <- df[begin:end, ]
        b <- cp_mean_single(
          df=subx,
          bootstrap_iter=bootstrap_iter,
          alpha=alpha,
          replace=replace)
        b$cpindex[1] <- b$cpindex[1] + begin - 1
        if (b$cpfound){
          cpct_now <- cpct_now + 1
          bhold <- rbind(bhold, b)
        }
      }
      begin <- cps_now[j]
    }
    # If no changepoints are found in current epoch, stop search
    if (cpct_now == cpct){
      break
    } else{
      # If this epoch identifies more changepoints than max allowable
      # Keep up to max by ranked p-value
      if (cpct_now > cpmax){
        bhold <- bhold[order(bhold$pvalue), ][c(1:(cpmax - nrow(a))), ]
      }
      a <- rbind(a, bhold)
      cpct <- nrow(a)
      epochct <- epochct + 1
    }
  }
} else{
  cat("\nNo signal found")
}


# # ORIGINAL CODE
# # -------------
#
# setwd("E:/Biosense/Distribute/R projects")
# library(foreign)
# flu<-read.dta("E:/Biosense/Distribute/tdistribr_20100327.dta")
# flu.nyc<-flu[flu$jno == 33.1,]
# age<-unique(flu.nyc$age_group)
# nage<-length(age)
# citype<-c("bca","perc")[2]
#
# ##Computer change points
# for (i in 1:nage) {
# 	flu.nyc.age<-flu.nyc[flu.nyc$age_group == age[i],]
# 	flu.nyc.age<-flu.nyc.age[order(flu.nyc.age$dateofvisit),]
# 	data<-flu.nyc.age
# 	a<-find_break(agegp=age[i],level=1,begin=1,rep=F,type=citype,x=data)
# 	for (level in 2:10) {
# 		a<-a[order(a$bpindex),]
# 		order<-c(1,a$bpindex,length(data$ratio)+1)
# 		begin<-order[1]
# 		jlevel<-min(2^(level-1),length(order)-1)
# 		for (j in 1: jlevel){
# 			end<-order[j+1]-1
# 			if (end > (begin+2)){
# 				subx<-data[begin:end,]
# 				if (a$conflevel[j] >= 95) {
# 					a1<-find_break(agegp=age[i],level=level,begin=begin,rep=F,
# 						type=citype,x=subx)
# 					a1$bpindex[1]<-a1$bpindex[1]+begin-1
# 					a<-rbind(a,a1)
# 					begin<-end+2
# 				}
# 				else if (a$conflevel[j] < 95){
# 					j<-j+1
# 					begin<-order[j+1]+1
# 				}
# 			}
#
# 			##2 means first split is (1,m) and the second is (m+2,n)
# 		}
# 	}
# 	##delete non-sig conflevel from dataframe a
# 	a<-a[-which(a$conflevel<95.0),]
# 	a<-a[order(a$level),]
# 	tablefilename<-paste("flu_jno_",33.1,".txt", sep="")
# 	write.table(a,file=tablefilename,append=(i>1),sep="\t",col.names=(i==1),row.names=FALSE)
# 	pngfilename<-paste("flu_jno_",33.1,"_agegp_",age[i],".png", sep="")
# 	png(file=pngfilename, bg="transparent",width=1024, height=768)
# 	time<-data$dateofvisit
# 	plot(as.Date(time,"%Y-%m-%d"),data$ratio,type="l",xaxt="n",
# 			ylab="Ratio",col="blue")
# 	title(main=paste("H1N1 flu change points at NYC for age group= ", age[i],sep=""))
# 	ticks.at <- seq(min(time), max(time), by = "months")
# 	## format the labels as abbreviated month names
# 	ticks.lab <- format(ticks.at, format = "%b")
# 	## indicator variable; is month January?
# 	m1 <- ticks.lab == "Jan"
# 	## plot small ticks and labels for months not Jan
# 	Axis(time, at = ticks.at[!m1], side = 1,
#       	labels = ticks.lab[!m1], las = 2, cex.axis = 0.7)
# 	## plot the default tick locations for years
# 	Axis(time, side = 1, las = 2)
# 	## add the box
# 	box()
# 	##Add breakpoint lines
# 	abline(v=a$bp)
# 	dev.off()
# }
#
# ##Self written Bootstrap function (not using library boot)
# find_break<-function(agegp,level,begin,rep,type,x) {
# 	x2<-cumsum(x$ratio-mean(x$ratio))
# 	bp<-which.max(abs(x2))+1 ##original change point
# 	diff<-max(x2)-min(x2)
# 	xlength<-length(x$ratio)
# 	n<-ifelse (xlength > 7, 1000, factorial(xlength))
# 	bpwo<-mat.or.vec(n,1)
# 	bpdiff<-mat.or.vec(n,1)
# 	for (n1 in 1: n){
# 		ratio<-sample(x$ratio,replace=rep)
# 		ratio2<-cumsum(ratio-mean(ratio))
# 		bpwo[n1]<-which.max(abs(ratio2))+1
# 		bpdiff[n1]<-max(ratio2)-min(ratio2)
# 	}
# 	conflevel<-sum(bpdiff <= diff)/n *100
#
# 	##Computer CI if the level is greater than 95%
# 	if (conflevel > 95) {
# 		bpwo<-sort(bpwo)
# 		bpquan<-c(0,0)
# 		if (type=="perc"){
# 			##percentile intervals
# 			alpha1<-0.025
# 			alpha2<-0.975
# 			bpquan<-quantile(bpwo,probs=c(0.025,0.975))
# 		}
# 		else if (type=="bca"){
# 			##BCA-CI
# 			z0<-qnorm(sum(bpwo < bp_index[1])/n)
# 			z0=ifelse(z0=="-Inf",-3,z0)
# 			z0=ifelse(z0=="Inf",3,z0)
# 			ratio<-x$ratio
# 			bp_jack<-mat.or.vec(xlength,1)
# 			for (n2 in 1:xlength){
# 				ratio2<-ratio[-n2]
# 				bp_jack[n2]<-which.max(cumsum(ratio2-mean(ratio2)))+1
# 			}
# 			L<-mean(bp_jack)-bp_jack
# 			alphahat<-sum(L^3)/(6*(sum(L^2))^1.5)
# 			##alphahat=ifelse(alphahat=="NaN",0,alphahat)
# 			alpha1<-max(0,pnorm(z0+(z0+qnorm(0.025))/(1-alphahat*(z0+qnorm(0.025)))))
# 			alpha2<-min(1,pnorm(z0+(z0+qnorm(0.975))/(1-alphahat*(z0+qnorm(0.975)))))
# 			##bpquan<-quantile(bpwo,probs=c(alpha1,alpha2))
# 			nquan=round(alpha1*n)
# 			nquan=ifelse(nquan==0,1,nquan)
# 			bpquan[1]<-bpwo[nquan]
# 			nquan=round(alpha2*n)
# 			nquan=ifelse(nquan==0,1,nquan)
# 			bpquan[2]<-bpwo[nquan]
# 		}
# 	}
# 	else {
# 		bpquan<-c(1,xlength)
# 	}
# 	lb<-x$dateofvisit[bpquan[1]]
# 	ub<-x$dateofvisit[bpquan[2]]
# 	from<-x$ratio[bpquan[1]]
# 	to<-x$ratio[bpquan[2]]
# 	a<-data.frame(age=agegp,level=level, bp=x$dateofvisit[bp], conflevel=conflevel,
# 		lb=lb,ub=ub,from=from,to=to,bpindex=bp)
# 	return(a)
# }
