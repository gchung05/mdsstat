# Seed a time series
df <- data.frame(time=seq(1:24),
                 event=c(rnorm(12, 35, 3),
                         rnorm(4, 45, 5),
                         rnorm(8, 53, 8)))
replace=F
boostrap_iter=1000
alpha=0.05
log2cpmax=9 # log2(9)=512 maximum changepoints findable
cpmax=100
min_seglen=4

# Single changepoint based on mean shift
cp_mean_single <- function(
  df,
  bootstrap_iter,
  alpha,
  replace
){
  # Initial changepoint estimate based on max absolute CUSUM difference
  x_ctr <- cumsum(df$event - mean(df$event))
  cp <- which.max(abs(x_ctr)) + 1
  diff <- max(x_ctr) - min(x_ctr)
  # Bootstrap null distribution
  n <- ifelse(length(df$event) > 7, boostrap_iter, factorial(length(df$event)))
  cpwo <- mat.or.vec(n, 1)
  cpdiff <- mat.or.vec(n, 1)
  for (n1 in 1:n){
    x_random <- sample(df$event, replace=replace)
    x_random_ctr <- cumsum(x_random - mean(x_random))
    cpwo[n1] <- which.max(abs(x_random_ctr)) + 1
    cpdiff[n1] <- max(x_random_ctr) - min(x_random_ctr)
  }
  conflevel <- sum(diff > cpdiff) / n

  # # UNIMPLEMENTED (function parameter: ci_type=c("perc", "bca"))
  # # (Xu, Kass-Hout et al.) Original estimation of time interval confidence
  # if (conflevel > (1 - alpha)){
  #   cpwo <- sort(cpwo)
  #   cpquan <- c(0, 0)
  #   alpha_lo <- alpha / 2
  #   alpha_hi <- 1 - (alpha_lo)
  #   if (ci_type[1] == "perc"){ # Percentile confidence intervals
  #     cpquan <- quantile(cpwo, probs=c(alpha_lo, alpha_hi))
  #   }
  #   else if (ci_type[1] == "bca"){ # BCA confidence intervals
  #     z0 <- qnorm(sum(cpwo < bp_index[1])/n) #bp_index is an external-scope variable
  #     z0=ifelse(z0=="-Inf", -3, z0)
  #     z0=ifelse(z0=="Inf", 3, z0)
  #     ratio <- df$event
  #     cp_jack <- mat.or.vec(length(df$event), 1)
  #     for (n2 in 1:length(df$event)){
  #       ratio2 <- ratio[-n2]
  #       cp_jack[n2] <- which.max(cumsum(ratio2 - mean(ratio2))) + 1
  #     }
  #     L <- mean(cp_jack) - cp_jack
  #     alphahat <- sum(L ^ 3) / (6 * (sum(L ^ 2)) ^ 1.5)
  #     ##alphahat=ifelse(alphahat=="NaN",0,alphahat)
  #     alpha_lo <- max(0, pnorm(z0 + (z0 + qnorm(0.025)) /
  #                                (1 - alphahat * (z0 + qnorm(alpha_lo)))))
  #     alpha_hi <- min(1, pnorm(z0 + (z0 + qnorm(0.975)) /
  #                                (1 - alphahat * (z0 + qnorm(alpha_hi)))))
  #     ##cpquan<-quantile(cpwo,probs=c(alpha_lo,alpha_hi))
  #     nquan = round(alpha_lo * n)
  #     nquan = ifelse(nquan == 0, 1, nquan)
  #     cpquan[1] <- cpwo[nquan]
  #     nquan = round(alpha_hi * n)
  #     nquan = ifelse(nquan == 0, 1, nquan)
  #     cpquan[2] <- bpwo[nquan]
  #   }
  #
  # } else{
  #   cpquan <- c(1, length(df$event))
  # }
  # lb <- df$time[cpquan[1]]
  # ub <- df$time[cpquan[2]]
  # from <- df$event[cpquan[1]]
  # to <- df$event[cpquan[2]]

  out <- data.frame(
    cp=df$time[cp],
    cpindex=cp,
    cpfound=(1 - conflevel) <= alpha,
    alpha=alpha,
    pvalue=(1 - conflevel),
    cpindex=cp)
  return(out)
}

# Initial changepoint estimate
cpct <- 0
a <- cp_mean_single(
  df=df,
  bootstrap_iter=bootstrap_iter,
  alpha=alpha,
  replace=replace)
if (a$cpfound){
  cpct <- cptct + 1
  cps_now <- c(1, a$cpindex, length(df$event) + 1)
  begin <- cps_now[1]
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
        newcps <- sss
      }
      a <- rbind(a, b)

    }
    begin <- cps_now[j]
  }
} else{
  # NO Signal yada yada yada
}


for (i in c(2:(log2cpmax + 1))){
  cps_now <- c(1, a$cpindex, length(df$event) + 1)
  begin <- cps_now[1]
  # Constrain the maximum logically possible # of segments to search over
  segments_now <- min(2 ^ (i - 1), length(cps_now) - 1)

}

# Constrain the maximum number of changepoints findable to 2 ^ log2cpmax
# or n slices of slices ...
for (i in c(2:(log2cpmax + 1))){
  cps_now <- c(1, a$cpindex, length(df$event) + 1)
  begin <- cps_now[1]
  # Constrain the maximum logically possible # of segments to search over
  segments_now <- min(2 ^ (i - 1), length(cps_now) - 1)
  for (j in c(1:segments_now)){
    end <- cps_now[j + 1] - 1
    if (end > (begin + min_seglen)){
      subx <- df[begin:end, ]
      if (a$cpfound[j]){ # This search strategy logic doesn't quite work.
        b <- cp_mean_single(
          df=subx,
          bootstrap_iter=bootstrap_iter,
          alpha=alpha,
          replace=replace)
        b$cpindex[1] <- b$cpindex[1] + begin - 1
        a <- rbind(a, b)
        begin <- end + cp_mindist
      }
      else{
        begin <- cps_now[j + 1] + 1
      }
    } else{
      begin <- cps_now[j + 1] + 1
    }
    # Changepoints must be spaced at least cp_mindist apart
  }
}






##delete non-sig conflevel from dataframe a
a<-a[-which(a$conflevel<95.0),]
a<-a[cps_now(a$level),]


setwd("E:/Biosense/Distribute/R projects")
library(foreign)
flu<-read.dta("E:/Biosense/Distribute/tdistribr_20100327.dta")
flu.nyc<-flu[flu$jno == 33.1,]
age<-unique(flu.nyc$age_group)
nage<-length(age)
citype<-c("bca","perc")[2]

##Computer change points
for (i in 1:nage) {
	flu.nyc.age<-flu.nyc[flu.nyc$age_group == age[i],]
	flu.nyc.age<-flu.nyc.age[order(flu.nyc.age$dateofvisit),]
	data<-flu.nyc.age
	a<-find_break(agegp=age[i],level=1,begin=1,rep=F,type=citype,x=data)
	for (level in 2:10) {
		a<-a[order(a$bpindex),]
		order<-c(1,a$bpindex,length(data$ratio)+1)
		begin<-order[1]
		jlevel<-min(2^(level-1),length(order)-1)
		for (j in 1: jlevel){
			end<-order[j+1]-1
			if (end > (begin+2)){
				subx<-data[begin:end,]
				if (a$conflevel[j] >= 95) {
					a1<-find_break(agegp=age[i],level=level,begin=begin,rep=F,
						type=citype,x=subx)
					a1$bpindex[1]<-a1$bpindex[1]+begin-1
					a<-rbind(a,a1)
					begin<-end+2
				}
				else if (a$conflevel[j] < 95){
					j<-j+1
					begin<-order[j+1]+1
				}
			}

			##2 means first split is (1,m) and the second is (m+2,n)
		}
	}
	##delete non-sig conflevel from dataframe a
	a<-a[-which(a$conflevel<95.0),]
	a<-a[order(a$level),]
	tablefilename<-paste("flu_jno_",33.1,".txt", sep="")
	write.table(a,file=tablefilename,append=(i>1),sep="\t",col.names=(i==1),row.names=FALSE)
	pngfilename<-paste("flu_jno_",33.1,"_agegp_",age[i],".png", sep="")
	png(file=pngfilename, bg="transparent",width=1024, height=768)
	time<-data$dateofvisit
	plot(as.Date(time,"%Y-%m-%d"),data$ratio,type="l",xaxt="n",
			ylab="Ratio",col="blue")
	title(main=paste("H1N1 flu change points at NYC for age group= ", age[i],sep=""))
	ticks.at <- seq(min(time), max(time), by = "months")
	## format the labels as abbreviated month names
	ticks.lab <- format(ticks.at, format = "%b")
	## indicator variable; is month January?
	m1 <- ticks.lab == "Jan"
	## plot small ticks and labels for months not Jan
	Axis(time, at = ticks.at[!m1], side = 1,
      	labels = ticks.lab[!m1], las = 2, cex.axis = 0.7)
	## plot the default tick locations for years
	Axis(time, side = 1, las = 2)
	## add the box
	box()
	##Add breakpoint lines
	abline(v=a$bp)
	dev.off()
}

##Self written Bootstrap function (not using library boot)
find_break<-function(agegp,level,begin,rep,type,x) {
	x2<-cumsum(x$ratio-mean(x$ratio))
	bp<-which.max(abs(x2))+1 ##original change point
	diff<-max(x2)-min(x2)
	xlength<-length(x$ratio)
	n<-ifelse (xlength > 7, 1000, factorial(xlength))
	bpwo<-mat.or.vec(n,1)
	bpdiff<-mat.or.vec(n,1)
	for (n1 in 1: n){
		ratio<-sample(x$ratio,replace=rep)
		ratio2<-cumsum(ratio-mean(ratio))
		bpwo[n1]<-which.max(abs(ratio2))+1
		bpdiff[n1]<-max(ratio2)-min(ratio2)
	}
	conflevel<-sum(bpdiff <= diff)/n *100

	##Computer CI if the level is greater than 95%
	if (conflevel > 95) {
		bpwo<-sort(bpwo)
		bpquan<-c(0,0)
		if (type=="perc"){
			##percentile intervals
			alpha1<-0.025
			alpha2<-0.975
			bpquan<-quantile(bpwo,probs=c(0.025,0.975))
		}
		else if (type=="bca"){
			##BCA-CI
			z0<-qnorm(sum(bpwo < bp_index[1])/n)
			z0=ifelse(z0=="-Inf",-3,z0)
			z0=ifelse(z0=="Inf",3,z0)
			ratio<-x$ratio
			bp_jack<-mat.or.vec(xlength,1)
			for (n2 in 1:xlength){
				ratio2<-ratio[-n2]
				bp_jack[n2]<-which.max(cumsum(ratio2-mean(ratio2)))+1
			}
			L<-mean(bp_jack)-bp_jack
			alphahat<-sum(L^3)/(6*(sum(L^2))^1.5)
			##alphahat=ifelse(alphahat=="NaN",0,alphahat)
			alpha1<-max(0,pnorm(z0+(z0+qnorm(0.025))/(1-alphahat*(z0+qnorm(0.025)))))
			alpha2<-min(1,pnorm(z0+(z0+qnorm(0.975))/(1-alphahat*(z0+qnorm(0.975)))))
			##bpquan<-quantile(bpwo,probs=c(alpha1,alpha2))
			nquan=round(alpha1*n)
			nquan=ifelse(nquan==0,1,nquan)
			bpquan[1]<-bpwo[nquan]
			nquan=round(alpha2*n)
			nquan=ifelse(nquan==0,1,nquan)
			bpquan[2]<-bpwo[nquan]
		}
	}
	else {
		bpquan<-c(1,xlength)
	}
	lb<-x$dateofvisit[bpquan[1]]
	ub<-x$dateofvisit[bpquan[2]]
	from<-x$ratio[bpquan[1]]
	to<-x$ratio[bpquan[2]]
	a<-data.frame(age=agegp,level=level, bp=x$dateofvisit[bp], conflevel=conflevel,
		lb=lb,ub=ub,from=from,to=to,bpindex=bp)
	return(a)
}



