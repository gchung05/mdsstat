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
#
#
#
