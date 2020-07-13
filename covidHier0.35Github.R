library(abind)
DATA=readRDS("covidHierData.rds")



#Some parameter name differences from paper:
#xi=dlt, c=ftw, initial testing tauI_t0=cvTl, final testing tauI_tf=tauI, tauA=tauEAf*tauI
#Mfact scales how much travel happens compared to reality; Mfact=1 is just movement from commuting data
#epsP allows voluntary distancing efficacy to differ from that of closures (epsP=eps in paper)
#Tg and Tl are the closure thresholds, n0 is fraction initially infected
#Msave is the travel matrix. Here Msave entries are numbers of commuters; msim3 converts them to proportions
pops=colSums(DATA$Msave)
parms=cbind(N=pops,Mfact=1,s=0.2,Tg=1,Tl=exp(-9.14),beta=6.41e-06,cvTl=0.045,tauI=0.3,tauEAf=0,dlt=0.27,epsP=0.62,eps=0.62,w=0.45,omg=44356,r=0.19,eta=0.8,ftw=4.785e-06,alpha=0.4,sig=0.4,rho=0.67,n0=0.0001,DATA$Msave[-50,]);

inCstart=325 # number of positive cases when province closed on Mar 17
inClen=75 #Duration of initial province closure
Resurge=50 #Duration of additional closures 
NP=ncol(parms)-nrow(parms)


#Defining model state space. Tn, Tk, Da, and Di are all untested, tested, asymptomatics, and infecteds respectively.
#Nt tracks cumulative # positive cases (including those recovered)
#C tracks # days since last closure, but in output msim converts all positive C entries into eps
Tn=c("A","SA","I","SI"); Tk=c("Ak","SAk","Ik","SIk"); Da=c("A","Ak","SA","SAk"); Di=c("I","Ik","SI","SIk"); Snm=c("S","E",Da,Di,"R","Nt","C");


#these functions handle all the state transitions. Input x is a matrix where 1st half of columns are source states and 2nd half are destination states
gpTransB=function(x,Prs,seed,nc=ncol(x)){ 
	xvp=cbind(as.vector(x[,1:(nc/2)]),as.vector(Prs))
	if(max(xvp[,1])==0) return(x); nz=(xvp[,1]*xvp[,2])>0; xvp[!nz,1]=0; 
	set.seed(seed); xvp[nz,1]=apply(matrix(xvp[nz,],ncol=2),1,function(y) rbinom(1,y[1],y[2])); 
	return(x+matrix(c(-xvp[,1],xvp[,1]),ncol=nc))
}
#gpTrans is a simplified version where one transition probability for all states. 
#If recovery=TRUE have individuals from 1st columns of x all transitioning into the last column
gpTrans=function(x,Pr,seed,Recovery=FALSE, nc=ncol(x)){
	xv=as.vector(x[,1:c(nc/2,nc-1)[1+Recovery]])
	if(max(xv)==0) return(x); set.seed(seed); xv[xv>0]=sapply(xv[xv>0], function(y) rbinom(1,y,Pr));
	if(Recovery){ Tr=matrix(xv,ncol=nc-1); return(x+cbind(-Tr,rowSums(Tr))); }; return(x+matrix(c(-xv,xv),ncol=nc));
}

#Transition probabilities for travel. Multinomial faster than many binomials. rs=TRUE returns just total # people going to each province
reshfl2=function(x,M,seed,rs=TRUE,L=ncol(M),xnm=diag(x)){
	set.seed(seed); if(max(x)>0) xnm[,x>0]=apply(matrix(rbind(x,M)[,x>0],nrow=L+1), 2, function(y) rmultinom(1,y[1],y[-1])); 
	if(rs) return(rowSums(xnm)); return(xnm); 
}

#Modifier of travel matrix as people sick and/or tested positive less likely to travel by a proportion pstay
Mstay=function(M,pstay,Mod=M*(-(diag(ncol(M))-1)*(1-pstay))) Mod+diag(1-colSums(Mod))

meansd=function(x,wtR=1+x*0,wts=wtR/sum(wtR)){ xm=weighted.mean(x,wts); return(c(xm,sqrt(sum(wts*(x-xm)^2)))); }

#Main function handling all within-day transitions over state space S
m3iter=function(S,parms,seed,Ns=parms[,"N"]){
	#First testing results come back
	S0k=S[,Tk]; S[,c(Tn,Tk)]=gpTransB(S[,c(Tn,Tk)],as.vector(parms[,"tauI"]%*%t(c(parms[1,c("tauEAf","tauEAf")],1,1))),seed+3); 
	#pos is local active case prevalence
	S[,"Nt"]=S[,"Nt"]+rowSums(S[,Tk]-S0k); pos=rowSums(S[,Tk])/Ns; 
	#Disease progression: pi is the fraction of people who never show symptoms (modeled implicitly)
	pi=0.2; S[,c(Di,"R")]=gpTrans(S[,c(Di,"R")],parms[1,"rho"],seed,TRUE); 
	S[,c(Da,"R")]=gpTrans(S[,c(Da,"R")],pi*prod(parms[1,c("sig","rho")]),seed,TRUE); 
	S[,c(Da,Di)]=gpTrans(S[,c(Da,Di)],parms[1,"sig"]*(1-pi),seed+1); 
	S[,c("E","A")]=gpTrans(S[,c("E","A")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+2); 
	S[,c("E","SA")]=gpTrans(S[,c("E","SA")],parms[1,"alpha"]*parms[1,"s"],seed+2.5);

	#Closure decisions: C tracks # days dC since last closure (ie units 1/eps rather than integers).
	#pastinC tracks whether initial closure (specified as Tg=Tl=0 or =1 in closeinappl) has expired
	rc=range(parms[,c("Tl","Tg")]); pastinC=rc[1]<1 & rc[2]>0; 
	dC=S[,"C"]/parms[1,"eps"]; decide=pastinC & (dC==0 | dC>Resurge-1); 
	Tg=parms[1,"Tg"]*(1 - (min(S[,"C"])>0 & any(!decide)));
	#if pastinC, then decide specifies which regions are reconsidering opening/closing (enacted closures last Resurge # days)
	#For multiscale approaches (0<Tl<1 and 0<Tg<1), global closures override local decisions
	closed=parms[1,"eps"]*(pos > (parms[,"Tl"]*decide - 10*((sum(S[,Tk])/sum(Ns))>Tg))); 
	if(pastinC | rc[2]==0) S[,"C"]=S[,"C"]*(!decide)+closed;
	
	#Make modified travel matrices for people feeling sick and/or tested positive, then implement travel
	M=Mc=parms[,-(1:NP)]; Mc=M[1:nrow(M),]*(1-closed); diag(Mc)=diag(Mc)+1-colSums(Mc); 
	McA=abind(Mc, Mstay(Mc,parms[1,"eta"]), Mstay(Mc,parms[1,"r"]), Mstay(Mc,1-(1-parms[1,"eta"])*(1-parms[1,"r"])), along=3);
	Ss=reshfl2(S[,"S"],Mc,seed+4,FALSE); Rearr=apply(rbind(seed+(5:14), c(1,1,1:2,1:2,3:4,3:4), S[,c("R","E",Da,Di)]), 2, function(x) reshfl2(x[-(1:2)],McA[,,x[2]],x[1]))
	#Population size-dependent transmission rate
	Pstar=rowSums(cbind(Ss,Rearr)); dlt=parms[1,"dlt"]; Pmod=(1-dlt)/Pstar + dlt/(1+parms[1,"ftw"]*Pstar);
	#Baseline infection probability. The last term with rowSums is how many people end up in each county for the day.
	Infect = parms[1,"beta"]*(parms[1,"w"]*(1-closed) + (1-parms[1,"w"])*(1-parms[1,"epsP"]*(1-exp(-parms[1,"omg"]*pos)))) * Pmod
	#Class-specific infection modifiers. Net infection Pr is then 1-Pr(avoid infection by anyone).
	modK=1-parms[1,"eta"]; modS=1/parms[1,"s"]-1; modA=cbind(Infect,modK*Infect,modS*Infect,modS*modK*Infect);
	Infects=1 - apply((1-cbind(modA,2*modA))^Rearr[,-(1:2)], 1, prod)
	#Implement infection and move susceptibles and newly exposeds back to home county
	S[,c("S","E")]=cbind(0,S[,"E"]) + colSums(gpTransB(cbind(Ss,0*Ss),round(Infects,10),seed+15))
	return(S)
}; FUN=m3iter

closein=function(threshs,tsNt,tm){ 
	if(max(tsNt)<inCstart) return(1+0*threshs) #Closures not yet started
	if((tm-which(tsNt>=inCstart)[1])<inClen) return(0*threshs) #Initial closure started
	return(threshs) #Initial closure over
}
#Implements decline in testing time from initial 10 days, and allows parameters (if specified in retest) to change after re-opening.
closeinappl=function(parms,TS,tm=dim(TS)[3],delayInit=10){
	Nta=colSums(t(t(TS[,"Nt",]))); parms[,c("Tl","Tg")]=closein(parms[,c("Tl","Tg")],Nta,tm);
	rc=range(parms[,c("Tl","Tg")]); if("retest"%in%ls(envir=.GlobalEnv) & rc[1]<1 & rc[2]>0) parms[,colnames(retest)]=retest;
	tauT=parms[,"tauI"]%*%t(c(parms[1,"cvTl"],1))
	parms[,"tauI"]=pmin(tauT[,1]-sum(Nta>50)*(tauT[,1]-tauT[,2])/55, tauT[,2])
	if(max(Nta>inCstart)==0) parms[,"omg"]=0; parms[,"cvTl"]=0; return(parms);
}
Resagg=function(TS,parms,plotgive=TRUE){  tmp=t(apply(TS,2:3,sum))/sum(parms[,"N"]); tmp2=cbind(tmp[,"R"],rowSums(tmp[,2:10])*1e2,rowSums(tmp[,Tk])*1000,tmp[,"C"]); if(plotgive) matplot(tmp2,type="l"); return(tmp2); }
#Function to implement simulations. InitInf sets which stages the initially sick people are in, InitInf=2 default is all sick initially exposed
msim3=function(FUN,parms,Trun=365,seed0=11,plotgive=FALSE,InitInf=c(2)){
	#assign initial infections
	L=nrow(parms); nI=length(InitInf); Ns=parms[,"N"]; 
	set.seed(seed0); Infs=NI0=apply(parms[,c("N","n0")], 1, function(x) rbinom(n=nI,x[1],prob=x[2]/nI)); if(nI>1){ Infs=t(Infs); NI0=rowSums(Infs); };
	
	TS=array(0,dim=c(L,length(Snm),1)); TS[,c(1,InitInf),1]=cbind(Ns-NI0,Infs); TS[,13,1]=0; colnames(TS)=Snm; 
	set.seed(seed0+2); Seeds=rnorm(Trun,1e6,1e6);
	Mn=round(parms[1,"Mfact"]*parms[,-(1:NP)]*(1-diag(L))); diag(Mn)=Ns-colSums(Mn); parms[,-(1:NP)]=Mn%*%diag(1/colSums(Mn));
	for(i in 2:Trun) TS=abind(TS,FUN(TS[,,i-1],closeinappl(parms,TS),seed=Seeds[i]),along=3); TS[,"C",]=parms[1,"eps"]*(TS[,"C",]>0);
	#Different levels of aggregation in model output. plotgive=3 or 3.5 give shortest output form (tracking only infections and costs) in integer format to reduce output size
	if(plotgive=="TS") return(TS); if(plotgive==TRUE) return(Resagg(TS,parms,plotgive==1));
	if(plotgive%in%c(3,3.5)){ TS[,"C",]=TS[,"C",]*matrix(Ns,nrow(parms),Trun); TSn=apply(TS,c(2,3),sum); out=matrix(as.integer(TSn),nrow(TSn),ncol(TSn)); if(plotgive==3.5) return(out[c(1,13),]); return(rbind(out,colMeans(TS[,"C",]>0))); }
	#In fitting also tracked mean and variance in proportion distancing (omgs) and proportion of infections in Toronto (fracTor)
	if(plotgive=="fit"){ omgs=t(apply(apply(TS[,Tk,],c(1,3),sum)*matrix(1/Ns,L,Trun), 2, function(x) meansd(1-exp(-parms[1,"omg"]*x),Ns)));
	propCits=t(TS[,"Nt",]); states=t(apply(TS,c(2,3),sum)); return(cbind(states,omgs,propCits)); }
}


#Running best-fit, local strategy model. Not changing any parameters upon re-opening, so retest does not exist. 
# Plotting function below.
x=(msim3(m3iter,parms,365,plotgive="TS",seed0=1,InitInf=c(2)));
SpatiotempPlot(log(apply(x[,2:10,],c(1,3),sum)),figtitle="Log # sick (E+A+I)",YaxN="County",XaxN="Day")
SpatiotempPlot(log(apply(x[,Tk,],c(1,3),sum)),figtitle="Log # sick and tested",YaxN="County",XaxN="Day")
SpatiotempPlot(x[,"C",],figtitle="Closures enacted (red)",YaxN="County",XaxN="Day")
SpatiotempPlot((diag(1/parms[,"N"])%*%apply(x[,2:10,],c(1,3),sum)),figtitle="Proportion of county pops sick",YaxN="County",XaxN="Day")



library(fields);
SpatiotempPlot=function (spacetime, Xax = 1:dim(spacetime)[2], Yax = 1:dim(spacetime)[1], 
    XaxN = "X", YaxN = "Y", figtitle = "Title", Zlim = c(0, max(spacetime, 
        na.rm = TRUE)), cont = NULL, cexAx = 1, contPlot = spacetime, 
    cexCont = 1.5 * cexAx, Conts = NULL, contSpan = 1, palette = 1) 
{
    require(fields)
    spacetime[is.na(spacetime)] = Zlim[1] - 1
    COL = rev(rainbow(1000, start = 0, end = 0.7))
    if (palette > 1) {
        require(pals)
        COL = list(parula(1000), head(tail(parula(1000), -50, -50)))[[palette - 1]]
    }
    if (length(Zlim) == 1) {
        Zlim = quantile(spacetime, c(Zlim, 1 - Zlim), na.rm = T)
        Rng = range(spacetime)
        if (Zlim[1] < Rng[1]) 
            Zlim[1] = Rng[1]
        if (Zlim[2] > Rng[2]) 
            Zlim[2] = Rng[2]
    }
    spacetime[which(is.na(spacetime), arr.ind = TRUE)] = max(Zlim) + 1
    image.plot(x = Xax, y = Yax, z = t(spacetime), zlim = Zlim, 
        xlab = XaxN, ylab = YaxN, cex.axis = cexAx, cex.lab = cexAx, 
        legend.cex = cexAx, main = figtitle, col = COL)
    box()
    if (!is.null(cont)) {
        if (abs(log(max(contPlot, na.rm = T), 10)) > 2) 
            options(scipen = -10)
        if (contSpan > 1) {
            smoo1 = t(apply(contPlot, 1, function(x) supsmu(1:ncol(contPlot), 
                x, span = contSpan/ncol(contPlot))$y))
            smoo2 = t(apply(smoo1, 1, function(x) supsmu(1:ncol(smoo1), 
                x, span = contSpan/ncol(smoo1))$y))
            contPlot = smoo2
        }
        if (is.null(Conts)) 
            contour(x = Xax, y = Yax, z = t(contPlot), add = T, 
                col = cont, lwd = 1.5 * cexCont, labcex = cexCont)
        if (!is.null(Conts)) 
            contour(x = Xax, y = Yax, z = t(contPlot), levels = Conts, 
                add = T, col = cont, lwd = 1.5 * cexCont, labcex = cexCont)
        if (abs(log(max(contPlot, na.rm = T), 10)) > 2) 
            options(scipen = 0)
    }
}





#MODEL FITTING:
# LLfunMassBeta=function(betaTry,ftw,blog=FALSE){
	# CONT=c(ftw,dlt=3,cvTl=0)[1:3]; if(blog) betaTry=exp(betaTry); 
	# inCstart=50; inClen=1000; pops=colSums(DATA$Msave); parms=cbind(N=pops,s=0.2,Tg=0.0000001,Tl=0.0001,beta=betaTry/mean(pops),cvTl=CONT[3],tauI=0.,tauEAf=0.5,dlt=CONT[2],eps=0,epsP=0,w=0.45,omg=1e4,r=0.2,eta=0.8,Mfact=1,ftw=CONT[1],alpha=0.4,sig=0.4,rho=0.67,n0=1e-5,DATA$Msave[-50,]); NP=ncol(parms)-nrow(parms);
	# x=1-msim3(m3iter,parms,2e2,plotgive=3.5,seed0=1,InitInf=c(2))[1,]/sum(pops); if(var(tail(x,20))>0.01) x=1-msim3(m3iter,parms,1e3,plotgive=3.5,seed0=1,InitInf=c(2))[1,]/sum(pops); 
	# return(abs(tail(x,1)-0.65))
# }
# library(doParallel); library(foreach); registerDoParallel(cores=6)
# storeSatInc.65dat=as.matrix(expand.grid(ftw=exp(seq(-17,0,len=16)),dlt=seq(0,0.85,len=10),cvTl=0)); 
# storeSatInc.65=foreach(i=1:nrow(storeSatInc.65dat),.combine=rbind,.packages=c("abind","nloptr"))%dopar%{ dt=storeSatInc.65dat[i,]; x=nloptr(opts=list("algorithm"="NLOPT_GN_DIRECT","xtol_rel"=1e-2,"maxeval"=40),blog=TRUE,ftw=dt,x0=10,eval_f=LLfunMassBeta,lb=-6,ub=14); 
  # if(x$objective>0.025) x=nloptr(opts=list("algorithm"="NLOPT_GN_DIRECT","xtol_rel"=1e-20,"maxeval"=200),blog=TRUE,ftw=dt,x0=10,eval_f=LLfunMassBeta,lb=-6,ub=14); c(dt,x$solution,x$objective); }
# storeSatIncTrm.65=storeSatInc.65[-which(storeSatInc.65[,2]==0),]
# stopImplicitCluster()


#Making a function pre-fitted to give the beta0I value that, for a given xi and c, leads to 65% infections without testing (ie, no mitigation)
splinefun2=function(dat,k=20,cats=1:3){
	require(mgcv); datDF=data.frame(dat); datnm=names(datDF);
    xgam=gam(as.formula(paste(datnm[1], "~", "s(", datnm[cats[2]], ",", datnm[cats[3]], ",k=", k, ")")), data = datDF)
	function(a,b){ newd=data.frame(a,b); names(newd)=names(xgam$var.summary);
	predict.gam(xgam, newdata = newd, type = "link"); }
}
storeSatIncTrm.65=DATA$storeSatIncTrm.65
beta.ftwSatInc.65Log=splinefun2(cbind(storeSatIncTrm.65[,4],log(storeSatIncTrm.65[,1]),storeSatIncTrm.65[,2]),k=100); beta.ftwSatInc.65=function(ftw,dlt) exp(beta.ftwSatInc.65Log(log(ftw),dlt));



library(nloptr); library(mgcv);
bincase=function(x,pdl=log(popd)){ pdl[3]=10; aggregate(x~as.numeric((pdl<5)+(pdl<5.5)+(pdl<6.3)+(pdl<6.4)+(pdl<7)+(pdl<9)),FUN=sum)[c(2:4,1,5:7),2]; }
popd=c(465.4,37.0,316.6,25.7,25.7,28.9,28.9,22.4,22.4,19.0,35.1,19.0,241.0,585.9,4149.5,1040.0,64.0,64.0,520.4,465.4,232.6,37.8,121.5,370.4,17.4,48.0,25.7,42.1,210.1,42.0,132.4,17.4,18.5,18.5,57.3,37.0,19.0,6.8,7.4,7.4,4.2,4.2,4.2,2.4,0.3,2.7,0.6,0.4,0.4)
casepbin=c(0.4,0.17,0.09,0.0656,0.085,0.11,0.08) #Empirical proportions of cases
LLfun=function(parmsTry,parmsFit=parmnames,extras=c(Trun=300,daystart=36,dayTfin=64,dayRfin=100,model=1),reps=5){
	print(parmsTry); parmsTry0=parmsTry;
	#Fitting these parameters on a log scale:
	parmsTry[parmsFit=="ftw"]=exp(parmsTry[parmsFit=="ftw"]); parmsTry[parmsFit=="omg"]=exp(parmsTry[parmsFit=="omg"]); 
	Beta=beta.ftwSatInc.65(parmsTry[parmsFit=="ftw"],parmsTry[parmsFit=="dlt"]); 
	parms=cbind(N=pops,Mfact=1,s=0.2,Tg=1,Tl=exp(-9.14),beta=Beta/mean(pops),cvTl=0.045,tauI=0.3,tauEAf=0,dlt=0.27,epsP=0.62,eps=0.62,w=0.45,omg=44356,r=0.19,eta=0.8,ftw=4.785e-06,alpha=0.4,sig=0.4,rho=0.67,n0=0.0001,DATA$Msave[-50,]);
	parmsTry=pmax(parmsTry,0); parms[,parmsFit]=matrix(parmsTry,nrow(parms),length(parmsFit),byrow=TRUE);
	if(!"epsP"%in%parmsFit) parms[,"epsP"]=parms[,"eps"]; if(is.na(reps)) return(parms);
	
	sim=msim3(m3iter,parms,extras["Trun"],plotgive="fit",seed0=1); for(i in 2:reps) sim=abind(sim, msim3(m3iter,parms,extras["Trun"],plotgive="fit",seed0=i*1e3), along=3);
	lenR=extras["dayRfin"]; FVN=7; simp=array(dim=c(lenR+1,5+FVN,0)); for(i in 1:reps){ 
		trg=sim[,,i][which(sim[,"Nt",i]>49)[1]+(0:lenR),]; if(nrow(trg)<lenR) next; 
		nts=t(apply(trg[,-(1:(length(Snm)+2))],1,bincase)); simp=abind(simp, cbind(trg[,"Nt"], rowSums(trg[,c(Da,Di,"R")])/trg[,"Nt"], trg[,length(Snm)+1:2], diag(0+1/rowSums(nts))%*%nts, rowSums(trg[,Tk])),  along=3)
	}; simp[,2,][simp[,2,]==Inf]=max(simp[,2,][simp[,2,]!=Inf]);
	if(length(dim(simp))<3) return(1e6);
	
	to=head(diff(head(DATA$npTrim2,-1)),-9); tp=apply(simp[1:(length(to)+1),1,],2,diff); 
	RV=meansd(simp[40:80,2,]); DV=rowMeans(simp[50,3:4,])+c(0,1e-6); sta=round(c(DV,NA,RV*c(1,1)),2); 
	FV=t(apply(simp[(length(to)+1)+(-25:1),5:(4+FVN),],2,meansd))
	
	if(plotgive){  cum=1+(max(to)>1e4); matplot(list(tp,log(tp))[[cum]],type="l",xlab="Day since 50th positive case",ylab="Daily reported cases",lwd=1.5,lty=1,col=1,main=paste0(sta,collapse="_"),ylim=rbind(range(to),log(range(to)))[cum,]); points(cbind(to,log(to))[,cum],pch=16,col=2); round(c(DV*(1:2),NA,RV*(1:2)),2); }
	if(plotgive==2){ plot(log(popd),log(1e5*sum(caseCt)*FV[,1]/pops),col=2,xlab="Log density",ylab="Log cases per 100,000"); points(log(popd),log(1e5*caseCt/pops),col=4); legend("bottomright",col=c(4,2),pch=16,c("Data","Model"),box.col=0); }

	Tlik=colMeans(apply(cbind(to,tp), 1, function(x) abs(dpois(x[1],lambda=x[-1],log=TRUE))))
	Tlik=mean(tail(Tlik,-30))*1; Rlik=abs(dnorm(DATA$testRat,RV[1],RV[2],log=TRUE)); Dlik=dnorm(DATA$Distd50,DV[1],DV[2]/2,log=TRUE); 
	Flik=5*sum(apply(cbind(FV,casepbin), 1, function(x) abs(dnorm(x[3],x[1],pmax(x[2],.01),log=TRUE))))
	pT=FV[1,1]/sum(FV[,1]); LL=sum(abs(Tlik),4*abs(Rlik),5*(RV[2]>12),0*(DV[1]>.6 | DV[1]<0.4),100*(pT>.4 | pT<0.2),abs(Dlik),abs(Flik)); 
	stats=round(c(parmsTry0,RV[1],DV[1],FV[c(1,2,6),1],Tlik,Rlik,Dlik,Flik,LL,min(c(1e8,REPORT[,ncol(REPORT)-1]),na.rm=TRUE)),2); assign("REPORT",rbind(REPORT,stats),.GlobalEnv); print(tail(REPORT,1));
	return(pmin(LL,1e8));
}



nrep=4; extras=c(Trun=260,daystart=36,dayTfin=64,dayRfin=120,model=3.5); inCstart=325; inClen=300; 
pfit=c(1:7)[-4]; fitnames=c("tauI","omg","eps","epsP","ftw","dlt","cvTl"); parmsFit=fitnames[pfit]; 
parmStr=c(0.4,12.5,0.68,0.6,-12,0.1,150/12e3); lowerBex=c(0.3,9.75,0.55,0.3,-15,0.1,100/12e3); upperBex=c(0.95,14,0.8,0.8,-9,0.85,600/12e3); 
REPORT=matrix(nrow=0,ncol=length(pfit)+11); colnames(REPORT)=c(fitnames[pfit],"RV","DV","FVt","FVp","FVs","Tlik","Rlik","Dlik","Flik","LL","pbLL");
x=nloptr(parmsFit=fitnames[pfit],opts=list("algorithm"="NLOPT_GN_ISRES","xtol_rel"=1e-5,"maxeval"=400,"print_level"=2), x0=parmStr[pfit],eval_f=LLfun,lb=lowerBex[pfit],ub=upperBex[pfit],reps=nrep,extras=extras); 
y=nloptr(parmsFit=fitnames[pfit],opts=list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1e-15,"maxeval"=1200,"print_level"=2), x0=x$solution,eval_f=LLfun,lb=lowerBex[pfit],ub=upperBex[pfit],reps=nrep,extras=extras); 


parms0=LLfun(parmsTry,parmsFit,extras,reps=NA) #Just for convinience to get the final parameter set from fitting
inCstart=325; parmsTry=c(0.3,10.7,0.615,-12.25,.27,540/12e3) #FINAL FIT
LLfun(parmsTry,parmsFit=parmsFit,extras=extras,reps=2)

