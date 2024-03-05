library("copula")
library("scatterplot3d") # Fancy 3D plain scatterplots
#library("ggplot2") 
#library("grid") 	     # Useful package to set ggplot plots one next to the other
library("CDVine")	     # For selecting an appropriate Copula Family and for Kendall's plot
#library("evd")    	     # For loading pgumbel, dgumbel ,rgumbel and also qgumbel
#library("ismev")	     # For gum.fit
library("VineCopula")
#library("graphics")		# For persp
library("pcaPP")		# For cor.fk()
library("fitdistrplus") # for fitdist()
#library("stats")
library("FAdist")
library("plot3D")  #for persp3d
library("rgl")	# for axis3d
library("asbio")
library("Deriv")      #for differentiation in Frank copula
library("visreg")
library("copBasic")
# install.packages("spacetime");install.packages("sp");install.packages("spcopula", repos="http://R-Forge.R-project.org",dependencies=TRUE)
library("spcopula")		#for kendallRP

library("survival")

remove(list=ls())     # for visreg2D
#GUDDU
P=read.csv("C:\\Users\\Musti Tanvir\\Downloads\\Guddudata.csv")[,3]
V=read.csv("C:\\Users\\Musti Tanvir\\Downloads\\Guddudata.csv")[,4]


cor = cor.test(P,V,method="pearson")
tau = cor.test(P,V,method="kendall",exact=FALSE)
rho = cor.test(P,V,method="spearman",exact=FALSE)
rbind(c("method", "estimate", "pvalue"),cbind( c("pearson", "kendall", "spearman"), as.numeric(c(cor$estimate,tau$estimate,cor$estimate)),as.numeric(c(cor$p.value,tau$p.value,cor$p.value))))
plot(P,V)
plot(pobs(P),pobs(V),col="blue", xlab= "u", ylab="v", cex.lab=1.5, pch=19,cex.axis=1.3)
summary(P)
summary(V)
plot(P,V)
chi.plot(P,V, cex.lab=1.4,cex.axis=1.3,pch=16, col="blue")
chiplot(cbind(P,V),which=1)
u=pobs(P); v=pobs(V)
BiCopKPlot(u,v,cex.axis=1.3,cex.lab=1.3)
pl=BiCopKPlot(u,v,PLOT=FALSE)
par(new=TRUE)
plot(pl$W.in,pl$Hi.sort,pch = 19,col="blue",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",cex.axis=1.3,cex.lab=1.3)
par(new=FALSE)

#
# kfuncCOP in copBasic is meeant for The Kendall (Distribution) Function of a Copula
#

##      COPULA 
#Independence test
library("FAdist")
m=pobs(cbind(P,V))  # Pseudo observations are the observations in the [0,1] interval.
empsamp=indepTestSim(length(P),p=2,N=10000)
i=indepTest(m,empsamp) # fail to accept the null hyp that P and V are independent
dependogram(i)
shapeP=-0.03659; scaleP=7074.7; locP=9137.0
shapeV=0.10064; scaleV=205780; locV=205660
paraP=c(shape=shapeP, scale=scaleP, loc=locP)
paraV=c(shape=shapeV, scale=scaleV, loc=locV)
distP=function(P) pgev(P,shape=shapeP, scale=scaleP, loc=locP)
distV = function(V) pgev(V,shape=shapeV, scale=scaleV, loc=locV)
cdfP=pgev(P,shape=shapeP, scale=scaleP, loc=locP)
cdfV=pgev(V,shape=shapeV, scale=scaleV, loc=locV)
rdistP=function(n) rgev(n,shape=shapeP, scale=scaleP, loc=locP)
rdistV = function(n) rgev(n,shape=shapeV, scale=scaleV, loc=locV)
ddistP=function(P) dgev(P,shape=shapeP, scale=scaleP, location=locP)
ddistV = function(V) dgev(V,shape=shapeV, scale=scaleV, loc=locV)
u=cdfP; v=cdfV;
# Upper tail coefficient estimator CFG
SUMM=0
for (i in 1:length(P))
{	SUMM=SUMM+log(sqrt(log(1/u[i])*log(1/v[i]))/log((1/max(u[i],v[i]))^2))
}
lambda_CFG=2-2*exp(1/length(P) * SUMM)

############################### GUMBEL COPULA
###############################################################

tau=cor.fk(P,V)
theta=1/(1-tau)
tail=2-2^(1/theta)
BiCopPar2TailDep(4, theta)  #family number 4 for gumbel
cdfP=distP(P)
cdfV=distV(V)

g=gumbelCopula(theta,dim=2)
r=rCopula(100000,g)
pbs=pobs(cbind(P,V))
plot(r, col="lightblue",pch=19, xlab="u", ylab="v", cex.lab=1.4)
points(pbs,pch=19,cex=0.7,col="blue")
r=cbind(qgev(r[,1], shape=-0.03659, scale=7074.7, loc=9137.0),qgev(r[,2],shape=0.10064, scale=205780, loc=205660))
par(mar=c(4,5,1,1))
plot(r,col="lightblue",pch=19,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), cex.lab=1.2,yaxt="n")
axis(2, at=seq(0,4)*10^6, labels=sprintf("%d",seq(0,4)*10^6))
points(P,V,pch=19,cex=0.7,col="blue")
rCDF_PVgum = exp(-((-log(u))^theta+(-log(v))^theta)^(1/theta))
CDF_f=function(x,y) exp(-((-log(distP(x)))^theta+(-log(distV(y)))^theta)^(1/theta))
par(mar=c(1,3.2,1,1))
Pseq=seq(min(P), max(P), length.out=30)
Vseq=seq(min(V),max(V), length.out=30)
cdfPV2=outer(Pseq, Vseq, function(x,y) CDF_f(x,y))
#colors
z=cdfPV2
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,1)
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp_lab_f <- function(top, bot, pmat, space, pos = c(-1, -1))
{
  coords_3d <- list(top = top, bot = bot)
  coords_2d <- lapply(coords_3d, function(v, pmat) unlist(trans3d(v[1], v[2], v[3], pmat)), pmat = pmat)
  coords_2d$mid <- (coords_2d$top + coords_2d$bot)/2  # mid is calculated from 2d-coordinates
  # coords_2d$mid <- unlist(trans3d(((top + bot)/2)[1], ((top + bot)/2)[2], ((top + bot)/2)[3], pmat)) if use mid in 3d-coordinates
  tb_diff <- coords_2d$top - coords_2d$bot
  angle <- 180/pi * atan2(tb_diff[2], tb_diff[1])
  names(angle) <- "angle"
  center <-  coords_2d$mid + sqrt(space^2 / sum(tb_diff^2)) * rev(abs(tb_diff)) * pos
  out <- list(angle = angle, center = as.data.frame(t(center)))
  return(out)
}
pmat <- persp(Pseq, Vseq, cdfPV2, theta = -40, phi = 20,expand=0.5)

persp(x=Pseq,y=Vseq,z=cdfPV2, col=fcol,xlab="", ylab="", zlab="\ncdf for Gumbel",box=TRUE,ticktype = "detailed",expand=0.5, zlim=c(0,1), phi=15, theta=-50, lphi=15, ltheta=-50,shade=0.7,cex.lab=1.2)
x_top <- c(max(Pseq), max(Vseq), min(cdfPV2))
x_bot <- c(min(Pseq), max(Vseq), min(cdfPV2))
y_top <- c(max(Pseq), max(Vseq), min(cdfPV2))
y_bot <- c(max(Pseq), min(Vseq), min(cdfPV2))
z_top <- c(max(Pseq), min(Vseq), max(cdfPV2))
z_bot <- c(max(Pseq), min(Vseq), min(cdfPV2))
xlab_param <- persp_lab_f(x_top, x_bot, pmat, 0.1, pos = c(1, -1))
ylab_param <- persp_lab_f(y_top, y_bot, pmat, 0.1)
zlab_param <- persp_lab_f(z_top, z_bot, pmat, 0.1)
text(xlab_param$center-c(-0.45,0.12), srt = xlab_param$angle+28, expression(P (m^3/s)),cex=1.2)
text(ylab_param$center-c(0.41,0.15), srt = ylab_param$angle+175, expression(V (day-m^3/s)),cex=1.2)

#persp3D(x=Pseq,y=Vseq,z=cdfPV2,box=TRUE,ticktype = "detailed",expand=0.5, zlim=c(0,1), phi=15, theta=-50,shade=0.2,xlab="\nP",ylab="\n\nV", zlab="\nProbability")
#axis3d('y', nticks=9, labels=seq(15000,500000, length=9))

# QQ PLOT FOR GUMBEL COPULA (KENDALL PLOT)
K= function(z) {z-log(z)*z/theta}

	#obtaining zi
z=0*P;
Ps=sort(P)
Vs=rep(0,length(P))
for (i in 1:length(P))
{	Vs[i]=V[which(P==Ps[i])]
}
for (i in 1:length(P))
{	z[i]=0;
	for (k in 1:length(P))
	{ 	if ((Ps[i]>Ps[k])&(Vs[i]>Vs[k])) z[i]=z[i]+1;
	}
	z[i]=z[i]/(length(P)-1);
}
Kn=0*P
for (i in 1:length(P))
{ 	Kn[i]= sum(z[i]>=z)/length(P)
}
Kp_gh=K(z)
par(mar=c(4,5,1,1))
plot(Kn,Kp_gh,pch=19, col="blue", xaxt="n", yaxt="n",xlab=expression(K[n]),ylab=expression(K[theta*n]),cex.lab=1.5,cex=1.2)
abline(0,1,lwd=2, col="red")
axis(2,cex=1.2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
axis(1,cex=1.2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
		#another approach from https://www.r-bloggers.com/kendalls-function-for-copulas/
i=rep(1:length(P),each=length(P))
j=rep(1:length(P),length(P))
S=((P[i]>P[j])&(V[i]>V[j]))
Z=tapply(S,i,sum)/(length(P)-1)
#plot(ecdf(Z))
plot(sort(Z),(1:length(P))/length(P),type="s",col="red",xlab="",ylab="")
rh=0.001
phi=function(t){(-log(t))^(theta)}
dphi=function(t){(phi(t+h)-phi(t-h))/2/h}
Kg=Vectorize(k)
par(new=TRUE)
Zs=sort(Z)
KK=Kg(sort(Z))
lines(Zs[which(KK!='NaN')],KK[which(KK!='NaN')],col="blue")
par(new=FALSE)


# GOODNESS-OF-FIT TEST
n=rep(0,length(P))
C_obs=rep(0,length(P))
Ps=sort(P)
Vs=rep(0,length(P))
for (i in 1:length(P))
{	Vs[i]=V[which(P==Ps[i])]
}
for (i in 1:length(P))
{	n[i]=0;
	for (k in 1:length(P))
	{ 	if ((Ps[i]<=Ps[k])&(Vs[i]<=Vs[k])) n[i]=n[i]+1;
	}
	C_obs[i]=(n[i]-0.44)/(length(P)+0.12);
}
CDF_f=function(x,y) exp(-((-log(distP(x)))^theta+(-log(distV(y)))^theta)^(1/theta))
C_comp=CDF_f(Ps,Vs)
k=7  #no of params
RMSE_gh=sqrt((1/(length(P)-k))*sum((C_obs-C_comp)^2))
AIC_gh=length(P)*log(RMSE_gh^2)+2*k

#   Sn statistic with p vaue
u=pobs(P)
v=pobs(V)
Cn=0*P
for (i in 1:length(P))
{	Cn[i]=sum(u[i]>=u & v[i]>=v)/length(P)
}
CDF_fuv=function(x,y) exp(-((-log(x))^theta+(-log(y))^theta)^(1/theta))
Cthetan=CDF_fuv(u,v)
Sn=sum((Cn-Cthetan)^2)
g=gumbelCopula(theta,dim=2)
gofCopula(g, cbind(P,V), method = "Sn",simulation = "pb",estim.method="itau")
plot(Cn,Cthetan,pch=19,xlim=c(0,1),ylim=c(0,1))
abline(0,1,lty=3,lwd=3,col="blue")

# RETURN PERIOD
TORgum=1/(1-CDF_PVgum)
Pseq2=seq(min(P), max(P)+5*min(P), length.out=60)
Vseq2=seq(2.5*10^5, max(V)+8*min(V), length.out=60)
TOR_gum=outer(Pseq2, Vseq2, function(x,y) 1/(1-(exp(-((-log(distP(x)))^theta+(-log(distV(y)))^theta)^(1/theta)))))
#colors
z=TOR_gum
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq2,y=Vseq2,z=TOR_gum, xlab="\nP", ylab="\n\nV", zlab="Time (years)",col=fcol,box=TRUE,expand=0.5,ticktype = "detailed", zlim=c(0,max(z)), phi=20, theta=-40)


#OR RETURN PERIOD
Pseq2=seq(0, max(P)+15*min(P), length.out=60)
Vseq2=seq(0, max(V)+35*min(V), length.out=60)
TOR_gum=outer(Pseq2, Vseq2, function(x,y) 1/(1-(exp(-((-log(distP(x)))^theta+(-log(distV(y)))^theta)^(1/theta)))))
g=gumbelCopula(theta,dim=2)
r=rCopula(100000,g)
pbs=pobs(cbind(P,V))
r=cbind(qgev(r[,1], shape=-0.03659, scale=7074.7, loc=9137.0),qgev(r[,2],shape=0.10064, scale=205780, loc=205660))
Pr=r[,1]
Vr=r[,2]
par(mar=c(4,5,1,1))
contour(x=Pseq2,y=Vseq2,z=TOR_gum,xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),cex.lab=1.2,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), levels=c(2,5,10, 15, seq(20,140, by=20),200),col="red",labcex=.9)
points(Pr,Vr,col="lightblue",pch=19)
points(P,V,pch=20,cex=1.5)
par(new=TRUE)
contour(x=Pseq2,y=Vseq2,z=TOR_gum,xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),lwd=2,xlab="",ylab="", levels=c(2,5,10, 15, seq(20,160, by=20),200,300,400,500),col="red",labcex=.9)

# AND RETURN PERIOD
T_AND_gum=1/(1-dist(P(P))-dist(V(V))+CDF_PVgum)

#Pseq2=seq(0, max(P)+15*min(P), length.out=60)
#Vseq2=seq(0, max(V)+50*min(V), length.out=60)
TAND_gum=outer(Pseq2, Vseq2, function(x,y) 1/(1-distP(x)-distV(y)+CDF_f(x,y)))
g=gumbelCopula(theta,dim=2)
r=rCopula(100000,g)
pbs=pobs(cbind(P,V))
rr=cbind(qgev(r[,1], shape=-0.03659, scale=7074.7, location=9137.0),qgev(r[,2],shape=0.10064, scale=205780, location=205660))
Pr=rr[,1]
Vr=rr[,2]
par(mar=c(5,5,1,1))
contour(x=Pseq2,y=Vseq2,z=TAND_gum,xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),cex.lab=1.2,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), levels=c(2,5,10, 15, seq(20,140, by=20),200),col="red",labcex=.9)
points(Pr,Vr,col="lightblue",pch=19)
points(P,V,pch=20,cex=1.5)
par(new=TRUE)
rcontour(x=Pseq2,y=Vseq2,z=TAND_gum,xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),lwd=2,xlab="",ylab="", levels=c(2,5,10, 15, seq(20,160, by=20),200,300,400,500),col="red",labcex=.9)


#SECONDARY

kendallRP(copula=g,cl=c(0.90,0.99,0.999,0.9999))
plot(seq(0.01,1,by=0.01),kendallRP(copula=g,cl=seq(0.01,1,by=0.01)),log="y", xlab=expression(t), ylab="kendall return period (years)",type="l",col="red")

Pseq2=seq(0, 80000, length.out=80)
Vseq2=seq(0, 80.5*10^6, length.out=80)
CDF_f=function(x,y) exp(-((-log(distP(x)))^theta+(-log(distV(y)))^theta)^(1/theta))
Tsec_gum=outer(Pseq2, Vseq2, function(x,y) 1/(1-K(CDF_f(x,y))))
g=gumbelCopula(theta,dim=2)
r=rCopula(100000,g)
pbs=pobs(cbind(P,V))
#rr=cbind(qgev(r[,1], shape=-0.03659, scale=7074.7, location=9137.0),qgev(r[,2],shape=0.10064, scale=205780, location=205660))
#Pr=rr[,1]
#Vr=rr[,2]
#par(mar=c(5,5,1,1))
contour(x=Pseq2,y=Vseq2,z=Tsec_gum,xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),cex.lab=1.2,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), levels=c(2,5,10, 15, seq(20,140, by=20),200),col="red",labcex=.9)
#points(Pr,Vr,col="lightblue",pch=19)
#points(P,V,pch=20,cex=1.5)
#par(new=TRUE)
contour(x=Pseq2,y=Vseq2,z=(1-1/Tsec_gum),xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),lwd=2,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), levels=c(0.99),col="red",labcex=.9)
max(z)
surCOP(1-u,1-v, cop=GHcop,par=theta)
Rp=1/surCOP(1-u,1-v, cop=GHcop,par=theta,dim=2)

contour(x=Pseq2,y=Vseq2,z=(1-1/Rp),xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),lwd=2,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), levels=c(0.9),col="red",labcex=.9)
max(z)
#CONDITIONAL RETURN PERIOD U>u, V<=v

Pseq2=seq(0, max(P)+25*min(P), length.out=60)
Vseq2=seq(0, max(V)+65*min(V), length.out=60)
Tcond_1=outer(Pseq2, Vseq2, function(x,y) 1/(1-CDF_f(x,y)/distP(x)))
g=gumbelCopula(theta,dim=2)
r=rCopula(100000,g)
pbs=pobs(cbind(P,V))
r=cbind(qgev(r[,1], shape=-0.03659, scale=7074.7, loc=9137.0),qgev(r[,2],shape=0.10064, scale=205780, loc=205660))
Pr=r[,1]
Vr=r[,2]
par(mar=c(5,5,1,1))
l=c(2,5,10, 15, seq(20,80, by=20), seq(100,500, by=100))
contour(x=Pseq2,y=Vseq2,z=Tcond_1,xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),cex.lab=1.2,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), levels=l,col="red",labcex=.9)
points(Pr,Vr,col="lightgrey",pch=19)
points(P,V,pch=19,cex=0.8)
par(new=TRUE)
contour(x=Pseq2,y=Vseq2,z=Tcond_1,xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),lwd=2,xlab="",ylab="", levels=l,col="red",labcex=.9)
#FINAL
Pseq3=c(1000,5000,10000,15000,20000)
#Vseq2=seq(0, max(V)+65*min(V), length.out=8)
Tcond_1f= function(x,y) 1/(1-CDF_f(x,y)/distP(x))
Tcond1=outer(Pseq3,Vseq2,Tcond_1f)
ltyy=c(2,3,4,5,6)
coll=c("blue","green","black","magenta", "orange", "purple")
i=1;
Vtemp1=Vseq2[which(Tcond1[i,]<500)]
plot(Vtemp1,Tcond1[i,][which(Tcond1[i,]<500)],log="y", xlab=expression(V (day-m^3/s)), ylab="Time (years)",type="l",xlim=c(0,max(Vtemp1)),ylim=c(1, 500),col=coll[i],lty=ltyy[i],lwd=2,cex.lab=1.2)
v1=Vseq2[60]
t1=Tcond_1[i,60]
p2=Pseq3[i]
text(max(Vseq2)+51000,t1, "P=",adj=c(1,0.5),col="red")
text(max(Vseq2)+51000,t1, Pseq3[i],adj=c(0,0.5),col="red")
for (i in 2:length(Pseq3))
{	par(new=TRUE)
	Vtemp=Vseq2[which(Tcond1[i,]<500)]
	plot(Vtemp,Tcond1[i,][which(Tcond1[i,]<500)],log="y",col=coll[i],lty=ltyy[i], xlab="", ylab="",type="l",xlim=c(0,max(Vtemp1)),ylim=c(1, 500),xaxt="n",yaxt="n",lwd=2)
	v1=Vseq2[60]
	t1=Tcond_1[i,60]
	p2=Pseq3[i]
	text(max(Vseq2)+51000,t1, "P=",adj=c(1,0.5),col="red")
	text(max(Vseq2)+51000,t1, Pseq3[i],adj=c(0,0.5),col="red")	
}

for( i in 1:length(Pseq3))
{	legend(0,500/(exp((i-1)/3)),legend=bquote("P=" ~.(Pseq3[i]) ~ m^3/s),bty="n",lty=ltyy[i],col=coll[i],lwd=2,cex=1.1,merge=TRUE)
}
par(new=FALSE)	
#
Vseq3=c(10^5,seq(4*10^5, 1000000, by=1.5*10^5))
Tcond_1f= function(x,y) 1/(1-CDF_f(x,y)/distP(x))
Tcond11=outer(Pseq2,Vseq3,Tcond_1f)
ltyy=c(2,3,4,5,6,2)
coll=c("blue","green","black","magenta", "orange", "purple")
i=1;
Ptemp=Pseq2[which(Tcond11[,i]<500)]
plot(Pseq2,Tcond11[,i],log="y", xlab=expression(P (m^3/s)), ylab="Time (years)",type="l",ylim=c(1,500),xlim=c(min(Pseq2),max(Pseq2)),col=coll[i],lty=ltyy[i],lwd=2,cex.lab=1.2)
v1=Vseq2[60]
t1=Tcond_1[i,60]
p2=Pseq3[i]
text(max(Vseq2)+51000,t1, "P=",adj=c(1,0.5),col="red")
text(max(Vseq2)+51000,t1, Pseq3[i],adj=c(0,0.5),col="red")
for (i in 2:length(Vseq3))
{	par(new=TRUE)
	plot(Pseq2,Tcond11[,i],log="y",col=coll[i],lty=ltyy[i], xlab="", ylab="",type="l",ylim=c(1,500),xlim=c(min(Pseq2),max(Pseq2)),xaxt="n",yaxt="n",lwd=2)
	v1=Vseq2[60]
	t1=Tcond_1[i,60]
	p2=Pseq3[i]
	text(max(Vseq2)+51000,t1, "P=",adj=c(1,0.5),col="red")
	text(max(Vseq2)+51000,t1, Pseq3[i],adj=c(0,0.5),col="red")	
}
par(new=FALSE)	
Vtext=sprintf("%d",Vseq3)
for( i in 1:length(Vseq3))
{	legend(35000,680/(exp((i-1)/3)),legend=bquote("V=" ~.(Vtext[i]) ~ day-m^3/s),bty="n",lty=ltyy[i],col=coll[i],lwd=2,cex=1.1,merge=TRUE)
}


#CONDITIONAL RETURN PERIOD U>u, V>v
Pseq2=seq(0, max(P)+25*min(P), length.out=60)
Vseq2=seq(0, max(V)+65*min(V), length.out=60)
Tcond_2=outer(Pseq2, Vseq2, function(x,y) (1/(1-distP(x)))*(1/(1-distP(x)-distV(y)+CDF_f(x,y))))
g=gumbelCopula(theta,dim=2)
r=rCopula(100000,g)
pbs=pobs(cbind(P,V))=cbind(qgev(r[,1], shape=-0.03659, scale=7074.7, location=9137.0),qgev(r[,2],shape=0.10064, scale=205780, location=205660))
Pr=r[,1]
Vr=r[,2]
par(mar=c(5,5,1,1))
l=c(2,5,10, 15, seq(20,80, by=20), seq(100,500, by=100))
contour(x=Pseq2,y=Vseq2,z=Tcond_1,xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),cex.lab=1.2,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), levels=l,col="red",labcex=.9)
points(Pr,Vr,col="lightgrey",pch=19)
points(P,V,pch=19,cex=0.8)
par(new=TRUE)
contour(x=Pseq2,y=Vseq2,z=Tcond_1,xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),lwd=2,xlab="",ylab="", levels=l,col="red",labcex=.9)

#FINAL V on x-axis
Pseq3=c(seq(1000,27000,by=5000))
Vseq2=seq(0, max(V)+65*min(V), length.out=30)
rTcond_2f= function(x,y) (1/(1-distP(x)))*(1/(1-distP(x)-distV(y)+CDF_f(x,y)))
Tcond2=outer(Pseq3,Vseq2,Tcond_2f)
ltyy=c(2,3,4,5,6,2)
coll=c("blue","green","black","magenta", "orange", "purple", "yellow")
i=1;
Vtemp1=Vseq2[Tcond2[i,]<500]
plot(Vseq2,Tcond2[i,],log="y", xlab=expression(V (day-m^3/s)), ylab="Time (years)",type="l",xlim=c(0,max(Vtemp1)),ylim=c(1, 500),col=coll[i],lty=ltyy[i],lwd=2,cex.lab=1.2)
for (i in 2:length(Pseq3))
{	par(new=TRUE)
	Vtemp=Vseq2[which(Tcond2[i,]<500)]
	plot(Vseq2,Tcond2[i,],log="y",col=coll[i],lty=ltyy[i], xlab="", ylab="",type="l",xlim=c(0,max(Vtemp1)),ylim=c(1, 500),xaxt="n",yaxt="n",lwd=2)
	v1=Vseq2[60]
	t1=Tcond_1[i,60]
	p2=Pseq3[i]
	text(max(Vseq2)+51000,t1, "P=",adj=c(1,0.5),col="red")
	text(max(Vseq2)+51000,t1, Pseq3[i],adj=c(0,0.5),col="red")	
}

ltext=bquote("P=" ~.(Pseq3[2]) ~ m^3/s)
#btyy=c(rep("n",times=(length(Pseq3)-1)),"o")
for( i in 1:length(Pseq3))
{	legend(1.15*10^6,10/(exp((i-1)/3)),legend=bquote("P=" ~.(Pseq3[i]) ~ m^3/s),bty="n",lty=ltyy[i],col=coll[i],lwd=2,cex=1.1,merge=TRUE)
}
par(new=FALSE)	

#FINAl P on x-axis
Vseq3=c(10^5,seq(4*10^5, 1700000, by=3*10^5))
Pseq2=seq(0,30000, length.out=12)
#Vseq2=seq(0, max(V)+65*min(V), length.out=8)
Tcond_2f= function(x,y) (1/(1-distP(x)))*(1/(1-distP(x)-distV(y)+CDF_f(x,y)))
Tcond2=outer(Pseq2,Vseq3,Tcond_2f)
ltyy=c(2,3,4,5,6,2)
coll=c("blue","green","black","magenta", "orange", "purple")
i=1;
Ptemp1=Pseq2[which(Tcond2[,i]<500)]
plot(Pseq2,Tcond2[,i],log="y", xlab=expression(P (m^3/s)), ylab="Time (years)",type="l",xlim=c(0,max(Ptemp1)),ylim=c(1, 500),col=coll[i],lty=ltyy[i],lwd=2,cex.lab=1.2)
for (i in 2:length(Vseq3))
{	par(new=TRUE)
	plot(Pseq2,Tcond2[,i],log="y",col=coll[i],lty=ltyy[i], xlab="", ylab="",type="l",xlim=c(0,max(Ptemp1)),ylim=c(1, 500),xaxt="n",yaxt="n",lwd=2)
}

Vtext=sprintf("%d",Vseq3)
#btyy=c(rep("n",times=(length(Pseq3)-1)),"o")
for( i in 1:length(Pseq3))
{	legend(14500,7.5/(exp((i-1)/3)),legend=bquote("V=" ~.(Vtext[i]) ~ day-m^3/s),bty="n",lty=ltyy[i],col=coll[i],cex=0.95,lwd=2,merge=TRUE)
}
par(new=FALSE)	

## Return Period calcs
medP=median(P)
medV=median(V)
mxP=max(P)
mxV=max(V)
TOR=function(x,y) 1/(1-CDF_f(x,y))
TAND=function(x,y) 1/(1-distP(x)-distV(y)+CDF_f(x,y))
Tsec=function(x,y) 1/(1-K(CDF_f(x,y)))
Tc1=function(x,y) 1/(1-CDF_f(x,y)/distP(x))
Tc2=function(x,y) (1/(1-distP(x)))*(1/(1-distP(x)-distV(y)+CDF_f(x,y)))
TORmed=TOR(medP,medV)
TORmax=TOR(mxP,mxV)
TANDmed=TAND(medP,medV)
TANDmax=TAND(mxP,mxV)
#Tsmed=Tsec(medP,medV)  #wrong method
#Tsmax=Tsec(mxP,mxV)
Tc1med=Tc1(medP,medV)
Tc1max=Tc1(mxP,mxV)
Tc2med=Tc2(medP,medV)
Tc2max=Tc2(mxP,mxV)
cbind(c(rep(c("TOR","TAND","Tc1","Tc2"),each=2)),round(c(TORmed,TORmax,TANDmed,TANDmax,Tc1med,Tc1max,Tc2med,Tc2max),3),rep(c("median","max"),4))



# AMH COPULA
# =============

# admissible range for tau is [-0.1817, 1/3]

tau=cor.fk(P,V)
f=function(theta) ((3*theta-2)/theta)-(2/3 * (1-1/theta)^2)*log(1-theta)-tau
theta = uniroot(f, c(0.1, 0.9), tol = 0.0001)$root
u=cdfP
v=cdfV
CDF_PV_AMH = u*v/(1-theta*(1-u)*(1-v))


g=amhCopula(theta,dim=2)
r=rCopula(100000,g)
pbs=pobs(cbind(P,V))
plot(r, col="lightblue",pch=19)
points(pbs,pch=19)
r=cbind(qgev(r[,1], shape=-0.03659, scale=7074.7, location=9137.0),qgev(r[,2],shape=0.10064, scale=205780, location=205660))
plot(rr,col="lightblue",pch=19)
points(P,V,pch=19,cex=0.7,col="blue")



Pseq=seq(0, max(P), length.out=60)
Vseq=seq(0, max(V), length.out=60)
cdfPV3=outer(Pseq, Vseq, function(x,y) distP(x)*distV(y)/(1-theta*(1-distP(x))*(1-distV(y))))
#colors
z=cdfPV3
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq,y=Vseq,z=cdfPV3, xlab="P", ylab="V", zlab="cdf for AMH",col=fcol,box=TRUE,expand=0.5,ticktype = "detailed", zlim=c(0,1), phi=20, theta=-40)

# QQ PLOT FOR AMH COPULA (KENDALL PLOT)
phi=function(z) log((1-theta*(1-z))/z)
dphi=function(z) (theta-1)/(z*(1-theta*(1-z)))
K= function(z) {z- phi(z)/dphi(z)}

	#obtaining zi
z=0*P;
Ps=sort(P)
Vs=rep(0,length(P))
for (i in 1:length(P))
{	Vs[i]=V[which(P==Ps[i])]
}
for (i in 1:length(P))
{	z[i]=0;
	for (k in 1:length(P))
	{ 	if ((Ps[i]>Ps[k])&(Vs[i]>Vs[k])) z[i]=z[i]+1;
	}
	z[i]=z[i]/(length(P)-1);
}
for (i in 1:length(P))
{ 	Kn[i]= sum(z[i]>=z)/length(P)
}
Kp=K(z)
plot(Kn,Kp	,pch=19, col="blue", xaxt="n", yaxt="n")
lines(Kn,Kn)
axis(2,cex=0.9,at=seq(0,1,0.2),labels=seq(0,1,0.2))
axis(1,cex=0.9,at=seq(0,1,0.2),labels=seq(0,1,0.2))
		#another approach from https://www.r-bloggers.com/kendalls-function-for-copulas/
i=rep(1:length(P),each=length(P))
j=rep(1:length(P),length(P))
S=((P[i]>P[j])&(V[i]>V[j]))
Z=tapply(S,i,sum)/(length(P)-1)
#plot(ecdf(Z))
plot(sort(Z),(1:length(P))/length(P),type="s",col="red",xlab="",ylab="")
h=0.001
phi=function(t) {log((1-theta*(1-t))/t)} 
dphi=function(t){(phi(t+h)-phi(t-h))/2/h}
kk=function(t){t-phi(t)/dphi(t)}
Kg=Vectorize(kk)
par(new=TRUE)
lines(sort(Z),Kg(sort(Z)),col="blue")
par(new=FALSE)



# GOODNESS-OF-FIT TEST
CDF_f=function(x,y) distP(x)*distV(y)/(1-theta*(1-distP(x))*(1-distV(y)))
C_comp=CDF_f(Ps,Vs)
MSE_amh=(1/(length(P)-5))*sum((C_obs-C_comp)^2)
RMSE_amh=sqrt(MSE_amh)
AIC_amh=length(P)*log(MSE_amh)+2*7

# RETURN PERIOD
TORamh=1/(1-CDF_PV_AMH)
Pseq2=seq(min(P), max(P)+150000, length.out=60)
Vseq2=seq(min(V), max(V)+250000, length.out=60)
TOR_amh=outer(Pseq2, Vseq2, function(x,y) 1/(1-CDF_f(x,y)))
#colors
z=TOR_amh
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq2,y=Vseq2,z=TOR_amh, xlab="P", ylab="V", zlab="Time (years)",col=fcol,box=TRUE,expand=0.5,ticktype = "detailed",  zlim=c(0,100), phi=20, theta=-40)

Pseq3=seq(4500, 16500, by=2000)
Vseq3=seq(20000,  1020000, length.out=60)
TOR_amh2=outer(Pseq3, Vseq3, function(x,y) 1/(1-CDF_f(x,y)))
for (i in 1:length(Pseq3))
{	
	plot(Vseq3,TOR_amh2[i,],col="blue", xlab="V", ylab="Time (years)",log="y",type="l",xlim=c(min(Vseq3),  max(Vseq3)),ylim=c( min(TOR_amh2), max(TOR_amh2)))
	v1=Vseq3[57]
	t1=TOR_amh2[i,57]
	p2=Pseq3[i]
	text(v1,t1, "P=",adj=c(1,0.5),col="red")
	text(v1,t1, Pseq3[i],adj=c(0,0.5),col="red")
	par(new=TRUE)
}
par(new=FALSE)	

# CookJohnson COPULA (a.k.a. Clayton)
# =======================================

tau=cor.fk(P,V)
theta = (2*tau)/(1-tau)
u=cdfP
v=cdfV
CDF_PV_CJ = (u^(-theta)+v^(-theta)-1)^(-1/theta)


g=claytonCopula(theta,dim=2)
r=rCopula(100000,g)
pbs=pobs(cbind(P,V))
plot(r, col="lightblue",pch=19, xlab="u", ylab="v", cex.lab=1.4)
points(pbs,pch=19,cex=0.7,col="blue")
rr=cbind(qgev(r[,1], shape=-0.03659, scale=7074.7, location=9137.0),qgev(r[,2],shape=0.10064, scale=205780, location=205660))
par(mar=c(4,5,1,1))
plot(rr,col="lightblue",pch=19,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), cex.lab=1.2,yaxt="n")
axis(2, at=seq(0,4)*10^6, labels=sprintf("%d",seq(0,4)*10^6))
points(P,V,pch=19,cex=0.7,col="blue")



CDF_f= function(x,y) (distP(x)^(-theta)+distV(y)^(-theta)-1)^(-1/theta)
Pseq=seq(min(P), max(P), length.out=50)
Vseq=seq(min(V), max(V), length.out=50)
cdfPV_CJ = outer(Pseq, Vseq, function(x,y) CDF_f(x,y))
#colors
z=cdfPV_CJ
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq,y=Vseq,z=cdfPV_CJ, xlab="P", ylab="V", zlab="cdf for clayton",col=fcol,box=TRUE,expand=0.5,ticktype = "detailed" ,zlim=c(0,1), phi=20, theta=-40)


# QQ PLOT FOR CLAYTON COPULA (KENDALL PLOT)
phi =function (z) (1/theta)*(z^(-theta)-1)
dphi=Deriv(phi,"z")
K= function(z) {z-phi(z)/dphi(z)}

	#obtaining zi
z=0*P;
Ps=sort(P)
Vs=rep(0,length(P))
for (i in 1:length(P))
{	Vs[i]=V[which(P==Ps[i])]
}
for (i in 1:length(P))
{	z[i]=0;
	for (k in 1:length(P))
	{ 	if ((Ps[i]>Ps[k])&(Vs[i]>Vs[k])) z[i]=z[i]+1;
	}
	z[i]=z[i]/(length(P)-1);
}
Kn=0*P
for (i in 1:length(P))
{ 	Kn[i]= sum(z[i]>=z)/length(P)
}
Kp_cj=K(z)
plot(Kn,Kp_cj,pch=19, col="blue", xaxt="n", yaxt="n",xlab=expression(K[n]),ylab=expression(K[theta*n]),cex.lab=1.5,cex=1.2)
abline(0,1,lwd=2, col="red")
axis(2,cex=1.2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
axis(1,cex=1.2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
		#another approach from https://www.r-bloggers.com/kendalls-function-for-copulas/
i=rep(1:length(P),each=length(P))
j=rep(1:length(P),length(P))
S=((P[i]>P[j])&(V[i]>V[j]))
Z=tapply(S,i,sum)/(length(P)-1)
#plot(ecdf(Z))
plot(sort(Z),(1:length(P))/length(P),type="s",col="red",xlab="",ylab="")
h=0.001
phi=function(t){(1/theta)*(t^(-theta)-1)}
dphi=Deriv(phi,"t")
kk=function(t){t-phi(t)/dphi(t)}
Kg=Vectorize(kk)
Zs=sort(Z)
KK=Kg(sort(Z))
lines(Zs[which(KK!='NaN')],KK[which(KK!='NaN')],col="blue")


# GOODNESS-OF-FIT TEST
C_comp=CDF_f(Ps,Vs)
k=7    #no of parameters
MSE_cj=(1/(length(P)-k))*sum((C_obs-C_comp)^2)
RMSE_cj=sqrt(MSE_cj)
AIC_cj=length(P)*log(MSE_cj)+2*k

CDF_fuv= function(x,y) (x^(-theta)+y^(-theta)-1)^(-1/theta)
Cthetan=CDF_fuv(u,v)
Sn=sum((Cn-Cthetan)^2)
g=claytonCopula(theta,dim=2)
gofCopula(g, cbind(P,V), method = "Sn",simulation = "pb",estim.method="itau")



# RETURN PERIOD

TORcj=1/(1-CDF_PV_CJ)
Pseq2=seq(min(P), max(P), length.out=50)
Vseq2=seq(min(V), max(V), length.out=50)
TOR_cj=outer(Pseq2, Vseq2, function(x,y) 1/(1-CDF_f(x,y)))
#colors
z=TOR_cj
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq2,y=Vseq2,z=TOR_cj, xlab="P", ylab="V", zlab="Time (years)", col=fcol,box=TRUE,expand=0.5,ticktype = "detailed" ,xlim=c(min(Pseq2),max(Pseq2)), ylim=c(min(Vseq2),max(Vseq2)),zlim=c(min(TOR_cj),max(TOR_cj)), phi=20, theta=-40)

Pseq3=seq(4500, 13500, by=1500)
Vseq3=seq(20000,  1020000, length.out=60)
TOR_cj2=outer(Pseq3, Vseq3, function(x,y) 1/(1-CDF_f(x,y)))
for (i in 1:length(Pseq3))
{	
	plot(Vseq3,TOR_cj2[i,],cex.axis=1.05,col="darkmagenta", xaxt="n", yaxt="n",xlab="", ylab="",log="y",type="l",xlim=c(min(Vseq3),  max(Vseq3+170000)),ylim=c( min(TOR_cj2), max(TOR_cj2)))
	v1=Vseq3[60]
	t1=TOR_cj2[i,60]
	p2=Pseq3[i]
	text(v1,t1,cex=1.10, "P=",adj=c(-0.2,0.5),col="goldenrod")
	text(v1,t1,cex=1.10, Pseq3[i],adj=c(-0.6,0.5),col="goldenrod")
	par(new=TRUE)
}
T=c(1,2,5,10,20,50,100,500)
axis(2,cex=0.9,at=T,labels=T)
Vseq4=seq(0, 1000000,length=6)
axis(1,cex=0.9,at=Vseq4,labels=sprintf("%.0f",Vseq4))
title(main=list("OR-RETURN PERIOD",col='midnightblue'), xlab=list("V (cusec)",col='midnightblue'),ylab=list("Time (years)",col='midnightblue'),col='midnightblue')
par(new=FALSE)




# FRANK COPULA
# =============

# admissible range for tau is [-1, 1]

tau=cor.fk(P,V)
tau
f1=function(a) a/(exp(a)-1)
f=function(theta) 1-(4/theta)*(1-(1/theta)*integrate(f1,0,theta)$value)-tau
theta = uniroot(f, c(0.1, 100), tol = 0.0001)$root
BiCopPar2TailDep(5, theta)  #family number 5 for FRANK
u=cdfP
v=cdfV
CDF_PV_FR = function (P,V) -1/theta * log(1+ (exp(-distP(P)*theta)-1)*(exp(-distV(V)*theta)-1)/(exp(-theta)-1))

g=frankCopula(theta,dim=2)
r=rCopula(100000,g)
pbs=pobs(cbind(P,V))
plot(r, col="lightblue",pch=19, xlab="u", ylab="v", cex.lab=1.4)
points(pbs,pch=19,cex=0.7,col="blue")
r=cbind(qgev(r[,1], shape=-0.03659, scale=7074.7, loc=9137.0),qgev(r[,2],shape=0.10064, scale=205780, loc=205660))
par(mar=c(4,5,1,1))
plot(rr,col="lightblue",pch=19,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), cex.lab=1.2,yaxt="n")
axis(2, at=seq(0,7)*10^6, labels=sprintf("%d",seq(0,7)*10^6))
points(P,V,pch=19,cex=0.7,col="blue")

Pseq=seq(0, max(P), length.out=60)
Vseq=seq(0, max(V), length.out=60)
cdfPV3=outer(Pseq, Vseq, function(x,y) CDF_PV_FR(x,y))
#colors
z=cdfPV3
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq,y=Vseq,z=cdfPV3, xlab="P", ylab="V", zlab="cdf for Frank copula",col=fcol,box=TRUE,expand=0.5,ticktype = "detailed", zlim=c(0,1), phi=20, theta=-40)

# QQ PLOT FOR FRANK COPULA (KENDALL PLOT)
phi=function(z) -log((exp(-theta*z)-1)/(exp(-theta)-1))
dphi=Deriv(phi,"z")
K= function(z) {z- phi(z)/dphi(z)}

	#obtaining zi
z=0*P;
Ps=sort(P)
Vs=rep(0,length(P))
for (i in 1:length(P))
{	Vs[i]=V[which(P==Ps[i])]
}
for (i in 1:length(P))
{	z[i]=0;
	for (k in 1:length(P))
	{ 	if ((Ps[i]>Ps[k])&(Vs[i]>Vs[k])) z[i]=z[i]+1;
	}
	z[i]=z[i]/(length(P)-1);
}
Kn=0*P
for (i in 1:length(P))
{ 	Kn[i]= sum(z[i]>=z)/length(P)
}
Kp_fr=K(z)
plot(Kn,Kp_fr,pch=19, col="blue", xaxt="n", yaxt="n",xlab=expression(K[n]),ylab=expression(K[theta*n]),cex.lab=1.5,cex=1.2)
abline(0,1,lwd=2, col="red")
axis(2,cex=1.2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
axis(1,cex=1.2,at=seq(0,1,0.2),labels=seq(0,1,0.2))

		#another approach from https://www.r-bloggers.com/kendalls-function-for-copulas/
i=rep(1:length(P),each=length(P))
j=rep(1:length(P),length(P))
S=((P[i]>P[j])&(V[i]>V[j]))
Z=tapply(S,i,sum)/(length(P)-1)
#plot(ecdf(Z))
plot(sort(Z),(1:length(P))/length(P),type="s",col="red",xlab="",ylab="")
h=0.001
phi=function(z) -log((exp(-theta*z)-1)/(exp(-theta)-1))
dphi=Deriv(phi,"z")
kk=function(t){t-phi(t)/dphi(t)}
Kg=Vectorize(kk)
par(new=TRUE)
lines(sort(Z),Kg(sort(Z)),col="blue", lty=4, lwd=2)
par(new=FALSE)



# GOODNESS-OF-FIT TEST
C_comp=CDF_PV_FR(Ps,Vs)
k=7   #No of parameters
MSE_fr=(1/(length(P)-k))*sum((C_obs-C_comp)^2)
RMSE_fr=sqrt(MSE_fr)
AIC_fr=length(P)*log(MSE_fr)+2*k

CDF_fuv = function (x,y) -1/theta * log(1+ (exp(-x*theta)-1)*(exp(-y*theta)-1)/(exp(-theta)-1))

Cthetan=CDF_fuv(u,v)
Sn=sum((Cn-Cthetan)^2)

g=frankCopula(theta,dim=2)
gofCopula(g, cbind(P,V), method = "Sn",simulation = "pb",estim.method="itau")




# RETURN PERIOD
TOR_FR=1/(1- CDF_PV_FR)
Pseq2=seq(min(P), max(P)+150000, length.out=60)
Vseq2=seq(min(V), max(V)+250000, length.out=60)
TOR_FR=outer(Pseq2, Vseq2, function(x,y) 1/(1- Cthetan))
#colors
z=TOR_FR
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq2,y=Vseq2,z=TOR_FR, xlab="P", ylab="V", zlab="Time (years)",col=fcol,box=TRUE,expand=0.5,ticktype = "detailed",  zlim=c(0,100), phi=20, theta=-40)

Pseq3=seq(4500, 16500, by=2000)
Vseq3=seq(20000,  1020000, length.out=60)
TORFR=outer(Pseq3, Vseq3, function(x,y) 1/(1- Cthetan))
for (i in 1:length(Pseq3))
{	
	plot(Vseq3,TOR_amh2[i,],col="blue", xlab="V", ylab="Time (years)",log="y",type="l",xlim=c(min(Vseq3),  max(Vseq3)),ylim=c( min(TOR_amh2), max(TOR_amh2)))
	v1=Vseq3[57]
	t1=TOR_amh2[i,57]
	p2=Pseq3[i]
	text(v1,t1, "P=",adj=c(1,0.5),col="red")
	text(v1,t1, Pseq3[i],adj=c(0,0.5),col="red")
	par(new=TRUE)
}
par(new=FALSE)
medP=median(P)
medV=median(V)
mxP=max(P)
mxV=max(V)
TOR=function(x,y) 1/(1- Cthetan)
TAND=function(x,y) 1/(1-distP(x)-distV(y)+ Cthetan)
Tsec=function(x,y) 1/(1-K(CDF_f(x,y)))
Tc1=function(x,y) 1/(1- Cthetan/distP(x))
Tc2=function(x,y) (1/(1-distP(x)))*(1/(1-distP(x)-distV(y)+ Cthetan))
TORmed=TOR(medP,medV)
TORmax=TOR(mxP,mxV)
TANDmed=TAND(medP,medV)
TANDmax=TAND(mxP,mxV)
#Tsmed=Tsec(medP,medV)  #wrong method
#Tsmax=Tsec(mxP,mxV)
Tc1med=Tc1(medP,medV)
Tc1max=Tc1(mxP,mxV)
Tc2med=Tc2(medP,medV)
Tc2max=Tc2(mxP,mxV)
cbind(c(rep(c("TOR","TAND","Tc1","Tc2"),each=2)),round(c(TORmed,TORmax,TANDmed,TANDmax,Tc1med,Tc1max,Tc2med,Tc2max),3),rep(c("median","max"),4))


	

# Joe COPULA
# =============

# admissible range for tau is [0, 1]

tau=cor.fk(P,V)
#qphi=function(theta) {(-log(1-(1-t)^theta))/(-(theta * (1-t)^(theta - 1)/(1 - (1-t)^theta)))}
#f=function(theta) {1 + 4* integrate(qphi, 0 , 1)$value -tau}
#theta = uniroot(f, c(0.1, 100), tol = 0.0001)$root
theta= copJoe@iTau(tau)
BiCopPar2TailDep(6, theta)  #family number 5 for FRANK
u=cdfP
v=cdfV
CDF_PV_JO = function (P,V) 1-( (1-distP(P))^theta  +  (1-distV(V))^theta   -  (1-distP(P))^theta  *  (1-distV(V))^theta)^(1/theta)

g=joeCopula(theta,dim=2)
r=rCopula(200000,g)
pbs=pobs(cbind(P,V))
plot(r, col="lightblue",pch=19, xlab="u", ylab="v", cex.lab=1.4)
points(pbs,pch=19,cex=0.7,col="blue")
r=cbind(qgev(r[,1], shape=-0.03659, scale=7074.7, loc=9137.0),qgev(r[,2],shape=0.10064, scale=205780, loc=205660))
par(mar=c(4,5,1,1))
plot(rr,col="lightblue",pch=19,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), cex.lab=1.2,yaxt="n")
axis(2, at=seq(0,7)*10^6, labels=sprintf("%d",seq(0,7)*10^6))
points(P,V,pch=19,cex=0.7,col="blue")



Pseq=seq(0, max(P), length.out=60)
Vseq=seq(0, max(V), length.out=60)
cdfPV3=outer(Pseq, Vseq, function(x,y) CDF_PV_JO(x,y))
#colors
rz=cdfPV3
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq,y=Vseq,z=cdfPV3, xlab="P", ylab="V", zlab="cdf for Joe Copula",col=fcol,box=TRUE,expand=0.5,ticktype = "detailed", zlim=c(0,1), phi=20, theta=-40)

# QQ PLOT FOR FRANK COPULA (KENDALL PLOT)
phi=function(t) (-log(1-(1-t)^theta))
dphi=Deriv(phi,"t")
K= function(z) {z- phi(z)/dphi(z)}

	#obtaining zi
z=0*P;
Ps=sort(P)
Vs=rep(0,length(P))
for (i in 1:length(P))
{	Vs[i]=V[which(P==Ps[i])] 
}
for (i in 1:length(P))
{	z[i]=0;
	for (k in 1:length(P))
	{ 	if ((Ps[i]>Ps[k])&(Vs[i]>Vs[k])) z[i]=z[i]+1;
	}
	z[i]=z[i]/(length(P)-1);
}
Kn=0*P
for (i in 1:length(P))
{ 	Kn[i]= sum(z[i]>=z)/length(P)
}
Kp_fr=K(z)
plot(Kn,Kp_fr,pch=19, col="blue", xaxt="n", yaxt="n",xlab=expression(K[n]),ylab=expression(K[theta*n]),cex.lab=1.5,cex=1.2)
abline(0,1,lwd=2, col="red")
axis(2,cex=1.2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
axis(1,cex=1.2,at=seq(0,1,0.2),labels=seq(0,1,0.2))
r
		#another approach from https://www.r-bloggers.com/kendalls-function-for-copulas/
i=rep(1:length(P),each=length(P))
j=rep(1:length(P),length(P))
S=((P[i]>P[j])&(V[i]>V[j]))
Z=tapply(S,i,sum)/(length(P)-1)
#plot(ecdf(Z))
plot(sort(Z),(1:length(P))/length(P),type="s",col="red",xlab="",ylab="")
h=0.001
phi=function(t) (-log(1-(1-t)^theta))
dphi=Deriv(phi,"t")
k=function(t){t-phi(t)/dphi(t)}
Kg=Vectorize(k)
par(new=TRUE)
lines(sort(Z),Kg(sort(Z)),col="blue", lty=4, lwd=2)
par(new=FALSE)



# GOODNESS-OF-FIT TEST
C_comp=CDF_PV_JO(Ps,Vs)
k=7    #No of parameters 3*2 in univariate dists (GEV) and 1 for copula
MSE_jo=(1/(length(P)-k))*sum((C_obs-C_comp)^2)
RMSE_jo=sqrt(MSE_jo)
AIC_jo=length(P)*log(MSE_jo)+2*k
CDF_PV_JO_xy = function (x,y) 1-( (1-x)^theta  +  (1-y)^theta   -  (1-x)^theta  *  (1-y)^theta)^(1/theta)
Cthetan=CDF_PV_JO_xy (u,v)
Sn=sum((Cn-Cthetan)^2)
g=joeCopula(theta,dim=2)
gofCopula(g, cbind(P,V), method = "Sn",simulation = "pb",estim.method="itau")

# RETURN PERIOD
TORamh=1/(1-CDF_PV_AMH)
Pseq2=seq(min(P), max(P)+150000, length.out=60)
Vseq2=seq(min(V), max(V)+250000, length.out=60)
TOR_amh=outer(Pseq2, Vseq2, function(x,y) 1/(1-CDF_f(x,y)))
#colors
z=TOR_amh
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq2,y=Vseq2,z=TOR_amh, xlab="P", ylab="V", zlab="Time (years)",col=fcol,box=TRUE,expand=0.5,ticktype = "detailed",  zlim=c(0,100), phi=20, theta=-40)

Pseq3=seq(4500, 16500, by=2000)
Vseq3=seq(20000,  1020000, length.out=60)
TOR_amh2=outer(Pseq3, Vseq3, function(x,y) 1/(1-CDF_f(x,y)))
for (i in 1:length(Pseq3))
{	
	plot(Vseq3,TOR_amh2[i,],col="blue", xlab="V", ylab="Time (years)",log="y",type="l",xlim=c(min(Vseq3),  max(Vseq3)),ylim=c( min(TOR_amh2), max(TOR_amh2)))
	v1=Vseq3[57]
	t1=TOR_amh2[i,57]
	p2=Pseq3[i]
	text(v1,t1, "P=",adj=c(1,0.5),col="red")
	text(v1,t1, Pseq3[i],adj=c(0,0.5),col="red")
	par(new=TRUE)
}
par(new=FALSE)	

### COMPARISON

bind(c("Copula","AIC","RMSE"),cbind(c("Clayton","Frank","Gumbel","Joe"),c(AIC_cj,AIC_fr,AIC_gh,AIC_jo),c(RMSE_cj,RMSE_fr,RMSE_gh,RMSE_jo)))


#for gumbel

theta=copGumbel@iTau(tau)
CDF_fuv= function(x,y) (x^(-theta)+y^(-theta)-1)^(-1/theta)
useq=seq(0,1, length.out=2000)
vseq=seq(0,1, length.out=2000)
u=pobs(P); v=pobs(V)
Cnfun=function (x,y) C.n(cbind(x,y), cbind(P,V))
CDF_emp=outer (useq,vseq,Cnfun)
CDF_th=outer (useq,vseq,CDF_fuv)
contour(x=useq,y=vseq,z=CDF_th, zlim=c(0,max(z)),cex.lab=1.5,xlab="u",ylab="v", levels=seq(0.1,0.9,by=0.1),col="red",lty=3,lwd=3,drawlabels=FALSE,cex.axis=1.2)
contour(x=useq,y=vseq,z=CDF_emp, zlim=c(0,max(z)),cex.lab=1.2,xlab="",ylab="", levels=seq(0.1,0.9,by=0.1),add=TRUE,xaxt="n",yaxt="n")


# for CLayton
theta=copClayton@iTau(tau)
CDF_fuv= function(x,y) (x^(-theta)+y^(-theta)-1)^(-1/theta)
CDF_th=outer (useq,vseq,CDF_fuv)
contour(x=useq,y=vseq,z=CDF_th, zlim=c(0,max(z)),cex.lab=1.5,xlab="u",ylab="v", levels=seq(0.1,0.9,by=0.1),col="red",lty=3,lwd=3,drawlabels=FALSE,cex.axis=1.2)
contour(x=useq,y=vseq,z=CDF_emp, zlim=c(0,max(z)),cex.lab=1.2,xlab="",ylab="", levels=seq(0.1,0.9,by=0.1), Rp)


# for frank
theta=copFrank@iTau(tau)
CDF_fuv = function (x,y) -1/theta * log(1+ (exp(-x*theta)-1)*(exp(-y*theta)-1)/(exp(-theta)-1))
CDF_th=outer (useq,vseq,CDF_fuv)
contour(x=useq,y=vseq,z=CDF_th, zlim=c(0,max(z)),cex.lab=1.5,xlab="u",ylab="v", levels=seq(0.1,0.9,by=0.1),col="red",lty=3,lwd=3,drawlabels=FALSE,cex.axis=1.2)
contour(x=useq,y=vseq,z=CDF_emp, zlim=c(0,max(z)),cex.lab=1.2,xlab="",ylab="", levels=seq(0.1,0.9,by=0.1),add=TRUE,xaxt="n",yaxt="n")

# for joe
theta=copJoe@iTau(tau)
CDF_fuv = function (x,y) 1-( (1-x)^theta  +  (1-y)^theta   -  (1-x)^theta  *  (1-y)^theta)^(1/theta)
CDF_fuv = function (x,y) -1/theta * log(1+ (exp(-x*theta)-1)*(exp(-y*theta)-1)/(exp(-theta)-1))
CDF_th=outer (useq,vseq,CDF_fuv)
contour(x=useq,y=vseq,z=CDF_th, zlim=c(0,max(z)),cex.lab=1.5,xlab="u",ylab="v", levels=seq(0.1,0.9,by=0.1),col="red",lty=3,lwd=3,drawlabels=FALSE,cex.axis=1.2)
contour(x=useq,y=vseq,z=CDF_emp, zlim=c(0,max(z)),cex.lab=1.2,xlab="",ylab="", levels=seq(0.1,0.9,by=0.1),add=TRUE,xaxt="n",yaxt="n")







#####################################################################
################      COPULA ON MATLAB DATA
#####################################################################

#data=read.csv("E:\\KU\\Flood and Copula\\SN\\Terbela_MATLAB.csv")

Pmat=read.csv("c:\\user\\Saba\\desktop\\saba\\Flood and Copula\\Terbela_MATLAB.csv")[,4]
Vmat=read.csv("c:\\user\\Saba\\Desktop\\saba\\Flood and Copula\\Terbela_MATLAB.csv")[,5]

#Independence test
m=pobs(cbind(Pmat,Vmat))  # Pseudo observations are the observations in the [0,1] interval.
empsamp=indepTestSim(length(Pmat),p=2,N=10000)
indepTest(m,empsamp) # fail to accept the null hyp that P and V are independent

# UNIVARIATE GUMBEL
dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q,a,b) exp(-exp((a-q)/b))
qgumbel <- function(p,a,b) a-b*log(-log(p))

scaP=sqrt(6)*sd(Pmat)/pi
locP= mean(Pmat)-0.5777*scaP

gumP<-fitdist(Pmat,"gumbel",start=list(a=mean(Pmat),b=sd(Pmat)))
bP=as.numeric(gumP$estimate[2])
aP=as.numeric(gumP$estimate[1])
cdfcomp(gumP)
denscomp(gumP)
ppcomp(gumP)
qqcomp(gumP)
cdfP=pgumbel(Pmat,aP,bP)

scaV=sqrt(6)*sd(Vmat)/pi
locV= mean(Vmat)-0.5777*scaV

gumV <- fitdist(Vmat,"gumbel",start=list(a=mean(Vmat),b=sd(Vmat)))
summary(gumV)
aV=as.numeric(gumV$estimate[1])
bV=as.numeric(gumV$estimate[2])
cdfcomp(gumV)
denscomp(gumV)
ppcomp(gumV)
qqcomp(gumV)
cdfV=pgumbel(Vmat,aV,bV)
plot(Vmat,cdfV)

# BIVARIATE GUMBEL
# ================


rho=cor(Pmat,Vmat)
theta=2*(1-cos(pi*sqrt(rho/6)))
CDF_PV = cdfP * cdfV * exp(-theta*(1/log(cdfP) + 1/log(cdfV))^(-1))

Pseq=seq(min(Pmat), max(Pmat)+min(Pmat), length.out=50)
Vseq=seq(min(Vmat), max(Vmat)+min(Vmat), length.out=50)
cdfPV1=outer(Pseq, Vseq, function(x,y) pgumbel(x, aP, bP) * pgumbel(y, aV, bV) * exp(-theta*(1/log(pgumbel(x, aP, bP)) + 1/log(pgumbel(y, aV, bV)))^(-1)))
#colors
z=cdfPV1
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,1)
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq,y=Vseq,z=cdfPV1, col=fcol,box=TRUE,ticktype = "detailed",expand=0.5, zlim=c(0,1), phi=20, theta=-40)

##      COPULA 
#Independence test
m=pobs(cbind(Pmat,Vmat))  # Pseudo observations are the observations in the [0,1] interval.
empsamp=indepTestSim(length(Pmat),p=2,N=10000)
indepTest(m,empsamp) # fail to accept the null hyp that P and V are independent

# GUMBEL COPULA
# =============

tau=cor.fk(Pmat,Vmat)
theta=1/(1-tau)
u=cdfP
v=cdfV
CDF_PVgum = exp(-((-log(u))^theta+(-log(v))^theta)^(1/theta))

Pseq=seq(min(Pmat)-1000, max(Pmat)+1000, length.out=50)
Vseq=seq(min(Vmat)-3000, max(Vmat)+3000, length.out=50)
cdfPV2=outer(Pseq, Vseq, function(x,y) exp(-((-log(pgumbel(x, aP, bP)))^theta+(-log(pgumbel(y, aV, bV)))^theta)^(1/theta)))
#colors
z=cdfPV2
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,1)
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq,y=Vseq,z=cdfPV2,col=fcol,box=TRUE,ticktype = "detailed",expand=0.5, zlim=c(0,1), phi=20, theta=-40)

TORgum=1/(1-CDF_PVgum)
Pseq2=seq(min(Pmat), max(Pmat)+5*min(Pmat), length.out=60)
Vseq2=seq(min(Vmat), max(Vmat)+5*min(Vmat), length.out=60)
TOR_gum=outer(Pseq2, Vseq2, function(x,y) 1/(1-(exp(-((-log(pgumbel(x, aP, bP)))^theta+(-log(pgumbel(y, aV, bV)))^theta)^(1/theta)))))
#colors
z=TOR_gum
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq2,y=Vseq2,z=TOR_gum,col=fcol,box=TRUE,expand=0.5,ticktype = "detailed", zlim=c(0,max(z)), phi=20, theta=-40)


# AMH COPULA
# =============

tau=cor.fk(Pmat,Vmat)
f=function(theta) ((3*theta-2)/theta)-(2/3 * (1-1/theta)^2)*log(1-theta)-tau
theta = uniroot(f, c(0.2, 0.8), tol = 0.0001)$root
u=cdfP
v=cdfV
CDF_PV_AMH = u*v/(1-theta*(1-u)*(1-v))

Pseq=seq(min(Pmat), max(Pmat), length.out=60)
Vseq=seq(min(Vmat), max(Vmat), length.out=60)
cdfPV3=outer(Pseq, Vseq, function(x,y) pgumbel(x, aP, bP)*pgumbel(y, aV, bV)/(1-theta*(1-pgumbel(x, aP, bP))*(1-pgumbel(y, aV, bV))))
#colors
z=cdfPV3
nrz <- nrow(z)
rncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq,y=Vseq,z=cdfPV3,col=fcol,box=TRUE,expand=0.5,ticktype = "detailed", zlim=c(0,max(z)), phi=20, theta=-40)

TORamh=1/(1-CDF_PV_AMH)
Pseq2=seq(min(Pmat), max(Pmat)+150000, length.out=60)
Vseq2=seq(min(Vmat), max(Vmat)+250000, length.out=60)
TOR_amh=outer(Pseq2, Vseq2, function(x,y) 1/(1-(pgumbel(x, aP, bP)*pgumbel(y, aV, bV)/(1-theta*(1-pgumbel(x, aP, bP))*(1-pgumbel(y, aV, bV))))))
#colors
z=TOR_amh
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
rzlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq2,y=Vseq2,z=TOR_amh,col=fcol,box=TRUE,expand=0.5,ticktype = "detailed",  zlim=c(0,max(z)), phi=20, theta=-40)


# CookJohnson COPULA (a.k.a. Clayton???)
# ==================

tau=cor.fk(Pmat,Vmat)
theta = (2*tau)/(1-tau)
u=cdfP
v=cdfV
CDF_PV_CJ = (u^(-theta)+v^(-theta)-1)^(-1/theta)

Pseq=seq(min(Pmat), max(Pmat), length.out=50)
Vseq=seq(min(Vmat), max(Vmat), length.out=50)
cdfPV_CJ = outer(Pseq, Vseq, function(x,y) ((pgumbel(x, aP, bP))^(-theta)+(pgumbel(y, aV, bV))^(-theta)-1)^(-1/theta))
#colors
z=cdfPV_CJ
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq,y=Vseq,z=cdfPV_CJ,col=fcol,box=TRUE,expand=0.5,ticktype = "detailed" ,zlim=c(0,1), phi=20, theta=-40)

TORcj=1/(1-CDF_PV_CJ)
Pseq2=seq(min(Pmat), max(Pmat), length.out=50)
Vseq2=seq(min(Vmat), max(Vmat), length.out=50)
TOR_cj=outer(Pseq2, Vseq2, function(x,y) 1/(1-((pgumbel(x, aP, bP))^(-theta)+(pgumbel(y, aV, bV))^(-theta)-1)^(-1/theta)))
#colors
z=TOR_cj
nrz <- nrow(z)
ncz <- ncol(z)
ncol=100
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
"cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
couleurs <- tail(jet.colors(1.2*ncol),ncol)
zlim=c(0,max(z))
fcol <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]
#colors
persp(x=Pseq2,y=Vseq2,z=TOR_cj, col=fcol,box=TRUE,expand=0.5,ticktype = "detailed" ,xlim=c(min(Pseq2),max(Pseq2)), ylim=c(min(Vseq2),max(Vseq2)),zlim=c(min(TOR_cj),max(TOR_cj)), phi=20, theta=-40)




#GUMBEL COPULA BY MAXIMUM LIKELIHOOD FUNCTION
model <- gumbelCopula(dim = 2)
mymvd = mvdc(model, margins = c("gumbel","gumbel"),paramMargins=list(list(loc=mean(Pmat),scale=sd(Pmat)),list(loc=mean(Vmat),scale=sd(Vmat))))
fitmvd=fitMvdc(m, mymvd,c(0,1,0,1,theta))
print(fitmvd)

# GUMBEL COPULA BY KENDALL METHOD

tau=cor.fk(P,V)
theta=1/(1-tau)
CDF_PV_kendall = (exp(-((-log(cdfP)))^theta +(-log(cdfV))^theta)^(1/theta))
#function(x,y) exp(-((-log(distP(x)))^theta+(-log(distV(y)))^theta)^(1/theta))

plot(criticalLevel(getKendallDistr(g), KRP = seq(1,1000,by=10), mu = 1, g,seq(1,1000,by=10), type="l"))
surCOP(1-u,1-v, cop=GHcop,par=theta)
Rp=1/surCOP(1-u,1-v, cop=GHcop,par=theta,dim=2)
plot(Pseq,Vseq,Rp=seq(1,1000,by=10),mu=1,type="l")
 surfuncCOP(1-u,1-v, cop=GHcop,par=theta,cor=tau)
Pseq2=seq(min(P), max(P)+5*min(P), length.out=60)
Vseq2=seq(2.5*10^5, max(V)+8*min(V), length.out=60)
#Pseq2=seq(0, max(P)+65*min(P), length.out=60)
#Vseq2=seq(0, max(V)+65*min(V), length.out=60)
contour(x=Pseq2,y=Vseq2,z=(1-1/Rp),xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),lwd=2,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), levels=c(0.1),col="red",labcex=.9)
 contour(x=Pseq,y=Vseq,z=Rp,zlim=c(0,max(z)),cex.lab=1.5,xlab="u",ylab="v", levels=seq(0.1,0.9,by=0.1),col="red",lty=3,lwd=3,drawlabels=FALSE,cex.axis=1.2)
contour(x=Pseq2,y=Vseq2,z=(1-1/Rp),xlim=c(0,max(Pseq2)),ylim=c(0,max(Vseq2)), zlim=c(0,max(z)),lwd=2,xlab=expression(P (m^3/s)),ylab=expression(V (day-m^3/s)), levels=c(0.99),col="red",labcex=.9)
Rp
max(z)

plot(criticalLevel(getKendallDistr(gum_hug), KRP = seq(1,1000,by=10), mu = 1, gum_hug,seq(1,1000,by=10), type="l", lty=1))
kendallRP(copula=g,cl=c(0.90,0.99,0.999,0.9999))
plot(seq(0.01,1,by=0.01),kendallRP(copula=g,cl=seq(0.01,1,by=0.01)),type="l",lty=1)
