#-------------------------------------------------------------------------
# Simulation of equations of terminal investment models.
#-------------------------------------------------------------------------

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
# Everything below is from September 2014 and includes different hypothetical species life history curves. These life history curves focus on different survival probabilities over time, and different fecundities over time. 

#In constructing this curves, I am deliberately avoiding creating situations in which expected lifetime fecundity is divergent with age -- In other words, I don't want a situation in which a species increases its survival probability to one and then never drops, because when integrating from zero to infinite age, offspring production also becomes ``infinite'' (i.e., the expected number of offspring produced at age infinity would be infinite -- more properly, the series is divergent). Instead, survival probability cannot terminate to equal one for sufficiently old individuals -- this allows us to not have an arbitrary terminal age, and to not have to have unrealistic expected fecundities because cumulative survival probability is always decreasing.
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX



# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
# Five different types of functions will be considered: Note that the FUNCTION.Cu ones look weird, but what this does is return the cumulative probability of surviving to age x given any entered ages (i.e., you can enter a vector of ages, or just one age, and get the cumulative survival probaiblity of each). It does this by first expanding so that the cumulative probability is calculated for 1 to the maximum age entered, but then returns only the ages asked for.
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX

# XXX TODO: Scale terminal investment. Scale maximum fecundity to one if you live to 20. Want terminal investment to be a percentage of age-specific fecundity. E.g., double the amount of offspring produced at age x to terminally invest. Will look at zero to 100% (essentially double). Can do the same for RV0 -- go from zero to 100%.

#-------------------------------------------------------------------------
# Plotting program -- makes a 3d plot from information calculated later
#-------------------------------------------------------------------------

rm(list=ls());

library(scatterplot3d);

TIplot <- function(dat,xax,yax){
	if(sum(dat)>0){
		Titable <- NULL;
		for(i in 1:dim(dat)[1]){
			for(j in 1:dim(dat)[2]){
				if(sum(dat[i,j,] > 0)){
					ages    <- which(dat[i,j,] > 0);	
					Fecu    <- rep(x=xax[i],times=length(ages));
					RRvO    <- rep(x=yax[j],times=length(ages));
					add     <- cbind(Fecu,RRvO,ages);
					Titable <- rbind(Titable,add); 
				}
			}
		}
		scatterplot3d(x=Titable[,1],y=Titable[,2],z=Titable[,3],
			pch=16,cex.symbols=1.0,zlim=c(0,20),xlim=c(xax[1],xax[length(xax)]),
            ylim=c(yax[1],yax[length(yax)]),
			xlab=expression(paste("TI in fecundity (",Delta,m[x],")")),
			ylab=expression(paste("TI in offspring reproductive value (",Delta,RV[0],")")),
			zlab="Age",color="black",y.margin.add=0.5,cex.lab=1.25);
	}else{
		scatterplot3d(x=0,y=0,z=0,type="n",xlim=c(0,1),ylim=c(0,1),
			pch=16,cex.symbols=1.0,zlim=c(0,20),
			xlab=expression(paste("TI in fecundity (",Delta,m[x],")")),
			ylab=expression(paste("TI in offspring reproductive value (",Delta,RV[0],")")),
			zlab="Age",color="black",y.margin.add=0.25,cex.lab=1.25,cex.axis=1.25);
	}
}


# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
# XXX                                                                 XXX
# XXX    Below 1 through 5 sets fecundity and survival curves         XXX
# XXX                                                                 XXX
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX


#-------------------------------------------------------------------------
# 1) Constant: Fecundity = {2,8}, Survival probability = {0.2,0.8}; pretty straightforward.
#-------------------------------------------------------------------------
# Low fecundities and survival probabilities:
constantF.Pr.lo <- function(age) rep(2,length(age)); # Fecundity at age
constantP.Pr.lo <- function(age) rep(0.2,length(age)); # Probability of survival age

constantP.Cu.lo <- function(age){ # Cumulative probability of surviving to age.
	all    <- 1:max(age);
	ret    <- rep(0,length(all));
	tot    <- 1; # Probability of being at age zero (must be one, starting point).
	for(i in 1:length(all)){
		ret[i] <- tot * constantP.Pr.lo(i);
		tot    <- ret[i];
	}
	return(ret[age]);
}

# High fecundities and survival probabilities:
constantF.Pr.hi <- function(age) rep(8,length(age)); # Fecundity at age
constantP.Pr.hi <- function(age) rep(0.8,length(age)); # Probability of survival age

constantP.Cu.hi <- function(age){ # Cumulative probability of surviving to age.
	all    <- 1:max(age);
	ret    <- rep(0,length(all));
	tot    <- 1; # Probability of being at age zero (must be one, starting point).
	for(i in 1:length(all)){
		ret[i] <- tot * constantP.Pr.hi(i);
		tot    <- ret[i];
	}
	return(ret[age]);
}
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# 2) Exponential slope up: These slopes increase, but not linearly with age (see above for the reasoning behind this -- don't want to accidentally create a divergent series). Instead, both fecundity and survival probability increase with age, but with slightly diminishing returns so that at high age classes, the fecundity is higher next year, as is survival probability, but the increase is not as good as in the previous year. Survival probability plateaus at ca 0.692 and 0.99 for low and high probability, respectively.
#-------------------------------------------------------------------------
# Low:
exponeupF.Pr.lo <- function(age){
	ret <- rep(0,length(age));
	tot <- 0;	
	for(i in 1:length(age)){
		ret[i] <- 1 * sum(1/((5/4)^((1:i)-1))) - 1;
	}
	return(ret);
}

exponeupP.Pr.lo <- function(age) {
	ret <- rep(0,length(age));
	tot <- 0;
	for(i in 1:length(age)){
		ret[i] <- (198/1000) * sum(1/((4/3)^((1:age[i])-1))) - 0.1;
	}
	return(ret);
}

exponeupP.Cu.lo <- function(age){
	all <- 1:max(age);
	ret <- rep(0,length(all));
	tot <- 1;	
	for(i in 1:length(all)){
		ret[i] <- tot * exponeupP.Pr.lo(i);
		tot    <- ret[i];
	}
	return(ret[age]);
}

# High:
exponeupF.Pr.hi <- function(age){
	ret <- rep(0,length(age));
	tot <- 0;	
	for(i in 1:length(age)){
		ret[i] <- 2 * sum(1/((10/9)^((1:i)-1)))
	}
	return(ret);
}

exponeupP.Pr.hi <- function(age) {
	ret <- rep(0,length(age));
	tot <- 0;
	for(i in 1:length(age)){
		ret[i] <- (198/1000) * sum(1/((5/4)^((1:age[i])-1)))
	}
	return(ret);
}

exponeupP.Cu.hi <- function(age){
	all <- 1:max(age);
	ret <- rep(0,length(all));
	tot <- 1;	
	for(i in 1:length(all)){
		ret[i] <- tot * exponeupP.Pr.hi(i);
		tot    <- ret[i];
	}
	return(ret[age]);
}
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# 3) Exponential slopes down: These slopes decrease exponentially (see above for reasoning). This allows for one continuous function in which fecundity and survival probability decreases constantly, asymptoting toward zero.
#-------------------------------------------------------------------------
# Low: 
exponednF.Pr.lo <- function(age){
	ret <- rep(0,length(age));
	for(i in 1:length(age)){
		ret[i] <- (15) * exp(-(10/40)*age[i]);
	}
	return(ret);
}

exponednP.Pr.lo <- function(age){
	ret <- rep(0,length(age));
	for(i in 1:length(age)){
		ret[i] <- (0.5) * exp(-(10/40)*age[i]);
	}
	return(ret);
}

exponednP.Cu.lo <- function(age){
	all <- 1:max(age);
	ret <- rep(0,length(all));
	tot <- 1;	
	for(i in 1:length(all)){
		ret[i] <- tot * exponednP.Pr.lo(i);
		tot    <- ret[i];
	}
	return(ret[age]);
}

# High:
exponednF.Pr.hi <- function(age){
	ret <- rep(0,length(age));
	for(i in 1:length(age)){
		ret[i] <- (30) * exp(-(5/40)*age[i]);
	}
	return(ret);
}

exponednP.Pr.hi <- function(age){
	ret <- rep(0,length(age));
	for(i in 1:length(age)){
		ret[i] <- (1) * exp(-(5/40)*age[i]);
	}
	return(ret);
}

exponednP.Cu.hi <- function(age){
	all <- 1:max(age);
	ret <- rep(0,length(all));
	tot <- 1;	
	for(i in 1:length(all)){
		ret[i] <- tot * exponednP.Pr.hi(i);
		tot    <- ret[i];
	}
	return(ret[age]);
}

#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# 4) Humped shaped distribution: Fecundity = 5 * (1/sqrt(2*pi))e^-(age - 10)^2. Survival probability = (1/sqrt(2*pi))e^-(age - 10)^2;
#-------------------------------------------------------------------------
# Low:
humpedshF.Pr.lo <- function(age,sd=5) 10 * (1/sd*sqrt(2*pi))*exp(-1*((age - 10)^2/(2*sd*sd)));
humpedshP.Pr.lo <- function(age,sd=5) (1/sd*sqrt(2*pi))*exp(-1*((age - 10)^2/(2*sd*sd)));

humpedshP.Cu.lo <- function(age,sd=5){
	all <- 1:max(age);
	ret <- rep(0,length(all));
	tot <- 1;	
	for(i in 1:length(all)){
		ret[i] <- tot * humpedshP.Pr.lo(i,sd);
		tot    <- ret[i];
	}
	return(ret[age]);
}

# High:
humpedshF.Pr.hi <- function(age,sd=5) 20 * (1/sd*sqrt(2*pi))*exp(-1*((age - 10)^2/(2*sd*sd)));
humpedshP.Pr.hi <- function(age,sd=5) 2 * (1/sd*sqrt(2*pi))*exp(-1*((age - 10)^2/(2*sd*sd)));

humpedshP.Cu.hi <- function(age,sd=5){
	all <- 1:max(age);
	ret <- rep(0,length(all));
	tot <- 1;	
	for(i in 1:length(all)){
		ret[i] <- tot * humpedshP.Pr.hi(i,sd);
		tot    <- ret[i];
	}
	return(ret[age]);
}
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# 5) U-shaped distribution: Fecundity = 5 * (1/sqrt(2*pi))e^-(age - 10)^2. Survival probability = (1/sqrt(2*pi))e^-(age - 10)^2;
#-------------------------------------------------------------------------
# Low
UshF.Pr.lo <- function(age,sd=5) 10 - 18*(1/sd*sqrt(2*pi))*exp(-1*((age - 10)^2/(2*sd*sd)));
UshP.Pr.lo <- function(age,sd=5) 0.6 - (1/sd*sqrt(2*pi))*exp(-1*((age - 10)^2/(2*sd*sd)));

UshP.Cu.lo <- function(age,sd=5){
	all <- 1:max(age);
	ret <- rep(0,length(all));
	tot <- 1;	
	for(i in 1:length(all)){
		ret[i] <- tot * UshP.Pr.lo(i,sd);
		tot    <- ret[i];
	}
	return(ret[age]);
}

# High:
UshF.Pr.hi <- function(age,sd=5) 15 - 18*(1/sd*sqrt(2*pi))*exp(-1*((age - 10)^2/(2*sd*sd)));
UshP.Pr.hi <- function(age,sd=5) 1 - (1/sd*sqrt(2*pi))*exp(-1*((age - 10)^2/(2*sd*sd)));

UshP.Cu.hi <- function(age,sd=5){
	all <- 1:max(age);
	ret <- rep(0,length(all));
	tot <- 1;	
	for(i in 1:length(all)){
		ret[i] <- tot * UshP.Pr.hi(i,sd);
		tot    <- ret[i];
	}
	return(ret[age]);
}


# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
# XXX                                                                 XXX
# XXX    Plots showing these different curves for ages from 1 to 20   XXX
# XXX                                                                 XXX
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX


x <- 1:20; # This is the number of years included, can increase beyond 20.

# FECUNDITIES:
#-------------------------------------------------------------------------
# Constant fecundity with age
# Low:
y <- constantF.Pr.lo(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
# High:
y <- constantF.Pr.hi(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
#-------------------------------------------------------------------------
# Exponential up fecundity with age
# Low:
y <- exponeupF.Pr.lo(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
# High:
y <- exponeupF.Pr.hi(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
#-------------------------------------------------------------------------
# Exponential down fecundity with age
# Low:
y <- exponednF.Pr.lo(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
# High:
y <- exponednF.Pr.hi(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
#-------------------------------------------------------------------------
# Humped-shaped fecundity with age
# Low:
y <- humpedshF.Pr.lo(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
# High:
y <- humpedshF.Pr.hi(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
#-------------------------------------------------------------------------
# U-shaped fecundity with age
# Low:
y <- UshF.Pr.lo(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
# High:
y <- UshF.Pr.hi(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
#-------------------------------------------------------------------------


# SURVIVAL PROBABILITIES:
#-------------------------------------------------------------------------
# Constant survival with age
# Low:
y <- constantP.Pr.lo(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
# High:
y <- constantP.Pr.hi(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
#-------------------------------------------------------------------------
# Exponential up survival with age
# Low:
y <- exponeupP.Pr.lo(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
# High:
y <- exponeupP.Pr.hi(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
#-------------------------------------------------------------------------
# Exponential down survival with age
# Low:
y <- exponednP.Pr.lo(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
# High:
y <- exponednP.Pr.hi(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
#-------------------------------------------------------------------------
# Humped-shaped survival with age
# Low:
y <- humpedshP.Pr.lo(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
# High:
y <- humpedshP.Pr.hi(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
#-------------------------------------------------------------------------
# U-shaped survival with age
# Low:
y <- UshP.Pr.lo(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
# High:
y <- UshP.Pr.hi(x);
par(mar=c(5,5,1,1),lwd=2);
plot(x=x,y=y,type="l",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),cex.lab=2,cex.axis=1.5);
polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80");
#-------------------------------------------------------------------------


# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
# XXX                                                                 XXX
# XXX    Comparing combinations of survival and fecundity             XXX
# XXX                                                                 XXX
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
# Here we compare different combinations of survival probability distributions and lifetime fecundity distributions to arrive at terminal investment values.
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX

# -----------------------------------------------------------------------
# Function to predict terminal investment (left-hand side of 13):
lhs <- function(Dmx, DRV0, mx, RV0, RVx){
	ti <- Dmx + (DRV0 * (mx + Dmx)) / (2 + RV0);
	if(ti > RVx){
		return(1);
	}else{
		return(0);
	}
}
# -----------------------------------------------------------------------

FecFun <- constantF.Pr.lo;
SurFun <- constantP.Cu.lo;

# -----------------------------------------------------------------------
# Function plot the space of terminal investment (requires lhs() and TIplot() [above]).
# yr: How many years to consider.
# yF: Shows maximum plotted (percentage) for increase in fecundity (default 100).
# yR: Shows maximum plotted for increase in offspring reproductive value.
TIspace <- function(FecFun,SurFun,yr=20,yF=1,yR=1){
	resvals <- rep(x=0, times=yr+1);
	termi   <- 10000;
	for(i in 1:(yr+1)){
		resvals[i] <- sum((SurFun(i:termi)/SurFun(i))*FecFun(i:termi));
	}
	Fpar  <- seq(from=0, to=yF, length.out=100);
	Rpar  <- seq(from=0, to=yR, length.out=100);
	Space <- array(data=0,dim=c(100,100,yr));
	for(age in 1:dim(Space)[3]){
		for(i in 1:dim(Space)[1]){
			for(j in 1:dim(Space)[2]){
				CHECK <- lhs( Dmx  =  Fpar[i] * FecFun(i),   # TI in fecundity 
						      DRV0 =  Rpar[j],                   # TI in RV of offspring
						      mx   =  FecFun(age),               # Age fecundity
						      RV0  =  resvals[1],                # Offspring RV 
						      RVx  =  resvals[age+1]             # RV at age x
						);			
				Space[i,j,age] <- CHECK;
			}
		}
		print(age);
	}
	TIplot(dat=Space,xax=Fpar,yax=Rpar);
	return(Space);
}
# ----------------------------------------------------------------------------

#       ----------------------------------------------------------------------
# 		LIST OF FECUNDITY FUNCTIONS:   	\	LIST OF SURVIVAL FUNCTIONS:
#		----------------------------	\	----------------------------------
# 		constantF.Pr.lo					\	constantP.Cu.lo	
# 		constantF.Pr.hi					\	constantP.Cu.hi
# 		exponeupF.Pr.lo					\	exponeupP.Cu.lo	
# 		exponeupF.Pr.hi					\	exponeupP.Cu.hi	
#		exponednF.Pr.lo					\	exponednP.Cu.lo	
#		exponednF.Pr.hi					\	exponednP.Cu.hi	
#		humpedshF.Pr.lo					\	humpedshP.Cu.lo	
#		humpedshF.Pr.hi					\	humpedshP.Cu.hi
#		UshF.Pr.lo						\	UshP.Cu.lo		
#		UshF.Pr.hi						\	UshP.Cu.hi		
#       ----------------------------------------------------------------------

Space <- TIspace(humpedshF.Pr.lo, SurFun=constantP.Cu.lo);

# ----------------------------------------------------------------------------


SP <- Space;
for(age in 2:(dim(Space)[3]-1)){
	for(i in 2:(dim(Space)[1]-1)){
		for(j in 2:(dim(Space)[2]-1)){
			if(Space[i,j,age]==1){
				if( Space[i-1,j,age]==1 && 
					Space[i+1,j,age]==1	&&
					Space[i,j+1,age]==1	&&
					Space[i,j-1,age]==1	&&
					Space[i,j,age+1]==1	&&
					Space[i,j,age-1]==1
				){
					SP[i,j,age] <- 0;
				}else{
					SP[i,j,age] <- age;
				}
			}
		}
	}
}

HI <- matrix(data=0,nrow=dim(SP)[1],ncol=dim(SP)[2]);
LO <- matrix(data=0,nrow=dim(SP)[1],ncol=dim(SP)[2]);
for(i in 1:(dim(SP)[1])){
	for(j in 1:(dim(SP)[2])){
		if(sum(SP[i,j,]>0)==0){
			break;
		}else{
			dat     <- SP[i,j,];
			dat     <- dat[dat>0];
			HI[i,j] <- max(SP[i,j,]);
			LO[i,j] <- min(SP[i,j,]);
		}
	}
}


persp(x = HI, theta = 45, phi = 30, r = 50,xlim=c(0,1),ylim=c(0,1),zlim=c(0,20));

# 
# library(plot3D);
# library(animation); 
# 
# 
# 
# surf3D(
# 
# 
# 
# xx <- seq(from=0, to=2, length.out=100);
# yy <- seq(from=0, to=1, length.out=100);
# 
# 
# TIplot(dat=SP,xax=xx,yax=yy);
# 
# 
# isosurf3Drgl(x=HI, colkey = FALSE, shade = 0.5,
#         box = FALSE, theta = 60)


#-------------------------------------------------------------------------
# BD code from the summary, and for making figures
#-------------------------------------------------------------------------

lhs <- function(Dmx, DRV0, mx, RV0, RVx){
    ti <- Dmx + (DRV0 * (mx + Dmx)) / (2 + RV0);
    if(ti > RVx){
        return(1);
    }else{
        return(0);
    }
}

TIspace <- function(FecFun,SurFun,yr=20,yF=1,yR=1){
    resvals <- rep(x=0, times=yr+1);
    termi   <- 1000;
    for(i in 1:(yr+1)){
        resvals[i] <- sum((SurFun(i:termi)/SurFun(i))*FecFun(i:termi));
    }
    Fpar  <- seq(from=0, to=yF, length.out=100);
    Rpar  <- seq(from=0, to=yR, length.out=100);
    Space <- array(data=0,dim=c(100,100,yr));
    for(age in 1:dim(Space)[3]){
        for(i in 1:dim(Space)[1]){
            for(j in 1:dim(Space)[2]){
                CHECK <- lhs( Dmx  =  Fpar[i] * FecFun(i),   # TI in fecundity 
                              DRV0 =  Rpar[j],                   # TI in RV of offspring
                              mx   =  FecFun(age),               # Age fecundity
                              RV0  =  resvals[1],                # Offspring RV 
                              RVx  =  resvals[age+1]             # RV at age x
                );			
                Space[i,j,age] <- CHECK;
            }
        }
    }
    return(Space);
}

Ffun  <- list(constantF.Pr.lo, constantF.Pr.hi, exponeupF.Pr.lo, exponeupF.Pr.hi,
              exponednF.Pr.lo, exponednF.Pr.hi, humpedshF.Pr.lo, humpedshF.Pr.hi,
              UshF.Pr.lo, UshF.Pr.hi);
Sfun  <- list(constantP.Cu.lo, constantP.Cu.hi, exponeupP.Cu.lo, exponeupP.Cu.hi,
              exponednP.Cu.lo, exponednP.Cu.hi, humpedshP.Cu.lo, humpedshP.Cu.hi,
              UshP.Cu.lo, UshP.Cu.hi);
Ffunl <-    c("constantF.Pr.lo", "constantF.Pr.hi", "exponeupF.Pr.lo", "exponeupF.Pr.hi",
              "exponednF.Pr.lo", "exponednF.Pr.hi", "humpedshF.Pr.lo", "humpedshF.Pr.hi",
              "UshF.Pr.lo", "UshF.Pr.hi");
Sfunl <-    c("constantP.Cu.lo", "constantP.Cu.hi", "exponeupP.Cu.lo", "exponeupP.Cu.hi",
              "exponednP.Cu.lo", "exponednP.Cu.hi", "humpedshP.Cu.lo", "humpedshP.Cu.hi",
              "UshP.Cu.lo", "UshP.Cu.hi");

rmat <- NULL;
for(i in 1:length(Ffun)){ 
    for(j in 1:length(Sfun)){
        chk  <- TIspace(FecFun=Ffun[[i]], SurFun=Sfun[[j]]);
        tot  <- sum(chk) / length(chk);
        res  <- c(Ffunl[i], Sfunl[j], tot);
        rmat <- rbind(rmat,res);
    }
}
colnames(rmat) <- c("Fecundity", "Survival", "TI sum");
ord  <- order(x=as.numeric(rmat[,3]),decreasing=TRUE);
print(rmat[ord,]);





## A figure showing all of the curves
setEPS();
postscript("AllCurves.eps",height=12,width=8);
x <- 1:20;
# ---------------------------------------------------
y <- constantF.Pr.lo(x);
par(mfrow=c(5,2),mar=c(0.25,0.5,0.25,1),lwd=3,oma=c(5,5,1,5));
plot(x=x,y=y,type="l",xlab="",ylab="",lwd=3,ylim=c(0,30),cex.lab=2,yaxt="n",xaxt="n");
axis(side=2,at=c(0,5,10,15,20,25),cex.axis=2);
y <- constantF.Pr.hi(x);
points(x=x,y=y,type="l",lwd=3,lty="dashed");
y <- constantP.Pr.lo(x);
text(x=19.5,y=28.5,labels="A",cex=2.5);
plot(x=x,y=y,type="l",xlab="",ylab="",lwd=3,ylim=c(0,1),cex.lab=2,yaxt="n",xaxt="n");
axis(side=4,at=c(0,0.2,0.4,0.6,0.8),cex.axis=2);
y <- constantP.Pr.hi(x);
points(x=x,y=y,type="l",lwd=3,lty="dashed");
text(x=19.5,y=28.5/30,labels="B",cex=2.5);
# ---------------------------------------------------
y <- exponeupF.Pr.lo(x);
plot(x=x,y=y,type="l",xlab="",ylab="",lwd=3,ylim=c(0,30),cex.lab=2,yaxt="n",xaxt="n");
axis(side=2,at=c(0,5,10,15,20,25),cex.axis=2);
y <- exponeupF.Pr.hi(x);
points(x=x,y=y,type="l",lwd=3,lty="dashed");
y <- exponeupP.Pr.lo(x);
text(x=19.5,y=28.5,labels="C",cex=2.5);
plot(x=x,y=y,type="l",xlab="",ylab="",lwd=3,ylim=c(0,1),cex.lab=2,yaxt="n",xaxt="n");
axis(side=4,at=c(0,0.2,0.4,0.6,0.8),cex.axis=2);
y <- exponeupP.Pr.hi(1:19);
points(x=1:19,y=y,type="l",lwd=3,lty="dashed");
text(x=19.5,y=28.5/30,labels="D",cex=2.5);
# ---------------------------------------------------
y <- exponednF.Pr.lo(x);
plot(x=x,y=y,type="l",xlab="",ylab="",lwd=3,ylim=c(0,30),cex.lab=2,yaxt="n",xaxt="n");
axis(side=2,at=c(0,5,10,15,20,25),cex.axis=2);
y <- exponednF.Pr.hi(x);
points(x=x,y=y,type="l",lwd=3,lty="dashed");
y <- exponednP.Pr.lo(x);
text(x=19.5,y=28.5,labels="E",cex=2.5);
plot(x=x,y=y,type="l",xlab="",ylab="",lwd=3,ylim=c(0,1),cex.lab=2,yaxt="n",xaxt="n");
axis(side=4,at=c(0,0.2,0.4,0.6,0.8),cex.axis=2);
y <- exponednP.Pr.hi(x);
points(x=x,y=y,type="l",lwd=3,lty="dashed");
text(x=19.5,y=28.5/30,labels="F",cex=2.5);
# ---------------------------------------------------
y <- humpedshF.Pr.lo(x);
plot(x=x,y=y,type="l",xlab="",ylab="",lwd=3,ylim=c(0,30),cex.lab=2,yaxt="n",xaxt="n");
axis(side=2,at=c(0,5,10,15,20,25),cex.axis=2);
y <- humpedshF.Pr.hi(x);
points(x=x,y=y,type="l",lwd=3,lty="dashed");
y <- humpedshP.Pr.lo(x);
text(x=19.5,y=28.5,labels="G",cex=2.5);
plot(x=x,y=y,type="l",xlab="",ylab="",lwd=3,ylim=c(0,1),cex.lab=2,yaxt="n",xaxt="n");
axis(side=4,at=c(0,0.2,0.4,0.6,0.8),cex.axis=2);
y <- humpedshP.Pr.hi(x);
points(x=x,y=y,type="l",lwd=3,lty="dashed");
text(x=19.5,y=28.5/30,labels="H",cex=2.5);
# ---------------------------------------------------
y <- UshF.Pr.lo(x);
plot(x=x,y=y,type="l",xlab="",ylab="",lwd=2,ylim=c(0,30),cex.lab=2,yaxt="n",xaxt="n");
axis(side=2,at=c(0,5,10,15,20,25),cex.axis=2);
axis(side=1,at=c(0,5,10,15,20),cex.axis=2);
y <- UshF.Pr.hi(x);
points(x=x,y=y,type="l",lwd=3,lty="dashed");
y <- UshP.Pr.lo(x);
text(x=19.5,y=28.5,labels="I",cex=2.5);
plot(x=x,y=y,type="l",xlab="",ylab="",lwd=3,ylim=c(0,1),cex.lab=2,yaxt="n",xaxt="n");
axis(side=4,at=c(0,0.2,0.4,0.6,0.8),cex.axis=2);
y <- UshP.Pr.hi(1:19);
points(x=1:19,y=y,type="l",lwd=3,lty="dashed");
text(x=19.5,y=28.5/30,labels="J",cex=2.5);
axis(side=1,at=c(0,5,10,15,20),cex.axis=2);

mtext(expression(paste("Age")),
      outer=TRUE,side=1,line=3.5,cex=1.75);

mtext(expression(paste("Fecundity")),
      outer=TRUE,side=2,line=2.5,cex=1.75);

mtext(expression(paste("Survival probability")),
      outer=TRUE,side=4,line=3,cex=1.75);

dev.off();


#-------------------------------------------------------------------------
# BD modified old code for making individual 3D plots
#-------------------------------------------------------------------------


# -----------------------------------------------------------------------
# Function plot the space of terminal investment (requires lhs() and TIplot() [above]).
# yr: How many years to consider.
# yF: Shows maximum plotted for increase in fecundity.
# yR: Shows maximum plotted for increase in offspring reproductive value.
TIspace <- function(FecFun,SurFun,PrDead=0.99){
    yr       <- sum(SurFun(1:1000)>(1-PrDead)) + 1;
    resvals1 <- rep(x=0, times=yr+1);
    resvals2 <- rep(x=0, times=yr+1);
    termi    <- 10000;
    RVal     <- sum((SurFun(1:termi))*FecFun(1:termi));
    for(i in 1:(yr+1)){
        # Reproductive value is resvals1:
        resvals1[i] <- sum((SurFun((i+0):termi)/SurFun(i))*FecFun((i+0):termi));
        # Residual reproductive value is resvals2:
        resvals2[i] <- sum((SurFun((i+1):termi)/SurFun(i))*FecFun((i+1):termi));
    }
    Fpar  <- seq(from=0, to=RVal, length.out=100);   # From 0 to 100% TI lifetime fec.
    Rpar  <- seq(from=0, to=RVal, length.out=100);   # From 0 to 100% TI offspring RV.
    Space <- array(data=0,dim=c(100,100,yr));
    for(age in 1:dim(Space)[3]){
        for(i in 1:dim(Space)[1]){
            for(j in 1:dim(Space)[2]){
                CHECK <- lhs( Dmx  =  Fpar[i] * FecFun(i),  # TI in fecundity 
                              DRV0 =  Rpar[j],              # TI in RV of offspring
                              mx   =  FecFun(age),          # Age fecundity
                              RV0  =  resvals1[1],          # Offspring RV 
                              RVx  =  resvals2[age]         # RV at age x
                );			
                Space[i,j,age] <- CHECK;
            }
        }
        print(age);
    }
    x <- seq(from=0,to=1,length=100);
    y <- seq(from=0,to=1,length=100);
    return(list(x,y,yr,Space)); # XXX modified output to use isosurf3D function XXX
}
# -----------------------------------------------------------------------

#loop to create all the plots
fecFunNames <- c("Constant Low Fecundity",
                 "Constant High Fecundity",
                 "Exponentially Increasing Low Fecundity",	
                 "Exponentially Increasing High Fecundity",	
                 "Exponentially Decreasing Low Fecundity",
                 "Exponentially Decreasing High Fecundity",	
                 "Hump-shaped Low Fecundity",
                 "Hump-shaped High Fecundity",
                 "U-shaped Low Fecundity",
                 "U-shaped High Fecundity")	
surFunNames <- c("Constant Low Survival Pr.",
                 "Constant High Survival Pr.",
                 "Exponentially Increasing Low Survival Pr.",
                 "Exponentially Increasing High Survival Pr.",
                 "Exponentially Decreasing Low Survival Pr.",
                 "Exponentially Decreasing High Survival Pr.",
                 "Hump-shaped Low Survival Pr.",
                 "Hump-shaped High Survival Pr.",
                 "U-shaped Low Survival Pr.",
                 "U-shaped High Survival Pr.")

fecFunlist <- list(constantF.Pr.lo,
                   constantF.Pr.hi,
                   exponeupF.Pr.lo,	
                   exponeupF.Pr.hi,	
                   exponednF.Pr.lo,
                   exponednF.Pr.hi,	
                   humpedshF.Pr.lo,
                   humpedshF.Pr.hi,
                   UshF.Pr.lo,
                   UshF.Pr.hi)	
surFunlist <- c(constantP.Cu.lo,
                constantP.Cu.hi,
                exponeupP.Cu.lo,
                exponeupP.Cu.hi,
                exponednP.Cu.lo,
                exponednP.Cu.hi,
                humpedshP.Cu.lo,
                humpedshP.Cu.hi,
                UshP.Cu.lo,
                UshP.Cu.hi)

SFnPrlist <- c(constantP.Pr.lo,
               constantP.Pr.hi,
               exponeupP.Pr.lo,
               exponeupP.Pr.hi,
               exponednP.Pr.lo,
               exponednP.Pr.hi,
               humpedshP.Pr.lo,
               humpedshP.Pr.hi,
               UshP.Pr.lo,
               UshP.Pr.hi)


library(plot3D);

# modify isosurf3D for what we want 
isosurf <- function(Space,FFn="Fecundity Function",SFn="Survival Function",SFnPr=NA){
    close.screen(all.screens=TRUE);
    opar <- par();
    #retrieve the fecundity function used by matching function name
    fecFun <- fecFunlist[[match(FFn,fecFunNames)]];
    #retrieve the survival function used by matching function name
    surFun <- surFunlist[[match(SFn,surFunNames)]];
    # matrix to split plotting screen
    fig.mat <- matrix(c(0,0.66,0,1,0.66,1,0.5,1,0.66,1,0,0.5),ncol=4,byrow=TRUE);
    par(split.screen(fig.mat)) # split the graph in 3
    screen(1,new=TRUE) #plot on screen 1
    par(mar=c(5,5,1,1));
    if(sum(Space[[4]])==0){
        plot(x=1:10,y=1:10,type="n",xaxt="n",yaxt="n",xlab="",ylab="");
        text(x=5.5,y=5.5,"No TI predicted",cex=1);
        titl <- paste(FFn,"\n",SFn)
        text(x=5.5,y=8.5,titl,cex=1);
        #return(NULL); removed since all the function should be evaluated
    }else{
        x <- Space[[1]] #added else argument for the function to work
        y <- Space[[2]]
        z <- (1:Space[[3]])
        colvar=Space[[4]]
        titl <- paste(FFn,"\n",SFn) # Commenting this out with main="" below.
        pmat <- isosurf3D(x, y, z, colvar=colvar, main="", # <- change to add title.
                          xlab="", ylab="", zlab="",ticktype="detailed", bty = "g",
                          level=0.5, theta=30, cex.axis=0.6);
        	text(trans3d(x=mean(c(min(x),max(x))), y=-0.2, z=min(z),pmat=pmat),
        		lab=expression(paste("TI in fecundity (",Delta,m[x],")")),
        		srt=-30,cex=1)	
        	text(trans3d(1.2, mean(c(min(y),max(y))), min(z), pmat=pmat),
        		lab=expression(paste("TI in offspring RV (",Delta,RV[0],")")),
        		srt=65,cex=1)
        	text(trans3d(-0.20, -0.1, mean(c(min(z),max(z))), pmat=pmat), lab="Age", srt=120, 
        	     xpd=TRUE,cex=1.25);
    }
    screen(2,new=TRUE) #plot on screen 1
    x <- 1:20
    y <- fecFun(x)
    plot(x=x,y=y,type="n",xlab="Age",ylab="Fecundity",lwd=2,ylim=c(0,30),
         cex.lab=1,cex.axis=0.8, mgp=c(2,1,0));
    polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80",lwd=1.5);
    SfPrFun <- SFnPr;
    screen(3,new=TRUE) #plot on screen 1
    y <- SfPrFun(x)
    plot(x=x,y=y,type="n",xlab="Age",ylab="Survival",lwd=2,ylim=c(0,1),
         cex.axis=0.8, cex.lab=1, mgp=c(2,1,0));
    polygon(x=c(x,rev(x)),y=c(rep(0,length(y)),rev(y)),col="grey80",lwd=1.5);
    par <- opar
}


Space <- TIspace( FecFun = fecFunlist[[i]], SurFun   = surFunlist[[j]]);
isosurf(Space,FFn=fecFunNames[i],SFn=surFunNames[j],SFnPr=SFnPrlist[[j]]);

## Below prints everything (takes a long time);

pdf(file=paste("Plots3Dsum.pdf",sep=""),onefile=TRUE);		
for(i in 1:10){
    print(fecFunNames[i]);
    for(j in 1:10){
        print(surFunNames[j]);
        Space <- TIspace( FecFun   = fecFunlist[[i]], 
                          SurFun   = surFunlist[[j]]
        );
        isosurf(Space,FFn=fecFunNames[i],SFn=surFunNames[j],SFnPr=SFnPrlist[[j]])
    }
}

dev.off();










