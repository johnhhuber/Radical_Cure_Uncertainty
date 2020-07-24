
OUTPUT <- read.table("Output/model_output.txt")

N_mosq = 2

mosq_species = c("faruati", "punctulatus", "koliensis")

mosq_comp = c( "S_M_dom", "E_M_dom", "I_M_dom",
		   "S_M_occ", "E_M_occ", "I_M_occ")

colnames(OUTPUT) <- c("time",
                      "S", "I_PCR", "I_LM", "D", "T", "P",
                      mosq_comp,
                      "N_pop", "PvPR_PCR", "PvPR_LM", "Pv_clin", 
                      "PvHR", "PvHR_batches", 
                      "new_PCR", "new_LM", "new_D", "new_BS", "new_PQ", "new_TQ",
                      "G6PD_tests", 
                      "PQ_dose_eff","TQ_dose_eff",
                      "EIR_dom", "EIR_occ", "LLIN_cov", "IRS_cov", "IVM_cov",
                      "CQ_treat", "PQ_treat", "TQ_treat",
                      "PQ_overtreat", "PQ_overtreat_9m", "TQ_overtreat", "TQ_overtreat_9m",
                      "pregnant", "cases_M_O16", "cases_M_U16", "cases_F_O16", "cases_F_U16", "cases_preg",
			    "A_par", "A_clin")


par(ask=TRUE)

#############################################
#############################################
##          ##                             ##
##  PLOT 1  ##  Human dynamics             ##
##          ##                             ##
#############################################
#############################################

par(mfrow=c(2,3))

plot(x=OUTPUT$time/365, y=OUTPUT$S, type='l',
xlab="time (years)", ylab="number", main="S")

plot(x=OUTPUT$time/365, y=OUTPUT$I_PCR, type='l',
xlab="time (years)", ylab="number", main="I_PCR")

plot(x=OUTPUT$time/365, y=OUTPUT$I_LM, type='l',
xlab="time (years)", ylab="number", main="I_LM")

plot(x=OUTPUT$time/365, y=OUTPUT$D, type='l',
xlab="time (years)", ylab="number", main="D")

plot(x=OUTPUT$time/365, y=OUTPUT$T, type='l',
xlab="time (years)", ylab="number", main="T")

plot(x=OUTPUT$time/365, y=OUTPUT$P, type='l',
xlab="time (years)", ylab="number", main="P")


#############################################
#############################################
##          ##                             ##
##  PLOT 2  ##  Mosquito dynamics          ##
##          ##                             ##
#############################################
#############################################

par(mfrow=c(2,3))


plot(x=OUTPUT$time/365, y=OUTPUT$S_M_dom, type='l',
xlab="time (years)", ylab="number", main="S_M_dom")

plot(x=OUTPUT$time/365, y=OUTPUT$E_M_dom, type='l',
xlab="time (years)", ylab="number", main="E_M_dom")

plot(x=OUTPUT$time/365, y=OUTPUT$I_M_dom, type='l',
xlab="time (years)", ylab="number", main="I_M_dom")


plot(x=OUTPUT$time/365, y=OUTPUT$S_M_occ, type='l',
xlab="time (years)", ylab="number", main="S_M_occ")

plot(x=OUTPUT$time/365, y=OUTPUT$E_M_occ, type='l',
xlab="time (years)", ylab="number", main="E_M_occ")

plot(x=OUTPUT$time/365, y=OUTPUT$I_M_occ, type='l',
xlab="time (years)", ylab="number", main="I_M_occ")






#############################################
#############################################
##                                         ##
##  PROCESS TREATMENT OUTPUT               ##
##                                         ##
#############################################
#############################################

INT_cov <- read.table("Intervention_coverage.txt")


MDA0_t  <- INT_cov[1,1+which(INT_cov[which(INT_cov[,1]=="MDA0_cover"),-1]>0)]
MDA1_t  <- INT_cov[1,1+which(INT_cov[which(INT_cov[,1]=="MDA1_cover"),-1]>0)]
MDA2_t  <- INT_cov[1,1+which(INT_cov[which(INT_cov[,1]=="MDA2_cover"),-1]>0)]
MSAT0_t <- INT_cov[1,1+which(INT_cov[which(INT_cov[,1]=="MSAT0_cover"),-1]>0)]
MSAT1_t <- INT_cov[1,1+which(INT_cov[which(INT_cov[,1]=="MSAT1_cover"),-1]>0)]
MSAT2_t <- INT_cov[1,1+which(INT_cov[which(INT_cov[,1]=="MSAT2_cover"),-1]>0)]
STAT1_t <- INT_cov[1,1+which(INT_cov[which(INT_cov[,1]=="STAT1_cover"),-1]>0)]
STAT2_t <- INT_cov[1,1+which(INT_cov[which(INT_cov[,1]=="STAT2_cover"),-1]>0)]

MDA0_CQ_cov  <- rep(NA, length(MDA0_t))

MDA1_CQ_cov  <- rep(NA, length(MDA1_t))
MDA1_PQ_cov  <- rep(NA, length(MDA1_t))

MDA2_CQ_cov  <- rep(NA, length(MDA2_t))
MDA2_PQ_cov  <- rep(NA, length(MDA2_t))
MDA2_TQ_cov  <- rep(NA, length(MDA2_t))

MSAT0_CQ_cov <- rep(NA, length(MSAT0_t))

MSAT1_CQ_cov <- rep(NA, length(MSAT1_t))
MSAT1_PQ_cov <- rep(NA, length(MSAT1_t))

MSAT2_CQ_cov <- rep(NA, length(MSAT2_t))
MSAT2_PQ_cov <- rep(NA, length(MSAT2_t))
MSAT2_TQ_cov <- rep(NA, length(MSAT2_t))

STAT1_CQ_cov <- rep(NA, length(STAT1_t))
STAT1_PQ_cov <- rep(NA, length(STAT1_t))

STAT2_CQ_cov <- rep(NA, length(STAT2_t))
STAT2_PQ_cov <- rep(NA, length(STAT2_t))
STAT2_TQ_cov <- rep(NA, length(STAT2_t))


if( length(MDA0_t) > 0 )
{
	for(i in 1:length(MDA0_t))
	{
		index_t <- which(OUTPUT$time/365 == MDA0_t[i])	

		MDA0_CQ_cov[i] <- OUTPUT$CQ_treat[index_t]/OUTPUT$N_pop[1]

		OUTPUT$CQ_treat[index_t] <- mean( OUTPUT$CQ_treat[c(index_t-1,index_t+1)] )

		OUTPUT$new_BS[index_t] <- mean( OUTPUT$new_BS[c(index_t-1,index_t+1)] )

		OUTPUT$cases_M_O16[index_t] <- mean( OUTPUT$cases_M_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_M_U16[index_t] <- mean( OUTPUT$cases_M_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_O16[index_t] <- mean( OUTPUT$cases_F_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_U16[index_t] <- mean( OUTPUT$cases_F_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_preg[index_t]  <- mean( OUTPUT$cases_preg[c(index_t-1,index_t+1)] )
	}
}


if( length(MDA1_t) > 0 )
{
	for(i in 1:length(MDA1_t))
	{
		index_t <- which(OUTPUT$time/365 == MDA1_t[i])	

		MDA1_CQ_cov[i] <- OUTPUT$CQ_treat[index_t]/OUTPUT$N_pop[1]
		MDA1_PQ_cov[i] <- OUTPUT$PQ_treat[index_t]/OUTPUT$N_pop[1]

		OUTPUT$CQ_treat[index_t] <- mean( OUTPUT$CQ_treat[c(index_t-1,index_t+1)] )
		OUTPUT$PQ_treat[index_t] <- mean( OUTPUT$PQ_treat[c(index_t-1,index_t+1)] )

		OUTPUT$new_BS[index_t] <- mean( OUTPUT$new_BS[c(index_t-1,index_t+1)] )

		OUTPUT$cases_M_O16[index_t] <- mean( OUTPUT$cases_M_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_M_U16[index_t] <- mean( OUTPUT$cases_M_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_O16[index_t] <- mean( OUTPUT$cases_F_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_U16[index_t] <- mean( OUTPUT$cases_F_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_preg[index_t]  <- mean( OUTPUT$cases_preg[c(index_t-1,index_t+1)] )
	}
}


if( length(MDA2_t) > 0 )
{
	for(i in 1:length(MDA2_t))
	{
		index_t <- which(OUTPUT$time/365 == MDA2_t[i])	

		MDA2_CQ_cov[i] <- OUTPUT$CQ_treat[index_t]/OUTPUT$N_pop[1]
		MDA2_PQ_cov[i] <- OUTPUT$PQ_treat[index_t]/OUTPUT$N_pop[1]
		MDA2_TQ_cov[i] <- OUTPUT$TQ_treat[index_t]/OUTPUT$N_pop[1]

		OUTPUT$CQ_treat[index_t] <- mean( OUTPUT$CQ_treat[c(index_t-1,index_t+1)] )
		OUTPUT$PQ_treat[index_t] <- mean( OUTPUT$PQ_treat[c(index_t-1,index_t+1)] )
		OUTPUT$TQ_treat[index_t] <- mean( OUTPUT$TQ_treat[c(index_t-1,index_t+1)] )

		OUTPUT$new_BS[index_t] <- mean( OUTPUT$new_BS[c(index_t-1,index_t+1)] )

		OUTPUT$cases_M_O16[index_t] <- mean( OUTPUT$cases_M_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_M_U16[index_t] <- mean( OUTPUT$cases_M_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_O16[index_t] <- mean( OUTPUT$cases_F_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_U16[index_t] <- mean( OUTPUT$cases_F_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_preg[index_t]  <- mean( OUTPUT$cases_preg[c(index_t-1,index_t+1)] )
	}
}


if( length(MSAT0_t) > 0 )
{
	for(i in 1:length(MSAT0_t))
	{
		index_t <- which(OUTPUT$time/365 == MSAT0_t[i])	

		MSAT0_CQ_cov[i] <- OUTPUT$CQ_treat[index_t]/OUTPUT$N_pop[1]

		OUTPUT$CQ_treat[index_t] <- mean( OUTPUT$CQ_treat[c(index_t-1,index_t+1)] )

		OUTPUT$new_BS[index_t] <- mean( OUTPUT$new_BS[c(index_t-1,index_t+1)] )

		OUTPUT$cases_M_O16[index_t] <- mean( OUTPUT$cases_M_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_M_U16[index_t] <- mean( OUTPUT$cases_M_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_O16[index_t] <- mean( OUTPUT$cases_F_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_U16[index_t] <- mean( OUTPUT$cases_F_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_preg[index_t]  <- mean( OUTPUT$cases_preg[c(index_t-1,index_t+1)] )
	}
}


if( length(MSAT1_t) > 0 )
{
	for(i in 1:length(MSAT1_t))
	{
		index_t <- which(OUTPUT$time/365 == MSAT1_t[i])	

		MSAT1_CQ_cov[i] <- OUTPUT$CQ_treat[index_t]/OUTPUT$N_pop[1]
		MSAT1_PQ_cov[i] <- OUTPUT$PQ_treat[index_t]/OUTPUT$N_pop[1]

		OUTPUT$CQ_treat[index_t] <- mean( OUTPUT$CQ_treat[c(index_t-1,index_t+1)] )
		OUTPUT$PQ_treat[index_t] <- mean( OUTPUT$PQ_treat[c(index_t-1,index_t+1)] )

		OUTPUT$new_BS[index_t] <- mean( OUTPUT$new_BS[c(index_t-1,index_t+1)] )

		OUTPUT$cases_M_O16[index_t] <- mean( OUTPUT$cases_M_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_M_U16[index_t] <- mean( OUTPUT$cases_M_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_O16[index_t] <- mean( OUTPUT$cases_F_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_U16[index_t] <- mean( OUTPUT$cases_F_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_preg[index_t]  <- mean( OUTPUT$cases_preg[c(index_t-1,index_t+1)] )
	}
}


if( length(MSAT2_t) > 0 )
{
	for(i in 1:length(MSAT2_t))
	{
		index_t <- which(OUTPUT$time/365 == MSAT2_t[i])	

		MSAT2_CQ_cov[i] <- OUTPUT$CQ_treat[index_t]/OUTPUT$N_pop[1]
		MSAT2_PQ_cov[i] <- OUTPUT$PQ_treat[index_t]/OUTPUT$N_pop[1]
		MSAT2_TQ_cov[i] <- OUTPUT$TQ_treat[index_t]/OUTPUT$N_pop[1]

		OUTPUT$CQ_treat[index_t] <- mean( OUTPUT$CQ_treat[c(index_t-1,index_t+1)] )
		OUTPUT$PQ_treat[index_t] <- mean( OUTPUT$PQ_treat[c(index_t-1,index_t+1)] )
		OUTPUT$TQ_treat[index_t] <- mean( OUTPUT$TQ_treat[c(index_t-1,index_t+1)] )

		OUTPUT$new_BS[index_t] <- mean( OUTPUT$new_BS[c(index_t-1,index_t+1)] )

		OUTPUT$cases_M_O16[index_t] <- mean( OUTPUT$cases_M_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_M_U16[index_t] <- mean( OUTPUT$cases_M_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_O16[index_t] <- mean( OUTPUT$cases_F_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_U16[index_t] <- mean( OUTPUT$cases_F_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_preg[index_t]  <- mean( OUTPUT$cases_preg[c(index_t-1,index_t+1)] )
	}
}


if( length(STAT1_t) > 0 )
{
	for(i in 1:length(STAT1_t))
	{
		index_t <- which(OUTPUT$time/365 == STAT1_t[i])	

		STAT1_CQ_cov[i] <- OUTPUT$CQ_treat[index_t]/OUTPUT$N_pop[1]
		STAT1_PQ_cov[i] <- OUTPUT$PQ_treat[index_t]/OUTPUT$N_pop[1]

		OUTPUT$CQ_treat[index_t] <- mean( OUTPUT$CQ_treat[c(index_t-1,index_t+1)] )
		OUTPUT$PQ_treat[index_t] <- mean( OUTPUT$PQ_treat[c(index_t-1,index_t+1)] )

		OUTPUT$new_BS[index_t] <- mean( OUTPUT$new_BS[c(index_t-1,index_t+1)] )

		OUTPUT$cases_M_O16[index_t] <- mean( OUTPUT$cases_M_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_M_U16[index_t] <- mean( OUTPUT$cases_M_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_O16[index_t] <- mean( OUTPUT$cases_F_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_U16[index_t] <- mean( OUTPUT$cases_F_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_preg[index_t]  <- mean( OUTPUT$cases_preg[c(index_t-1,index_t+1)] )
	}
}


if( length(STAT2_t) > 0 )
{
	for(i in 1:length(STAT2_t))
	{
		index_t <- which(OUTPUT$time/365 == STAT2_t[i])	

		STAT2_CQ_cov[i] <- OUTPUT$CQ_treat[index_t]/OUTPUT$N_pop[1]
		STAT2_CQ_cov[i] <- OUTPUT$PQ_treat[index_t]/OUTPUT$N_pop[1]
		STAT2_CQ_cov[i] <- OUTPUT$TQ_treat[index_t]/OUTPUT$N_pop[1]

		OUTPUT$CQ_treat[index_t] <- mean( OUTPUT$CQ_treat[c(index_t-1,index_t+1)] )
		OUTPUT$PQ_treat[index_t] <- mean( OUTPUT$PQ_treat[c(index_t-1,index_t+1)] )
		OUTPUT$TQ_treat[index_t] <- mean( OUTPUT$TQ_treat[c(index_t-1,index_t+1)] )

		OUTPUT$new_BS[index_t] <- mean( OUTPUT$new_BS[c(index_t-1,index_t+1)] )

		OUTPUT$cases_M_O16[index_t] <- mean( OUTPUT$cases_M_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_M_U16[index_t] <- mean( OUTPUT$cases_M_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_O16[index_t] <- mean( OUTPUT$cases_F_O16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_F_U16[index_t] <- mean( OUTPUT$cases_F_U16[c(index_t-1,index_t+1)] )
		OUTPUT$cases_preg[index_t]  <- mean( OUTPUT$cases_preg[c(index_t-1,index_t+1)] )
	}
}






N_window = 50

N_tt = length(OUTPUT$CQ_treat)


CQ_smooth = rep(NA, N_tt )
PQ_smooth = rep(NA, N_tt )
TQ_smooth = rep(NA, N_tt )
new_BS_smooth = rep(NA, N_tt )

for(i in 1:N_tt)
{
	CQ_smooth[i] = 365*mean( OUTPUT$CQ_treat[max(1,i-N_window):min(i+N_window,N_tt)]/OUTPUT$N_pop[max(1,i-N_window):min(i+N_window,N_tt)])
	PQ_smooth[i] = 365*mean( OUTPUT$PQ_treat[max(1,i-N_window):min(i+N_window,N_tt)]/OUTPUT$N_pop[max(1,i-N_window):min(i+N_window,N_tt)])
	TQ_smooth[i] = 365*mean( OUTPUT$TQ_treat[max(1,i-N_window):min(i+N_window,N_tt)]/OUTPUT$N_pop[max(1,i-N_window):min(i+N_window,N_tt)])
	new_BS_smooth[i] = mean( OUTPUT$new_BS[max(1,i-N_window):min(i+N_window,N_tt)] )
}


cases_M_U16_smooth = rep(NA, N_tt )
cases_M_O16_smooth = rep(NA, N_tt )
cases_F_U16_smooth = rep(NA, N_tt )
cases_F_O16_smooth = rep(NA, N_tt )
cases_preg_smooth = rep(NA, N_tt )


for(i in 1:N_tt)
{
	cases_M_U16_smooth[i] = mean( OUTPUT$cases_M_U16[max(1,i-N_window):min(i+N_window,N_tt)] )
	cases_M_O16_smooth[i] = mean( OUTPUT$cases_M_O16[max(1,i-N_window):min(i+N_window,N_tt)] )
	cases_F_U16_smooth[i] = mean( OUTPUT$cases_F_U16[max(1,i-N_window):min(i+N_window,N_tt)] )
	cases_F_O16_smooth[i] = mean( OUTPUT$cases_F_O16[max(1,i-N_window):min(i+N_window,N_tt)] )
	cases_preg_smooth[i]  = mean( OUTPUT$cases_preg[max(1,i-N_window):min(i+N_window,N_tt)] )
}



#############################################
#############################################
##          ##                             ##
##  PLOT 3  ##  Extra check                ##
##          ##                             ##
#############################################
#############################################



par(mfrow=c(1,3))


#############################
## PANEL 1: Pregnancy

plot(x=OUTPUT$time/365, y=OUTPUT$pregnant/OUTPUT$N_pop, 
type='l', col="black",
ylim=c(0,max(OUTPUT$pregnant/OUTPUT$N_pop)),
xlab="time (years)", ylab="proportion", 
main="Pregnancy")


#############################
## PANEL 2: A_par

plot(x=OUTPUT$time/365, y=OUTPUT$A_par, 
type='l', col="black",
ylim=c(0,max(OUTPUT$A_par)),
xlab="time (years)", ylab="proportion", 
main="Parasite immunity")



#############################
## PANEL 3: A_clin

plot(x=OUTPUT$time/365, y=OUTPUT$A_clin, 
type='l', col="black",
ylim=c(0,max(OUTPUT$A_clin)),
xlab="time (years)", ylab="proportion", 
main="Clinical immunity")




#############################################
#############################################
##          ##                             ##
##  PLOT 3  ##  Prevalence                 ##
##          ##                             ##
#############################################
#############################################


# par(mfrow=c(3,3))
par(mfrow=c(2,4))


#############################
## PANEL 1: Prevalence

plot(x=OUTPUT$time/365, y=OUTPUT$PvPR_LM/OUTPUT$N_pop, 
type='l', col="black",
ylim=c(0,max(OUTPUT$PvHR/OUTPUT$N_pop)),
xlab="time (years)", ylab="prevalence", main="Prevalence")

points(x=OUTPUT$time/365, y=OUTPUT$PvPR_PCR/OUTPUT$N_pop, 
type='l', col="red")

points(x=OUTPUT$time/365, y=OUTPUT$PvHR/OUTPUT$N_pop, 
type='l', col="gold")


legend(x="topright", bty='n', cex=1,
legend=c( "hypnozoites", "PCR", "LM"),
fill=c( "gold",  "red",  "black" ),
border=c( "gold",  "red",  "black" ) )


#############################
## PANEL 2: Clinical incidence

plot(x=OUTPUT$time/365, y=new_BS_smooth*365*1000/OUTPUT$N_pop, 
type='l',
xlab="time (years)", ylab="detected cases per year", main="Annual Parasite Index (API) per 1000")


#############################
## PANEL 3: EIR

plot(x=OUTPUT$time/365, y=365*OUTPUT$EIR_dom, type='l', 
ylim=c(0, max( c(365*OUTPUT$EIR_dom, 365*OUTPUT$EIR_occ) )),
#ylim=c(0,50),
xlab="time (years)", ylab="EIR (ibppy)", main="EIR")

points(x=OUTPUT$time/365, y=365*OUTPUT$EIR_occ, 
type='l', col="forestgreen")




#############################
## PANEL 4: Hypnozoite broods

plot(x=OUTPUT$time/365, y=OUTPUT$PvHR_batch/OUTPUT$N_pop, type='l', ylim=c(0,max(OUTPUT$PvHR_batch/OUTPUT$N_pop)),
xlab="time (years)", ylab="hypnozoite batches", 
main="Mean number of hypnozoite batches per person")



#############################
## PANEL 5: Age and gender stratified cases

plot(x=OUTPUT$time/365, y=cases_M_U16_smooth, 
type='l', col="royalblue",
ylim=1.2*c(0, max( c(cases_M_U16_smooth, cases_M_O16_smooth, cases_F_U16_smooth, cases_F_O16_smooth) ) ),
xlab="time (years)", ylab="daily cases", 
main="Age and gender stratified cases")

points(x=OUTPUT$time/365, y=cases_M_O16_smooth, 
type='l', col="royalblue", lty="dashed")

points(x=OUTPUT$time/365, y=cases_F_O16_smooth, 
type='l', col="green", lty="dashed")

points(x=OUTPUT$time/365, y=cases_F_U16_smooth, 
type='l', col="green", lty="dashed")

points(x=OUTPUT$time/365, y=cases_preg_smooth, 
type='l', col="forestgreen", lty="dashed")





legend(x="topright", bty='n', cex=1,
legend=c( "males", "females", "pregnant"),
fill=c( "royalblue", "green", "forestgreen" ),
border=c( "royalblue", "green", "forestgreen" ) )



#############################
## PANEL 6: Case management treatment coverage


plot(OUTPUT$time/365, y=CQ_smooth, type='l', col="blue",
ylim=c(0,1),
xlab="time (years)", ylab="treatments per person per year", main="Case management")


points(x=OUTPUT$time/365, y=PQ_smooth, type='l', col="magenta")

points(x=OUTPUT$time/365, y=TQ_smooth, type='l', col="orange")





legend(x="topright", bty='n', cex=1,
legend=c( "blood-stage", "primaquine", "tafenoquine"),
fill=c( "blue", "magenta", "orange" ),
border=c( "blue", "magenta", "orange" ) )


#############################
## PANEL 7: Primaquine and tafenoquine overtreatment


plot(100, y=100, 
ylim=c(0,1), xlim=range(OUTPUT$time/365),
xlab="time (years)", ylab="proportion of population", 
main="8-aminoquinoline overtreatment")

points(x=OUTPUT$time/365, y=OUTPUT$TQ_overtreat_9m/OUTPUT$N_pop + OUTPUT$PQ_overtreat_9m/OUTPUT$N_pop, 
type='l', lwd=2, col="magenta")

points(x=OUTPUT$time/365, y=OUTPUT$TQ_overtreat_9m/OUTPUT$N_pop, 
type='l', lwd=2, col="orange")


legend(x="topright", bty='n', cex=1,
legend=c("primaquine", "tafenoquine"),
fill=c(  "magenta", "orange" ),
border=c(  "magenta", "orange" ) )




#############################
## PANEL 8: Treatment doses

plot(x=OUTPUT$time/365, y=cumsum(OUTPUT$new_BS),
     type='l', col="blue",
     ylim=c(0,max(cumsum(OUTPUT$new_BS))),
     xlab="time (years)", ylab="Cumulative doses", main="Treatment doses and effective doses")

points(x=OUTPUT$time/365, y=cumsum(OUTPUT$new_PQ),
       type='l', col="black")

points(x=OUTPUT$time/365, y=cumsum(OUTPUT$PQ_dose_eff),
       type='l', col="magenta")

points(x=OUTPUT$time/365, y=cumsum(OUTPUT$new_TQ),
       type='l', col="red")

points(x=OUTPUT$time/365, y=cumsum(OUTPUT$TQ_dose_ef),
       type='l', col="orange")

legend(x="topleft", bty='n', cex=1,
       legend=c( "BS_doses","PQ doses administered", "PQ doses effective",
                 "TQ doses administered", "TQ doses effective"),
       fill=c("blue","black","magenta","red","orange"),
       border=c("blue","black","magenta","red","orange") )







#############################
## PANEL 9: G6PD testing
#
#plot(x=OUTPUT$time/365, y=cumsum(OUTPUT$G6PD_tests),
#     type='l', col="black",
#     ylim=c(0,max(cumsum(OUTPUT$G6PD_tests))),
#     xlab="time (years)", ylab="Cumulative G6PD tests", main="G6PD testing")
#
#legend(x="topleft", bty='n', cex=1,
#       legend=c( "G6PD tests administered"),
#       fill=c("black"),
#       border=c("black") )




