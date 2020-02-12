rm(list=ls())


###############################################################
###############################################################
##     ##  Case management of clinical episodes with 
##  0  ##  chloraquine (or another blood-stage drug).
##     ##  Includes information for coverage, efficacy and
#########  duration of prophylaxis.
 
CM0_years      <-  c( )        ## Start of case management regimen
CM0_cover      <-  c( )        ## Coverage of chloroquine (accounts for proportion seeking care)
CM0_CQ_eff     <-  c( )        ## Efficacy of chloroquine drugs 
CM0_CQ_proph   <-  c( )        ## Duration of chloroquine prophylaxis


###############################################################
###############################################################
##     ##  Case management of clinical episodes with 
##  1  ##  a hypnozoiticidal drug (primaquine).
##     ##  Includes information for coverage, efficacy and
#########  duration of prophylaxis.

CM1_years           <- c( 2015 )     ## Start of case management regimen
CM1_cover           <- c( 0.95 )     ## Coverage of chloroquine (accounts for proportion seeking care)
CM1_CQ_eff          <- c( 0.899 )    ## chloroquine drug efficacy (on its own)
CM1_CQ_eff_wPQ      <- c( 0.946 )    ## chloroquine drug efficacy (co-administered with PQ)
CM1_CQ_proph        <- c( 28 )       ## duration of prophylaxis of chloroquine 
CM1_PQ_eff          <- c( 0.686 )    ## primaquine efficacy
CM1_PQ_proph        <- c( 8 )        ## duration of primaquine prophylaxis (number of days + 1)
CM1_PQ_adhere       <- c( 0.667  )   ## adherence to full primaquine regimen
CM1_PQ_lowage       <- c( 180 )      ## youngest age for primaquine treatment regimen
CM1_PQ_G6PD_risk    <- c( 1 )        ## Risk in G6PD-deficient individuals (0 = no risk; 1 = risk)
CM1_PQ_CYP2D6_risk  <- c( 1 )        ## Risk of not working in low CYP2D6 metabolizers  (0 = no risk; 1 = risk, i.e. doesn't work)
CM1_PQ_preg_risk    <- c( 1 )        ## Risk in pregnant women  (0 = no risk; 1 = risk)
CM1_G6PD_test       <- c( 1 )	       ## Is G6PD testing implemented in case management (0 = no; 1 = yes)


###############################################################
###############################################################
##     ##  Case management of clinical episodes with 
##  2  ##  a hypnozoiticidal drug (tafenoquine).
##     ##  Includes information for coverage, efficacy and
#########  duration of prophylaxis.

CM2_years           <- c( )     ## Start of case management regimen
CM2_cover           <- c( )     ## Coverage of chloroquine (accounts for proportion seeking care)
CM2_CQ_eff          <- c( )     ## chloroquine efficacy (on its own)
CM2_CQ_eff_wPQ      <- c( )     ## chloroquine efficacy (co-administered with PQ)
CM2_CQ_proph        <- c( )     ## duration of prophylaxis of chloroquine (co-administered with PQ)
CM2_PQ_eff          <- c( )     ## primaquine efficacy
CM2_PQ_proph        <- c( )     ## duration of primaquine prophylaxis (number of days + 1)
CM2_PQ_adhere       <- c( )     ## adherence to full primaquine regimen
CM2_PQ_lowage       <- c( )     ## youngest age for primaquine treatment regimen
CM2_PQ_G6PD_risk    <- c( )     ## Risk in G6PD-deficient individuals (0 = no risk; 1 = risk)
CM2_PQ_CYP2D6_risk  <- c( )     ## Risk of not working in low CYP2D6 metabolizers  (0 = no risk; 1 = risk, i.e. doesn't work)
CM2_PQ_preg_risk    <- c( )     ## Risk in pregnant women  (0 = no risk; 1 = risk)
CM2_TQ_eff          <- c( )     ## tafenoquine efficacy
CM2_TQ_proph        <- c( )     ## duration of tafenoquine prophylaxis (number of days + 1)
CM2_TQ_adhere       <- c( )     ## adherence to full tafenoquine regimen
CM2_TQ_lowage       <- c( )     ## youngest age for tafenoquine treatment regimen
CM2_TQ_G6PD_risk    <- c( )     ## Risk in G6PD-deficient individuals (0 = no risk; 1 = risk)
CM2_TQ_CYP2D6_risk  <- c( )     ## Risk of not working in low CYP2D6 metabolizers  (0 = no risk; 1 = risk, i.e. doesn't work)
CM2_TQ_preg_risk    <- c( )     ## Risk in pregnant women  (0 = no risk; 1 = risk)
CM2_G6PD_test       <- c( )	  ## Is G6PD testing implemented in case management (0 = no; 1 = yes)


###############################################################
###############################################################
##     ##  LLINs
##  3  ##  At each timepoint, x% of the population receive a new net.
##     ##  If they have an existing net it is replaced.
#########

LLIN_years <- c( ) ## c(2008, 2011, 2015, 2018, 2021, 2024, 2027, 2030, 2033, 2036, 2039, 2042)
LLIN_cover <- c( ) ## c( 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)


###############################################################
###############################################################
##     ##  IRS
##  4  ##  At each timepoint, x% of the population is protected by spraying
##     ##  Note the model does not have household structure and protection
#########  status is on the individual level.

IRS_years <- c( )
IRS_cover <- c( )


###############################################################
###############################################################
##     ##  Mass drug administration with chloroquine (or another blood-stage drug).
##  5  ##  Includes information for coverage, efficacy and
##     ##  duration of prophylaxis.
#########

MDA0_years     <- c( )      ## Time of MDA programme
MDA0_cover     <- c( )      ## Coverage of blood-stage drugs        
MDA0_CQ_eff    <- c( )      ## Efficacy of blood-stage drugs
MDA0_CQ_proph  <- c( )      ## Duration of blood-stage prophylaxis 


###############################################################
###############################################################
##     ##  Mass drug administration with chloroquine
##  6  ##  and a hypnozoiticidal drug (primaquine).
##     ##  Includes information for coverage, efficacy and
#########  duration of prophylaxis.

MDA1_years           <- c( )     ## Time of MDA programme
MDA1_cover           <- c( )     ## Coverage of blood-stage drugs        
MDA1_CQ_eff          <- c( )     ## Efficacy of chloroquine (on its own)
MDA1_CQ_eff_wPQ      <- c( )     ## Efficacy of chloroquine (co-administered with PQ)
MDA1_CQ_proph        <- c( )     ## Duration of chloroquine prophylaxis
MDA1_PQ_eff          <- c( )     ## Efficacy of primaquine
MDA1_PQ_proph        <- c( )     ## Duration of primaquine prophylaxis
MDA1_PQ_adhere       <- c( )     ## Primaquine adherence
MDA1_PQ_lowage       <- c( )     ## Lower age limit for primaquine treatment (in days)
MDA1_PQ_G6PD_risk    <- c( )     ## Risk in G6PD - deficient individuals
MDA1_PQ_CYP2D6_risk  <- c( )     ## Risk of not working in low CYP2D6 metabolizers
MDA1_PQ_preg_risk    <- c( )     ## Risk in pregnant women
MDA1_G6PD_test       <- c( )     ## G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


###############################################################
###############################################################
##     ##  Mass drug administration with chloroquine
##  7  ##  and a hypnozoiticidal drug (tafenoquine).
##     ##  Includes information for coverage, efficacy and
#########  duration of prophylaxis.

MDA2_years           <- c( )     ## Time of MDA programme
MDA2_cover           <- c( )     ## Coverage of chloroquine        
MDA2_CQ_eff          <- c( )     ## Efficacy of chloroquine (on its own)
MDA2_CQ_eff_wPQ      <- c( )     ## Efficacy of chloroquine (co-administered with PQ)
MDA2_CQ_proph        <- c( )     ## Duration of chloroquine prophylaxis
MDA2_PQ_eff          <- c( )     ## Efficacy of primaquine
MDA2_PQ_proph        <- c( )     ## Duration of primaquine prophylaxis
MDA2_PQ_adhere       <- c( )     ## Primaquine adherence
MDA2_PQ_lowage       <- c( )     ## Lower age limit for primaquine treatment (in days)
MDA2_PQ_G6PD_risk    <- c( )     ## Risk in G6PD - deficient individuals
MDA2_PQ_CYP2D6_risk  <- c( )     ## Risk of not working in low CYP2D6 metabolizers
MDA2_PQ_preg_risk    <- c( )     ## Risk in pregnant women
MDA2_TQ_eff          <- c( )     ## Efficacy of tafenoquine
MDA2_TQ_proph        <- c( )     ## Duration of tafenoquine prophylaxis
MDA2_TQ_adhere       <- c( )     ## Tafenoquine adherence
MDA2_TQ_lowage       <- c( )     ## Lower age limit for tafenoquine treatment (in days)
MDA2_TQ_G6PD_risk    <- c( )     ## Risk in G6PD - deficient individuals
MDA2_TQ_CYP2D6_risk  <- c( )     ## Risk of not working in low CYP2D6 metabolizers
MDA2_TQ_preg_risk    <- c( )     ## Risk in pregnant women
MDA2_G6PD_test       <- c( )     ## G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


###############################################################
###############################################################
##     ##  Mass drug administration with chloroquine.
##  8  ##  Includes information for coverage, efficacy and
##     ##  duration of prophylaxis.
#########

MSAT0_years       <-  c( )      ## Time of MDA programme
MSAT0_cover       <-  c( )      ## Coverage of blood-stage drugs  
MSAT0_RDT_PCR     <-  c( )      ## What diagnostic tool: 1 for RDT (=LM); 2 for PCR
MSAT0_sens        <-  c( )      ## Sensitivity of diagnostic tool
MSAT0_CQ_eff      <-  c( )      ## Efficacy of blood-stage drugs 
MSAT0_CQ_proph    <-  c( )      ## Duration of blood-stage prophylaxis


###############################################################
###############################################################
##     ##  Mass drug administration with chloroquine
##  9  ##  and primaquine.
##     ##  Includes information for coverage, efficacy and
#########  duration of prophylaxis.

MSAT1_years           <-  c( )      ## Time of MDA programme
MSAT1_cover           <-  c( )      ## Coverage of chloroquine  
MSAT1_RDT_PCR         <-  c( )      ## What diagnostic tool: 1 for RDT (=LM); 2 for PCR
MSAT1_sens            <-  c( )      ## Sensitivity of diagnostic tool
MSAT1_CQ_eff          <-  c( )      ## Efficacy of chloroquine (on its own) 
MSAT1_CQ_eff_wPQ      <-  c( )      ## Efficacy of chloroquine (co-administered with PQ) 
MSAT1_CQ_proph        <-  c( )      ## Duration of chloroquine prophylaxis
MSAT1_PQ_eff          <-  c( )      ## Efficacy of primaquine
MSAT1_PQ_proph        <-  c( )      ## Duration of primaquine prophylaxis
MSAT1_PQ_adhere       <-  c( )      ## Primaquine adherence
MSAT1_PQ_lowage       <-  c( )      ## Lower age limit for primaquine treatment (in days)
MSAT1_PQ_G6PD_risk    <-  c( )      ## Risk in G6PD - deficient individuals
MSAT1_PQ_CYP2D6_risk  <-  c( )      ## Risk of not working in low CYP2D6 metabolizers
MSAT1_PQ_preg_risk    <-  c( )      ## Risk in pregnant women
MSAT1_G6PD_test       <-  c( )      ## G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


###############################################################
###############################################################
##      ##  Mass drug administration with chloroquine
##  10  ##  and tafenoquine.
##      ##  Includes information for coverage, efficacy and
##########  duration of prophylaxis.

MSAT2_years           <-  c( )      ## Time of MDA programme
MSAT2_cover           <-  c( )      ## Coverage  
MSAT2_RDT_PCR         <-  c( )      ## What diagnostic tool: 1 for RDT (=LM); 2 for PCR
MSAT2_sens            <-  c( )      ## Sensitivity of diagnostic tool
MSAT2_CQ_eff          <-  c( )      ## Efficacy of chloroquine (on its own) 
MSAT2_CQ_eff_wPQ      <-  c( )      ## Efficacy of chloroquine (co-administered with PQ)
MSAT2_CQ_proph        <-  c( )      ## Duration of chloroquine prophylaxis
MSAT2_PQ_eff          <-  c( )      ## Efficacy of primaquine
MSAT2_PQ_proph        <-  c( )      ## Duration of primaquine prophylaxis
MSAT2_PQ_adhere       <-  c( )      ## Primaquine adherence
MSAT2_PQ_lowage       <-  c( )      ## Lower age limit for primaquine treatment (in days)
MSAT2_PQ_G6PD_risk    <-  c( )      ## Risk in G6PD - deficient individuals
MSAT2_PQ_CYP2D6_risk  <-  c( )      ## Risk of not working in low CYP2D6 metabolizers
MSAT2_PQ_preg_risk    <-  c( )      ## Risk in pregnant women
MSAT2_TQ_eff          <-  c( )      ## Efficacy of tafenoquine
MSAT2_TQ_proph        <-  c( )      ## Duration of tafenoquine prophylaxis
MSAT2_TQ_adhere       <-  c( )      ## Tafenoquine adherence
MSAT2_TQ_lowage       <-  c( )      ## Lower age limit for tafenoquine treatment (in days)
MSAT2_TQ_G6PD_risk    <-  c( )      ## Risk in G6PD - deficient individuals
MSAT2_TQ_CYP2D6_risk  <-  c( )      ## Risk of not working in low CYP2D6 metabolizers
MSAT2_TQ_preg_risk    <-  c( )      ## Risk in pregnant women
MSAT2_G6PD_test       <-  c( )      ## G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


###############################################################
###############################################################
##      ##  Serological testing and treatment with chloroquine
##  11  ##  and a hypnozoiticidal drug (primaquine).
##      ##  Includes information for coverage, efficacy and
##########  duration of prophylaxis

STAT1_years           <- c( 2020 )       ## Time of MDA programme
STAT1_cover           <- c( 0.8 )        ## Coverage
STAT1_sens            <- c( 0.8 )        ## Sensitivity of diagnostic tool
STAT1_spec            <- c( 0.8 )        ## Specificity of diagnostic tool
STAT1_RDT_PCR         <- c( 1 )          ## What diagnostic tool (extra parasitological testing)): 0 = nothing; 1 for RDT(= LM); 2 for PCR
STAT1_CQ_eff          <- c( 0.899 )      ## Efficacy of chloroquine (on its own)
STAT1_CQ_eff_wPQ      <- c( 0.946 )      ## Efficacy of chloroquine (co-administered with PQ)
STAT1_CQ_proph        <- c( 28 )         ## Duration of chloroquine prophylaxis
STAT1_PQ_eff          <- c( 0.686 )      ## Efficacy of primaquine
STAT1_PQ_proph        <- c( 8 )          ## Duration of primaquine prophylaxis
STAT1_PQ_adhere       <- c( 0.667 )      ## Primaquine adherence
STAT1_PQ_lowage       <- c(  180 )       ## Lower age limit for primaquine treatment (in days)
STAT1_PQ_G6PD_risk    <- c( 1 )          ## Risk in G6PD - deficient individuals
STAT1_PQ_CYP2D6_risk  <- c( 1 )          ## Risk of not working in low CYP2D6 metabolizers
STAT1_PQ_preg_risk    <- c( 1 )          ## Risk in pregnant women
STAT1_G6PD_test       <- c( 1 )          ## G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


###############################################################
###############################################################
##      ##  Serologicla testing and treatment with a chloroquine
##  12  ##  and a hypnozoiticidal drug (tafenoquine).
##      ##  Includes information for coverage, efficacy and
##########  duration of prophylaxis

STAT2_years           <- c( )        ## Time of MDA programme
STAT2_cover           <- c( )        ## Coverage of blood
STAT2_sens            <- c( )        ## Sensitivity of diagnostic tool
STAT2_spec            <- c( )        ## Specificity of diagnostic tool
STAT2_RDT_PCR         <- c( )        ## What diagnostic tool (extra parasitological testing)): 0 = nothing; 1 for RDT(= LM); 2 for PCR
STAT2_CQ_eff          <- c( )        ## Efficacy of chloroquine (on its own)
STAT2_CQ_eff_wPQ      <- c( )        ## Efficacy of chloroquine (co-administered with PQ)
STAT2_CQ_proph        <- c( )        ## Duration of chloroquine prophylaxis
STAT2_PQ_eff          <- c( )        ## Efficacy of primaquine
STAT2_PQ_proph        <- c( )        ## Duration of primaquine prophylaxis
STAT2_PQ_adhere       <- c( )        ## Primaquine adherence
STAT2_PQ_lowage       <- c( )        ## Lower age limit for primaquine treatment (in days)
STAT2_PQ_G6PD_risk    <- c( )        ## Risk in G6PD - deficient individuals
STAT2_PQ_CYP2D6_risk  <- c( )        ## Risk of not working in low CYP2D6 metabolizers
STAT2_PQ_preg_risk    <- c( )        ## Risk in pregnant women
STAT2_TQ_eff          <- c( )        ## Efficacy of tafenoquine
STAT2_TQ_proph        <- c( )        ## Duration of tafenoquine prophylaxis
STAT2_TQ_adhere       <- c( )        ## Tafenoquine adherence
STAT2_TQ_lowage       <- c( )        ## Lower age limit for tafenoquine treatment (in days)
STAT2_TQ_G6PD_risk    <- c( )        ## Risk in G6PD - deficient individuals
STAT2_TQ_CYP2D6_risk  <- c( )        ## Risk of not working in low CYP2D6 metabolizers
STAT2_TQ_preg_risk    <- c( )        ## Risk in pregnant women
STAT2_G6PD_test       <- c( )        ## G6PD testing for PQ eligible cases  (0 = no test; 1 = test)


###############################################################
###############################################################
##      ##  Ivermectin.
##  13  ##  This is a simplified model, based loosely on 
##      ##  Slater et al JID 2014. I'll need to fully match
##########  things up with Hannah's model.
 
IVM_years     <- c( )     ## Time of ivermectin
IVM_cover     <- c( )     ## Coverage of ivermectin
d_IVM_0       <- c( )     ## Proportion of mosquitoes killed at baseline 
IVM_half_life <- c( )     ## half-life of ivermectin (days)



years <- sort(unique( c( CM0_years, CM1_years, CM2_years,
                         LLIN_years, IRS_years, 
                         MDA0_years, MDA1_years, MDA2_years, 
                         MSAT0_years, MSAT1_years, MSAT2_years,
                         STAT1_years, STAT2_years,
				 IVM_years ) ))


INT_cov <- matrix(-1, nrow=length(years), ncol=151)
colnames(INT_cov) <- c( "years", 
                        "CM0_cover", "CM0_CQ_eff", "CM0_CQ_proph",  
				"CM1_cover", "CM1_CQ_eff", "CM1_CQ_eff_wPQ", "CM1_CQ_proph", 
				"CM1_PQ_eff", "CM1_PQ_proph", "CM1_PQ_adhere", "CM1_PQ_lowage", "CM1_PQ_G6PD_risk", "CM1_PQ_CYP2D6_risk", "CM1_PQ_preg_risk", "CM1_G6PD_test",      
				"CM2_cover", "CM2_CQ_eff", "CM2_CQ_eff_wPQ", "CM2_CQ_proph",
				"CM2_PQ_eff", "CM2_PQ_proph", "CM2_PQ_adhere", "CM2_PQ_lowage", "CM2_PQ_G6PD_risk", "CM2_PQ_CYP2D6_risk", "CM2_PQ_preg_risk",
				"CM2_TQ_eff", "CM2_TQ_proph", "CM2_TQ_adhere", "CM2_TQ_lowage", "CM2_TQ_G6PD_risk", "CM2_TQ_CYP2D6_risk", "CM2_TQ_preg_risk", "CM2_G6PD_test",     
				"LLIN_cover", 
                        "IRS_cover",
				"MDA0_cover", "MDA0_CQ_eff", "MDA0_CQ_proph",
				"MDA1_cover", "MDA1_CQ_eff", "MDA1_CQ_eff_wPQ", "MDA1_CQ_proph",
				"MDA1_PQ_eff", "MDA1_PQ_proph", "MDA1_PQ_adhere", "MDA1_PQ_lowage", "MDA1_PQ_G6PD_risk", "MDA1_PQ_CYP2D6_risk", "MDA1_PQ_preg_risk", "MDA1_G6PD_test",
				"MDA2_cover", "MDA2_CQ_eff", "MDA2_CQ_eff_wPQ", "MDA2_CQ_proph",
				"MDA2_PQ_eff", "MDA2_PQ_proph", "MDA2_PQ_adhere", "MDA2_PQ_lowage", "MDA2_PQ_G6PD_risk", "MDA2_PQ_CYP2D6_risk", "MDA2_PQ_preg_risk",
				"MDA2_TQ_eff", "MDA2_TQ_proph", "MDA2_TQ_adhere", "MDA2_TQ_lowage", "MDA2_TQ_G6PD_risk", "MDA2_TQ_CYP2D6_risk", "MDA2_TQ_preg_risk", "MDA2_G6PD_test",
				"MSAT0_cover", "MSAT0_RDT_PCR", "MSAT0_sens", "MSAT0_CQ_eff", "MSAT0_CQ_proph",
				"MSAT1_cover", "MSAT1_RDT_PCR", "MSAT1_sens", "MSAT1_CQ_eff", "MSAT1_CQ_eff_wPQ", "MSAT1_CQ_proph",
				"MSAT1_PQ_eff", "MSAT1_PQ_proph", "MSAT1_PQ_adhere", "MSAT1_PQ_lowage", "MSAT1_PQ_G6PD_risk", "MSAT1_PQ_CYP2D6_risk", "MSAT1_PQ_preg_risk", "MSAT1_G6PD_test",
				"MSAT2_cover", "MSAT2_RDT_PCR", "MSAT2_sens", "MSAT2_CQ_eff", "MSAT2_CQ_eff_wPQ", "MSAT2_CQ_proph",
				"MSAT2_PQ_eff", "MSAT2_PQ_proph", "MSAT2_PQ_adhere", "MSAT2_PQ_lowage", "MSAT2_PQ_G6PD_risk", "MSAT2_PQ_CYP2D6_risk", "MSAT2_PQ_preg_risk",
				"MSAT2_TQ_eff", "MSAT2_TQ_proph", "MSAT2_TQ_adhere", "MSAT2_TQ_lowage", "MSAT2_TQ_G6PD_risk", "MSAT2_TQ_CYP2D6_risk", "MSAT2_TQ_preg_risk", "MSAT2_G6PD_test",
				"STAT1_cover", "STAT1_sens", "STAT1_spec", "STAT1_RDT_PCR", "STAT1_CQ_eff", "STAT1_CQ_eff_wPQ", "STAT1_CQ_proph",
				"STAT1_PQ_eff", "STAT1_PQ_proph", "STAT1_PQ_adhere", "STAT1_PQ_lowage", "STAT1_PQ_G6PD_risk", "STAT1_PQ_CYP2D6_risk", "STAT1_PQ_preg_risk", "STAT1_G6PD_test",
				"STAT2_cover", "STAT2_sens", "STAT2_spec", "STAT2_RDT_PCR", "STAT2_CQ_eff", "STAT2_CQ_eff_wPQ", "STAT2_CQ_proph",
				"STAT2_PQ_eff", "STAT2_PQ_proph", "STAT2_PQ_adhere", "STAT2_PQ_lowage", "STAT2_PQ_G6PD_risk", "STAT2_PQ_CYP2D6_risk", "STAT2_PQ_preg_risk",
				"STAT2_TQ_eff", "STAT2_TQ_proph", "STAT2_TQ_adhere", "STAT2_TQ_lowage", "STAT2_TQ_G6PD_risk", "STAT2_TQ_CYP2D6_risk", "STAT2_TQ_preg_risk", "STAT2_G6PD_test",
 	      		"IVM_cover", "d_IVM_0", "IVM_half_life" )


INT_cov[,1] = years

for(i in 1:nrow(INT_cov))
{
	if( INT_cov[i,1] %in% CM0_years )
	{
		INT_cov[i,2] = CM0_cover[which(CM0_years==INT_cov[i,1])]
		INT_cov[i,3] = CM0_CQ_eff[which(CM0_years==INT_cov[i,1])]
		INT_cov[i,4] = CM0_CQ_proph[which(CM0_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% CM1_years )
	{
		INT_cov[i,5]  = CM1_cover[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,6]  = CM1_CQ_eff[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,7]  = CM1_CQ_eff_wPQ[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,8]  = CM1_CQ_proph[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,9]  = CM1_PQ_eff[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,10] = CM1_PQ_proph[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,11] = CM1_PQ_adhere[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,12] = CM1_PQ_lowage[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,13] = CM1_PQ_G6PD_risk[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,14] = CM1_PQ_CYP2D6_risk[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,15] = CM1_PQ_preg_risk[which(CM1_years==INT_cov[i,1])]
		INT_cov[i,16] = CM1_G6PD_test[which(CM1_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% CM2_years )
	{
		INT_cov[i,17] = CM2_cover[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,18] = CM2_CQ_eff[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,19] = CM2_CQ_eff_wPQ[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,20] = CM2_CQ_proph[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,21] = CM2_PQ_eff[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,22] = CM2_PQ_proph[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,23] = CM2_PQ_adhere[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,24] = CM2_PQ_lowage[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,25] = CM2_PQ_G6PD_risk[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,26] = CM2_PQ_CYP2D6_risk[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,27] = CM2_PQ_preg_risk[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,28] = CM2_TQ_eff[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,29] = CM2_TQ_proph[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,30] = CM2_TQ_adhere[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,31] = CM2_TQ_lowage[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,32] = CM2_TQ_G6PD_risk[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,33] = CM2_TQ_CYP2D6_risk[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,34] = CM2_TQ_preg_risk[which(CM2_years==INT_cov[i,1])]
		INT_cov[i,35] = CM2_G6PD_test[which(CM2_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% LLIN_years )
	{
		INT_cov[i,36] = LLIN_cover[which(LLIN_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% IRS_years )
	{
		INT_cov[i,37] = IRS_cover[which(IRS_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% MDA0_years )
	{
		INT_cov[i,38] = MDA0_cover[which(MDA0_years==INT_cov[i,1])]
		INT_cov[i,39] = MDA0_CQ_eff[which(MDA0_years==INT_cov[i,1])]
		INT_cov[i,40] = MDA0_CQ_proph[which(MDA0_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% MDA1_years )
	{
		INT_cov[i,41] = MDA1_cover[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,42] = MDA1_CQ_eff[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,43] = MDA1_CQ_eff_wPQ[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,44] = MDA1_CQ_proph[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,45] = MDA1_PQ_eff[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,46] = MDA1_PQ_proph[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,47] = MDA1_PQ_adhere[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,48] = MDA1_PQ_lowage[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,49] = MDA1_PQ_G6PD_risk[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,50] = MDA1_PQ_CYP2D6_risk[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,51] = MDA1_PQ_preg_risk[which(MDA1_years==INT_cov[i,1])]
		INT_cov[i,52] = MDA1_G6PD_test[which(MDA1_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% MDA2_years )
	{
		INT_cov[i,53] = MDA2_cover[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,54] = MDA2_CQ_eff[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,55] = MDA2_CQ_eff_wPQ[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,56] = MDA2_CQ_proph[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,57] = MDA2_PQ_eff[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,58] = MDA2_PQ_proph[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,59] = MDA2_PQ_adhere[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,60] = MDA2_PQ_lowage[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,61] = MDA2_PQ_G6PD_risk[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,62] = MDA2_PQ_CYP2D6_risk[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,63] = MDA2_PQ_preg_risk[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,64] = MDA2_TQ_eff[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,65] = MDA2_TQ_proph[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,66] = MDA2_TQ_adhere[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,67] = MDA2_TQ_lowage[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,68] = MDA2_TQ_G6PD_risk[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,69] = MDA2_TQ_CYP2D6_risk[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,70] = MDA2_TQ_preg_risk[which(MDA2_years==INT_cov[i,1])]
		INT_cov[i,71] = MDA2_G6PD_test[which(MDA2_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% MSAT0_years )
	{
		INT_cov[i,72] = MSAT0_cover[which(MSAT0_years==INT_cov[i,1])]
		INT_cov[i,73] = MSAT0_RDT_PCR[which(MSAT0_years==INT_cov[i,1])]
		INT_cov[i,74] = MSAT0_sens[which(MSAT0_years==INT_cov[i,1])]
		INT_cov[i,75] = MSAT0_CQ_eff[which(MSAT0_years==INT_cov[i,1])]		
		INT_cov[i,76] = MSAT0_CQ_proph[which(MSAT0_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% MSAT1_years )
	{
		INT_cov[i,77] = MSAT1_cover[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,78] = MSAT1_RDT_PCR[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,79] = MSAT1_sens[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,80] = MSAT1_CQ_eff[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,81] = MSAT1_CQ_eff_wPQ[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,82] = MSAT1_CQ_proph[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,83] = MSAT1_PQ_eff[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,84] = MSAT1_PQ_proph[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,85] = MSAT1_PQ_adhere[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,86] = MSAT1_PQ_lowage[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,87] = MSAT1_PQ_G6PD_risk[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,88] = MSAT1_PQ_CYP2D6_risk[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,89] = MSAT1_PQ_preg_risk[which(MSAT1_years==INT_cov[i,1])]
		INT_cov[i,90] = MSAT1_G6PD_test[which(MSAT1_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% MSAT2_years )
	{
		INT_cov[i,91]  = MSAT2_cover[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,92]  = MSAT2_RDT_PCR[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,93]  = MSAT2_sens[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,94]  = MSAT2_CQ_eff[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,95]  = MSAT2_CQ_eff_wPQ[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,96]  = MSAT2_CQ_proph[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,97]  = MSAT2_PQ_eff[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,98]  = MSAT2_PQ_proph[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,99] = MSAT2_PQ_adhere[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,100] = MSAT2_PQ_lowage[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,101] = MSAT2_PQ_G6PD_risk[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,102] = MSAT2_PQ_CYP2D6_risk[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,103] = MSAT2_PQ_preg_risk[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,104] = MSAT2_TQ_eff[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,105] = MSAT2_TQ_proph[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,106] = MSAT2_TQ_adhere[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,107] = MSAT2_TQ_lowage[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,108] = MSAT2_TQ_G6PD_risk[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,109] = MSAT2_TQ_CYP2D6_risk[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,110] = MSAT2_TQ_preg_risk[which(MSAT2_years==INT_cov[i,1])]
		INT_cov[i,111] = MSAT2_G6PD_test[which(MSAT2_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% STAT1_years )
	{
		INT_cov[i,112] = STAT1_cover[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,113] = STAT1_sens[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,114] = STAT1_spec[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,115] = STAT1_RDT_PCR[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,116] = STAT1_CQ_eff[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,117] = STAT1_CQ_eff_wPQ[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,118] = STAT1_CQ_proph[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,119] = STAT1_PQ_eff[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,120] = STAT1_PQ_proph[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,121] = STAT1_PQ_adhere[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,122] = STAT1_PQ_lowage[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,123] = STAT1_PQ_G6PD_risk[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,124] = STAT1_PQ_CYP2D6_risk[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,125] = STAT1_PQ_preg_risk[which(STAT1_years==INT_cov[i,1])]
		INT_cov[i,126] = STAT1_G6PD_test[which(STAT1_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% STAT2_years )
	{
		INT_cov[i,127] = STAT2_cover[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,128] = STAT2_sens[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,129] = STAT2_spec[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,130] = STAT2_RDT_PCR[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,131] = STAT2_CQ_eff[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,132] = STAT2_CQ_eff_wPQ[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,133] = STAT2_CQ_proph[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,134] = STAT2_PQ_eff[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,135] = STAT2_PQ_proph[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,136] = STAT2_PQ_adhere[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,137] = STAT2_PQ_lowage[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,138] = STAT2_PQ_G6PD_risk[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,139] = STAT2_PQ_CYP2D6_risk[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,140] = STAT2_PQ_preg_risk[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,141] = STAT2_TQ_eff[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,142] = STAT2_TQ_proph[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,143] = STAT2_TQ_adhere[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,144] = STAT2_TQ_lowage[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,145] = STAT2_TQ_G6PD_risk[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,146] = STAT2_TQ_CYP2D6_risk[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,147] = STAT2_TQ_preg_risk[which(STAT2_years==INT_cov[i,1])]
		INT_cov[i,148] = STAT2_G6PD_test[which(STAT2_years==INT_cov[i,1])]
	}

	if( INT_cov[i,1] %in% IVM_years )
	{
		INT_cov[i,149]  = IVM_cover[which(IVM_years==INT_cov[i,1])]
		INT_cov[i,150] = d_IVM_0[which(IVM_years==INT_cov[i,1])]
		INT_cov[i,151] = IVM_half_life[which(IVM_years==INT_cov[i,1])]
	}
}
	


INT_cov <- rbind( INT_cov, rep(-1,ncol(INT_cov)) )
INT_cov <- t(INT_cov)

write.table(INT_cov, file="intervention_coverage.txt", col.names=FALSE)



