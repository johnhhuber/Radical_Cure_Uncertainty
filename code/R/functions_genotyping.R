# install necessary functions 
require(icenReg)

# function to reshape the data into format for calculating efficacy at varying genotyping accuracies
createDataFrame <- function(output_participant, # output for the participant information 
                            output_recurrent, # output for the recurrent infections 
                            output_trial, # output for the trial endpoints 
                            genotyping_sens_spec) # different combinations of sensitivity and specificity of genotyping 
{
  # create data frame 
  n_participants <- nrow(output_participant)
  # create data frame 
  n_participants <- nrow(output_participant)
  df <- data.frame(participant_id = rep(NA, n_participants),
                   trial_arm = rep(NA, n_participants),
                   days_enrolled = rep(NA, n_participants))
  
  # add in columns for the different sensitivity and specificity combinations
  df[paste('num_relapse_PCR_', 'SENS_', genotyping_sens_spec[,'sensitivity'] * 100, '_SPEC_', genotyping_sens_spec[,'specificity'] * 100, sep = '')] <- NA
  df[paste('num_relapse_LM_', 'SENS_', genotyping_sens_spec[,'sensitivity'] * 100, '_SPEC_', genotyping_sens_spec[,'specificity'] * 100, sep = '')] <- NA
  df[paste('num_relapse_D_', 'SENS_', genotyping_sens_spec[,'sensitivity'] * 100, '_SPEC_', genotyping_sens_spec[,'specificity'] * 100, sep = '')] <- NA
  
  df[paste('time_LM_', 'SENS_', genotyping_sens_spec[,'sensitivity'] * 100, '_SPEC_', genotyping_sens_spec[,'specificity'] * 100, sep = '')] <- NA
  
  # add in particpant information 
  df$participant_id <- output_participant$Participant_ID
  df$trial_arm <- output_participant$Trial_Arm
  df$days_enrolled <- output_participant$Dropout_Date - output_participant$Enrollment_Date
  
  # loop through the different different sensitivity and specificity combos and calculate the number of genotyped relapse for each participant
  for(ss in 1:nrow(genotyping_sens_spec))
  {
    # simulate genotyping process 
    genotyped_relapses <- sapply(df$participant_id, function(x){countGenotypedRelapses(output_recurrent = output_recurrent,
                                                                output_participant = output_participant,
                                                                participant_id = x,
                                                                sensitivity = genotyping_sens_spec[ss, 'sensitivity'],
                                                                specificity = genotyping_sens_spec[ss, 'specificity'])})
    
    # PCR-based detection 
    col_pcr <- paste('num_relapse_PCR_', 'SENS_', genotyping_sens_spec[ss,'sensitivity'] * 100, '_SPEC_', genotyping_sens_spec[ss,'specificity'] * 100, sep = '')
    df[,col_pcr] <- genotyped_relapses['PCR',]
    
    # LM-based detection
    col_pcr <- paste('num_relapse_LM_', 'SENS_', genotyping_sens_spec[ss,'sensitivity'] * 100, '_SPEC_', genotyping_sens_spec[ss,'specificity'] * 100, sep = '')
    df[,col_pcr] <- genotyped_relapses['LM',]
    
    # D-based detection
    col_pcr <- paste('num_relapse_D_', 'SENS_', genotyping_sens_spec[ss,'sensitivity'] * 100, '_SPEC_', genotyping_sens_spec[ss,'specificity'] * 100, sep = '')
    df[,col_pcr] <- genotyped_relapses['D',]
    
    # get the time of first genotyped recurrent infection 
    time_genotyped_recurrent <- sapply(df$participant_id, function(x){getTimeFirstGenotyped(output_recurrent = output_recurrent,
                                                                                            output_participant = output_participant,
                                                                                            output_trial = output_trial,
                                                                                            participant_id = x,
                                                                                            sensitivity = genotyping_sens_spec[ss, 'sensitivity'],
                                                                                            specificity = genotyping_sens_spec[ss, 'specificity'])})
    
    col_time <- paste('time_LM_', 'SENS_', genotyping_sens_spec[ss,'sensitivity'] * 100, '_SPEC_', genotyping_sens_spec[ss,'specificity'] * 100, sep = '')
    df[,col_time] <- time_genotyped_recurrent
  }
  
  # return data frame
  return(df)
}

# function to get the time of the first recurrent infection that is genotyped as a true treatment failure 
getTimeFirstGenotyped <- function(output_recurrent,
                                  output_participant,
                                  output_trial,
                                  participant_id,
                                  sensitivity,
                                  specificity,
                                  days_left_censored = 32,
                                  detection_type = 'LM')
{
  # subset to not include the period of left censoring
  enrollment_date <- output_participant$Enrollment_Date[output_participant$Participant_ID == participant_id]
  recurrent_censored <- subset(output_recurrent, Participant_ID == participant_id & ((Infection_Date - enrollment_date) > days_left_censored))
  recurrent_all <- subset(output_recurrent, Participant_ID == participant_id)
  followup <- subset(output_trial, Participant_ID == participant_id & Days_Since_Enrollment > days_left_censored)
  
  # count the number of true treatment failure relapses and the number of other recurrent infections
  is_true_relapse_all <- grepl('RELAPSE_.*_PRE_ENROLLMENT', recurrent_all$Infection_Type)
  is_true_relapse_censored <- is_true_relapse_all[which(recurrent_all$Infection_Date %in% recurrent_censored$Infection_Date)]
  
  # simulate the number of genotyped relapses and other recurrent infections 
  is_genotyped_relapse_all <- rbinom(n = length(is_true_relapse_all), size = 1, prob = ifelse(is_true_relapse_all, sensitivity, 1 - specificity))
  is_genotyped_relapse_censored <- is_genotyped_relapse_all[which(recurrent_all$Infection_Date %in% recurrent_censored$Infection_Date)]
  
  # get the timing of the first recurrent infection
  day_first_recurrent <- max(followup$Days_Since_Enrollment) + 1
  if(detection_type == 'LM')
  {
    if(any(followup$LM_Infection))
    {
      index_recurrent_infection <- tail(which(recurrent_all$Infection_Date - enrollment_date <= min(followup$Days_Since_Enrollment[followup$LM_Infection == 1])), n = 1)
      if(is_genotyped_relapse_all[index_recurrent_infection])
      {
        day_first_recurrent <- min(day_first_recurrent, min(followup$Days_Since_Enrollment[followup$LM_Infection == 1]))
      }
    }
    if(nrow(recurrent_censored) > 0)
    {
      if(any((grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)) & is_genotyped_relapse_censored))
      {
        day_first_recurrent <- min(day_first_recurrent, min(recurrent_censored$Infection_Date[(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)) & is_genotyped_relapse_censored] - enrollment_date))
      }
    }
  }
  
  if(detection_type == 'PCR')
  {
    if(any(followup$PCR_Infection))
    {
      index_recurrent_infection <- tail(which(recurrent_all$Infection_Date - enrollment_date <= min(followup$Days_Since_Enrollment[followup$PCR_Infection == 1])), n = 1)
      if(is_genotyped_relapse_all[index_recurrent_infection])
      {
        day_first_recurrent <- min(day_first_recurrent, min(followup$Days_Since_Enrollment[followup$PCR_Infection == 1]))
      }
    }
    if(nrow(recurrent_censored) > 0)
    {
      if(any((grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)) & is_genotyped_relapse_censored))
      {
        day_first_recurrent <- min(day_first_recurrent, min(recurrent_censored$Infection_Date[(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)) & is_genotyped_relapse_censored] - enrollment_date))
      }
    }
  }
  
  # return output 
  return(day_first_recurrent)
}

# function to count the number of recurrent infections of a specific type 
countGenotypedRelapses <- function(output_recurrent,
                                   output_participant,
                                   participant_id, 
                                   sensitivity,
                                   specificity)
{
  # get enrollment date 
  enrollment_date <- output_participant$Enrollment_Date[output_participant$Participant_ID == participant_id]
  
  # subset data
  output_recurrent$Infection_Type <- as.character(output_recurrent$Infection_Type)
  recurrent <- subset(output_recurrent, Participant_ID == participant_id & Infection_Date > enrollment_date)

  # count the number of true treatment failure relapses and the number of other recurrent infections
  is_true_relapse <- grepl('RELAPSE_.*_PRE_ENROLLMENT', recurrent$Infection_Type)

  # simulate the number of genotyped relapses and other recurrent infections 
  is_genotyped_relapse <- rbinom(n = length(is_true_relapse), size = 1, prob = ifelse(is_true_relapse, sensitivity, 1 - specificity))
  
  # stratify by detection type 
  num_relapses_genotyped_PCR = sum(is_genotyped_relapse[grepl('_PCR', recurrent$Infection_Type)])
  num_relapses_genotyped_LM = sum(is_genotyped_relapse[grepl('_LM', recurrent$Infection_Type)])
  num_relapses_genotyped_D = sum(is_genotyped_relapse[grepl('_D', recurrent$Infection_Type)])
  
  # generate counts 
  num_relapses_genotyped <- c(PCR = num_relapses_genotyped_PCR + num_relapses_genotyped_LM + num_relapses_genotyped_D,
                              LM = num_relapses_genotyped_LM + num_relapses_genotyped_D,
                              D = num_relapses_genotyped_D)
  
  # return count
  return(num_relapses_genotyped)
}

# function to compute the incidence efficacy measures 
calcIncidEff <- function(df)
{
  # get the appropriate columns to calculate efficacy for 
  incid_cols <- grepl('num_', colnames(df))
  
  # calculate the number of recurrent events in treatment and placebo arms 
  num_events_treatment <- colSums(df[df$trial_arm == 'TREATMENT',incid_cols])
  num_events_placebo <- colSums(df[df$trial_arm == 'PLACEBO', incid_cols])
  
  # calculate the number of person days in treatment and placebo arms 
  person_days_treatment <- sum(df$days_enrolled[df$trial_arm == 'TREATMENT'])
  person_days_placebo <- sum(df$days_enrolled[df$trial_arm == 'PLACEBO'])
  
  # calculate eff and return 
  RR <- (num_events_treatment / person_days_treatment) / (num_events_placebo / person_days_placebo)
  eff <- 1 - RR
  return(eff)
}

# fucntion to compute the efficacy based on Cox Proportional Hazards model
calcCoxEff <- function(df,
					   days_followup,
					   days_left_censored = 32)
{
  # get the appropriate columns to calculate efficay for 
  time_cols <- which(grepl('time_', colnames(df)))
  
  # create storage vector for efficacy estimates 
  eff_cph <- rep(NA, length(time_cols))
  
  # calculate efficacy for each genotyping scenario 
  for(ee in 1:length(time_cols))
  {
    # data organization 
    df_subset <- subset(df, df[,time_cols[ee]] >= 0)
    df_subset$censor_recurrent_LM <- 0
    df_subset$censor_recurrent_LM[df_subset[,time_cols[ee]] == (df_subset$days_enrolled + 1)] <- 1
    df_subset[df$censor_recurrent_LM == 1, time_cols[ee]] = df_subset[df$censor_recurrent_LM == 1, time_cols[ee]] - 1 
    df_subset$status <- (!df_subset$censor_recurrent_LM) * 1 + 1
    df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
    
    # get the updated number of days of follow-up in light of left censoring 
  	days_followup <- c(0, days_followup[days_followup >= days_left_censored])

  	# specify the lower and upper bounds for each interval 
  	df_subset$time_upper_LM <- sapply(1:nrow(df_subset), function(i){ifelse(df_subset$censor_recurrent_LM[i], Inf, df_subset[i,time_cols[ee]])})
  	df_subset$time_lower_LM <- sapply(1:nrow(df_subset), function(i){ifelse(df_subset$censor_recurrent_LM[i], df_subset[i,time_cols[ee]], tail(days_followup[days_followup < df_subset[i,time_cols[ee]]], n = 1))})

    # fit cox proportional hazard
    #cph <- coxph(formula = Surv(df_subset[,time_cols[ee]], status) ~ trial_arm_num, data = df_subset)
    cph <- ic_sp(Surv(time = df_subset$time_lower_LM, time2 = df_subset$time_upper_LM, type = 'interval2') ~ trial_arm_num, data = df_subset, model = 'ph')
    eff_cph[ee] <- 1 - exp(as.numeric(cph$coefficients))
  }
  
  names(eff_cph) <- colnames(df)[time_cols]
  
  # return output
  return(eff_cph)
}
