# install necessary functions 
require(icenReg)

# function to count the number of recurrent infections of a specific type 
countRecurrentInfs <- function(output_recurrent,
                               output_participant,
                               participant_id, 
                               inf_types)
{
  # get enrollment date 
  enrollment_date <- output_participant$Enrollment_Date[output_participant$Participant_ID == participant_id]
  
  # subset data 
  output_recurrent$Infection_Type <- as.character(output_recurrent$Infection_Type)
  recurrent <- subset(output_recurrent, Participant_ID == participant_id & Infection_Type %in% inf_types & Infection_Date > enrollment_date)
  
  # return count
  return(nrow(recurrent))
}

# function to count the number of recurrent infections within the context of a trial 
countRecurrentInfsTrial <- function(output_recurrent,
                                    output_trial,
                                    output_participant,
                                    participant_id,
                                    days_left_censored = 32,
                                    detection_type = 'LM')
{
  # subset to not left censory
  enrollment_date <- output_participant$Enrollment_Date[output_participant$Participant_ID == participant_id]
  recurrent <- subset(output_recurrent, Participant_ID == participant_id & ((Infection_Date - enrollment_date) > days_left_censored))
  followup <- subset(output_trial, Participant_ID == participant_id & Days_Since_Enrollment > days_left_censored)

  # count the number of recurrent infections 
  if(detection_type == 'D')
  {
    num_recurrent <- sum(grepl('_D', recurrent$Infection_Type) | grepl('_T', recurrent$Infection_Type))
  }
  
  if(detection_type == 'LM')
  {
    # first count the number of times that they tested postive on followup 
    num_positive_followup <- sum(followup$LM_Infection)

    # calculate the number of symptomatic/treated infections independent of testing during follow-up. In doing so, we avoid double-counting
    num_symptomatic <-  sum(grepl('_D', recurrent$Infection_Type) | grepl('_T', recurrent$Infection_Type))
    num_overlap <- 0
    
    if(any(followup$LM_Infection))
    {
      dates_positive_followup <- followup$Days_Since_Enrollment[followup$LM_Infection == 1] + enrollment_date
      infection_type_followup <- sapply(dates_positive_followup, function(d){recurrent$Infection_Type[tail(which(recurrent$Infection_Date < d), n = 1)]})
      infection_type_followup <- unlist(lapply(infection_type_followup, function(l){ifelse(length(l) == 0, NA, l)}))
      
      num_overlap <- sum(grepl('_D', infection_type_followup) | grepl('_T', infection_type_followup))
    }
    
    # get the union of the different events 
    num_recurrent <- num_positive_followup + num_symptomatic - num_overlap
  }
  
  if(detection_type == 'PCR')
  {
    # first count the number of times that they tested postive on followup 
    num_positive_followup <- sum(followup$PCR_Infection)
    
    # calculate the number of symptomatic/treated infections independent of testing during follow-up. In doing so, we avoid double-counting
    num_symptomatic <-  sum(grepl('_D', recurrent$Infection_Type) | grepl('_T', recurrent$Infection_Type))
    num_overlap <- 0
    
    if(any(followup$PCR_Infection))
    {
      dates_positive_followup <- followup$Days_Since_Enrollment[followup$PCR_Infection == 1] + enrollment_date
      infection_type_followup <- sapply(dates_positive_followup, function(d){recurrent$Infection_Type[tail(which(recurrent$Infection_Date < d), n = 1)]})
      infection_type_followup <- unlist(lapply(infection_type_followup, function(l){ifelse(length(l) == 0, NA, l)}))
      
      num_overlap <- sum(grepl('_D', infection_type_followup) | grepl('_T', infection_type_followup))
    }
    
    # get the union of the different events 
    num_recurrent <- num_positive_followup + num_symptomatic - num_overlap
  }
  
  # return output 
  return(num_recurrent)
}

# function to determine the timing of the first recurrent infection 
getTimeFirstRecurrent <- function(output_recurrent,
                                  output_trial,
                                  output_participant,
                                  participant_id,
                                  days_left_censored = 32,
                                  detection_type = 'LM')
{
  
  # subset to not include the first 29 days 
  enrollment_date <- output_participant$Enrollment_Date[output_participant$Participant_ID == participant_id]
  recurrent <- subset(output_recurrent, Participant_ID == participant_id & ((Infection_Date - enrollment_date) > days_left_censored))
  followup <- subset(output_trial, Participant_ID == participant_id & Days_Since_Enrollment > days_left_censored)
  
  # get the timing of the first recurrent infection
  day_first_recurrent <- max(followup$Days_Since_Enrollment) + 1
  if(detection_type == 'LM')
  {
    if(any(followup$LM_Infection))
    {
      day_first_recurrent <- min(day_first_recurrent, min(followup$Days_Since_Enrollment[followup$LM_Infection == 1]))
    }
    if(any(grepl('_D', recurrent$Infection_Type) | grepl('_T', recurrent$Infection_Type)))
    {
      day_first_recurrent <- min(day_first_recurrent, min(recurrent$Infection_Date[grepl('_D', recurrent$Infection_Type) | grepl('_T', recurrent$Infection_Type)] - enrollment_date))
    }
  }
  
  if(detection_type == 'PCR')
  {
    if(any(followup$PCR_Infection))
    {
      day_first_recurrent <- min(day_first_recurrent, min(followup$Days_Since_Enrollment[followup$PCR_Infection == 1]))
    }
    if(any(grepl('_D', recurrent$Infection_Type) | grepl('_T', recurrent$Infection_Type)))
    {
      day_first_recurrent <- min(day_first_recurrent, min(recurrent$Infection_Date[grepl('_D', recurrent$Infection_Type) | grepl('_T', recurrent$Infection_Type)] - enrollment_date))
    }
  }

  if(detection_type == 'D')
  {
    if(any(grepl('_D', recurrent$Infection_Type) | grepl('_T', recurrent$Infection_Type)))
    {
      day_first_recurrent <- min(day_first_recurrent, min(recurrent$Infection_Date[grepl('_D', recurrent$Infection_Type) | grepl('_T', recurrent$Infection_Type)] - enrollment_date))
    }
  }
  
  # return output 
  return(day_first_recurrent)
}

# function to get the timing of the first relapse infection 
getTimeFirstRelapseAll <- function(output_recurrent,
                                   output_trial,
                                   output_participant,
                                   participant_id,
                                   days_left_censored = 5,
                                   detection_type = 'LM')
{
  # subset to not include the first 29 days 
  enrollment_date <- output_participant$Enrollment_Date[output_participant$Participant_ID == participant_id]
  recurrent <- subset(output_recurrent, Participant_ID == participant_id & Infection_Date > enrollment_date)
  followup <- subset(output_trial, Participant_ID == participant_id & Days_Since_Enrollment > days_left_censored)
  
  # get the timing of the first recurrent infection
  day_first_relapse <- max(followup$Days_Since_Enrollment) + 1
  if(detection_type == 'LM')
  {
    if(nrow(recurrent) > 0)
    {
      day_first_relapse <- min(day_first_relapse, min(recurrent$Infection_Date[(grepl('_LM', recurrent$Infection_Type) | grepl('_D', recurrent$Infection_Type) | grepl('_T', recurrent$Infection_Type)) & grepl('RELAPSE_.*_PRE_ENROLLMENT', recurrent$Infection_Type)] - enrollment_date))
    }
  }
  
  if(detection_type == 'PCR')
  {
    if(nrow(recurrent) > 0)
    {
      day_first_relapse <- min(day_first_relapse, min(recurrent$Infection_Date[grepl('RELAPSE_.*_PRE_ENROLLMENT', recurrent$Infection_Type)] - enrollment_date))
    }
  }
  
  # return output 
  return(day_first_relapse)
}

# function to determine the timing of the first relapse infection 
getTimeFirstRelapseTrial <- function(output_recurrent,
                                output_trial,
                                output_participant,
                                participant_id,
                                days_left_censored = 32,
                                detection_type = 'LM',
                                relapse_type = 'ANY')
{
  
  # subset to not include the first 29 days 
  enrollment_date <- output_participant$Enrollment_Date[output_participant$Participant_ID == participant_id]
  recurrent_censored <- subset(output_recurrent, Participant_ID == participant_id & ((Infection_Date - enrollment_date) > days_left_censored))
  recurrent <- subset(output_recurrent, Participant_ID == participant_id)
  followup <- subset(output_trial, Participant_ID == participant_id & Days_Since_Enrollment > days_left_censored)
  
  # get the timing of the first recurrent infection
  day_first_relapse <- max(followup$Days_Since_Enrollment) + 1
  if(detection_type == 'D')
  {
    if(any(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)))
    {
      if(relapse_type == 'ANY')
      {
        day_first_relapse <- min(day_first_relapse, min(recurrent_censored$Infection_Date[(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)) & grepl('RELAPSE', recurrent_censored$Infection_Type)] - enrollment_date))
      }
      if(relapse_type == 'PRE')
      {
        day_first_relapse <- min(day_first_relapse, min(recurrent_censored$Infection_Date[(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)) & grepl('RELAPSE_.*_PRE_ENROLLMENT', recurrent_censored$Infection_Type)] - enrollment_date))
      }
    }
  }

  if(detection_type == 'LM')
  {
    if(any(followup$LM_Infection))
    {
      min_day <- min(followup$Days_Since_Enrollment[followup$LM_Infection == 1])
      days_recurrent <- recurrent$Infection_Date - enrollment_date
      index_last_recurrent <- tail(which(days_recurrent <= min_day), n = 1)
      
      if(relapse_type == 'ANY')
      {
        if(grepl('RELAPSE', recurrent$Infection_Type[index_last_recurrent]))
        {
          day_first_relapse <- min(day_first_relapse, min_day)
        }
      }
      
      if(relapse_type == 'PRE')
      {
        if(grepl('RELAPSE_.*_PRE_ENROLLMENT', recurrent$Infection_Type[index_last_recurrent]))
        {
          day_first_relapse <- min(day_first_relapse, min_day)
        }
      }
    }
    if(any(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)))
    {
      if(relapse_type == 'ANY')
      {
        day_first_relapse <- min(day_first_relapse, min(recurrent_censored$Infection_Date[(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)) & grepl('RELAPSE', recurrent_censored$Infection_Type)] - enrollment_date))
      }
      if(relapse_type == 'PRE')
      {
        day_first_relapse <- min(day_first_relapse, min(recurrent_censored$Infection_Date[(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)) & grepl('RELAPSE_.*_PRE_ENROLLMENT', recurrent_censored$Infection_Type)] - enrollment_date))
      }
    }
  }
  
  if(detection_type == 'PCR')
  {
    if(any(followup$PCR_Infection))
    {
      min_day <- min(followup$Days_Since_Enrollment[followup$PCR_Infection == 1])
      days_recurrent <- recurrent$Infection_Date - enrollment_date
      index_last_recurrent <- tail(which(days_recurrent <= min_day), n = 1)
      
      if(relapse_type == 'ANY')
      {
        if(grepl('RELAPSE', recurrent$Infection_Type[index_last_recurrent]))
        {
          day_first_relapse <- min(day_first_relapse, min_day)
        }
      }
      
      if(relapse_type == 'PRE')
      {
        if(grepl('RELAPSE_.*_PRE_ENROLLMENT', recurrent$Infection_Type[index_last_recurrent]))
        {
          day_first_relapse <- min(day_first_relapse, min_day)
        }
      }
    }
    if(any(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)))
    {
      if(relapse_type == 'ANY')
      {
        day_first_relapse <- min(day_first_relapse, min(recurrent_censored$Infection_Date[(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)) & grepl('RELAPSE', recurrent_censored$Infection_Type)] - enrollment_date))
      }
      if(relapse_type == 'PRE')
      {
        day_first_relapse <- min(day_first_relapse, min(recurrent_censored$Infection_Date[(grepl('_D', recurrent_censored$Infection_Type) | grepl('_T', recurrent_censored$Infection_Type)) & grepl('RELAPSE_.*_PRE_ENROLLMENT', recurrent_censored$Infection_Type)] - enrollment_date))
      }
    }
  }
  
  # return output 
  return(day_first_relapse)
}


# function to reshape data into format for calculating efficacy 
createEffDataFrame <- function(output_participant, # output for the participant information
                               output_recurrent, # output for recurrent infections 
                               output_trial) # output for the trial endpoints
{
  # create data frame 
  n_participants <- nrow(output_participant)
  df <- data.frame(participant_id = rep(NA, n_participants),
                   trial_arm = rep(NA, n_participants),
                   days_enrolled = rep(NA, n_participants),
                   time_recurrent_LM = rep(NA, n_participants),
                   time_recurrent_PCR = rep(NA, n_participants),
                   time_recurrent_D = rep(NA, n_participants),
                   time_any_relapse_LM = rep(NA, n_participants),
                   time_any_relapse_PCR = rep(NA, n_participants),
                   time_relapse_LM_trial = rep(NA, n_participants),
                   time_relapse_PCR_trial = rep(NA, n_participants),
                   time_relapse_D_trial = rep(NA, n_participants),
                   time_relapse_LM_all = rep(NA, n_participants),
                   time_relapse_PCR_all = rep(NA, n_participants),
                   time_relapse_D_all = rep(NA, n_participants),
                   censor_recurrent_LM = rep(NA, n_participants),
                   censor_recurrent_PCR = rep(NA, n_participants),
                   censor_any_relapse_LM = rep(NA, n_participants),
                   censor_any_relapse_PCR = rep(NA, n_participants),
                   censor_relapse_LM_trial = rep(NA, n_participants),
                   censor_relapse_PCR_trial = rep(NA, n_participants),
                   censor_relapse_D_trial = rep(NA, n_participants),
                   censor_relapse_LM_all = rep(NA, n_participants),
                   censor_relapse_PCR_all = rep(NA, n_participants),
                   censor_relapse_D_all = rep(NA, n_participants),
                   num_any_PCR = rep(NA, n_participants),
                   num_any_LM = rep(NA, n_participants),
                   num_any_D = rep(NA, n_participants),
                   num_any_relapse_PCR = rep(NA, n_participants),
                   num_any_relapse_LM = rep(NA, n_participants),
                   num_any_relapse_D = rep(NA, n_participants),
                   num_relapse_PCR = rep(NA, n_participants),
                   num_relapse_LM = rep(NA, n_participants),
                   num_relapse_D = rep(NA, n_participants),
                   num_trial_PCR = rep(NA, n_participants),
                   num_trial_LM = rep(NA, n_participants),
                   num_trial_D = rep(NA, n_participants))
  
  # add in particpant information 
  df$participant_id <- output_participant$Participant_ID
  df$trial_arm <- output_participant$Trial_Arm
  df$days_enrolled <- output_participant$Dropout_Date - output_participant$Enrollment_Date
  
  # count number of each infection type 
  df$num_any_PCR <- sapply(df$participant_id, function(x){countRecurrentInfs(output_recurrent = output_recurrent,
                                                                             output_participant = output_participant,
                                                                             participant_id = x,
                                                                             inf_types = c('REINFECTION_PCR', 'REINFECTION_LM', 'REINFECTION_D', 'REINFECTION_T',
                                                                                           'RELAPSE_PCR_PRE_ENROLLMENT', 'RELAPSE_LM_PRE_ENROLLMENT', 'RELAPSE_D_PRE_ENROLLMENT', 'RELAPSE_T_PRE_ENROLLMENT',
                                                                                           'RELAPSE_PCR_POST_ENROLLMENT', 'RELAPSE_LM_POST_ENROLLMENT', 'RELAPSE_D_POST_ENROLLMENT', 'RELAPSE_T_POST_ENROLLMENT'))})
  
  df$num_any_LM <- sapply(df$participant_id, function(x){countRecurrentInfs(output_recurrent = output_recurrent,
                                                                            output_participant = output_participant,
                                                                            participant_id = x,
                                                                            inf_types = c('REINFECTION_LM', 'REINFECTION_D', 'REINFECTION_T',
                                                                                          'RELAPSE_LM_PRE_ENROLLMENT', 'RELAPSE_D_PRE_ENROLLMENT', 'RELAPSE_T_PRE_ENROLLMENT',
                                                                                          'RELAPSE_LM_POST_ENROLLMENT', 'RELAPSE_D_POST_ENROLLMENT', 'RELAPSE_T_POST_ENROLLMENT'))})
  
  df$num_any_D <- sapply(df$participant_id, function(x){countRecurrentInfs(output_recurrent = output_recurrent,
                                                                           output_participant = output_participant,
                                                                           participant_id = x,
                                                                           inf_types = c('REINFECTION_D', 'REINFECTION_T',
                                                                                          'RELAPSE_D_PRE_ENROLLMENT', 'RELAPSE_T_PRE_ENROLLMENT',
                                                                                         'RELAPSE_D_POST_ENROLLMENT', 'RELAPSE_T_POST_ENROLLMENT'))})
  
  df$num_any_relapse_PCR <- sapply(df$participant_id, function(x){countRecurrentInfs(output_recurrent = output_recurrent,
                                                                                 output_participant = output_participant,
                                                                                 participant_id = x,
                                                                                 inf_types = c('RELAPSE_PCR_PRE_ENROLLMENT', 'RELAPSE_LM_PRE_ENROLLMENT', 'RELAPSE_D_PRE_ENROLLMENT', 'RELAPSE_T_PRE_ENROLLMENT',
                                                                                               'RELAPSE_PCR_POST_ENROLLMENT', 'RELAPSE_LM_POST_ENROLLMENT', 'RELAPSE_D_POST_ENROLLMENT', 'RELAPSE_T_POST_ENROLLMENT'))})
  
  df$num_any_relapse_LM <- sapply(df$participant_id, function(x){countRecurrentInfs(output_recurrent = output_recurrent,
                                                                                output_participant = output_participant,
                                                                                participant_id = x,
                                                                                inf_types = c('RELAPSE_LM_PRE_ENROLLMENT', 'RELAPSE_D_PRE_ENROLLMENT', 'RELAPSE_T_PRE_ENROLLMENT',
                                                                                              'RELAPSE_LM_POST_ENROLLMENT', 'RELAPSE_D_POST_ENROLLMENT', 'RELAPSE_T_POST_ENROLLMENT'))})
  
  df$num_any_relapse_D <- sapply(df$participant_id, function(x){countRecurrentInfs(output_recurrent = output_recurrent,
                                                                               output_participant = output_participant,
                                                                               participant_id = x,
                                                                               inf_types = c('RELAPSE_D_PRE_ENROLLMENT', 'RELAPSE_T_PRE_ENROLLMENT',
                                                                                             'RELAPSE_D_POST_ENROLLMENT', 'RELAPSE_T_POST_ENROLLMENT'))})
  
  df$num_relapse_PCR <- sapply(df$participant_id, function(x){countRecurrentInfs(output_recurrent = output_recurrent,
                                                                                 output_participant = output_participant,
                                                                                 participant_id = x,
                                                                                 inf_types = c('RELAPSE_PCR_PRE_ENROLLMENT', 'RELAPSE_LM_PRE_ENROLLMENT', 'RELAPSE_D_PRE_ENROLLMENT', 'RELAPSE_T_PRE_ENROLLMENT'))})
  
  df$num_relapse_LM <- sapply(df$participant_id, function(x){countRecurrentInfs(output_recurrent = output_recurrent,
                                                                                output_participant = output_participant,
                                                                                participant_id = x,
                                                                                inf_types = c('RELAPSE_LM_PRE_ENROLLMENT', 'RELAPSE_D_PRE_ENROLLMENT', 'RELAPSE_T_PRE_ENROLLMENT'))})
  
  df$num_relapse_D <- sapply(df$participant_id, function(x){countRecurrentInfs(output_recurrent = output_recurrent,
                                                                               output_participant = output_participant,
                                                                               participant_id = x,
                                                                               inf_types = c('RELAPSE_D_PRE_ENROLLMENT', 'RELAPSE_T_PRE_ENROLLMENT'))})

  # get incidence rates in the context of a trial 
  df$num_trial_PCR <- sapply(df$participant_id, function(x){countRecurrentInfsTrial(output_recurrent = output_recurrent,
                                                                                    output_trial = output_trial,
                                                                                    output_participant = output_participant,
                                                                                    participant_id = x,
                                                                                    detection_type = 'PCR')})

  df$num_trial_LM <- sapply(df$participant_id, function(x){countRecurrentInfsTrial(output_recurrent = output_recurrent,
                                                                                    output_trial = output_trial,
                                                                                    output_participant = output_participant,
                                                                                    participant_id = x,
                                                                                    detection_type = 'LM')})

  df$num_trial_D <- sapply(df$participant_id, function(x){countRecurrentInfsTrial(output_recurrent = output_recurrent,
                                                                                    output_trial = output_trial,
                                                                                    output_participant = output_participant,
                                                                                    participant_id = x,
                                                                                    detection_type = 'D')})

  
  # add in information about when they had their first detected recurrent infection 
  df$time_recurrent_LM <- sapply(df$participant_id, function(x){getTimeFirstRecurrent(output_recurrent = output_recurrent,
                                                                                      output_trial = output_trial,
                                                                                      output_participant = output_participant,
                                                                                      participant_id = x,
                                                                                      detection_type = 'LM')})
  df$time_recurrent_PCR <- sapply(df$participant_id, function(x){getTimeFirstRecurrent(output_recurrent = output_recurrent,
                                                                                       output_trial = output_trial,
                                                                                       output_participant = output_participant,
                                                                                       participant_id = x,
                                                                                       detection_type = 'PCR')})

  df$time_recurrent_D <- sapply(df$participant_id, function(x){getTimeFirstRecurrent(output_recurrent = output_recurrent,
                                                                                     output_trial = output_trial,
                                                                                     output_participant = output_participant,
                                                                                     participant_id = x,
                                                                                     detection_type = 'D')})

  
  df$time_any_relapse_LM <- sapply(df$participant_id, function(x){getTimeFirstRelapseTrial(output_recurrent = output_recurrent,
                                                                                  output_trial = output_trial,
                                                                                  output_participant = output_participant,
                                                                                  participant_id = x,
                                                                                  detection_type = 'LM',
                                                                                  relapse_type = 'ANY')})
  
  df$time_any_relapse_PCR <- sapply(df$participant_id, function(x){getTimeFirstRelapseTrial(output_recurrent = output_recurrent,
                                                                                   output_trial = output_trial,
                                                                                   output_participant = output_participant,
                                                                                   participant_id = x,
                                                                                   detection_type = 'PCR',
                                                                                   relapse_type = 'ANY')})
  
  df$time_relapse_LM_trial <- sapply(df$participant_id, function(x){getTimeFirstRelapseTrial(output_recurrent = output_recurrent,
                                                                                      output_trial = output_trial,
                                                                                      output_participant = output_participant,
                                                                                      participant_id = x,
                                                                                      detection_type = 'LM',
                                                                                      relapse_type = 'PRE')})
  
  df$time_relapse_PCR_trial <- sapply(df$participant_id, function(x){getTimeFirstRelapseTrial(output_recurrent = output_recurrent,
                                                                                       output_trial = output_trial,
                                                                                       output_participant = output_participant,
                                                                                       participant_id = x,
                                                                                       detection_type = 'PCR',
                                                                                       relapse_type = 'PRE')})

  df$time_relapse_D_trial <- sapply(df$participant_id, function(x){getTimeFirstRelapseTrial(output_recurrent = output_recurrent,
                                                                                      output_trial = output_trial,
                                                                                      output_participant = output_participant,
                                                                                      participant_id = x,
                                                                                      detection_type = 'D',
                                                                                      relapse_type = 'PRE')})
  
  df$time_relapse_LM_all <- sapply(df$participant_id, function(x){getTimeFirstRelapseAll(output_recurrent = output_recurrent,
                                                                                         output_trial = output_trial,
                                                                                         output_participant = output_participant,
                                                                                         participant_id = x,
                                                                                         detection_type = 'LM')})
  
  df$time_relapse_PCR_all <- sapply(df$participant_id, function(x){getTimeFirstRelapseAll(output_recurrent = output_recurrent,
                                                                                         output_trial = output_trial,
                                                                                         output_participant = output_participant,
                                                                                         participant_id = x,
                                                                                         detection_type = 'PCR')})
  
  df$time_relapse_D_all <- sapply(df$participant_id, function(x){getTimeFirstRelapseAll(output_recurrent = output_recurrent,
                                                                                         output_trial = output_trial,
                                                                                         output_participant = output_participant,
                                                                                         participant_id = x,
                                                                                         detection_type = 'D')})
  # determine whether the individual was censored 
  df$censor_recurrent_LM <- df$censor_recurrent_PCR <- df$censor_recurrent_D <- df$censor_any_relapse_LM <- df$censor_any_relapse_PCR <- df$censor_relapse_LM_trial <- df$censor_relapse_PCR_trial <- df$censor_relapse_D_trial <- df$censor_relapse_LM_all <- df$censor_relapse_PCR_all <- df$censor_relapse_D_all <- 0
  

  df$censor_recurrent_LM[df$time_recurrent_LM == (df$days_enrolled + 1)] <- 1
  df$censor_recurrent_PCR[df$time_recurrent_PCR == (df$days_enrolled + 1)] <- 1
  df$censor_recurrent_D[df$time_recurrent_D == (df$days_enrolled + 1)] <- 1
  df$censor_any_relapse_LM[df$time_any_relapse_LM == (df$days_enrolled + 1)] <- 1
  df$censor_any_relapse_PCR[df$time_any_relapse_PCR == (df$days_enrolled + 1)] <- 1
  df$censor_relapse_LM_trial[df$time_relapse_LM_trial == (df$days_enrolled + 1)] <- 1
  df$censor_relapse_PCR_trial[df$time_relapse_PCR_trial == (df$days_enrolled + 1)] <- 1
  df$censor_relapse_D_trial[df$time_relapse_D_trial == (df$days_enrolled + 1)] <- 1
  df$censor_relapse_LM_all[df$time_relapse_LM_all == (df$days_enrolled + 1)] <- 1
  df$censor_relapse_PCR_all[df$time_relapse_PCR_all == (df$days_enrolled + 1)] <- 1
  df$censor_relapse_D_all[df$time_relapse_D_all == (df$days_enrolled + 1)] <- 1
  
  df$time_recurrent_LM[df$censor_recurrent_LM == 1] = df$time_recurrent_LM[df$censor_recurrent_LM == 1] - 1
  df$time_recurrent_PCR[df$censor_recurrent_PCR == 1] = df$time_recurrent_PCR[df$censor_recurrent_PCR == 1] - 1
  df$time_recurrent_D[df$censor_recurrent_D == 1] = df$time_recurrent_D[df$censor_recurrent_D == 1] - 1
  df$time_any_relapse_LM[df$censor_any_relapse_LM == 1] = df$time_any_relapse_LM[df$censor_any_relapse_LM == 1] - 1
  df$time_any_relapse_PCR[df$censor_any_relapse_PCR == 1] = df$time_any_relapse_PCR[df$censor_any_relapse_PCR == 1] - 1
  df$time_relapse_LM_trial[df$censor_relapse_LM_trial == 1] = df$time_relapse_LM_trial[df$censor_relapse_LM_trial == 1] - 1
  df$time_relapse_PCR_trial[df$censor_relapse_PCR_trial == 1] = df$time_relapse_PCR_trial[df$censor_relapse_PCR_trial == 1] - 1
  df$time_relapse_D_trial[df$censor_relapse_D_trial == 1] = df$time_relapse_D_trial[df$censor_relapse_D_trial == 1] - 1
  df$time_relapse_LM_all[df$censor_relapse_LM_all == 1] = df$time_relapse_LM_all[df$censor_relapse_LM_all == 1] - 1
  df$time_relapse_PCR_all[df$censor_relapse_PCR_all == 1] = df$time_relapse_PCR_all[df$censor_relapse_PCR_all == 1] - 1
  df$time_relapse_D_all[df$censor_relapse_D_all == 1] = df$time_relapse_D_all[df$censor_relapse_D_all == 1] - 1
  
  # return data frame 
  return(df)
}

# function to compute the efficacy measures 
calcCoxEff <- function(df,
                       days_followup,
                       days_left_censored = 32)
{
  ##### recurrent LM #####
  # data organizaton 
  df_subset <- subset(df, time_recurrent_LM >= 0)
  df_subset$status <- (!df_subset$censor_recurrent_LM) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # get the updated number of days of follow-up in light of left censoring 
  days_followup <- c(0, days_followup[days_followup >= days_left_censored])
  
  # specify the lower and upper bounds for each interval 
  df_subset$time_upper_LM <- ifelse(df_subset$censor_recurrent_LM, Inf, df_subset$time_recurrent_LM)
  df_subset$time_lower_LM <- sapply(1:nrow(df_subset), function(i){ifelse(df_subset$censor_recurrent_LM[i], df_subset$time_recurrent_LM[i], tail(days_followup[days_followup < df_subset$time_recurrent_LM[i]], n = 1))})
  
  # fit cox proportional hazard
  #cph <- coxph(formula = Surv(time_recurrent_LM, status) ~ trial_arm_num, data = df_subset)
  #cph <- survreg(formula = Surv(time = df_subset$time_lower_LM, time2 = df_subset$time_upper_LM, type = 'interval2') ~ trial_arm_num, data = df_subset, dist = 'exponential')
  cph <- ic_sp(Surv(time = df_subset$time_lower_LM, time2 = df_subset$time_upper_LM, type = 'interval2') ~ trial_arm_num, data = df_subset, model = 'ph')
  
  # calcuate efficacy and return
  cph_recurrent_LM <- 1 - exp(as.numeric(cph$coefficients))
  
  ##### recurrent PCR #####
  # data organizaton 
  df_subset <- subset(df, time_recurrent_PCR >= 0)
  df_subset$status <- (!df_subset$censor_recurrent_PCR) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # specify the lower and upper bounds for each interval 
  df_subset$time_upper_PCR <- ifelse(df_subset$censor_recurrent_PCR, Inf, df_subset$time_recurrent_PCR)
  df_subset$time_lower_PCR <- sapply(1:nrow(df_subset), function(i){ifelse(df_subset$censor_recurrent_PCR[i], df_subset$time_recurrent_PCR[i], tail(days_followup[days_followup < df_subset$time_recurrent_PCR[i]], n = 1))})
  
  # fit cox proportional hazard
  #cph <- coxph(formula = Surv(time_recurrent_PCR, status) ~ trial_arm_num, data = df_subset)
  #cph <- coxph(formula = Surv(time = time_lower_PCR, time2 = time_upper_PCR) ~ trial_arm_num, data = df_subset)
  cph <- ic_sp(Surv(time = df_subset$time_lower_PCR, time2 = df_subset$time_upper_PCR, type = 'interval2') ~ trial_arm_num, data = df_subset, model = 'ph')
  
  # calcuate efficacy and return
  cph_recurrent_PCR <- 1 - exp(as.numeric(cph$coefficients))

   ##### recurrent D #####
  # data organizaton 
  df_subset <- subset(df, time_recurrent_D >= 0)
  df_subset$status <- (!df_subset$censor_recurrent_D) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # specify the lower and upper bounds for each interval 
  df_subset$time_upper_D <- ifelse(df_subset$censor_recurrent_D, Inf, df_subset$time_recurrent_D)
  df_subset$time_lower_D <- sapply(1:nrow(df_subset), function(i){ifelse(df_subset$censor_recurrent_D[i], df_subset$time_recurrent_D[i], tail(days_followup[days_followup < df_subset$time_recurrent_D[i]], n = 1))})
  
  # fit cox proportional hazard
  cph <- ic_sp(Surv(time = df_subset$time_lower_D, time2 = df_subset$time_upper_D, type = 'interval2') ~ trial_arm_num, data = df_subset, model = 'ph')
  
  # calcuate efficacy and return
  cph_recurrent_D <- 1 - exp(as.numeric(cph$coefficients))
  
  ##### any relapse LM #####
  # data organizaton 
  df_subset <- subset(df, time_any_relapse_LM >= 0)
  df_subset$status <- (!df_subset$censor_any_relapse_LM) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # fit cox proportional hazard
  cph <- coxph(formula = Surv(time_any_relapse_LM, status) ~ trial_arm_num, data = df_subset)
  
  # calcuate efficacy and return
  cph_any_relapse_LM <- 1 - exp(as.numeric(cph$coefficients))
  
  ##### any relapse PCR #####
  # data organizaton 
  df_subset <- subset(df, time_any_relapse_PCR >= 0)
  df_subset$status <- (!df_subset$censor_any_relapse_PCR) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # fit cox proportional hazard
  cph <- coxph(formula = Surv(time_any_relapse_PCR, status) ~ trial_arm_num, data = df_subset)
  
  # calcuate efficacy and return
  cph_any_relapse_PCR <- 1 - exp(as.numeric(cph$coefficients))
  
  ##### relapse LM trial #####
  # data organizaton 
  df_subset <- subset(df, time_relapse_LM_trial >= 0)
  df_subset$status <- (!df_subset$censor_relapse_LM_trial) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # fit cox proportional hazard
  cph <- coxph(formula = Surv(time_relapse_LM_trial, status) ~ trial_arm_num, data = df_subset)
  
  # calcuate efficacy and return
  cph_relapse_LM_trial <- 1 - exp(as.numeric(cph$coefficients))
  
  ##### relapse PCR trial #####
  # data organizaton 
  df_subset <- subset(df, time_relapse_PCR_trial >= 0)
  df_subset$status <- (!df_subset$censor_relapse_PCR_trial) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # fit cox proportional hazard
  cph <- coxph(formula = Surv(time_relapse_PCR_trial, status) ~ trial_arm_num, data = df_subset)
  
  # calcuate efficacy and return
  cph_relapse_PCR_trial <- 1 - exp(as.numeric(cph$coefficients))

  ##### relapse D trial #####
  # data organizaton 
  df_subset <- subset(df, time_relapse_D_trial >= 0)
  df_subset$status <- (!df_subset$censor_relapse_D_trial) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # fit cox proportional hazard
  cph <- coxph(formula = Surv(time_relapse_D_trial, status) ~ trial_arm_num, data = df_subset)
  
  # calcuate efficacy and return
  cph_relapse_D_trial <- 1 - exp(as.numeric(cph$coefficients))
  
  ##### relapse LM all #####
  # data organizaton 
  df_subset <- subset(df, time_relapse_LM_all >= 0)
  df_subset$status <- (!df_subset$censor_relapse_LM_all) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # fit cox proportional hazard
  cph <- coxph(formula = Surv(time_relapse_LM_all, status) ~ trial_arm_num, data = df_subset)
  
  # calcuate efficacy and return
  cph_relapse_LM_all <- 1 - exp(as.numeric(cph$coefficients))
  
  ##### relapse PCR all #####
  # data organizaton 
  df_subset <- subset(df, time_relapse_PCR_all >= 0)
  df_subset$status <- (!df_subset$censor_relapse_PCR_all) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # fit cox proportional hazard
  cph <- coxph(formula = Surv(time_relapse_PCR_all, status) ~ trial_arm_num, data = df_subset)
  
  # calcuate efficacy and return
  cph_relapse_PCR_all <- 1 - exp(as.numeric(cph$coefficients))

    ##### relapse D all #####
  # data organizaton 
  df_subset <- subset(df, time_relapse_D_all >= 0)
  df_subset$status <- (!df_subset$censor_relapse_D_all) * 1 + 1
  df_subset$trial_arm_num <- ifelse(df_subset$trial_arm == 'TREATMENT', 2, 1)
  
  # fit cox proportional hazard
  cph <- coxph(formula = Surv(time_relapse_D_all, status) ~ trial_arm_num, data = df_subset)
  
  # calcuate efficacy and return
  cph_relapse_D_all <- 1 - exp(as.numeric(cph$coefficients))
  
  ### organize into efficacy vector 
  eff <- c(cph_recurrent_LM = cph_recurrent_LM,
           cph_recurrent_PCR = cph_recurrent_PCR,
           cph_recurrent_D = cph_recurrent_D,
           cph_any_relapse_LM = cph_any_relapse_LM,
           cph_any_relapse_PCR = cph_any_relapse_PCR,
           cph_relapse_LM_trial = cph_relapse_LM_trial,
           cph_relapse_PCR_trial = cph_relapse_PCR_trial,
           cph_relapse_D_trial = cph_relapse_D_trial,
           cph_relapse_LM_all = cph_relapse_LM_all,
           cph_relapse_PCR_all = cph_relapse_PCR_all,
           cph_relapse_D_all = cph_relapse_D_all)

  return(eff)
}


# function to compute the incidence efficacy measures 
calcIncidEff <- function(df)
{
  # data organization
  df <- subset(df, time_recurrent_LM >= 0 | time_recurrent_PCR >= 0 | time_any_relapse_LM >= 0 | time_any_relapse_PCR >= 0 | time_relapse_LM_trial >= 0 | time_relapse_PCR_trial >= 0 | time_relapse_LM_all >= 0 | time_relapse_PCR_all >= 0)
  
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

# function to compute the risk-based efficacy measure 
calcRiskEff <- function(df)
{
  # data organization
  df <- subset(df, time_recurrent_LM >= 0 | time_recurrent_PCR >= 0 | time_any_relapse_LM >= 0 | time_any_relapse_PCR >= 0 | time_relapse_LM_trial >= 0 | time_relapse_PCR_trial >= 0 | time_relapse_LM_all >= 0 | time_relapse_PCR_all >= 0)
  
  # get the risk efficacy for PCR-detectable recurrent infections 
  any_event_treatment <- sum(df$num_trial_PCR[df$trial_arm == "TREATMENT"] > 0)
  any_event_placebo <- sum(df$num_trial_PCR[df$trial_arm == "PLACEBO"] > 0)
  
  risk_trial_PCR = 1 - (any_event_treatment / any_event_placebo)

# get the risk efficacy for LM-detectable recurrent infections 
  any_event_treatment <- sum(df$num_trial_LM[df$trial_arm == "TREATMENT"] > 0)
  any_event_placebo <- sum(df$num_trial_LM[df$trial_arm == "PLACEBO"] > 0)
  
  risk_trial_LM = 1 - (any_event_treatment / any_event_placebo)

  # get the risk efficacy for PCR-detectable recurrent infections 
  any_event_treatment <- sum(df$num_trial_D[df$trial_arm == "TREATMENT"] > 0)
  any_event_placebo <- sum(df$num_trial_D[df$trial_arm == "PLACEBO"] > 0)
  
  risk_trial_D = 1 - (any_event_treatment / any_event_placebo)

  # get the risk efficacy for true PCR-detectabe treatment failures
  any_event_treatment <- sum(df$num_relapse_PCR[df$trial_arm == "TREATMENT"] > 0)
  any_event_placebo <- sum(df$num_relapse_PCR[df$trial_arm == "PLACEBO"] > 0)
  
  risk_relapse_PCR = 1 - (any_event_treatment / any_event_placebo)
  
  # get the risk efficacy for true LM-detectable treatment failures
  any_event_treatment <- sum(df$num_relapse_LM[df$trial_arm == "TREATMENT"] > 0)
  any_event_placebo <- sum(df$num_relapse_LM[df$trial_arm == "PLACEBO"] > 0)
  
  risk_relapse_LM = 1 - (any_event_treatment / any_event_placebo)
  
  # get the risk efficacy for true symptomatic treatment failures
  any_event_treatment <- sum(df$num_relapse_D[df$trial_arm == "TREATMENT"] > 0)
  any_event_placebo <- sum(df$num_relapse_D[df$trial_arm == "PLACEBO"] > 0)
  
  risk_relapse_D = 1 - (any_event_treatment / any_event_placebo)
  
  # organize into an efficacy vector and return 
  eff <- c(risk_relapse_PCR = risk_relapse_PCR,
           risk_relapse_LM = risk_relapse_LM,
           risk_relapse_D = risk_relapse_D,
           risk_trial_PCR = risk_trial_PCR,
           risk_trial_LM = risk_trial_LM,
           risk_trial_D = risk_trial_D)
  
  return(eff)
}
