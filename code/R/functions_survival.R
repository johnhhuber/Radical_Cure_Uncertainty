getTimeFirstEvent <- function(recurr, indiv, id, detection_type)
{
  # subset to individuals with that id 
  recurr <- recurr[-match(id, recurr$Participant_ID),]
  recurr <- subset(recurr, Participant_ID == id)
  
  # get timing of infection 
  if(nrow(recurr) > 0)
  {
    if(detection_type == 'PCR')
    {
      recurr <- subset(recurr, grepl('_PCR', Infection_Type) || grepl('_LM', Infection_Type) || grepl('_D', Infection_Type) || grepl('_T', Infection_Type))
    }
    if(detection_type == 'LM')
    {
      recurr <- subset(recurr, grepl('_LM', Infection_Type) || grepl('_D', Infection_Type) || grepl('_T', Infection_Type))
    }
    if(detection_type == 'D')
    {
      recurr <- subset(recurr, grepl('_D', Infection_Type) || grepl('_T', Infection_Type))
    }
    
    if(nrow(recurr) > 0)
    {
      time <- min(recurr$Infection_Date) - indiv$Enrollment_Date[indiv$Participant_ID == id]
    }else{
      time <- indiv$Dropout_Date[indiv$Participant_ID == id] - indiv$Enrollment_Date[indiv$Participant_ID == id] + 1 
    }
  }else{
    time <- indiv$Dropout_Date[indiv$Participant_ID == id] - indiv$Enrollment_Date[indiv$Participant_ID == id] + 1
  }
  return(time)
}

# generate survival curve 
genSurvCurv <- function(recurr, indiv, detection_type)
{
  # get timing of infections for placebo and treatment arms
  id_placebo <- indiv$Participant_ID[indiv$Trial_Arm == 'PLACEBO']
  id_treatment <- indiv$Participant_ID[indiv$Trial_Arm == 'TREATMENT']
  
  time_placebo <- sapply(id_placebo, getTimeFirstEvent, recurr = recurr, indiv = indiv, detection_type = detection_type)
  time_treatment <- sapply(id_treatment, getTimeFirstEvent, recurr = recurr, indiv = indiv, detection_type = detection_type)
  
  # remove -inf observations if any 
  time_placebo <- time_placebo[time_placebo > 0]
  time_treatment <- time_treatment[time_treatment > 0]
  
  # get the maximum time of followup 
  max_followup <- max(indiv$Dropout_Date - indiv$Enrollment_Date)
  times_followup <- 0:max_followup
  
  # construct survival curves 
  surv_placebo <- sapply(times_followup, function(t){sum(time_placebo > t) / length(time_placebo)})
  surv_treatment <- sapply(times_followup, function(t){sum(time_treatment > t) / length(time_treatment)})
  
  # store in dataframe and return 
  df <- data.frame(t = times_followup, surv_placebo = surv_placebo, surv_treatment = surv_treatment)
  return(df)
}

# function to write survival curve 
genSurvDf <- function(recurr, indiv)
{
  
  # PCR-detectable infections 
  df_PCR <- genSurvCurv(recurr = recurr, indiv = indiv, detection_type = 'PCR')
  df <- data.frame(t = df_PCR$t)
  df$surv_placebo_PCR <- df_PCR$surv_placebo
  df$surv_treatment_PCR <- df_PCR$surv_treatment
  
  # LM-detectable infections 
  df_LM <- genSurvCurv(recurr = recurr, indiv = indiv, detection_type = 'LM')
  df$surv_placebo_LM <- df_LM$surv_placebo
  df$surv_treatment_LM <- df_LM$surv_treatment
  
  
  # clinical infections 
  df_D <- genSurvCurv(recurr = recurr, indiv = indiv, detection_type = 'D')
  df$surv_placebo_D <- df_D$surv_placebo
  df$surv_treatment_D <- df_D$surv_treatment

  # return output
  return(df)
}
