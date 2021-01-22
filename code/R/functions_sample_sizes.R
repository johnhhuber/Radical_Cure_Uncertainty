# function to calculate sample size using the risk-based measure of efficacy
calcSampleSizeRisk <- function(alpha_level, # the probability of a type I error
                               power, # desired power of trial
                               prop_placebo, # proportion of individuals that experience event in placebo arm
                               efficacy)
{
  # get the z-scores associated with the alpha-level and the power 
  z1 = abs(qnorm(1 - (alpha_level / 2)))
  z2 = abs(qnorm(1 - power))
  
  # calculate the proportion in treatment arm
  prop_treatment <- efficacy * prop_placebo
  
  # get the average proportion between two arms
  prop_mean <- 0.5 * (prop_placebo + prop_treatment)
  
  # calculate the sample size in each trial arm
  n_treatment <- ((z1 + z2)^2) * 2 * prop_mean * (1 - prop_mean) / ((prop_placebo - prop_treatment) ^ 2)
  n_placebo <- n_treatment
  
  # return output 
  return(c(treatment = n_treatment,
           placebo = n_placebo))
}

# function to calculate the sample size using incidence-based measure of efficacy 
calcSampleSizeIncidence <- function(alpha_level, # probability of a type I error 
                                    power, # desired power of trial 
                                    incidence_rate_placebo, # incidence rate of event in placebo arm 
                                    efficacy) # efficacy of intervention 
{
  # get the z-scores associated with the alpha-level and the power 
  z1 = abs(qnorm(1 - (alpha_level / 2)))
  z2 = abs(qnorm(1 - power))
  
  # get the incidence rate of the treatment arm 
  incidence_rate_treatment <- incidence_rate_placebo * (1 - efficacy)
  
  # calculate the sample size in each trial arm 
  n_treatment <- (((z1 + z2)^2) * (incidence_rate_treatment + incidence_rate_placebo)) / ((incidence_rate_placebo - incidence_rate_treatment)^2)
  n_placebo <- n_treatment
  
  # return output 
  return(c(treatment = n_treatment,
           placebo = n_placebo))
}

# function to calculate the power of a trial for a given sample size calculated using incidence-based measures of efficacy
calcPowerIncidence <- function(alpha_level, # probability of a type I error
                               sample_size, # sample size of each trial arm
                               incidence_rate_placebo, # incidence rate of event in the placebo arm
                               efficacy) # efficacy of the intervention 
{
  # get the z-score associated with the alpha-level 
  z1 = abs(qnorm(1 - (alpha_level / 2)))
  
  # get the incidence rate of the treatment arm 
  incidence_rate_treatment <- incidence_rate_placebo * (1 - efficacy)
  
  # get the z-score for power 
  z2 <- sqrt(sample_size / (incidence_rate_placebo + incidence_rate_treatment)) * abs(incidence_rate_placebo - incidence_rate_treatment) - z1
  
  # get the power associated with the z-score 
  power <- pnorm(abs(z2))
  
  # return output
  return(power)
}

# function to calculate the sample size using incidence-based measure of efficacy and lower limit of desired efficacy
calcSampleSizeIncidenceLowerLim = function(alpha_level, # probability of a type I error
                                           power, # desired power of trial
                                           incidence_rate_placebo, # incidence rate of event in placebo arm
                                           efficacy, # efficacy of intervention 
                                           efficacy_lower) # lower limit of the desired efficacy
{
  # get the z-scores associated with the alpha-level and the power 
  z1 = abs(qnorm(1 - (alpha_level / 2)))
  z2 = abs(qnorm(1 - power))
  
  # get the relative risks 
  RR <- 1 - efficacy
  RR_lower <- 1 - efficacy_lower
  
  # get the incidence rate of the treatment arm 
  incidence_rate_treatment <- incidence_rate_placebo * (1 - efficacy)
  
  # calculate the sample size in each trial arm 
  n_treatment <- (((z1 + z2)^2) * ((1/incidence_rate_treatment) + (1/incidence_rate_placebo))) / ((log(RR / RR_lower))^2)
  n_placebo <- n_treatment
  
  # return output 
  return(c(treatment = n_treatment,
           placebo = n_placebo))
}

# function to calculate the sample size necessary for Cox calculation
calcSampleSizeCox = function(alpha_level, # probability of a type I error
                             power, # desired power of trial
                             allocation_ratio, # allocation ratio between treatment and control arms
                             incidence_rate_placebo, # incidence rate of event in placebo arm
                             efficacy) # efficacy of the intervention 
{
  # get the z-scores associated with the alpha-level and the power 
  z1 = abs(qnorm(1 - (alpha_level / 2)))
  z2 = abs(qnorm(1 - power))
  
  # get the incidence rate of the treatment arm 
  incidence_rate_treatment <- incidence_rate_placebo * (1 - efficacy)
  
  # calculate the number of events necessary for the trial 
  events = ((((allocation_ratio + 1)^2) / allocation_ratio) * ((z1 + z2)^2) / (log(1 - efficacy)^2))
  
  # calculate the sample size in each trial arm 
  n_placebo <- events / (allocation_ratio * incidence_rate_treatment + incidence_rate_placebo)
  n_treatment <- allocation_ratio * n_placebo
  
  # return output 
  return(c(treatment = n_treatment,
           placebo = n_placebo))
}

# function to calculate the power of 
calcPowerCox <- function(alpha_level, # probability of a type I error 
                         sample_size_placebo, # sample size in the placebo arm
                         allocation_ratio, # allocation ratio between treatment and control arms
                         incidence_rate_placebo, # incidence rate of events in the control arm 
                         efficacy) # efficacy of the intervention
{
  # get the z-score associated with the alpha-level 
  z1 = abs(qnorm(1 - (alpha_level / 2)))
  
  # get the incidence rate of the treatment arm 
  incidence_rate_treatment <- incidence_rate_placebo * (1 - efficacy)
  
  # calculate the total number of events needed to be observed in the trial 
  events <- sample_size_placebo * (allocation_ratio * incidence_rate_treatment + incidence_rate_placebo)
  
  # get the z-score for power 
  z2 = sqrt((events * allocation_ratio * (log(1 - efficacy)^2))/((allocation_ratio + 1)^2)) - z1
  
  # get the power associated with the z-score 
  power <- pnorm(abs(z2))
  
  # return output
  return(power)
}


calcPowerIncidence(alpha_level = 0.05,
                   sample_size = 1000,
                   incidence_rate_placebo = 5,
                   efficacy = 0.75)
