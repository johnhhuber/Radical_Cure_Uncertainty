#include "Trial.h"
#include "randlib.h"

using namespace std;

void Trial::update(Params& theta, Population& POP, int curr_day)
{
  // only update if we are at the point where the trial is ongoing
  if(curr_day >= recruitment_start_date)
  {
    // if the first date of the trial, open up to enrollment
    if(curr_day == recruitment_start_date)
    {
      is_enrollment_open = true;
    }

    // if enrollment is open, enroll participants
    if(is_enrollment_open)
    {
      enrollParticipants(POP, curr_day);
    }

    // update the status of the participants if necessary and collect the data
    updateParticipants(POP, theta, curr_day);
  }
}



void Trial::enrollParticipants(Population& POP, int curr_day)
{
  // std::vector<Individual> eligible_participants = POP->people;
  std::vector<Individual>::iterator itr_eligible = POP.people.begin();

  // loop through population and look to enroll participants
  while((treatment_arm.size() < treatment_arm_sample_size || placebo_arm.size() < placebo_arm_sample_size) && itr_eligible != POP.people.end())
  {
    // check if we can enroll the person in the trial
    // to be in the trial, need to have light-microscopy detectable infection and show up in the health system
    if(!(itr_eligible->enrolled_in_trial) && (itr_eligible->participant_ID < 0) && (itr_eligible->T_last_Symp_BS == 0))
    {
      std::string arm_allocation;
      // randomly allocate to treatment or placebo arms
      if(genunf(0.0, 1.0) < prob_treatment_arm_allocation)
      {
        if(treatment_arm.size() < treatment_arm_sample_size)
        {
          arm_allocation = "TREATMENT";
        }else{
          arm_allocation = "PLACEBO";
        }
      }
      else
      {
        if(placebo_arm.size() < placebo_arm_sample_size)
        {
          arm_allocation = "PLACEBO";
        }
        else
        {
          arm_allocation = "TREATMENT";
        }
      }

      num_enrolled++;

      // update individual information
      itr_eligible->enrolled_in_trial = true;
      itr_eligible->trial_arm = arm_allocation;
      itr_eligible->participant_ID = num_enrolled;
      itr_eligible->enrollment_date = curr_day;

      // add to appropriate trial arm
      if(arm_allocation == "TREATMENT")
      {
        // treatment_arm.push_back(&(*itr_eligible));
        treatment_arm.push_back(itr_eligible->participant_ID);
      }else{
        // placebo_arm.push_back(&(*itr_eligible));
        placebo_arm.push_back(itr_eligible->participant_ID);
      }

      // add participant to participant data structure
      participant_data.insert(std::pair<int, std::tuple<string, int, double>>(itr_eligible->participant_ID, make_tuple(itr_eligible->trial_arm, itr_eligible->enrollment_date, itr_eligible->zeta_het)));

      // create vector in map for recording trial data
      std::vector<std::tuple<int, bool>> individual_LM_recurrent_data;
      std::vector<std::tuple<int, string>> individual_all_recurrent_data;

      record_LM_recurrent_infections.insert(std::pair<int, std::vector<std::tuple<int, bool>>>(itr_eligible->participant_ID, individual_LM_recurrent_data));
      record_all_recurrent_infections.insert(std::pair<int, std::vector<std::tuple<int, string>>>(itr_eligible->participant_ID, individual_all_recurrent_data));
    }

    // increment iterator
    itr_eligible++;
  }

  // update the status of enrollment
  is_enrollment_open = (treatment_arm.size() < treatment_arm_sample_size || placebo_arm.size() < placebo_arm_sample_size);
}



void Trial::updateParticipants(Population& POP, Params& theta, int curr_day)
{
  std::vector<int>::iterator itr_treatment_arm;
  std::vector<int>::iterator itr_placebo_arm;

  int enrollment_date;
  int days_since_enrollment;

  // loop through the treatment arm and perform actions for each participant
  // this includes administering radical cure on the first day or recording
  // the status of a LM-detectable recurrent infection
  for(itr_treatment_arm = treatment_arm.begin(); itr_treatment_arm != treatment_arm.end(); itr_treatment_arm++)
  {
    std::vector<Individual>::iterator participant = std::find_if(POP.people.begin(), POP.people.end(), [&](const Individual& indiv){return indiv.participant_ID == *itr_treatment_arm;});
    // check that the participant is still enrolled in the trial
    if(participant->enrolled_in_trial)
    {
      enrollment_date = participant->enrollment_date;
      days_since_enrollment = curr_day - enrollment_date;

      // check whether date of followup
      bool is_followup_date = (std::find(trial_followup_dates.begin(), trial_followup_dates.end(), days_since_enrollment) != trial_followup_dates.end());
      if(is_followup_date)
      {
        // record whether or not the individual has LM-detectable infection
        bool has_LM_infection = (participant->I_LM || participant->I_D || participant->T);

        // record status of recurrent infection
        record_LM_recurrent_infections.find(participant->participant_ID)->second.push_back(make_tuple(days_since_enrollment, has_LM_infection));

        // adminster treatment if necessary
        administerTreatment(theta, &(*participant), days_since_enrollment);
      }

      // check whether the individual has a new recurrent infection and record if true
      bool has_new_infection = (participant->I_PCR_new || participant->I_LM_new || participant->I_D_new);
      if(has_new_infection)
      {
        string infection_type;
        // check if it is a symptomatic episode
        if(participant->I_D_new)
        {
          if(participant->CQ_treat)
          {
            infection_type = "T";
          }else{
            infection_type = "D";
          }
        }else{
          if(participant->I_LM_new)
          {
            infection_type = "LM";
          }else{
            infection_type = "PCR";
          }
        }

        // add to record of all recurrent infections
        record_all_recurrent_infections.find(participant->participant_ID)->second.push_back(make_tuple(curr_day, infection_type));
      }

      // if we have reached end of trial for this individual, they are no longer enrolled
      if(days_since_enrollment == trial_duration)
      {
        participant->enrolled_in_trial = false;
      }
    }
  }

  // loop through the placebo arm and perform actions for each participants
  // this involves recording the status of a LM-detectable recurrent infection
  for(itr_placebo_arm = placebo_arm.begin(); itr_placebo_arm != placebo_arm.end(); itr_placebo_arm++)
  {
    std::vector<Individual>::iterator participant = std::find_if(POP.people.begin(), POP.people.end(), [&](const Individual& indiv){return indiv.participant_ID == *itr_placebo_arm;});

    // check that the participant is still enrolled in the trial
    if(participant->enrolled_in_trial)
    {
      enrollment_date = participant->enrollment_date;
      days_since_enrollment = curr_day - enrollment_date;

      // check whether date of followup
      bool is_followup_date = (std::find(trial_followup_dates.begin(), trial_followup_dates.end(), days_since_enrollment) != trial_followup_dates.end());
      if(is_followup_date)
      {
        // record whether or not the individual has LM-detectable infection
        bool has_LM_infection = (participant->I_LM || participant->I_D || participant->T);

        // record status of recurrent infection
        record_LM_recurrent_infections.find(participant->participant_ID)->second.push_back(make_tuple(days_since_enrollment, has_LM_infection));

        // adminster treatment if necessary
        administerTreatment(theta, &(*participant), days_since_enrollment);
      }

      // check whether the individual has a new recurrent infection and record if true
      bool has_new_infection = (participant->I_PCR_new || participant->I_LM_new || participant->I_D_new);
      if(has_new_infection)
      {
        string infection_type;
        // check if it is a symptomatic episode
        if(participant->I_D_new)
        {
          if(participant->CQ_treat)
          {
            infection_type = "T";
          }else{
            infection_type = "D";
          }
        }else{
          if(participant->I_LM_new)
          {
            infection_type = "LM";
          }else{
            infection_type = "PCR";
          }
        }

        // add to record of all recurrent infections
        record_all_recurrent_infections.find(participant->participant_ID)->second.push_back(make_tuple(curr_day, infection_type));
      }
      // if we have reached end of trial for this individual, they are no longer enrolled
      if(days_since_enrollment == trial_duration)
      {
        participant->enrolled_in_trial = false;
      }
    }
  }
}



void Trial::administerTreatment(Params &theta, Individual *trial_participant, int days_since_enrollment)
{
  if(days_since_enrollment == 0)
  {
    if(trial_participant->trial_arm == "TREATMENT" && trial_participant->enrolled_in_trial)
    {
      // randomly draw whether individual receives effective PQ treatment
      if(genunf(0.0, 1.0) < trial_PQ_eff)
      {
        // account for effects of PQ on hypnozoites and prophylaxis
        trial_participant->Hyp = 0;
        trial_participant->AQ8_proph = 1;
        trial_participant->AQ8_proph_timer = theta.CM_PQ_proph;

        // clear developing hypnozoites
        for(int z = 0; z < trial_participant->lam_bite_track.size(); z++)
        {
          trial_participant->lam_bite_track[z] = 0.0;
        }
        for(int z = 0; z < trial_participant->lam_rel_track.size(); z++)
        {
          trial_participant->lam_rel_track[z] = 0.0;
        }

        // PQ can improve efficacy of CQ, so account for this effect in clearance
        // of BS infection
        if(trial_participant->I_LM)
        {
          if(genunf(0.0, 1.0) < theta.CM_CQ_eff_wPQ)
          {
            trial_participant->S     = 0;
            trial_participant->I_PCR = 0;
            trial_participant->I_LM  = 0;
            trial_participant->I_D   = 0;
            trial_participant->T     = 1;
            trial_participant->P     = 0;
          }
        }
      }
    }
  }
}


void Trial::writeParticipantData(Population& POP)
{
  // open output file stream
  ofstream file_out;
  file_out.open(output_file_participants);

  // write the header
  file_out << "Participant_ID," << "Trial_Arm," << "Enrollment_Date," << "Zeta_Het" << std::endl;

  // write out the data
  std::map<int, std::tuple<string, int, double>>::iterator itr_participants = participant_data.begin();
  for(; itr_participants != participant_data.end(); itr_participants++)
  {
    file_out << itr_participants->first << "," << std::get<0>(itr_participants->second) << "," << std::get<1>(itr_participants->second) << "," << std::get<2>(itr_participants->second) << std::endl;
  }

  // std::vector<int>::iterator itr_treatment = treatment_arm.begin();
  // for(; itr_treatment != treatment_arm.end(); itr_treatment++)
  // {
  //   std::vector<Individual>::iterator participant = std::find_if(POP.people.begin(), POP.people.end(), [&](const Individual& indiv){return indiv.participant_ID == *itr_treatment;});
  //   file_out << participant->participant_ID << "," << participant->trial_arm << "," << participant->enrollment_date << ',' << participant->zeta_het << std::endl;
  // }
  //
  // std::vector<int>::iterator itr_placebo = placebo_arm.begin();
  // for(; itr_placebo != placebo_arm.end(); itr_placebo++)
  // {
  //   std::vector<Individual>::iterator participant = std::find_if(POP.people.begin(), POP.people.end(), [&](const Individual& indiv){return indiv.participant_ID == *itr_placebo;});
  //   file_out << participant->participant_ID << "," << participant->trial_arm << "," << participant->enrollment_date << ',' << participant->zeta_het << std::endl;
  // }
  file_out.close();
}



void Trial::writeTrialOutcomes()
{
  // open output file stream
  ofstream file_out;
  file_out.open(output_file_trial);

  // write the header
  file_out << "Participant_ID," << "Days_Since_Enrollment," << "LM_Infection" << endl;

  // write out the data
  std::map<int, std::vector<std::tuple<int, bool>>>::iterator itr_trial = record_LM_recurrent_infections.begin();
  for(; itr_trial != record_LM_recurrent_infections.end(); itr_trial++)
  {
    std::vector<std::tuple<int, bool>>::iterator itr_followup_stat = itr_trial->second.begin();
    for(; itr_followup_stat != itr_trial->second.end(); itr_followup_stat++)
    {
      file_out << itr_trial->first << "," << std::get<0>(*itr_followup_stat) << ',' << std::get<1>(*itr_followup_stat) << std::endl;
    }
  }
  file_out.close();
}



void Trial::writeAllRecurrentInfs()
{
  // open output file stream
  ofstream file_out;
  file_out.open(output_file_recurrent_infs);

  // write the header
  file_out << "Participant_ID," << "Infection_Date," << "Infection_Type" << std::endl;

  // write out the data
  std::map<int, std::vector<std::tuple<int, string>>>::iterator itr_infs = record_all_recurrent_infections.begin();
  for(; itr_infs != record_all_recurrent_infections.end(); itr_infs++)
  {
    std::vector<std::tuple<int, string>>::iterator itr_individual = itr_infs->second.begin();
    for(; itr_individual != itr_infs->second.end(); itr_individual++)
    {
      file_out << itr_infs->first << "," << std::get<0>(*itr_individual) << "," << std::get<1>(*itr_individual) << std::endl;
    }
  }
  file_out.close();
}



void Trial::Initialize(std::string input_file, std::string output_File_Participants, std::string output_File_Recurrent, std::string output_File_Trial)
{
  // read in parameters
  readParamFile(input_file);

  // specify output files
  output_file_participants = output_File_Participants;
  output_file_trial = output_File_Trial;
  output_file_recurrent_infs = output_File_Recurrent;

  // initialize other necessary parameters
  num_enrolled = 0;
  prob_treatment_arm_allocation = (1.0 * treatment_arm_sample_size) / (1.0 * (treatment_arm_sample_size + placebo_arm_sample_size));
  is_enrollment_open = false;
}



void Trial::readParamFile(std::string input_file)
{
  cout << "Reading in trial parameter file............." << endl;
  cout << endl;

  string discard;

  std::ifstream file_in(input_file);

  if(file_in.fail())
  {
    cout << "Failure reading in data." << endl;
  }

  // read in data
  std::string line;
  while(getline(file_in, line))
  {
    std::string parameter_name;

    std::replace(line.begin(), line.end(), ',', ' ');
    stringstream ss(line);

    // read in parmeter name
    ss >> parameter_name;

    // assign the parameter value to the correct variable on the basis
    // of the parameter name
    if(parameter_name == "treatment_arm_sample_size")
    {
      ss >> treatment_arm_sample_size;
    }
    if(parameter_name == "placebo_arm_sample_size")
    {
      ss >> placebo_arm_sample_size;
    }
    if(parameter_name == "recruitment_start_date")
    {
      ss >> recruitment_start_date;
    }
    if(parameter_name == "trial_duration")
    {
      ss >> trial_duration;
    }
    if(parameter_name == "dropout_rate")
    {
      ss >> dropout_rate;
    }
    if(parameter_name == "trial_PQ_eff")
    {
      ss >> trial_PQ_eff;
    }
    if(parameter_name.find("followup_date_") != string::npos)
    {
      int followup_date;
      ss >> followup_date;
      trial_followup_dates.push_back(followup_date);
    }
  }
}
