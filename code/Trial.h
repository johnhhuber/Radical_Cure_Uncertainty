// Include guard
#ifndef PVIVAX_MODEL_TRIAL
#define PVIVAX_MODEL_TRIAL

#include "Individual.h"
#include "Population.h"
#include "Params.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <tuple>
#include <sstream>


/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// 8.1.1.  Define a class for clinical trials                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

class Trial
{
public:

  // constructor
  Trial(){};

  // destructor
  ~Trial(){};

  // parameters
  bool is_enrollment_open; // boolean specifying whether or not individuals can be enrolled in the trial
  int treatment_arm_sample_size; // sample size for the treatment arm
  int placebo_arm_sample_size; // sample size for the placebo arm

  int num_enrolled;

  double prob_treatment_arm_allocation; // probability of allocating an individual to the treatment arm

  std::vector<int> treatment_arm; // vector of participant ids enrolled in the treatment arm
  std::vector<int> placebo_arm; // vector of participant ids enrolled in the placebo arm

  int recruitment_start_date; // date at which recruitment for the trial begins

  int trial_duration; // duration of the trial (in days)
  std::vector<int> trial_followup_dates; // dates of follow-up for the trial

  double dropout_rate; // probability of an individual dropping out of trial

  double trial_PQ_eff; // efficacy of PQ treatment in the trial
  double trial_PQ_lowage; // lowest age of individuals that can be enrolled in the trial

  std::string output_file_participants;
  std::string output_file_trial;
  std::string output_file_recurrent_infs;

  std::map<int, std::vector<std::tuple<int, bool>>> record_LM_recurrent_infections;
  std::map<int, std::vector<std::tuple<int, string>>> record_all_recurrent_infections;

  std::map<int, std::tuple<string, int, int, double, double>> participant_data;

  // member functions
  void readParamFile(std::string);
  void Initialize(std::string, std::string, std::string, std::string);
  void update(Params&, Population&, int); // function to update the trial for each day
  void enrollParticipants(Population&, int); // function to enroll participants in the trial
  void updateParticipants(Population&, Params&, int); // function to update the status of the participants in the trial
  void administerTreatment(Params&, Individual *, int); // function to treatment individuals of the treatment and participant arm
  void writeParticipantData();
  void writeTrialOutcomes();
  void writeAllRecurrentInfs();
};

#endif
