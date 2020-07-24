/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///                                                                       ///
///  Individual-based Plasmodium vivax transmission model.                ///
///                                                                       ///
///  Dr Michael White                                                     ///
///  Institut Pasteur                                                     ///
///  michael.white@pasteur.fr                                             ///
///                                                                       ///
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

#include "Individual.h"

#include <cmath>
#include "randlib.h"
#include <iostream>

////////////////////////////////////////////////////////////
//                                                        //
//  Function declarations                                 //
//                                                        //
////////////////////////////////////////////////////////////

int CH_sample(double *xx, int nn);
double G6PD_SD_BioSensor(double G6PD_true);


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//        //                                                                               //
// 4.2.1. //  THE MODEL (within humans at least!)                                          //
//        //  Stochastic moves between compartments for each individual                    //
////////////  for a fixed time step                                                        //
////////////                                                                               //
////////////  TO DO: (i) lots of small stuff to speed things up                            //
////////////         (ii) might be able to save on some multiplications by not normalising //
////////////             vector of probabilities and doing it instead inside CH_sample     //
////////////         (iii) add new state for MDA prophylaxis                               //
////////////                                                                               //
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

void Individual::state_mover(Params& theta, double lam_bite)
{
	/////////////////////////////////////////////
	// Track force of infection - it is assumed
	// 8 AQ prophylaxis prevents new blood-stage infections

	if (AQ8_proph == 0)
	{
		lam_bite_track.push_back(lam_bite);
	}
	else
	{
		lam_bite_track.push_back(0.0);
	}

	lam_bite_track.erase(lam_bite_track.begin());

	lam_bite_lag = lam_bite_track[0];


	lam_rel_track.push_back(((double)Hyp)*theta.ff);
	lam_rel_track.erase(lam_rel_track.begin());

	lam_rel_lag = lam_rel_track[0];


	lam_H_lag = lam_bite_lag + lam_rel_lag;


	///////////////////////////////////
	// Indicators for new events

	I_PCR_new = 0;
	I_LM_new = 0;
	I_D_new  = 0;

	CQ_treat = 0;
	PQ_treat = 0;
	TQ_treat = 0;

	G6PD_test = 0;

	CQ_effective = 0;

	PQ_effective = 0;
	PQ_overtreat = 0;
	PQ_overtreat_9m = 0;

	TQ_effective = 0;
	TQ_overtreat = 0;
	TQ_overtreat_9m = 0;

	///////////////////////////////////
	// Probability of exiting states

	double S_out;
	double I_PCR_out;
	double I_LM_out;
	double I_D_out;
	double T_out;
	double P_out;


	///////////////////////////////////////////////////
	// Vector of probabilities for competing hazards
	//
	// These are stored in the parameter structure for convenience.      
	// It would perhaps be more natural to have them specific for each   
	// individual, but that would require storing them N_pop times.      

	double S_move[4];
	double I_PCR_move[5];
	double I_LM_move[4];
	double I_D_move[2];
	double T_move[2];
	double P_move[2];


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  0  //  S: Susceptible                                           //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (S == 1)
	{
		S_out = lam_H_lag;

		if (exp(-t_step * S_out) < genunf(0, 1))
		{
			//theta.r_PCR   = 1.0/( theta.d_PCR_min + (theta.d_PCR_max-theta.d_PCR_min)/( 1.0 + pow((A_par+A_par_mat)*theta.A_PCR_50pc_inv,theta.K_PCR) )); 
			theta.phi_LM = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min)/(1.0 + pow((A_par + A_par_mat)*theta.A_LM_50pc_inv, theta.K_LM));
			theta.phi_D  = theta.phi_D_min  + (theta.phi_D_max - theta.phi_D_min)/(1.0 + pow((A_clin + A_clin_mat)*theta.A_D_50pc_inv, theta.K_D));


			S_move[0] = (1.0 - theta.phi_LM);                                      // Move to I_PCR  //  lam_H_lag*(1.0-theta.phi_LM)/S_out; 
			S_move[1] = theta.phi_LM*(1.0 - theta.phi_D);                          // Move to I_LM   //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)/S_out;
			S_move[2] = theta.phi_LM*theta.phi_D*(1.0 - theta.CM_cover);           // Move to I_D    //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)*(1.0-theta.CM_cover)/S_out;
			S_move[3] = theta.phi_LM*theta.phi_D*theta.CM_cover;                   // Move to T      //  lam_H_lag*theta.phi_LM*(1.0-theta.phi_D)*theta.CM_cover/S_out;

			CH_move = CH_sample(S_move, 4);

			////////////////////////////////
			// S -> I_PCR

			if (CH_move == 0)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_PCR_new = 1;

				T_last_BS = 0.0;

				S = 0;
				I_PCR = 1;

				return;
			}

			////////////////////////////////
			// S -> I_LM

			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_LM_new = 1;
				I_PCR_new = 1;

				T_last_BS = 0.0;

				S = 0;
				I_LM = 1;

				return;
			}

			////////////////////////////////
			// S -> I_D

			if (CH_move == 2)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;

				T_last_BS = 0.0;

				S = 0;
				I_D = 1;

				return;
			}

			////////////////////////////////
			// S -> T

			if (CH_move == 3)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				/////////////////////////////////////////////////////////////////////
				// Case management: treatment of symptomatic cases

				case_management(theta);


				return;
			}
		}
	}


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  1  //  I_PCR: PCR detectable BS infections                      //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (I_PCR == 1)
	{
		theta.r_PCR = 1.0 / (theta.d_PCR_min + (theta.d_PCR_max - theta.d_PCR_min) / (1.0 + pow((A_par + A_par_mat)*theta.A_PCR_50pc_inv, theta.K_PCR)));

		I_PCR_out = lam_H_lag + theta.r_PCR;

		if (exp(-t_step * I_PCR_out) < genunf(0, 1))
		{
			theta.phi_LM = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min) / (1.0 + pow((A_par + A_par_mat)*theta.A_LM_50pc_inv, theta.K_LM));
			theta.phi_D = theta.phi_D_min + (theta.phi_D_max - theta.phi_D_min) / (1.0 + pow((A_clin + A_clin_mat)*theta.A_D_50pc_inv, theta.K_D));

			I_PCR_move[0] = theta.r_PCR / I_PCR_out;                                                           // Move to S 
			I_PCR_move[1] = lam_H_lag*(1 - theta.phi_LM) / I_PCR_out;                                          // Move to I_PCR
			I_PCR_move[2] = lam_H_lag*theta.phi_LM*(1.0 - theta.phi_D) / I_PCR_out;                            // Move to I_LM
			I_PCR_move[3] = lam_H_lag*theta.phi_LM*theta.phi_D*(1.0 - theta.CM_cover) / I_PCR_out;             // Move to I_D
			I_PCR_move[4] = lam_H_lag*theta.phi_LM*theta.phi_D*theta.CM_cover / I_PCR_out;                     // Move to T

			CH_move = CH_sample(I_PCR_move, 5);

			////////////////////////////////
			// I_PCR -> S

			if (CH_move == 0)
			{
				I_PCR = 0;
				S = 1;

				return;
			}

			////////////////////////////////
			// I_PCR -> I_PCR
			// 
			// This is a super-infection event and we assume boosting of immunity

			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				return;
			}

			////////////////////////////////
			// I_PCR -> I_LM

			if (CH_move == 2)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_PCR_new = 1;
				I_LM_new = 1;

				I_PCR = 0;
				I_LM = 1;

				return;
			}

			////////////////////////////////
			// I_PCR -> I_D

			if (CH_move == 3)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;

				I_PCR = 0;
				I_D = 1;

				return;
			}


			////////////////////////////////
			// I_PCR -> T

			if (CH_move == 4)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				/////////////////////////////////////////////////////////////////////
				// Case management: treatment of symptomatic cases

				case_management(theta);


				return;
			}

		}
	}


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  2  //  I_LM: LM detectable BS infection                         //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (I_LM == 1)
	{
		I_LM_out = lam_H_lag + theta.r_LM;

		if (exp(-t_step * I_LM_out) < genunf(0, 1))
		{
			//theta.r_PCR = 1.0/( theta.d_PCR_min + (theta.d_PCR_max-theta.d_PCR_min)/( 1.0 + pow((A_par+A_par_mat)*theta.A_PCR_50pc_inv,theta.K_PCR) ) ); 
			theta.phi_LM = theta.phi_LM_min + (theta.phi_LM_max - theta.phi_LM_min) / (1.0 + pow((A_par + A_par_mat)*theta.A_LM_50pc_inv, theta.K_LM));
			theta.phi_D = theta.phi_D_min + (theta.phi_D_max - theta.phi_D_min) / (1.0 + pow((A_clin + A_clin_mat)*theta.A_D_50pc_inv, theta.K_D));


			I_LM_move[0] = theta.r_LM / I_LM_out;                                                // Move to I_PCR
			I_LM_move[1] = lam_H_lag*(1.0 - theta.phi_D) / I_LM_out;                             // Move to I_LM
			I_LM_move[2] = lam_H_lag*theta.phi_D*(1.0 - theta.CM_cover) / I_LM_out;              // Move to I_D
			I_LM_move[3] = lam_H_lag*theta.phi_D*theta.CM_cover / I_LM_out;                      // Move to T

			CH_move = CH_sample(I_LM_move, 4);

			////////////////////////////////
			// I_LM -> I_PCR

			if (CH_move == 0)
			{
				I_LM = 0;
				I_PCR = 1;

				return;
			}


			////////////////////////////////
			// I_LM -> I_LM
			//
			// Super-infection boosts immunity

			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				return;
			}


			////////////////////////////////
			// I_LM -> I_D

			if (CH_move == 2)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				I_PCR_new = 1;
				I_LM_new = 1;
				I_D_new = 1;

				I_LM = 0;
				I_D = 1;

				return;
			}


			////////////////////////////////
			// I_LM -> T

			if (CH_move == 3)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				/////////////////////////////////////////////////////////////////////
				// Case management: treatment of symptomatic cases

				case_management(theta);

				return;
			}

		}
	}


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  3  //  I_D: Clinical disease                                    //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (I_D == 1)
	{
		I_D_out = lam_H_lag + theta.r_D;

		if (exp(-t_step * I_D_out) < genunf(0, 1))
		{
			I_D_move[0] = theta.r_D / I_D_out;              // Move to I_LM
			I_D_move[1] = lam_H_lag / I_D_out;              // Move to D

			CH_move = CH_sample(I_D_move, 2);

			////////////////////////////////
			// I_D -> I_LM

			if (CH_move == 0)
			{
				I_D = 0;
				I_LM = 1;

				return;
			}


			////////////////////////////////
			// I_D -> I_D

			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				return;
			}

		}
	}


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  4  //  T : Clinical episode under treatment                     //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (T == 1)
	{
		T_out = lam_H_lag + theta.r_T;

		if (exp(-t_step * T_out) < genunf(0, 1))
		{
			T_move[0] = theta.r_T / T_out;              // Move to P
			T_move[1] = lam_H_lag / T_out;              // Move to T

			CH_move = CH_sample(T_move, 2);

			////////////////////////////////
			// T -> P

			if (CH_move == 0)
			{
				T = 0;
				P = 1;

				return;
			}


			////////////////////////////////
			// T -> T
			//
			// Note that we don't allow immune boosting but we do
			// allow for the accumulation of new hypnozoites


			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				return;
			}

		}
	}


	//////////////////////////////////////////////////////////////////////
	//     //                                                           // 
	//  5  //  P: Treatment prophylaxis                                 //
	//     //                                                           //
	//////////////////////////////////////////////////////////////////////

	if (P == 1)
	{
		P_out = lam_H_lag + theta.r_P;

		if (exp(-t_step * P_out) < genunf(0, 1))
		{
			P_move[0] = theta.r_P / P_out;              // Move to S
			P_move[1] = lam_H_lag / P_out;              // Move to P

			CH_move = CH_sample(P_move, 2);

			////////////////////////////////
			// P -> S

			if (CH_move == 0)
			{
				P = 0;
				S = 1;

				return;
			}


			////////////////////////////////
			// P -> P
			//
			// Note that we don't allow immune boosting but we do
			// allow for the accumulation of new hypnozoites


			if (CH_move == 1)
			{
				if (lam_bite_lag / lam_H_lag > genunf(0.0, 1.0))
				{
					if (AQ8_proph == 0)
					{
						Hyp = Hyp + 1;
					}
				}

				if (A_par_boost == 1)
				{
					A_par = A_par + 1.0;
					A_par_timer = theta.u_par;
					A_par_boost = 0;
				}

				if (A_clin_boost == 1)
				{
					A_clin = A_clin + 1.0;
					A_clin_timer = theta.u_clin;
					A_clin_boost = 0;
				}

				return;
			}

		}
	}

}



////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//         //                                                                             //
// 4.2.2.  //  Ageing and immune boosting                                                 //
//         //                                                                             //  
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void Individual::ager(Params& theta)
{
	/////////////////////////
	// Ageing

	age = age + t_step;


	/////////////////////////
	// Loss of hypnozoites

	if (Hyp > K_max)
	{
		Hyp = K_max;
	}

	if (1.0 - exp(-t_step * theta.gamma_L*Hyp) > genunf(0, 1))
	{
		Hyp = Hyp - 1;
	}


	/////////////////////////
	// Immune decay

	A_par  = A_par*theta.A_par_decay;
	A_clin = A_clin*theta.A_clin_decay;


	///////////////////////////////////////////////////
	// Maternal immunity decays exponentially for the first year
	// and is then set to zero afterwards

	if (age < 365.0)
	{
		A_par_mat  = A_par_mat*theta.mat_decay;
		A_clin_mat = A_clin_mat*theta.mat_decay;
	}
	else {
		A_par_mat  = 0.0;
		A_clin_mat = 0.0;
	}


	/////////////////////////
	// Of child-bearing age? ~ 20 years
	// age in (18, 22) years

	if (gender == 2)
	{
		if ((age > theta.preg_age_low) && (age <= theta.preg_age_high))
		{
			preg_age = 1;
		}
		else {
			preg_age = 0;
		}

		if (pregnant == 1)
		{
			preg_timer = preg_timer + 1.0;
		}

		if (pregnant == 0)
		{
			if (preg_age == 1)
			{
				if (genunf(0.0, 1.0) < theta.Preg_daily)
				{
					pregnant   = 1;
					preg_timer = 0.0;
				}
			}
		}

		if (pregnant == 1)
		{
			if (preg_timer > 270.0)
			{
				pregnant = 0;
				preg_timer = 0.0;
			}
		}
	}


	//////////////////////////////////////////////////////
	// Switches for refractory period of immune boosting

	if (A_par_boost == 0)
	{
		A_par_timer = A_par_timer - t_step;

		if (A_par_timer < 0.0)
		{
			A_par_boost = 1;
		}
	}

	if (A_clin_boost == 0)
	{
		A_clin_timer = A_clin_timer - t_step;

		if (A_clin_timer < 0.0)
		{
			A_clin_boost = 1;
		}
	}


	//////////////////////////////////////////////////////
	// Switches for primaquine prophylaxis

	if (AQ8_proph == 1)
	{
		AQ8_proph_timer = AQ8_proph_timer - t_step;

		if (AQ8_proph_timer < 0.0)
		{
			AQ8_proph = 0;
		}
	}


	//////////////////////////////////////////////////////
	// Time since last blood-stage infection

	if (S == 1)
	{
		T_last_BS = T_last_BS + 1.0;
	}
}




////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//         //                                                                             //
// 4.2.3.  //  Case management: treatment of symptomatic cases                            //
//         //                                                                             //  
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void Individual::case_management(Params& theta)
{
	/////////////////////////////////////////////////////////////////////
	// Baseline treatment regimen: chloroquine

	if (theta.CM_regimen == 0)
	{
		CQ_treat = 1;

		I_PCR_new = 1;
		I_LM_new = 1;
		I_D_new = 1;

		G6PD_test = 0;

		T_last_BS = 0.0;

		/////////////////////////////////////////////////////////////////////
		// ACTION: administer CQ

		if (genunf(0.0, 1.0) < theta.CM_CQ_eff)
		{
			S     = 0;
			I_PCR = 0;
			I_LM  = 0;
			I_D   = 0;
			T     = 1;
			P     = 0;
		}else 
		{
			S     = 0;
			I_PCR = 0;
			I_LM  = 1;
			I_D   = 0;
			T     = 0;
			P     = 0;
		}

	}


	/////////////////////////////////////////////////////////////////////
	// Primaquine treatment regimen: blood-stage drugs and primaquine
	// This is based on Narimane's primaquine treatment pathway.

	if (theta.CM_regimen == 1)
	{
		///////////////////////////////////////////////////////////
		// Intitialise

		CQ_treat = 1;
		PQ_treat = 0;

		G6PD_test = 0;

		I_PCR_new = 1;
		I_LM_new  = 1;
		I_D_new   = 1;

		T_last_BS = 0.0;


		/////////////////////////////////////////////////////////////////////
		// Exclude PQ because of young age

		if (age > theta.CM_PQ_lowage)
		{
			PQ_treat = 1;
		}

		/////////////////////////////////////////////////////////////////////
		// Exclude PQ because of pregancy

		if ((theta.CM_PQ_preg_risk == 1) && (pregnant == 1))
		{
			PQ_treat = 0;
		}
		

		//////////////////////////////////////////////////////////////////////////////////////////
		// Is G6PD testing administered to those >6 months and not pregnant? If so, count test

		if ((PQ_treat == 1) && (theta.CM_G6PD_test == 1))
		{
			G6PD_test = 1;
		}


		/////////////////////////////////////////////////////////////////////
		// Exclude PQ because of G6PD deficiency if there is testing

		if( (G6PD_test == 1) && (theta.CM_PQ_G6PD_risk == 1) )
		{
			G6PD_read = G6PD_SD_BioSensor(G6PD_activity);

			if (G6PD_read < 3.0)
			{
				PQ_treat = 0;
			}
		}


		/////////////////////////////////////////////////////////////////////
		// Is PQ effective?

		PQ_effective = 0;

		if (genunf(0.0, 1.0) < theta.CM_PQ_eff)
		{
			PQ_effective = 1;
		}


		/////////////////////////////////////////////////////////////////////
		// Is PQ adhered to?

		if (genunf(0.0, 1.0) < (1 - theta.CM_PQ_adhere))
		{
			PQ_effective = 0;
		}


		/////////////////////////////////////////////////////////////////////
		// Is PQ metabolised?

		if ((theta.CM_PQ_CYP2D6_risk == 1) && (CYP2D6 == 1))
		{
			PQ_effective = 0;
		}


		/////////////////////////////////////////////////////////////////////////////////
		// ACTION: administer PQ 

		if ((PQ_treat == 1) && (PQ_effective == 1))
		{
			Hyp = 0;                                  // Hypnozoites cleared

			AQ8_proph = 1;                            // Put under prophylaxis
			AQ8_proph_timer = theta.CM_PQ_proph;      // Time for prophylaxis set

			// Developing liver hepatic stages killed

			for (int z = 0; z < lam_bite_track.size(); z++)
			{
				lam_bite_track[z] = 0.0;
			}
			for (int z = 0; z < lam_rel_track.size(); z++)
			{
				lam_rel_track[z] = 0.0;
			}
		}		

		/////////////////////////////////////////////////////////////////////////////////
		// ACTION: administer CQ 

		if (PQ_treat == 0)
		{
			if (genunf(0.0, 1.0) < theta.CM_CQ_eff)
			{
				S     = 0;
				I_PCR = 0;
				I_LM  = 0;
				I_D   = 0;
				T     = 1;
				P     = 0;
			}
			else
			{
				S     = 0;
				I_PCR = 0;
				I_LM  = 1;
				I_D   = 0;
				T     = 0;
				P     = 0;
			}
		}

		if (PQ_treat == 1)
		{
			if( genunf(0.0, 1.0) < theta.CM_CQ_eff_wPQ )
			{
				S     = 0;
				I_PCR = 0;
				I_LM  = 0;
				I_D   = 0;
				T     = 1;
				P     = 0;
			}
			else
			{
				S     = 0;
				I_PCR = 0;
				I_LM  = 1;
				I_D   = 0;
				T     = 0;
				P     = 0;
			}
		}

	}


	/////////////////////////////////////////////////////////////////////
	// Tafenoquine treatment regimen: blood-stage drugs and primaquine
	// This is based on Narimane's tafenoquine treatment pathway.

	if (theta.CM_regimen == 2)
	{
		///////////////////////////////////////////////////////////
		// Initialise

		CQ_treat = 1;
		PQ_treat = 0;
		TQ_treat = 0;

		I_PCR_new = 1;
		I_LM_new  = 1;
		I_D_new   = 1;

		G6PD_test = 0;

		T_last_BS = 0.0;


		/////////////////////////////////////////////////////////////////////
		// Primaquine administered to people 6 months to 16 years

		if ((age > theta.CM_PQ_lowage) && (age < theta.CM_TQ_lowage))
		{
			PQ_treat = 1;

			/////////////////////////////////////////////////////////////////////
			// Exclude PQ because of pregancy

			if ((theta.CM_PQ_preg_risk == 1) && (pregnant == 1))
			{
				PQ_treat = 0;
			}

			//////////////////////////////////////////////////////////////////////////////////////////
			// Is G6PD testing administered to those >6 months and not pregnant? If so, count test

			if ((PQ_treat == 1) && (theta.CM_G6PD_test == 1))
			{
				G6PD_test = 1;
			}

			/////////////////////////////////////////////////////////////////////
			// Exclude PQ because of G6PD deficiency if there is testing

			if ((G6PD_test == 1) && (theta.CM_PQ_G6PD_risk == 1))
			{
				G6PD_read = G6PD_SD_BioSensor(G6PD_activity);

				if (G6PD_read < 3.0)
				{
					PQ_treat = 0;
				}
			}


			/////////////////////////////////////////////////////////////////////
			// Is PQ effective?

			PQ_effective = 0;

			if (genunf(0.0, 1.0) < theta.CM_PQ_eff)
			{
				PQ_effective = 1;
			}

			/////////////////////////////////////////////////////////////////////
			// Is PQ adhered to?

			if (genunf(0.0, 1.0) < (1 - theta.CM_PQ_adhere))
			{
				PQ_effective = 0;
			}

			/////////////////////////////////////////////////////////////////////
			// Is PQ metabolised?

			if ((theta.CM_PQ_CYP2D6_risk == 1) && (CYP2D6 == 1))
			{
				PQ_effective = 0;
			}

			/////////////////////////////////////////////////////////////////////////////////
			// ACTION: administer PQ

			if ((PQ_treat == 1) && (PQ_effective == 1))
			{
				Hyp = 0;                                // Hypnozoites cleared

				AQ8_proph = 1;                          // Put under prophylaxis
				AQ8_proph_timer = theta.CM_PQ_proph;    // Timer for prophylaxis set

				// Developing liver hepatic stages killed

				for (int z = 0; z < lam_bite_track.size(); z++)
				{
					lam_bite_track[z] = 0.0;
				}
				for (int z = 0; z < lam_rel_track.size(); z++)
				{
					lam_rel_track[z] = 0.0;
				}
			}
		}


		/////////////////////////////////////////////////////////////////////
		// Tafenoquine administered to people older than 16 years

		if (age >= theta.CM_TQ_lowage)
		{
			TQ_treat = 1;
			PQ_treat = 0;

			/////////////////////////////////////////////////////////////////////
			// Exclude TQ because of pregancy

			if ((theta.CM_TQ_preg_risk == 1) && (pregnant == 1))
			{
				TQ_treat = 0;
			}

			//////////////////////////////////////////////////////////////////////////////////////////
			// Is G6PD testing is administered to those >16 years and not pregnant? If so, count test

			if ((TQ_treat == 1) && (theta.CM_G6PD_test == 1))
			{
				G6PD_test = 1;
			}

			/////////////////////////////////////////////////////////////////////
			// Exclude TQ because of G6PD deficiency

			if ((G6PD_test == 1) && (theta.CM_PQ_G6PD_risk == 1))
			{
				G6PD_read = G6PD_SD_BioSensor(G6PD_activity);

				if (G6PD_read < 7.0)
				{
					TQ_treat = 0;

					if (G6PD_read > 3.0)
					{
						PQ_treat = 1;
					}
				}
			}


			/////////////////////////////////////////////////////////////////////
			// Is TQ effective?

			TQ_effective = 0;

			if (TQ_treat == 1)
			{
				if (genunf(0.0, 1.0) < theta.CM_TQ_eff)
				{
					TQ_effective = 1;
				}

				/////////////////////////////////////////////////////////////////////
				// Is TQ adhered to?

				if (genunf(0.0, 1.0) < 1 - theta.CM_TQ_adhere)
				{
					TQ_effective = 0;
				}

				////////////////////////////////////////////////////////////////////
				// Is TQ metabolised?

				if ((theta.CM_TQ_CYP2D6_risk == 1) && (CYP2D6 == 1))
				{
					TQ_effective = 0;
				}
			}

			/////////////////////////////////////////////////////////////////////
			// Is PQ effective?

			PQ_effective = 0;

			if (PQ_treat == 1)
			{
				if (genunf(0.0, 1.0) < theta.CM_PQ_eff)
				{
					PQ_effective = 1;
				}

				/////////////////////////////////////////////////////////////////////
				// Is PQ adhered to?

				if (genunf(0.0, 1.0) < 1 - theta.CM_PQ_adhere)
				{
					PQ_effective = 0;
				}

				////////////////////////////////////////////////////////////////////
				// Is PQ metabolised?

				if ((theta.CM2_PQ_CYP2D6_risk == 1) && (CYP2D6 == 1))
				{
					PQ_effective = 0;
				}
			}

			/////////////////////////////////////////////////////////////////////////////////
			// ACTION: administer PQ or TQ

			if ((TQ_treat == 1) && (TQ_effective == 1))
			{
				Hyp = 0;                               // Hypnozoites cleared

				AQ8_proph = 1;                         // Put under prophylaxis
				AQ8_proph_timer = theta.CM_TQ_proph;   // Timer for prophylaxis set  

				// Developing liver hepatic stages killed

				for (int z = 0; z < lam_bite_track.size(); z++)
				{
					lam_bite_track[z] = 0.0;
				}
				for (int z = 0; z < lam_rel_track.size(); z++)
				{
					lam_rel_track[z] = 0.0;
				}
			}

			if ((PQ_treat == 1) && (PQ_effective == 1))
			{
				Hyp = 0;                                   // Hypnozoites cleared

				AQ8_proph = 1;                             // Put under prophylaxis
				AQ8_proph_timer = theta.CM_PQ_proph;       // Timer for prophylaxis set

				// Developing liver hepatic stages killed

				for (int z = 0; z < lam_bite_track.size(); z++)
				{
					lam_bite_track[z] = 0.0;
				}
				for (int z = 0; z < lam_rel_track.size(); z++)
				{
					lam_rel_track[z] = 0.0;
				}
			}
		}


		/////////////////////////////////////////////////////////////////////////////////
		// ACTION: administer CQ

		if ((PQ_treat == 0) && (TQ_treat == 0))
		{
			if (genunf(0.0, 1.0) < theta.CM_CQ_eff)
			{
				S     = 0;
				I_PCR = 0;
				I_LM  = 0;
				I_D   = 0;
				T     = 1;
				P     = 0;
			}
			else
			{
				S     = 0;
				I_PCR = 0;
				I_LM  = 1;
				I_D   = 0;
				T     = 0;
				P     = 0;
			}
		}

		if( PQ_treat == 1 )
		{
			if (genunf(0.0, 1.0) < theta.CM_CQ_eff_wPQ)
			{
				S     = 0;
				I_PCR = 0;
				I_LM  = 0;
				I_D   = 0;
				T     = 1;
				P     = 0;
			}
			else
			{
				S     = 0;
				I_PCR = 0;
				I_LM  = 1;
				I_D   = 0;
				T     = 0;
				P     = 0;
			}
		}

		if( TQ_treat == 1 )
		{
			S     = 0;
			I_PCR = 0;
			I_LM  = 0;
			I_D   = 0;
			T     = 1;
			P     = 0;
		}

	}

}



////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//         //                                                                             //
// 4.2.4.  //  Individual-level vector control updater                                    //
//         //                                                                             //
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void Individual::vec_con_updater(Params& theta)
{

	///////////////////////////////////////////
	// Is LLIN lost

	if (LLIN == 1)
	{
		if (theta.P_LLIN_loss > genunf(0, 1))
		{
			LLIN = 0;
		}
	}

	///////////////////////////////////////////
	// IRS turned off if sufficiently long ago (3 years???)

	if (IRS_age > 1095.0)
	{
		IRS = 0;
	}

	///////////////////////////////////////////
	// IVM turned off if sufficiently long ago (30 days???)

	if (IVM_age > 180.0)
	{
		IVM = 0;
	}



	if (LLIN == 1 && IRS == 0 && IVM == 0)
	{
		LLIN_age = LLIN_age + t_step;

		for (int v = 0; v < N_mosq; v++)
		{
			r_LLIN[v] = theta.r_LLIN_net[v] + (r_LLIN[v] - theta.r_LLIN_net[v])*theta.P_PYR_decay;
			d_LLIN[v] = d_LLIN[v] * theta.P_PYR_decay;
			s_LLIN[v] = 1.0 - r_LLIN[v] - d_LLIN[v];

			w_VC[v] = 1 - theta.PSI_bed[v] + theta.PSI_bed[v] * s_LLIN[v];
			y_VC[v] = w_VC[v];
			z_VC[v] = theta.PSI_bed[v] * r_LLIN[v];
		}
	}

	if (LLIN == 0 && IRS == 1 && IVM == 0)
	{
		IRS_age = IRS_age + t_step;

		for (int v = 0; v < N_mosq; v++)
		{
			r_IRS[v] = r_IRS[v] * theta.P_IRS_decay;
			d_IRS[v] = d_IRS[v] * theta.P_IRS_decay;
			s_IRS[v] = 1.0 - r_IRS[v] - d_IRS[v];

			w_VC[v] = 1 - theta.PSI_indoors[v] + theta.PSI_indoors[v] * (1.0 - r_IRS[v])*s_IRS[v];
			y_VC[v] = 1 - theta.PSI_indoors[v] + theta.PSI_indoors[v] * (1.0 - r_IRS[v]);
			z_VC[v] = theta.PSI_indoors[v] * r_IRS[v];
		}
	}

	if (LLIN == 1 && IRS == 1 && IVM == 0)
	{
		LLIN_age = LLIN_age + t_step;
		IRS_age = IRS_age + t_step;

		for (int v = 0; v < N_mosq; v++)
		{
			r_LLIN[v] = theta.r_LLIN_net[v] + (r_LLIN[v] - theta.r_LLIN_net[v])*theta.P_PYR_decay;
			d_LLIN[v] = d_LLIN[v] * theta.P_PYR_decay;
			s_LLIN[v] = 1.0 - r_LLIN[v] - d_LLIN[v];

			r_IRS[v] = r_IRS[v] * theta.P_IRS_decay;
			d_IRS[v] = d_IRS[v] * theta.P_IRS_decay;
			s_IRS[v] = 1.0 - r_IRS[v] - d_IRS[v];

			w_VC[v] = 1.0 - theta.PSI_indoors[v] + theta.PSI_bed[v] * (1.0 - r_IRS[v])*s_LLIN[v] * s_IRS[v] + (theta.PSI_indoors[v] - theta.PSI_bed[v])*(1.0 - r_IRS[v])*s_IRS[v];
			y_VC[v] = 1.0 - theta.PSI_indoors[v] + theta.PSI_bed[v] * (1.0 - r_IRS[v])*s_LLIN[v] + (theta.PSI_indoors[v] - theta.PSI_bed[v])*(1.0 - r_IRS[v]);
			z_VC[v] = theta.PSI_bed[v] * (1.0 - r_IRS[v])*r_LLIN[v] + theta.PSI_indoors[v] * r_IRS[v];
		}
	}

	if (LLIN == 0 && IRS == 0 && IVM == 0)
	{
		for (int v = 0; v < N_mosq; v++)
		{
			w_VC[v] = 1.0;
			y_VC[v] = 1.0;
			z_VC[v] = 0.0;
		}
	}


	if (LLIN == 1 && IRS == 0 && IVM == 1)
	{
		LLIN_age = LLIN_age + t_step;
		IVM_age  = IVM_age + t_step;


		for (int v = 0; v < N_mosq; v++)
		{
			r_LLIN[v] = theta.r_LLIN_net[v] + (r_LLIN[v] - theta.r_LLIN_net[v])*theta.P_PYR_decay;
			d_LLIN[v] = d_LLIN[v] * theta.P_PYR_decay;
			s_LLIN[v] = 1.0 - r_LLIN[v] - d_LLIN[v];

			d_IVM = d_IVM * theta.P_IVM_decay;
			s_IVM = 1.0 - d_IVM;

			w_VC[v] = s_IVM * (1 - theta.PSI_bed[v] + theta.PSI_bed[v] * s_LLIN[v]);
			y_VC[v] = w_VC[v];
			z_VC[v] = theta.PSI_bed[v] * r_LLIN[v];
		}
	}

	if (LLIN == 0 && IRS == 1 && IVM == 1)
	{
		IRS_age = IRS_age + t_step;
		IVM_age = IVM_age + t_step;

		for (int v = 0; v < N_mosq; v++)
		{
			r_IRS[v] = r_IRS[v] * theta.P_IRS_decay;
			d_IRS[v] = d_IRS[v] * theta.P_IRS_decay;
			s_IRS[v] = 1.0 - r_IRS[v] - d_IRS[v];

			d_IVM = d_IVM * theta.P_IVM_decay;
			s_IVM = 1.0 - d_IVM;

			w_VC[v] = s_IVM * (1 - theta.PSI_indoors[v] + theta.PSI_indoors[v] * (1.0 - r_IRS[v])*s_IRS[v]);
			y_VC[v] = 1 - theta.PSI_indoors[v] + theta.PSI_indoors[v] * (1.0 - r_IRS[v]);      // ivermectin will not stop the mosquito feeding
			z_VC[v] = theta.PSI_indoors[v] * r_IRS[v];
		}
	}

	if (LLIN == 1 && IRS == 1 && IVM == 1)
	{
		LLIN_age = LLIN_age + t_step;
		IRS_age  = IRS_age + t_step;
		IVM_age  = IVM_age + t_step;

		for (int v = 0; v < N_mosq; v++)
		{
			r_LLIN[v] = theta.r_LLIN_net[v] + (r_LLIN[v] - theta.r_LLIN_net[v])*theta.P_PYR_decay;
			d_LLIN[v] = d_LLIN[v] * theta.P_PYR_decay;
			s_LLIN[v] = 1.0 - r_LLIN[v] - d_LLIN[v];

			r_IRS[v] = r_IRS[v] * theta.P_IRS_decay;
			d_IRS[v] = d_IRS[v] * theta.P_IRS_decay;
			s_IRS[v] = 1.0 - r_IRS[v] - d_IRS[v];

			d_IVM = d_IVM * theta.P_IVM_decay;
			s_IVM = 1.0 - d_IVM;

			w_VC[v] = s_IVM * (1.0 - theta.PSI_indoors[v] + theta.PSI_bed[v] * (1.0 - r_IRS[v])*s_LLIN[v] * s_IRS[v] + (theta.PSI_indoors[v] - theta.PSI_bed[v])*(1.0 - r_IRS[v])*s_IRS[v]);
			y_VC[v] = 1.0 - theta.PSI_indoors[v] + theta.PSI_bed[v] * (1.0 - r_IRS[v])*s_LLIN[v] + (theta.PSI_indoors[v] - theta.PSI_bed[v])*(1.0 - r_IRS[v]);        // ivermectin will not stop the mosquito feeding
			z_VC[v] = theta.PSI_bed[v] * (1.0 - r_IRS[v])*r_LLIN[v] + theta.PSI_indoors[v] * r_IRS[v];
		}
	}

	if (LLIN == 0 && IRS == 0 && IVM == 1)
	{
		IVM_age = IVM_age + t_step;

		d_IVM = d_IVM * theta.P_IVM_decay;
		s_IVM = 1.0 - d_IVM;

		for (int v = 0; v < N_mosq; v++)
		{
			w_VC[v] = s_IVM;
			y_VC[v] = 1.0;
			z_VC[v] = 0.0;
		}
	}

}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  4.2.5. G6PD diagnosis                                                   //
//         Based on data from Pal et al AJTMH 2019 study of US venous blood //  
//         The values returned (1.5, 5.0, 10.0) are representative          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

double G6PD_SD_BioSensor(double G6PD_true)
{
	double Pal_30_70[3] = { 0.4230769, 0.3846154, 0.1923077};
	double G6PD_read;
	int CH_move;

	//////////////////////////////////////
	// CASE 1: True G6PD activity < 30%

	if (G6PD_true <= 3.0)
	{
		G6PD_read = 1.5;
	}


	//////////////////////////////////////
	// CASE 2: True G6PD activity > 30% & < 70%

	if ((G6PD_true > 3.0) && (G6PD_true <= 7.0))
	{
		CH_move = CH_sample(Pal_30_70, 3);

		if (CH_move == 0)
		{
			G6PD_read = 1.5;
		}
		if (CH_move == 1)
		{
			G6PD_read = 5.0;
		}
		if (CH_move == 2)
		{
			G6PD_read = 10.0;
		}
	}

	//////////////////////////////////////
	// CASE 3: True G6PD activity > 70%

	if (G6PD_true > 7.0)
	{
		if (genunf(0.0, 1.0) < 0.9703947)
		{
			G6PD_read = 10.0;
		}
		else
		{
			G6PD_read = 5.0;
		}
	}

	//////////////////////////////////////
	// Return reading

	return G6PD_read;

}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  4.2.6. Competing hazards sampler                                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

int CH_sample(double *xx, int nn)
{
	vector<double> xx_cum(nn);

	xx_cum[0] = xx[0];

	for (int k = 1; k < nn; k++)
	{
		xx_cum[k] = xx_cum[k - 1] + xx[k];
	}

	int index = 0;
	double unif = genunf(0, 1);

	if (unif < xx_cum[0])
	{
		return index;
	}

	for (int k = 1; k < nn; k++)
	{
		if ((unif > xx_cum[k - 1]) & (unif < xx_cum[k]))
		{
			index = k;

			return index;
		}
	}

	return index;
}
