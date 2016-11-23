#include "scalehmm.h"

// ============================================================
// Hidden Markov Model implemented with scaling strategy
// ============================================================

// Public =====================================================

// Constructor and Destructor ---------------------------------
ScaleHMM::ScaleHMM(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & distances, Rcpp::NumericVector startProbs_initial, Rcpp::NumericMatrix transProbs_initial, double transDist, Rcpp::DataFrame emissionParams_initial, int verbosity)
{
	if (verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	this->xvariate = UNIVARIATE;
	this->verbosity = verbosity;
	this->NDATA = obs.size();
	this->NSTATES = startProbs_initial.size();
	this->distances = distances;
	this->scalefactoralpha = Rcpp::NumericVector(this->NDATA);
	this->scalealpha = Rcpp::NumericMatrix(this->NDATA, this->NSTATES);
	this->scalebeta = Rcpp::NumericMatrix(this->NDATA, this->NSTATES);
	this->densities = Rcpp::NumericMatrix(this->NSTATES, this->NDATA);
	this->gamma = Rcpp::NumericMatrix(this->NSTATES, this->NDATA);
	this->sumgamma = Rcpp::NumericVector(this->NSTATES);
	this->sumxi = Rcpp::NumericMatrix(this->NSTATES, this->NSTATES);
	this->loglik = -INFINITY;
	this->dloglik = INFINITY;

	// Initial probabilities
	this->transProbs = Rcpp::clone(transProbs_initial);
	this->transExp = Rcpp::NumericVector(this->NDATA);
	for (int t=0; t<this->NDATA; t++)
	{
		this->transExp[t] = exp(- this->distances[t] / transDist);
		if (std::isnan(this->transExp[t]))
		{
			throw nan_detected;
		}
	}
	this->startProbs = Rcpp::clone(startProbs_initial);

	// Make log version of observations for faster updating
	this->logObs = Rcpp::NumericVector(this->NDATA);
	this->log1mObs = Rcpp::NumericVector(this->NDATA);
	double cutoff = 1e-10; // cutoff to avoid -Inf in log
	double log_cutoff = log(cutoff);
	double log_1m_cutoff = log(1-cutoff);
	for(int t=0; t<obs.size(); t++)
	{
		if (obs[t] <= cutoff)
		{
			this->logObs[t] = log_cutoff;
			this->log1mObs[t] = log_1m_cutoff;
		}
		else if (obs[t] >= 1-cutoff)
		{
			this->logObs[t] = log_1m_cutoff;
			this->log1mObs[t] = log_cutoff;
		}
		else
		{
			this->logObs[t] = log(obs[t]);
			this->log1mObs[t] = log(1-obs[t]);
		}
	}

	// Emission densities
	this->emissionParams = Rcpp::clone(emissionParams_initial);
	Rcpp::CharacterVector emissionTypes = this->emissionParams["type"];
	Rcpp::NumericVector a_s = Rcpp::as<Rcpp::NumericVector>(this->emissionParams["a"]);
	Rcpp::NumericVector b_s = Rcpp::as<Rcpp::NumericVector>(this->emissionParams["b"]);
	// UNmethylated
	Beta_mirror * d0 = new Beta_mirror(obs, this->logObs, this->log1mObs, a_s[0], b_s[0]);
	this->emissionDensities.push_back(d0);
	// Hemimethylated
	Beta_symmetric * d1 = new Beta_symmetric(obs, this->logObs, this->log1mObs, a_s[1], b_s[1]);
	this->emissionDensities.push_back(d1);
	// Methylated
	Beta_mirror * d2 = new Beta_mirror(obs, this->logObs, this->log1mObs, a_s[2], b_s[2]);
	this->emissionDensities.push_back(d2);
}

ScaleHMM::ScaleHMM(const Rcpp::IntegerVector & obs, const Rcpp::NumericVector & distances, Rcpp::NumericVector startProbs_initial, Rcpp::NumericMatrix transProbs_initial, double transDist, Rcpp::DataFrame emissionParams_initial, int verbosity)
{
	if (verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	this->xvariate = UNIVARIATE;
	this->verbosity = verbosity;
	this->NDATA = obs.size();
	this->NSTATES = startProbs_initial.size();
	this->distances = distances;
	this->scalefactoralpha = Rcpp::NumericVector(this->NDATA);
	this->scalealpha = Rcpp::NumericMatrix(this->NDATA, this->NSTATES);
	this->scalebeta = Rcpp::NumericMatrix(this->NDATA, this->NSTATES);
	this->densities = Rcpp::NumericMatrix(this->NSTATES, this->NDATA);
	this->gamma = Rcpp::NumericMatrix(this->NSTATES, this->NDATA);
	this->sumgamma = Rcpp::NumericVector(this->NSTATES);
	this->sumxi = Rcpp::NumericMatrix(this->NSTATES, this->NSTATES);
	this->loglik = -INFINITY;
	this->dloglik = INFINITY;

	// Initial probabilities
	this->transProbs = Rcpp::clone(transProbs_initial);
	this->transExp = Rcpp::NumericVector(this->NDATA);
	for (int t=0; t<this->NDATA; t++)
	{
		this->transExp[t] = exp(-this->distances[t] / transDist);
		if (std::isnan(this->transExp[t]))
		{
			throw nan_detected;
		}
	}
	this->startProbs = Rcpp::clone(startProbs_initial);

	// Make vector of positions of unique observations for faster updating
	this->obs_unique = Rcpp::unique(obs);
	this->obs_unique.sort();
	Rcpp::IntegerVector uobsind_per_obs = Rcpp::IntegerVector(this->obs_unique[this->obs_unique.size()-1]+1);
	int i = 0;
	for (int j=0; j<uobsind_per_obs.size(); j++)
	{
		if (this->obs_unique[i] == j)
		{
			uobsind_per_obs[j] = i;
			i += 1;
		}
		if (this->verbosity>=5) Rprintf("uobsind_per_obs[j=%d] = %d\n", j, uobsind_per_obs[j]);
	}
	// Get the index for each observation
	this->uobsind_per_t = Rcpp::IntegerVector(obs.size());
	for (int t=0; t<obs.size(); t++)
	{
		this->uobsind_per_t[t] = uobsind_per_obs[obs[t]];
	}

	// Emission densities
	this->emissionParams = Rcpp::clone(emissionParams_initial);
	Rcpp::CharacterVector emissionTypes = this->emissionParams["type"];
	Rcpp::NumericVector sizes = Rcpp::as<Rcpp::NumericVector>(this->emissionParams["size"]);
	Rcpp::NumericVector probs = Rcpp::as<Rcpp::NumericVector>(this->emissionParams["prob"]);
	for (int i=0; i<this->NSTATES; i++)
	{
		std::string dtype = Rcpp::as<std::string>(emissionTypes[i]);
		if (dtype.compare("delta") == 0)
		{
			// Zero Inflation
			ZeroInflation * d = new ZeroInflation(obs);
			this->emissionDensities.push_back(d);
		}
		else if (dtype.compare("dnbinom") == 0)
		{
			// Negative Binomial
			NegativeBinomial * d = new NegativeBinomial(obs, this->obs_unique, this->uobsind_per_t, sizes[i], probs[i]);
			this->emissionDensities.push_back(d);
		}
	}

}

ScaleHMM::ScaleHMM(const Rcpp::IntegerMatrix & multi_obs, const Rcpp::NumericVector & distances, Rcpp::NumericVector startProbs_initial, Rcpp::NumericMatrix transProbs_initial, double transDist, Rcpp::List emissionParamsList, int verbosity, const Rcpp::List & cor_mat_inv, const Rcpp::NumericVector & determinant, const Rcpp::DataFrame & statedef)
{
	if (verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	this->xvariate = MULTIVARIATE;
	this->verbosity = verbosity;
	this->NDATA = multi_obs.nrow();
	this->NSTATES = startProbs_initial.size();
	this->NMOD = multi_obs.ncol();
	this->distances = distances;
	this->transProbs = Rcpp::NumericMatrix(this->NSTATES, this->NSTATES);
	this->scalefactoralpha = Rcpp::NumericVector(this->NDATA);
	this->scalealpha = Rcpp::NumericMatrix(this->NDATA, this->NSTATES);
	this->scalebeta = Rcpp::NumericMatrix(this->NDATA, this->NSTATES);
	this->densities = Rcpp::NumericMatrix(this->NSTATES, this->NDATA);
// 	this->tdensities = Rcpp::NumericMatrix(this->NDATA, this->NSTATES);
	this->startProbs = Rcpp::NumericVector(this->NSTATES);
	this->gamma = Rcpp::NumericMatrix(this->NSTATES, this->NDATA);
	this->sumgamma = Rcpp::NumericVector(this->NSTATES);
	this->sumxi = Rcpp::NumericMatrix(this->NSTATES, this->NSTATES);
	this->loglik = -INFINITY;
	this->dloglik = INFINITY;

	// Initial probabilities
	this->transProbs = Rcpp::clone(transProbs_initial);
	this->transExp = Rcpp::NumericVector(this->NDATA);
	for (int t=0; t<this->NDATA; t++)
	{
		this->transExp[t] = exp(-this->distances[t] / transDist);
		if (std::isnan(this->transExp[t]))
		{
			throw nan_detected;
		}
	}
	this->startProbs = Rcpp::clone(startProbs_initial);

	// Emission densities
	int ndistr = emissionParamsList.size();
	Rcpp::IntegerVector state_def;
	this->emissionParamsList = emissionParamsList;
	for (int i=0; i<this->NSTATES; i++)
	{
		// Get row of state definition
		for (int j=0; j<ndistr; j++)
		{
			Rcpp::IntegerVector jcol = statedef[j];
			state_def[j] = jcol[i];
		}
		MVCopulaApproximation * MVdens = new MVCopulaApproximation(multi_obs, state_def, this->emissionParamsList, Rcpp::as<Rcpp::NumericMatrix>(cor_mat_inv[i]), determinant[i]);
		this->emissionDensities.push_back(MVdens);
	}

}

ScaleHMM::~ScaleHMM()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	for (int i=0; i<this->emissionDensities.size(); i++)
	{
		delete this->emissionDensities[i];
	}
}

// Methods ----------------------------------------------------
Rcpp::List ScaleHMM::viterbi(double eps, double maxiter, double maxtime)
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	// measuring the time
	this->baumWelchStartTime_sec = time(NULL);

	if (this->xvariate == UNIVARIATE)
	{
		this->print_uni_iteration(0);
	}
	else if (this->xvariate == MULTIVARIATE)
	{
		this->print_multi_iteration(0);
	}
	R_CheckUserInterrupt();

	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
		
	try { this->calc_densities(); } catch(...) { throw; }
	R_CheckUserInterrupt();

	try { this->forward(); } catch(...) { throw; }
	R_CheckUserInterrupt();

	try { this->backward(); } catch(...) { throw; }
	R_CheckUserInterrupt();

	this->calc_loglikelihood();
	if(std::isnan(this->loglik))
	{
		throw nan_detected;
	}

	this->calc_sumxi();
	R_CheckUserInterrupt();

	this->calc_sumgamma();
	R_CheckUserInterrupt();

	if (this->xvariate == UNIVARIATE)
	{
		this->print_uni_iteration(1);
	}
	else if (this->xvariate == MULTIVARIATE)
	{
		this->print_multi_iteration(1);
	}

	/* Aggregate information to return */

	// Calculate posterior weights
	Rcpp::NumericVector weights = this->calc_weights();

	// Maximize over posteriors
	Rcpp::IntegerVector imax(this->NDATA);
	double maxgamma;
	for (int t=0; t<this->NDATA; t++)
	{
		maxgamma = -1;
		for (int i=0; i<this->NSTATES; i++)
		{
			if (this->gamma(i,t) > maxgamma)
			{
				maxgamma = this->gamma(i,t);
				imax[t] = i;
			}
		}
	}

	// Convergence information
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	Rcpp::List convergenceInfo = Rcpp::List::create(Rcpp::Named("logliks") = this->loglik,
																					 				Rcpp::Named("time") = this->baumWelchTime_real);

	// Collect results in list
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("convergenceInfo") = convergenceInfo,
																				 Rcpp::Named("transProbs") = this->transProbs,
																				 Rcpp::Named("startProbs") = this->startProbs,
																				 Rcpp::Named("weights") = weights,
																				 Rcpp::Named("posteriors") = this->gamma,
																				 Rcpp::Named("states") = imax,
																				 Rcpp::Named("densities") = this->densities);

	// Emission parameters
	if (this->xvariate == UNIVARIATE)
	{
		result.push_back(this->emissionParams, "emissionParams");
	}

	return result;
										 
}


Rcpp::List ScaleHMM::baumWelch(double eps, double maxiter, double maxtime)
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);

	// Parallelization settings
// 	#ifdef _OPENMP
// 	omp_set_nested(1);
// 	#endif
	
	// measuring the time
	this->baumWelchStartTime_sec = time(NULL);

	if (this->xvariate == UNIVARIATE)
	{
		this->print_uni_iteration(0);
	}
	else if (this->xvariate == MULTIVARIATE)
	{
		this->print_multi_iteration(0);
		if (this->verbosity>=1) Rprintf("HMM: Precomputing densities ...\n");
		try { this->calc_densities(); } catch(...) { throw; }
		this->print_multi_iteration(0);
		// Print densities
// 		int bs = 100;
// 		char buffer [100];
// 		int cx;
// 		for (int t=0; t<this->NDATA; t++)
// 		{
// 			cx = 0;
// 			cx += snprintf(buffer+cx, bs-cx, "t=%d\t", t);
// 			for (int i=0; i<this->NSTATES; i++)
// 			{
// 				cx += snprintf(buffer+cx, bs-cx, "%.6f\t", this->densities(i,t));
// 			}
// 		}
				
	}

	R_CheckUserInterrupt();
	// Do the Baum-Welch
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	int iteration = 0;
	std::vector<double> logliks;
	logliks.push_back(-INFINITY);
	while (((this->baumWelchTime_real < maxtime) or (maxtime < 0)) and ((iteration < maxiter) or (maxiter < 0)))
	{

		iteration++;
		
		if (this->xvariate == UNIVARIATE)
		{
			try { this->calc_densities(); } catch(...) { throw; }
			R_CheckUserInterrupt();
		}

		try { this->forward(); } catch(...) { throw; }
		R_CheckUserInterrupt();

		try { this->backward(); } catch(...) { throw; }
		R_CheckUserInterrupt();

		this->calc_loglikelihood();
		logliks.push_back(this->loglik);
		if(std::isnan(logliks[iteration]))
		{
			throw nan_detected;
		}
		this->dloglik = logliks[iteration] - logliks[iteration-1];

		this->calc_sumxi();
		R_CheckUserInterrupt();

		this->calc_sumgamma();
		R_CheckUserInterrupt();

		// Print information about current iteration
		if (this->xvariate == UNIVARIATE)
		{
			this->print_uni_iteration(iteration);
		}
		else if (this->xvariate == MULTIVARIATE)
		{
			this->print_multi_iteration(iteration);
		}

		// Check convergence
		if(this->dloglik < eps) //it has converged
		{
			if (this->verbosity>=1) Rprintf("HMM: Convergence reached!\n");
			break;
		} else {// not converged
			this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
			if (iteration == maxiter)
			{
				if (this->verbosity>=1) Rprintf("HMM: Maximum number of iterations reached!\n");
			}
			else if ((this->baumWelchTime_real >= maxtime) and (maxtime >= 0))
			{
				if (this->verbosity>=1) Rprintf("HMM: Exceeded maximum time!\n");
			}
		}
		
		// Updating initial probabilities startProbs and transition matrix transProbs
		if (this->verbosity>=2) Rprintf("  updating transition matrix\n");
		for (int i=0; i<this->NSTATES; i++)
		{
			this->startProbs[i] = this->gamma(i,0);
			if (this->verbosity>=4) Rprintf("  startProbs[%d] = %g\n", i, startProbs[i]);
			if (this->sumgamma[i] == 0)
			{
// 				Rprintf("Not reestimating transProbs(%d,x) because sumgamma[%d] = 0\n", i, i);
			}
			else
			{
				for (int j=0; j<this->NSTATES; j++)
				{
					this->transProbs(i,j) = this->sumxi(i,j) / this->sumgamma[i];
					if (std::isnan(this->transProbs(i,j)))
					{
						throw nan_detected;
					}
				}
			}
		}

		if (this->xvariate == UNIVARIATE)
		{
			// Update the parameters of the distribution
			if (this->verbosity>=2) Rprintf("  updating distribution parameters\n");
			if (this->emissionDensities[0]->get_name() == BETA_MIRROR)
			{
				const int rows1[] = {0,2};
				this->emissionDensities[0]->update(this->gamma, rows1);
				this->emissionDensities[2]->set_a(this->emissionDensities[0]->get_b());
				this->emissionDensities[2]->set_b(this->emissionDensities[0]->get_a());
				const int rows2[] = {1};
				this->emissionDensities[1]->update(this->gamma, rows2);
			}
			else if ((this->emissionDensities[0]->get_name() == NEGATIVE_BINOMIAL) | (this->emissionDensities[0]->get_name() == ZERO_INFLATION))
			{ 
				for (int i=0; i<this->NSTATES; i++)
				{
					const int rows[] = {i};
					this->emissionDensities[i]->update(this->gamma, rows);
				}
			}
			R_CheckUserInterrupt();
		}

	} /* main loop end */
    
	/* Aggregate information to return */

	// Calculate posterior weights
	Rcpp::NumericVector weights = this->calc_weights();

	// Maximize over posteriors
	Rcpp::IntegerVector imax(this->NDATA);
	double maxgamma;
	for (int t=0; t<this->NDATA; t++)
	{
		maxgamma = -1;
		for (int i=0; i<this->NSTATES; i++)
		{
			if (this->gamma(i,t) > maxgamma)
			{
				maxgamma = this->gamma(i,t);
				imax[t] = i;
			}
		}
	}

	// Convergence information
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	Rcpp::List convergenceInfo = Rcpp::List::create(Rcpp::Named("logliks") = logliks,
																					 				Rcpp::Named("time") = this->baumWelchTime_real);

	// Collect results in list
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("convergenceInfo") = convergenceInfo,
																				 Rcpp::Named("transProbs") = this->transProbs,
																				 Rcpp::Named("startProbs") = this->startProbs,
																				 Rcpp::Named("weights") = weights,
																				 Rcpp::Named("posteriors") = this->gamma,
																				 Rcpp::Named("states") = imax,
																				 Rcpp::Named("densities") = this->densities);

	// Emission parameters
	if (this->xvariate == UNIVARIATE)
	{
		// Emission densities
		Rcpp::CharacterVector emissionTypes = this->emissionParams["type"];
		if (this->emissionDensities[0]->get_name() == BETA_MIRROR)
		{
			Rcpp::NumericVector a_s = this->emissionParams["a"];
			Rcpp::NumericVector b_s = this->emissionParams["b"];
			for (int irow=0; irow<this->NSTATES; irow++)
			{
				std::string dtype = Rcpp::as<std::string>(emissionTypes[irow]);
		 		if (dtype.compare("dbeta") == 0)
		 		{
		 			a_s[irow] = this->emissionDensities[irow]->get_a();
		 			b_s[irow] = this->emissionDensities[irow]->get_b();
		 		}
			}
		}
		else if ((this->emissionDensities[0]->get_name() == NEGATIVE_BINOMIAL) | (this->emissionDensities[0]->get_name() == ZERO_INFLATION))
		{ 
			Rcpp::NumericVector sizes = this->emissionParams["size"];
			Rcpp::NumericVector probs = this->emissionParams["prob"];
			Rcpp::NumericVector means = this->emissionParams["mu"];
			Rcpp::NumericVector vars = this->emissionParams["var"];
			for (int irow=0; irow<this->NSTATES; irow++)
			{
				std::string dtype = Rcpp::as<std::string>(emissionTypes[irow]);
				if (dtype.compare("dnbinom") == 0)
				{
					sizes[irow] = this->emissionDensities[irow]->get_size();
					probs[irow] = this->emissionDensities[irow]->get_prob();
					means[irow] = this->emissionDensities[irow]->get_mean();
					vars[irow] = this->emissionDensities[irow]->get_variance();
				}
			}
		}
		result.push_back(this->emissionParams, "emissionParams");
	}

	return result;
										 
}


Rcpp::NumericVector ScaleHMM::calc_weights()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	Rcpp::NumericVector weights(this->NSTATES);
	#pragma omp parallel for
	for (int i=0; i<this->NSTATES; i++)
	{
		// Calculate weights by summing posteriors (gamma)
// 		double sum_over_gammas_per_state = 0;
// 		for (int t=0; t<this->NDATA; t++)
// 		{
// 			sum_over_gammas_per_state += this->gamma(i,t);
// 		}
// 		weights[i] = sum_over_gammas_per_state / this->NDATA;

		// Weights by using sumgamma !Be careful if you swap states somewhere!
		weights[i] = ( this->sumgamma[i] + this->gamma(i,this->NDATA-1) ) / this->NDATA;
	}
	return(weights);
}

void ScaleHMM::calc_weights(Rcpp::NumericVector & weights)
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	#pragma omp parallel for
	for (int i=0; i<this->NSTATES; i++)
	{
		// Calculate weights by summing posteriors (gamma)
// 		double sum_over_gammas_per_state = 0;
// 		for (int t=0; t<this->NDATA; t++)
// 		{
// 			sum_over_gammas_per_state += this->gamma(i,t);
// 		}
// 		weights[i] = sum_over_gammas_per_state / this->NDATA;

		// Weights by using sumgamma !Be careful if you swap states somewhere!
		weights[i] = ( this->sumgamma[i] + this->gamma(i,this->NDATA-1) ) / this->NDATA;
	}
}

// Getters and Setters ----------------------------------------
void ScaleHMM::get_posteriors(Rcpp::NumericMatrix & post)
{
	if (this->verbosity>=3) Rprintf("%s\n", __PRETTY_FUNCTION__);
	for (int i=0; i<this->NSTATES; i++)
	{
		for (int t=0; t<this->NDATA; t++)
		{
			post(i,t) = this->gamma(i,t);
		}
	}
}

double ScaleHMM::get_posterior(int i, int t)
{
	if (this->verbosity>=3) Rprintf("%s\n", __PRETTY_FUNCTION__);
	return(this->gamma(i,t));
}

double ScaleHMM::get_density(int i, int t)
{
	if (this->verbosity>=3) Rprintf("%s\n", __PRETTY_FUNCTION__);
	return(this->densities(i,t));
}

double ScaleHMM::get_startProbs(int i)
{
	if (this->verbosity>=3) Rprintf("%s\n", __PRETTY_FUNCTION__);
	return( this->startProbs[i] );
}

double ScaleHMM::get_transProbs(int i, int j)
{
	if (this->verbosity>=3) Rprintf("%s\n", __PRETTY_FUNCTION__);
	return( this->transProbs(i,j) );
}

double ScaleHMM::get_loglik()
{
	if (this->verbosity>=3) Rprintf("%s\n", __PRETTY_FUNCTION__);
	return( this->loglik );
}

// Private ====================================================
// Methods ----------------------------------------------------
void ScaleHMM::forward()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
// 	clock_t time = clock(), dtime;

// 	if (this->xvariate==UNIVARIATE)
// 	{

		std::vector<double> alpha(this->NSTATES);
		// Initialization
		this->scalefactoralpha[0] = 0.0;
		for (int i=0; i<this->NSTATES; i++)
		{
			alpha[i] = this->startProbs[i] * this->densities(i,0);
			this->scalefactoralpha[0] += alpha[i];
		}
		for (int i=0; i<this->NSTATES; i++)
		{
			this->scalealpha(0,i) = alpha[i] / this->scalefactoralpha[0];
		}
		// Induction
		double transProbDistance, transProbHelp;
		for (int t=1; t<this->NDATA; t++)
		{
			transProbHelp = 1.0/this->NSTATES * ( 1.0 - this->transExp[t] );
			this->scalefactoralpha[t] = 0.0;
			for (int i=0; i<this->NSTATES; i++)
			{
				double helpsum = 0.0;
				for (int j=0; j<this->NSTATES; j++)
				{
					if (this->distances[t] > 0)
					{
						transProbDistance = this->transProbs(j,i) * this->transExp[t] + transProbHelp;
					}
					else
					{
						transProbDistance = this->transProbs(j,i);
					}
					helpsum += this->scalealpha(t-1,j) * transProbDistance;
// 					helpsum += this->scalealpha(t-1,j) * this->transProbs(j,i);
				}
				alpha[i] = helpsum * this->densities(i,t);
				this->scalefactoralpha[t] += alpha[i];
			}
			for (int i=0; i<this->NSTATES; i++)
			{
				this->scalealpha(t,i) = alpha[i] / this->scalefactoralpha[t];
				if(std::isnan(this->scalealpha(t,i)))
				{
					if (this->verbosity>=4) Rprintf("scalealpha(t=%d,i=%d) = %g, alpha[i=%d] = %g\n", t, i, scalealpha(t,i), i, alpha[i]);
					if (this->verbosity>=4) Rprintf("scalefactoralpha[t=%d] = %g, scalefactoralpha[t-1=%d] = %g\n", t, scalefactoralpha[t], t-1, scalefactoralpha[t-1]);
					if (this->verbosity>=4) Rprintf("densities(i=%d,t=%d) = %g, startProbs[i=%d] = %g\n", i, t, densities(i,t), i, startProbs[i]);
					for (int j=0; j<this->NSTATES; j++)
					{
						if (this->verbosity>=4) Rprintf("  transProbs(j=%d,i=%d) = %g, transExp[t=%d] = %g, startProbs[j=%d] = %g\n", j, i, transProbs(j,i), t, transExp[t], j, startProbs[j]);
						if (this->verbosity>=4) Rprintf("  densities(j=%d,t=%d) = %g, densities(j=%d,t-1=%d) = %g\n", j, t, densities(j,t), j, t-1, densities(j,t-1));
						if (this->verbosity>=4) Rprintf("  scalealpha(t-1=%d,j=%d) = %g\n", t-1, j, scalealpha(t-1,j));
					}
					throw nan_detected;
				}
			}
		}

// 	}
// 	else if (this->xvariate==MULTIVARIATE)
// 	{
// 
// 		std::vector<double> alpha(this->NSTATES);
// 		// Initialization
// 		this->scalefactoralpha[0] = 0.0;
// 		for (int i=0; i<this->NSTATES; i++)
// 		{
// 			alpha(i] = this->startProbs[i] * this->tdensities[0,i);
// 			this->scalefactoralpha[0] += alpha[i];
// 		}
// 		for (int i=0; i<this->NSTATES; i++)
// 		{
// 			this->scalealpha(0,i) = alpha[i] / this->scalefactoralpha[0];
// 		}
// 		// Induction
// 		for (int t=1; t<this->NDATA; t++)
// 		{
// 			this->scalefactoralpha[t] = 0.0;
// 			for (int i=0; i<this->NSTATES; i++)
// 			{
// 				double helpsum = 0.0;
// 				for (int j=0; j<this->NSTATES; j++)
// 				{
// 					helpsum += this->scalealpha(t-1,j) * this->transProbs(j,i);
// 				}
// 				alpha(i] = helpsum * this->tdensities[t,i);
// 				this->scalefactoralpha[t] += alpha[i];
// 			}
// 			for (int i=0; i<this->NSTATES; i++)
// 			{
// 				this->scalealpha(t,i) = alpha[i] / this->scalefactoralpha[t];
// 				if(std::isnan(this->scalealpha(t,i)))
// 				{
// 					for (int j=0; j<this->NSTATES; j++)
// 					{
// 					}
// 					throw nan_detected;
// 				}
// 			}
// 		}
// 
// 	}

// 	dtime = clock() - time;
}

void ScaleHMM::backward()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
// 	clock_t time = clock(), dtime;

// 	if (this->xvariate==UNIVARIATE)
// 	{

		std::vector<double> beta(this->NSTATES);
		// Initialization
		for (int i=0; i<this->NSTATES; i++)
		{
			beta[i] = 1.0;
		}
		for (int i=0; i<this->NSTATES; i++)
		{
			this->scalebeta(NDATA-1,i) = beta[i] / this->scalefactoralpha[NDATA-1];
		}
		// Induction
		double transProbDistance, transProbHelp;
		for (int t=this->NDATA-2; t>=0; t--)
		{
			transProbHelp = 1.0/this->NSTATES * ( 1.0 - this->transExp[t+1] );
			for (int i=0; i<this->NSTATES; i++)
			{
				beta[i] = 0.0;
				for(int j=0; j<this->NSTATES; j++)
				{
					if (this->distances[t+1] > 0)
					{
						transProbDistance = this->transProbs(i,j) * this->transExp[t+1] + transProbHelp;
					}
					else
					{
						transProbDistance = this->transProbs(i,j);
					}
					beta[i] += transProbDistance * this->densities(j,t+1) * this->scalebeta(t+1,j);
// 					beta[i] += this->transProbs(i,j) * this->densities(j,t+1) * this->scalebeta(t+1,j);
				}
			}
			for (int i=0; i<this->NSTATES; i++)
			{
				this->scalebeta(t,i) = beta[i] / this->scalefactoralpha[t];
				if (std::isnan(this->scalebeta(t,i)))
				{
					for (int j=0; j<this->NSTATES; j++)
					{
					}
					throw nan_detected;
				}
			}
		}

// 	}
// 	else if (this->xvariate==MULTIVARIATE)
// 	{
// 
// 		std::vector<double> beta(this->NSTATES);
// 		// Initialization
// 		for (int i=0; i<this->NSTATES; i++)
// 		{
// 			beta[i] = 1.0;
// 		}
// 		for (int i=0; i<this->NSTATES; i++)
// 		{
// 			this->scalebeta(NDATA-1,i) = beta[i] / this->scalefactoralpha[NDATA-1];
// 		}
// 		// Induction
// 		for (int t=this->NDATA-2; t>=0; t--)
// 		{
// 			for (int i=0; i<this->NSTATES; i++)
// 			{
// 				beta[i] = 0.0;
// 				for(int j=0; j<this->NSTATES; j++)
// 				{
// 					beta(i] += this->transProbs[i,j) * this->tdensities(t+1,j) * this->scalebeta(t+1,j);
// 				}
// 			}
// 			for (int i=0; i<this->NSTATES; i++)
// 			{
// 				this->scalebeta(t,i) = beta[i] / this->scalefactoralpha[t];
// 				if (std::isnan(this->scalebeta(t,i)))
// 				{
// 					for (int j=0; j<this->NSTATES; j++)
// 					{
// 					}
// 					throw nan_detected;
// 				}
// 			}
// 		}
// 
// 	}

// 	dtime = clock() - time;
}

void ScaleHMM::calc_sumgamma()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
// 	clock_t time = clock(), dtime;

	// Initialize the sumgamma
	for (int i=0; i<this->NSTATES; i++)
	{
		this->sumgamma[i] = 0.0;
	}

	// Compute the gammas (posteriors) and sumgamma
	#pragma omp parallel for
	for (int i=0; i<this->NSTATES; i++)
	{
		for (int t=0; t<this->NDATA; t++)
		{
			this->gamma(i,t) = this->scalealpha(t,i) * this->scalebeta(t,i) * this->scalefactoralpha[t];
			this->sumgamma[i] += this->gamma(i,t);
		}
	}
	// Subtract the last value because sumgamma goes only until NDATA-1 and we computed until NDATA to get also loggamma at NDATA
	for (int i=0; i<this->NSTATES; i++)
	{
		this->sumgamma[i] -= this->gamma(i,NDATA-1);
	}

	if (this->verbosity>=6)
	{
		for (int t=0; t<this->NDATA; t++)
		{
			for (int i=0; i<this->NSTATES; i++)
			{
				Rprintf("gamma(i=%d,t=%d) = %g, scalealpha(t=%d,i=%d) = %g, scalebeta(t=%d,i=%d) = %g, scalefactoralpha[t=%d] = %g, densities(i=%d,t=%d) = %g\n", i, t, gamma(i,t), t, i, scalealpha(t,i), t, i, scalebeta(t,i), t, scalefactoralpha[t], i, t, densities(i,t));
			}
		}
	}
// 	dtime = clock() - time;
}

void ScaleHMM::calc_sumxi()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
// 	clock_t time = clock(), dtime;

	double xi;
	// Initialize the sumxi
	for (int i=0; i<this->NSTATES; i++)
	{
		for (int j=0; j<this->NSTATES; j++)
		{
			this->sumxi(i,j) = 0.0;
		}
	}	

// 	if (this->xvariate==UNIVARIATE)
// 	{

		double transProbDistance, transProbHelp;

		for (int t=0; t<this->NDATA-1; t++)
		{
			transProbHelp = 1.0/this->NSTATES * ( 1.0 - this->transExp[t+1] );
// 		#pragma omp parallel for
			for (int i=0; i<this->NSTATES; i++)
			{
				for (int j=0; j<this->NSTATES; j++)
				{
					if (this->distances[t+1] > 0)
					{
						transProbDistance = this->transProbs(i,j) * this->transExp[t+1] + transProbHelp;
					}
					else
					{
						transProbDistance = this->transProbs(i,j);
					}
					xi = this->scalealpha(t,i) * transProbDistance * this->densities(j,t+1) * this->scalebeta(t+1,j);
// 					xi = this->scalealpha(t,i) * this->transProbs(i,j) * this->densities(j,t+1) * this->scalebeta(t+1,j);
					this->sumxi(i,j) += xi;
				}
			}
		}

// 	}
// 	else if (this->xvariate==MULTIVARIATE)
// 	{
// 
// 		#pragma omp parallel for
// 		for (int i=0; i<this->NSTATES; i++)
// 		{
// 			for (int t=0; t<this->NDATA-1; t++)
// 			{
// 				for (int j=0; j<this->NSTATES; j++)
// 				{
// 					xi = this->scalealpha(t,i) * this->transProbs(i,j) * this->tdensities(t+1,j) * this->scalebeta(t+1,j);
// 					this->sumxi(i,j) += xi;
// 				}
// 			}
// 		}
// 
// 	}

// 	dtime = clock() - time;
}

void ScaleHMM::calc_loglikelihood()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
// 	clock_t time = clock(), dtime;

	this->loglik = 0.0;
	for (int t=0; t<this->NDATA; t++)
	{
		this->loglik += log(this->scalefactoralpha[t]);
	}

// 	dtime = clock() - time;
}

void ScaleHMM::calc_densities()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
//	clock_t time = clock(), dtime;
	// Errors thrown inside a #pragma must be handled inside the thread
	std::vector<bool> nan_encountered(this->NSTATES);
	#pragma omp parallel for
	for (int i=0; i<this->NSTATES; i++)
	{
		try
		{
			Rcpp::NumericMatrix::Row row = this->densities(i,Rcpp::_);
			this->emissionDensities[i]->calc_densities(row);
		}
		catch(std::exception& e)
		{
			if (strcmp(e.what(),"nan detected")==0) { nan_encountered[i]=true; }
			else { throw; }
		}
	}
	for (int i=0; i<this->NSTATES; i++)
	{
		if (nan_encountered[i]==true)
		{
			throw nan_detected;
		}
	}

	// Check if the density for all states is close to zero and correct to prevent NaNs
// 	double zero_cutoff = 1.18e-37; // 32-bit precision is 1.18e-38
	double zero_cutoff = 2.23e-307; // 64-bit precision is 2.23e-308
	std::vector<double> temp(this->NSTATES);
	// t=0
	for (int i=0; i<this->NSTATES; i++)
	{
		temp[i] = this->densities(i,0);
	}
	if (*std::max_element(temp.begin(), temp.end()) < zero_cutoff)
	{
		for (int i=0; i<this->NSTATES; i++)
		{
			this->densities(i,0) = zero_cutoff;
		}
	}
	// t>0
	for (int t=1; t<this->NDATA; t++)
	{
		for (int i=0; i<this->NSTATES; i++)
		{
			temp[i] = this->densities(i,t);
		}
		if (*std::max_element(temp.begin(), temp.end()) < zero_cutoff)
		{
			for (int i=0; i<this->NSTATES; i++)
			{
				this->densities(i,t) = this->densities(i,t-1);
			}
		}
	}

// 	dtime = clock() - time;
}

void ScaleHMM::print_uni_iteration(int iteration)
{
	if (this->verbosity>=1)
	{
		this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
		int bs = 106;
		char buffer [106];
		if (iteration % 20 == 0)
		{
			snprintf(buffer, bs, "%10s%20s%20s%15s", "Iteration", "log(P)", "dlog(P)", "Time in sec");
			Rprintf("%s\n", buffer);
		}
		if (iteration == 0)
		{
			snprintf(buffer, bs, "%10s%20s%20s%*d", "0", "-inf", "-", 15, this->baumWelchTime_real);
		}
		else if (iteration == 1)
		{
			snprintf(buffer, bs, "%*d%*f%20s%*d", 10, iteration, 20, this->loglik, "inf", 15, this->baumWelchTime_real);
		}
		else
		{
			snprintf(buffer, bs, "%*d%*f%*f%*d", 10, iteration, 20, this->loglik, 20, this->dloglik, 15, this->baumWelchTime_real);
		}
		Rprintf("%s\n", buffer);

		// Flush Rprintf statements to R console
		R_FlushConsole();
	}
}

void ScaleHMM::print_multi_iteration(int iteration)
{
	if (this->verbosity>=1)
	{
		this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
		int bs = 86;
		char buffer [86];
		if (iteration % 20 == 0)
		{
			snprintf(buffer, bs, "%10s%20s%20s%15s", "Iteration", "log(P)", "dlog(P)", "Time in sec");
			Rprintf("%s\n", buffer);
		}
		if (iteration == 0)
		{
			snprintf(buffer, bs, "%10s%20s%20s%*d", "0", "-inf", "-", 15, this->baumWelchTime_real);
		}
		else if (iteration == 1)
		{
			snprintf(buffer, bs, "%*d%*f%20s%*d", 10, iteration, 20, this->loglik, "inf", 15, this->baumWelchTime_real);
		}
		else
		{
			snprintf(buffer, bs, "%*d%*f%*f%*d", 10, iteration, 20, this->loglik, 20, this->dloglik, 15, this->baumWelchTime_real);
		}
		Rprintf("%s\n", buffer);

		// Flush Rprintf statements to R console
		R_FlushConsole();
	}
}

