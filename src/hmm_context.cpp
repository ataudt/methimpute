#include "hmm_context.h"

// ============================================================
// Hidden Markov Model implemented with scaling strategy
// ============================================================

// Public =====================================================

// Constructor and Destructor ---------------------------------
HMM_context::HMM_context(const Rcpp::IntegerVector & obs_total, const Rcpp::IntegerVector & obs_meth, const Rcpp::IntegerVector & context, const Rcpp::IntegerVector & transitionContext, const Rcpp::NumericVector & distances, Rcpp::NumericVector startProbs_initial, Rcpp::List transProbs_initial, Rcpp::NumericVector transDist, Rcpp::List emissionParams_initial, int min_obs, int verbosity, int update_procedure)
{
	if (verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	this->xvariate = UNIVARIATE;
	this->verbosity = verbosity;
	this->NDATA = obs_total.size();
	this->NSTATES = startProbs_initial.size();
	this->UPDATE_PROCEDURE = update_procedure;
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
	this->transDist = transDist;
	this->transitionContext = transitionContext;

	// Initial probabilities
	this->transProbsList = Rcpp::clone(transProbs_initial);
	this->transExp = Rcpp::NumericVector(this->NDATA);
	for (int t=1; t<this->NDATA; t++)
	{
		if (this->distances[t] == INFINITY)
		{
			this->transExp[t] = 0;
		}
		else
		{
			this->transExp[t] = exp(-this->distances[t] / transDist[transitionContext[t]]);
		}
		if (std::isnan(this->transExp[t]))
		{
			if (this->verbosity>=4) Rprintf("transExp[t=%d] = %g, distances[t] = %g, transitionContext[t] = %d, transDist[%d] = %g\n", t, transExp[t], distances[t], transitionContext[t], transitionContext[t], transDist[transitionContext[t]]);
			throw nan_detected;
		}
	}
	this->startProbs = Rcpp::clone(startProbs_initial);

	// Emission densities
	this->emissionParamsList = Rcpp::clone(emissionParams_initial);
	for (int i=0; i<this->NSTATES; i++)
	{
		Rcpp::DataFrame probsDF = Rcpp::as<Rcpp::DataFrame>(this->emissionParamsList[i]);
		Rcpp::NumericVector probs = probsDF["prob"];
		// Binomial test
		BinomialTestContext * d = new BinomialTestContext(obs_total, obs_meth, context, probs, min_obs, this->verbosity);
		this->emissionDensities.push_back(d);
	}

}

HMM_context::~HMM_context()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	for (int i=0; i<this->emissionDensities.size(); i++)
	{
		delete this->emissionDensities[i];
	}
}

// Methods ----------------------------------------------------
Rcpp::List HMM_context::forward_backward(double eps, double maxiter, double maxtime)
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

// 	this->calc_sumxi();
// 	R_CheckUserInterrupt();

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
																									Rcpp::Named("dloglik") = this->dloglik,
																					 				Rcpp::Named("time") = this->baumWelchTime_real);

	// Collect results in list
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("convergenceInfo") = convergenceInfo,
																				 Rcpp::Named("transProbs") = this->transProbsList,
																				 Rcpp::Named("transDist") = this->transDist,
																				 Rcpp::Named("startProbs") = this->startProbs,
																				 Rcpp::Named("weights") = weights,
																				 Rcpp::Named("posteriors") = this->gamma,
																				 Rcpp::Named("states") = imax,
																				 Rcpp::Named("densities") = this->densities);

	// Emission parameters
	if (this->xvariate == UNIVARIATE)
	{
		// Emission densities
		if (this->emissionDensities[0]->get_name() == BETA_MIRROR)
		{
			result.push_back(this->emissionParams, "emissionParams");
		}
		else if ((this->emissionDensities[0]->get_name() == NEGATIVE_BINOMIAL) | (this->emissionDensities[0]->get_name() == ZERO_INFLATION))
		{ 
			result.push_back(this->emissionParams, "emissionParams");
		}
		else if (this->emissionDensities[0]->get_name() == BINOMIAL_TEST)
		{ 
			result.push_back(this->emissionParams, "emissionParams");
		}
		else if (this->emissionDensities[0]->get_name() == BINOMIAL_TEST_CONTEXT)
		{ 
			result.push_back(this->emissionParamsList, "emissionParams");
		}
	}

	return result;
										 
}


Rcpp::List HMM_context::baumWelch(double eps, double maxiter, double maxtime)
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

	// Store information about parameters
	Rcpp::NumericVector probs0, probs1;
	if (this->xvariate == UNIVARIATE)
	{
		if (this->NSTATES == 2)
		{
			Rcpp::NumericVector probs0Current = this->emissionDensities[0]->get_probs();
			Rcpp::NumericVector probs1Current = this->emissionDensities[1]->get_probs();
			for (int c=0; c<probs0Current.size(); c++)
			{
				probs0.push_back(probs0Current[c]);
				probs1.push_back(probs1Current[c]);
			}
		}
		else if (this->NSTATES == 3)
		{
			Rcpp::NumericVector probs0Current = this->emissionDensities[0]->get_probs();
			Rcpp::NumericVector probs1Current = this->emissionDensities[2]->get_probs();
			for (int c=0; c<probs0Current.size(); c++)
			{
				probs0.push_back(probs0Current[c]);
				probs1.push_back(probs1Current[c]);
			}
		}
	}
					


	R_CheckUserInterrupt();
	// Do the Baum-Welch
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	int iteration = 0;
	Rcpp::NumericVector logliks, dlogliks;
	logliks.push_back(-INFINITY);
	dlogliks.push_back(INFINITY);
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
		dlogliks.push_back(this->dloglik);
		// Store information about parameters
		if (this->xvariate == UNIVARIATE)
		{
			if (this->NSTATES == 2)
			{
				Rcpp::NumericVector probs0Current = this->emissionDensities[0]->get_probs();
				Rcpp::NumericVector probs1Current = this->emissionDensities[1]->get_probs();
				for (int c=0; c<probs0Current.size(); c++)
				{
					probs0.push_back(probs0Current[c]);
					probs1.push_back(probs1Current[c]);
				}
			}
			else if (this->NSTATES == 3)
			{
				Rcpp::NumericVector probs0Current = this->emissionDensities[0]->get_probs();
				Rcpp::NumericVector probs1Current = this->emissionDensities[2]->get_probs();
				for (int c=0; c<probs0Current.size(); c++)
				{
					probs0.push_back(probs0Current[c]);
					probs1.push_back(probs1Current[c]);
				}
			}
		}


// 		this->calc_sumxi();
// 		R_CheckUserInterrupt();

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
		if(fabs(this->dloglik) < eps) //it has converged
		{
			if (this->verbosity>=1) Rprintf("HMM: Convergence reached!\n");
			break;
		} else {// not converged
			this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
			if (iteration == maxiter)
			{
				if (this->verbosity>=1) Rprintf("HMM: Maximum number of iterations reached!\n");
				break;
			}
			else if ((this->baumWelchTime_real >= maxtime) and (maxtime >= 0))
			{
				if (this->verbosity>=1) Rprintf("HMM: Exceeded maximum time!\n");
				break;
			}
		}
		
		// Updating initial probabilities startProbs and transition matrix transProbs
		this->update_startProbs();
		this->update_transProbs();
// 		this->update_transDist();

		// Update the parameters of the distribution
		if (this->verbosity>=2) Rprintf("  updating distribution parameters\n");
		if (this->emissionDensities[0]->get_name() == BINOMIAL_TEST)
		{ 
			if (this->NSTATES == 2)
			{
				for (int i=0; i<this->NSTATES; i++)
				{
					int rows[] = {i};
					this->emissionDensities[i]->update(this->gamma, rows);
					if (this->verbosity>=4) Rprintf("  emissionDensities[%d]: prob = %g\n", i, emissionDensities[i]->get_prob());
				}
			} else if (this->NSTATES == 3) { // epi-heterozygosity
					double r = this->emissionDensities[0]->get_prob();
					double p = this->emissionDensities[2]->get_prob();
					const int rows1[] = {0,1};
					this->emissionDensities[0]->update_constrained(this->gamma, rows1, p);
					if (this->verbosity>=4) Rprintf("  emissionDensities[%d]: prob = %g\n", 0, emissionDensities[0]->get_prob());
					const int rows2[] = {2,1};
					this->emissionDensities[2]->update_constrained(this->gamma, rows2, r);
					if (this->verbosity>=4) Rprintf("  emissionDensities[%d]: prob = %g\n", 2, emissionDensities[2]->get_prob());
					this->emissionDensities[1]->set_prob(0.5*(this->emissionDensities[0]->get_prob() + this->emissionDensities[2]->get_prob()));
					if (this->verbosity>=4) Rprintf("  emissionDensities[%d]: prob = %g\n", 1, emissionDensities[1]->get_prob());
			}
		}
		else if (this->emissionDensities[0]->get_name() == BINOMIAL_TEST_CONTEXT)
		{ 
			if (this->NSTATES == 2 | this->UPDATE_PROCEDURE == 1)
			{
				for (int i=0; i<this->NSTATES; i++)
				{
					int rows[] = {i};
					this->emissionDensities[i]->update(this->gamma, rows);
				}
			} else if (this->NSTATES == 3) { // epi-heterozygosity
					Rcpp::NumericVector rs = this->emissionDensities[0]->get_probs();
					Rcpp::NumericVector ps = this->emissionDensities[2]->get_probs();
					const int rows1[] = {0,1};
					this->emissionDensities[0]->update_constrained_context(this->gamma, rows1, ps);
					const int rows2[] = {2,1};
					this->emissionDensities[2]->update_constrained_context(this->gamma, rows2, rs);
					// Update heterozyous state
					rs = this->emissionDensities[0]->get_probs();
					ps = this->emissionDensities[2]->get_probs();
					Rcpp::NumericVector qs = Rcpp::NumericVector(rs.size());
					for (int c=0; c<rs.size(); c++)
					{
						qs[c] = 0.5*(rs[c] + ps[c]);
					}
					this->emissionDensities[1]->set_probs(qs);
			}
		}
		R_CheckUserInterrupt();

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

	// Store information about parameters
	Rcpp::List parameterInfo = Rcpp::List::create(Rcpp::Named("probsUN") = probs0,
																								Rcpp::Named("probsM") = probs1);

	// Convergence information
	this->baumWelchTime_real = difftime(time(NULL),this->baumWelchStartTime_sec);
	Rcpp::List convergenceInfo = Rcpp::List::create(Rcpp::Named("logliks") = logliks,
																									Rcpp::Named("dlogliks") = dlogliks,
																									Rcpp::Named("parameterInfo") = parameterInfo,
																									Rcpp::Named("time") = this->baumWelchTime_real);

	// Collect results in list
	Rcpp::List result = Rcpp::List::create(Rcpp::Named("convergenceInfo") = convergenceInfo,
																				 Rcpp::Named("transProbs") = this->transProbsList,
																				 Rcpp::Named("transDist") = this->transDist,
																				 Rcpp::Named("startProbs") = this->startProbs,
																				 Rcpp::Named("weights") = weights,
																				 Rcpp::Named("posteriors") = this->gamma,
																				 Rcpp::Named("states") = imax,
																				 Rcpp::Named("densities") = this->densities);

	// Emission parameters
	if (this->xvariate == UNIVARIATE)
	{
		// Emission densities
		if (this->emissionDensities[0]->get_name() == BETA_MIRROR)
		{
			Rcpp::CharacterVector emissionTypes = this->emissionParams["type"];
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
			result.push_back(this->emissionParams, "emissionParams");
		}
		else if ((this->emissionDensities[0]->get_name() == NEGATIVE_BINOMIAL) | (this->emissionDensities[0]->get_name() == ZERO_INFLATION))
		{ 
			Rcpp::CharacterVector emissionTypes = this->emissionParams["type"];
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
			result.push_back(this->emissionParams, "emissionParams");
		}
		else if (this->emissionDensities[0]->get_name() == BINOMIAL_TEST)
		{ 
			Rcpp::CharacterVector emissionTypes = this->emissionParams["type"];
			Rcpp::NumericVector probs = this->emissionParams["prob"];
			for (int irow=0; irow<this->NSTATES; irow++)
			{
				std::string dtype = Rcpp::as<std::string>(emissionTypes[irow]);
				if (dtype.compare("dbinom") == 0)
				{
					probs[irow] = this->emissionDensities[irow]->get_prob();
				}
			}
			result.push_back(this->emissionParams, "emissionParams");
		}
		else if (this->emissionDensities[0]->get_name() == BINOMIAL_TEST_CONTEXT)
		{ 
// 			for (int i=0; i<this->NSTATES; i++)
// 			{
// 				Rcpp::DataFrame probsDF = Rcpp::as<Rcpp::DataFrame>(this->emissionParamsList[i]);
// 				Rcpp::NumericVector probs = probsDF["prob"];
// 				Rcpp::NumericVector ps = this->emissionDensities[i]->get_probs();
// 				for (int c=0; c<ps.size(); c++)
// 				{
// 					probs[c] = ps[c];
// 				}
// 			}
			result.push_back(this->emissionParamsList, "emissionParams");
		}
	}

	return result;
										 
}


Rcpp::NumericVector HMM_context::calc_weights()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	Rcpp::NumericVector weights(this->NSTATES);
	#pragma omp parallel for
	for (int i=0; i<this->NSTATES; i++)
	{
		// Calculate weights by summing posteriors (gamma)
		double sum_over_gammas_per_state = 0;
		for (int t=0; t<this->NDATA; t++)
		{
			sum_over_gammas_per_state += this->gamma(i,t);
		}
		weights[i] = sum_over_gammas_per_state / this->NDATA;

// 		// Weights by using sumgamma !Be careful if you swap states somewhere!
// 		weights[i] = ( this->sumgamma[i] + this->gamma(i,this->NDATA-1) ) / this->NDATA;
	}
	return(weights);
}

void HMM_context::calc_weights(Rcpp::NumericVector & weights)
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	#pragma omp parallel for
	for (int i=0; i<this->NSTATES; i++)
	{
		// Calculate weights by summing posteriors (gamma)
		double sum_over_gammas_per_state = 0;
		for (int t=0; t<this->NDATA; t++)
		{
			sum_over_gammas_per_state += this->gamma(i,t);
		}
		weights[i] = sum_over_gammas_per_state / this->NDATA;

// 		// Weights by using sumgamma !Be careful if you swap states somewhere!
// 		weights[i] = ( this->sumgamma[i] + this->gamma(i,this->NDATA-1) ) / this->NDATA;
	}
}

// Getters and Setters ----------------------------------------
void HMM_context::get_posteriors(Rcpp::NumericMatrix & post)
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

double HMM_context::get_posterior(int i, int t)
{
	if (this->verbosity>=3) Rprintf("%s\n", __PRETTY_FUNCTION__);
	return(this->gamma(i,t));
}

double HMM_context::get_density(int i, int t)
{
	if (this->verbosity>=3) Rprintf("%s\n", __PRETTY_FUNCTION__);
	return(this->densities(i,t));
}

double HMM_context::get_startProbs(int i)
{
	if (this->verbosity>=3) Rprintf("%s\n", __PRETTY_FUNCTION__);
	return( this->startProbs[i] );
}

double HMM_context::get_loglik()
{
	if (this->verbosity>=3) Rprintf("%s\n", __PRETTY_FUNCTION__);
	return( this->loglik );
}

// Private ====================================================
// Methods ----------------------------------------------------
void HMM_context::forward()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	std::vector<double> alpha(this->NSTATES);
	Rcpp::NumericMatrix transProbs;
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
		transProbs = Rcpp::as<Rcpp::NumericMatrix>(this->transProbsList[this->transitionContext[t]]);
		transProbHelp = 1.0/this->NSTATES * ( 1.0 - this->transExp[t] );
		this->scalefactoralpha[t] = 0.0;
		for (int i=0; i<this->NSTATES; i++)
		{
			double helpsum = 0.0;
			for (int j=0; j<this->NSTATES; j++)
			{
				if (this->distances[t] > 0)
				{
					transProbDistance = transProbs(j,i) * this->transExp[t] + transProbHelp;
				}
				else
				{
					transProbDistance = transProbs(j,i);
				}
				helpsum += this->scalealpha(t-1,j) * transProbDistance;
// 					helpsum += this->scalealpha(t-1,j) * transProbs(j,i);
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

}

void HMM_context::backward()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	std::vector<double> beta(this->NSTATES);
	Rcpp::NumericMatrix transProbs;
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
		transProbs = Rcpp::as<Rcpp::NumericMatrix>(this->transProbsList[this->transitionContext[t+1]]);
		transProbHelp = 1.0/this->NSTATES * ( 1.0 - this->transExp[t+1] );
		for (int i=0; i<this->NSTATES; i++)
		{
			beta[i] = 0.0;
			for(int j=0; j<this->NSTATES; j++)
			{
				if (this->distances[t+1] > 0)
				{
					transProbDistance = transProbs(i,j) * this->transExp[t+1] + transProbHelp;
				}
				else
				{
					transProbDistance = transProbs(i,j);
				}
				beta[i] += transProbDistance * this->densities(j,t+1) * this->scalebeta(t+1,j);
// 					beta[i] += transProbs(i,j) * this->densities(j,t+1) * this->scalebeta(t+1,j);
			}
		}
		for (int i=0; i<this->NSTATES; i++)
		{
			this->scalebeta(t,i) = beta[i] / this->scalefactoralpha[t];
			if (std::isnan(this->scalebeta(t,i)))
			{
				for (int j=0; j<this->NSTATES; j++)
				{
// 						Rprintf("scalebeta(t+1=%d,j=%d) = %f, densities(j=%d,t+1=%d) = %f\n", t+1, j, scalebeta(t+1,j), j, t+1, densities(j,t+1));
				}
// 					Rprintf("scalebeta(t=%d,i=%d) = %f, beta[i=%d] = %f, scalefactoralpha[t=%d] = %f\n", t, i, scalebeta(t,i), i, beta[i], t, scalefactoralpha[t]);
				throw nan_detected;
			}
		}
	}

}

void HMM_context::calc_sumgamma()
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
		}
		for (int t=0; t<this->NDATA-1; t++)
		{
			this->sumgamma[i] += this->gamma(i,t);
		}
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

void HMM_context::calc_sumxi()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	double xi;
	Rcpp::NumericMatrix transProbs;
	// Initialize the sumxi
	for (int i=0; i<this->NSTATES; i++)
	{
		for (int j=0; j<this->NSTATES; j++)
		{
			this->sumxi(i,j) = 0.0;
		}
	}	

	double transProbDistance, transProbHelp;

	for (int t=0; t<this->NDATA-1; t++)
	{
		transProbs = Rcpp::as<Rcpp::NumericMatrix>(this->transProbsList[this->transitionContext[t+1]]);
		transProbHelp = 1.0/this->NSTATES * ( 1.0 - this->transExp[t+1] );
// 		#pragma omp parallel for
		for (int i=0; i<this->NSTATES; i++)
		{
			for (int j=0; j<this->NSTATES; j++)
			{
				if (this->distances[t+1] > 0)
				{
					transProbDistance = transProbs(i,j) * this->transExp[t+1] + transProbHelp;
				}
				else
				{
					transProbDistance = transProbs(i,j);
				}
				xi = this->scalealpha(t,i) * transProbDistance * this->densities(j,t+1) * this->scalebeta(t+1,j);
// 					xi = this->scalealpha(t,i) * transProbs(i,j) * this->densities(j,t+1) * this->scalebeta(t+1,j);
				this->sumxi(i,j) += xi;
			}
		}
	}

}

void HMM_context::calc_loglikelihood()
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

void HMM_context::calc_densities()
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
// 	double zero_cutoff = 1.18e-27; // 32-bit precision is 1.18e-38
	double zero_cutoff = 2.23e-207; // 64-bit precision is 2.23e-308
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
// 			Rprintf("corrected at t=%d\n", t);
			for (int i=0; i<this->NSTATES; i++)
			{
				this->densities(i,t) = this->densities(i,t-1);
			}
		}
	}

// 	dtime = clock() - time;
}

void HMM_context::update_startProbs()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	for (int i=0; i<this->NSTATES; i++)
	{
		this->startProbs[i] = this->gamma(i,0);
		if (this->verbosity>=4) Rprintf("  startProbs[%d] = %g\n", i, startProbs[i]);
	}
}

void HMM_context::update_transProbs()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);

	Rcpp::NumericMatrix transProbs;
	Rcpp::NumericMatrix transProbsCurrent;

	for (int tc=0; tc<this->transProbsList.size(); tc++)
	{
		transProbs = Rcpp::as<Rcpp::NumericMatrix>(this->transProbsList[tc]);
		transProbsCurrent = Rcpp::clone(transProbs);

		#pragma omp parallel for
		for (int i=0; i<this->NSTATES; i++)
		{
			double dist_f; // Correction factor due to distance dependency
			double xir; // By transProbs reduced xi
			double denominator = 0.0;
			std::vector<double> numerators(this->NSTATES);
			for (int j=0; j<this->NSTATES; j++)
			{
				numerators[j] = 0.0;
				for (int t=0; t<this->NDATA-1; t++)
				{
					if (this->transitionContext[t+1] == tc)
					{
						// Calculate xir
						dist_f = this->transExp[t+1] * transProbsCurrent(i,j);
						xir = this->scalealpha(t,i) * this->densities(j,t+1) * this->scalebeta(t+1,j);

						// Numerator
						numerators[j] += xir * dist_f;
					}
				}
			}
			
			// Denominator
			denominator = 0.0;
			for (int j=0; j<this->NSTATES; j++)
			{
				denominator += numerators[j];
			}

			// Update
			for (int j=0; j<this->NSTATES; j++)
			{
				transProbs(i,j) = numerators[j] / denominator;
				// Check for nan
				if (std::isnan(transProbs(i,j)))
				{
					throw nan_detected;
				}
			}

		}
	}
}

void HMM_context::update_transDist()
{
	if (this->verbosity>=2) Rprintf("%s\n", __PRETTY_FUNCTION__);
	double xir;
	double eps = 1e-4;
	double kmax = 1000;
	double F, dFdX, FdivM; // F = dL/dX
	double A, A1, A2;
	double Ahelp;
	double x;
	Rcpp::NumericMatrix transProbs;
	// Update with Newton Method
	for (int tc=0; tc<this->transProbsList.size(); tc++)
	{
		transProbs = Rcpp::as<Rcpp::NumericMatrix>(this->transProbsList[tc]);
		x = this->transDist[tc];

		for (int k=0; k<kmax; k++)
		{
			F = dFdX = 0.0;
			#pragma omp parallel for
			for (int i=0; i<this->NSTATES; i++)
			{
				for (int j=0; j<this->NSTATES; j++)
				{
					for(int t=0; t<this->NDATA-1; t++)
					{
						if (this->transitionContext[t+1] == tc)
						{
							xir = this->scalealpha(t,i) * this->densities(j,t+1) * this->scalebeta(t+1,j);
							Ahelp = (transProbs(i,j) - 1.0/this->NSTATES) * exp(-this->distances[t+1]/x);
							A = Ahelp + 1.0/this->NSTATES;
							A1 = Ahelp * this->distances[t+1] * pow(x,-2);
							A2 = A1 * ( this->distances[t+1] * pow(x,-2) - 2 / x );
							F += xir * A1;
							dFdX += xir * ( A2 - A1*A1 / A);
						}
					}
				}
			}
			FdivM = F/dFdX;
			if (FdivM < x)
			{
				x = x-FdivM;
			}
			else if (FdivM >= x)
			{
				x = x/2.0;
			}
			if(fabs(F)<eps)
			{
				break;
			}
			if (k == kmax-1)
			{
				Rprintf("WARNING: %s: Increase kmax to achieve convergence.\n", __PRETTY_FUNCTION__);
			}
		}
		this->transDist[tc] = x;

	}

	// Update transExp
	for (int t=1; t<this->NDATA; t++)
	{
		if (this->distances[t] == INFINITY)
		{
			this->transExp[t] = 0;
		}
		else
		{
			this->transExp[t] = exp(-this->distances[t] / transDist[transitionContext[t]]);
		}
		if (std::isnan(this->transExp[t]))
		{
			if (this->verbosity>=4) Rprintf("transExp[t=%d] = %g, distances[t] = %g, transitionContext[t] = %d, transDist[%d] = %g\n", t, transExp[t], distances[t], transitionContext[t], transitionContext[t], transDist[transitionContext[t]]);
			throw nan_detected;
		}
	}

}

void HMM_context::print_uni_iteration(int iteration)
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

void HMM_context::print_multi_iteration(int iteration)
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

