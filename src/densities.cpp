#include "densities.h"

// ============================================================
// Zero-inflated Negative Binomial
// ============================================================

// Constructor and Destructor ---------------------------------
ZiNB::ZiNB(const Rcpp::IntegerVector & obs, double size, double prob, double w, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs = obs;
	this->prob = prob;
	this->size = size;
	this->w = w;
	this->lxfactorials = NULL;
	this->max_obs = Rcpp::max(this->obs);
	this->lxfactorials = Rcpp::NumericVector(max_obs+1);
	this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of Calloc
	this->lxfactorials[1] = 0.0;
	for (int j=2; j<=max_obs; j++)
	{
		this->lxfactorials[j] = this->lxfactorials[j-1] + log(j);
	}

	// Make vector of positions of unique observations for faster updating
	this->obs_unique = Rcpp::unique(obs);
	this->obs_unique.sort();
	Rcpp::IntegerVector uobsind_per_obs = Rcpp::IntegerVector(this->obs_unique[this->obs_unique.size()-1]+1);
	int i = 0;
	for (int j=0; j<uobsind_per_obs.size(); j++)
	{
		if (this->obs_unique[i] = j)
		{
			uobsind_per_obs[j] = i;
			i += 1;
		}
	}
	// Get the index for each observation
	this->uobsind_per_t = Rcpp::IntegerVector(this->obs.size());
	for (int t=0; t<this->obs.size(); t++)
	{
		this->uobsind_per_t[t] = uobsind_per_obs[this->obs[t]];
	}

}

ZiNB::ZiNB(const Rcpp::IntegerVector & obs, const Rcpp::IntegerVector & obs_unique, const Rcpp::IntegerVector & uobsind_per_t, double size, double prob, double w, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs = obs;
	this->obs_unique = obs_unique;
	this->uobsind_per_t = uobsind_per_t;
	this->prob = prob;
	this->size = size;
	this->w = w;
	this->lxfactorials = NULL;
	this->max_obs = Rcpp::max(this->obs);
	this->lxfactorials = Rcpp::NumericVector(max_obs+1);
	this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of Calloc
	this->lxfactorials[1] = 0.0;
	for (int j=2; j<=max_obs; j++)
	{
		this->lxfactorials[j] = this->lxfactorials[j-1] + log(j);
	}

}

ZiNB::~ZiNB()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
}

// Methods ----------------------------------------------------
void ZiNB::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	double obs_j;
	lGammaR = lgamma(this->size);
	// Select strategy for computing gammas, redundant since obs_unique.size() always < obs.size()
	if (this->obs_unique.size() <= this->obs.size())
	{
		std::vector<double> logdens_per_uobs(this->obs_unique.size());
		for (int j=0; j<=this->obs_unique.size(); j++)
		{
			obs_j = this->obs_unique[j];
			if (obs_j == 0)
			{
				logdens_per_uobs[j] = log( this->w + (1-this->w) * exp( lgamma(this->size + 0) - lGammaR - lxfactorials[0] + this->size * logp + 0 * log1minusp ) );
			}
			else
			{
				logdens_per_uobs[j] = log(1-this->w) + lgamma(this->size + obs_j) - lGammaR - lxfactorials[obs_j] + this->size * logp + obs_j * log1minusp;
			}
		}
		for (int t=0; t<this->obs.size(); t++)
		{
			logdens[t] = logdens_per_uobs[this->uobsind_per_t[t]];
			if (std::isnan(logdens[t]))
			{
				throw nan_detected;
			}
		}
	}
	else
	{
		for (int t=0; t<this->obs.size(); t++)
		{
			lGammaRplusX = lgamma(this->size + this->obs[t]);
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			if (obs[t] == 0)
			{
				logdens[t] = log( this->w + (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp ) );
			}
			else
			{
				logdens[t] = log(1-this->w) + lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp;
			}
			if (std::isnan(logdens[t]))
			{
				throw nan_detected;
			}
		}
	}
}

void ZiNB::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	double obs_j;
	lGammaR = lgamma(this->size);
	// Select strategy for computing gammas, redundant since obs_unique.size() always < obs.size()
	if (this->obs_unique.size() <= this->obs.size())
	{
		std::vector<double> dens_per_uobs(this->obs_unique.size());
		for (int j=0; j<=this->obs_unique.size(); j++)
		{
			obs_j = this->obs_unique[j];
			if (obs_j == 0)
			{
				dens_per_uobs[j] = this->w + (1-this->w) * exp( lgamma(this->size + 0) - lGammaR - lxfactorials[0] + this->size * logp + 0 * log1minusp );
			}
			else
			{
				dens_per_uobs[j] = (1-this->w) * exp( lgamma(this->size + obs_j) - lGammaR - lxfactorials[obs_j] + this->size * logp + obs_j * log1minusp );
			}
		}
		for (int t=0; t<this->obs.size(); t++)
		{
			dens[t] = dens_per_uobs[this->uobsind_per_t[t]];
			if (std::isnan(dens[t]))
			{
				throw nan_detected;
			}
		}
	}
	else
	{
		for (int t=0; t<this->obs.size(); t++)
		{
			lGammaRplusX = lgamma(this->size + this->obs[t]);
			lxfactorial = this->lxfactorials[(int) this->obs[t]];
			if (obs[t] == 0)
			{
				dens[t] = this->w + (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
			}
			else
			{
				dens[t] = (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
			}
			if (std::isnan(dens[t]))
			{
				throw nan_detected;
			}
		}
	}
}

void ZiNB::calc_CDFs(Rcpp::NumericMatrix::Row & CDF)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR;
	lGammaR=lgamma(this->size);
	std::vector<double> precomputed_CDF(this->max_obs+1);
	double dens;

	// Calculate for j=0
	precomputed_CDF[0] = this->w + (1-this->w) * exp( lgamma(this->size) - lGammaR - this->lxfactorials[0] + this->size * logp );
	// Calculate for j>0
	for (int j=1; j<=this->max_obs; j++)
	{
		dens = (1-this->w) * exp( lgamma(this->size + j) - lGammaR - this->lxfactorials[j] + this->size * logp + j * log1minusp );
		if (std::isnan(dens))
		{
			throw nan_detected;
		}
		precomputed_CDF[j] = precomputed_CDF[j-1] + dens;
		if (precomputed_CDF[j] >= 1)
		{
			precomputed_CDF[j] = precomputed_CDF[j-1]; 
		}
	}
	for (int t=0; t<this->obs.size(); t++)
	{
		CDF[t] = precomputed_CDF[(int)obs[t]];
		if (std::isnan(CDF[t]))
		{
			throw nan_detected;
		}
	}
}

void ZiNB::calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR;
	lGammaR=lgamma(this->size);
	std::vector<double> precomputed_logCDF(this->max_obs+1);
	double logdens;

	// Calculate for j=0
	precomputed_logCDF[0] = log( this->w + (1-this->w) * exp( lgamma(this->size) - lGammaR - this->lxfactorials[0] + this->size * logp ) );
	// Calculate for j>0
	for (int j=1; j<=this->max_obs; j++)
	{
		logdens = log(1-this->w) + lgamma(this->size + j) - lGammaR - this->lxfactorials[j] + this->size * logp + j * log1minusp;
		if (std::isnan(logdens))
		{
			throw nan_detected;
		}
		precomputed_logCDF[j] = log( exp(precomputed_logCDF[j-1]) + exp(logdens) );
		if (precomputed_logCDF[j] >= 0)
		{
			precomputed_logCDF[j] = precomputed_logCDF[j-1]; 
		}
	}
	for (int t=0; t<this->obs.size(); t++)
	{
		logCDF[t] = precomputed_logCDF[(int)obs[t]];
		if (std::isnan(logCDF[t]))
		{
			throw nan_detected;
		}
	}
}

double ZiNB::getLogDensityAt(int x)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	double logdens;
	// Calculate variance
	double mean = 0, variance = 0;
	for(int t=0; t<this->obs.size(); t++)
	{
		mean += obs[t];
	}
	mean = mean / this->obs.size();
	for(int t=0; t<this->obs.size(); t++)
	{
		variance += pow(obs[t] - mean, 2);
	}
	variance = variance / this->obs.size();
	// Calculate logdensity
	lGammaR = lgamma(this->size);
	lGammaRplusX = lgamma(this->size + x);
	lxfactorial = this->lxfactorials[x];
	if (x == 0)
	{
		logdens = log( this->w + (1-this->w) * exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + x * log1minusp ) );
	}
	else
	{
		logdens = log(1-this->w) + lGammaRplusX - lGammaR - lxfactorial + this->size * logp + x * log1minusp;
	}
	if (std::isnan(logdens))
	{
		throw nan_detected;
	}
	
	return(logdens);
}

// Getter and Setter ------------------------------------------
double ZiNB::get_mean()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return (1-this->w)*this->size*(1-this->prob)/this->prob;
}

double ZiNB::get_variance()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return (1-this->w)*this->size*(1-this->prob)/this->prob/this->prob; //TODO: Is this correct?
}

DensityName ZiNB::get_name()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(ZERO_INFLATED_NEGATIVE_BINOMIAL);
}

double ZiNB::get_size()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
  return(this->size);
}

double ZiNB::get_prob()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
  return(this->prob);
}

double ZiNB::get_w()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(this->w);
}


// ============================================================
// Binomial test
// ============================================================

// Constructor and Destructor ---------------------------------
BinomialTest::BinomialTest() { }

BinomialTest::BinomialTest(const Rcpp::IntegerVector & obs_total, const Rcpp::IntegerVector & obs_test, double prob, int min_obs, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs_total = obs_total;
	this->obs_test = obs_test;
	this->prob = prob;
	this->min_obs = min_obs;

}

BinomialTest::~BinomialTest()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
}

// Methods ----------------------------------------------------
void BinomialTest::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logdens_min = log(1.0/this->min_obs);
	double logprob = log(this->prob);
	for (int t=0; t<this->obs_total.size(); t++)
	{
		if (this->obs_total[t] < this->min_obs)
		{
			logdens[t] = logdens_min;
		}
		else
		{
			logdens[t] = R::dbinom(this->obs_test[t], this->obs_total[t], this->prob, true);
// 			logdens[t] = this->lxfactorials[this->obs_total[t]] - this->lxfactorials[this->obs_test[t]] - this->lxfactorials[this->obs_total[t] - this->obs_test[t]] + this->obs_test[t]*logprob + (this->obs_total[t] - this->obs_test[t]) * logprob; // precision problems
		}
		if (std::isnan(logdens[t]))
		{
			throw nan_detected;
		}
	}
} 

void BinomialTest::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logprob = log(this->prob);
	double dens_min = 1.0/this->min_obs;
	for (int t=0; t<this->obs_total.size(); t++)
	{
		if (this->obs_total[t] < this->min_obs)
		{
			dens[t] = dens_min;
		}
		else
		{
			dens[t] = R::dbinom(this->obs_test[t], this->obs_total[t], this->prob, false);
// 			dens[t] = exp( this->lxfactorials[this->obs_total[t]] - this->lxfactorials[this->obs_test[t]] - this->lxfactorials[this->obs_total[t] - this->obs_test[t]] + this->obs_test[t]*logprob + (this->obs_total[t] - this->obs_test[t]) * logprob ); // precision problems
		}
		if (std::isnan(dens[t]))
		{
			throw nan_detected;
		}
	}
} 

void BinomialTest::update(const Rcpp::NumericMatrix & weights, const int * rows)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	// Update prob (p)
	double numerator, denominator;
	numerator = denominator = 0.0;
	for (int t=0; t<this->obs_total.size(); t++)
	{
		if (this->obs_total[t] >= this->min_obs)
		{
			numerator += weights(rows[0],t) * (this->obs_test[t]);
			denominator += weights(rows[0],t) * (this->obs_total[t]);
		}
	}
	this->prob = numerator/denominator; // Update this->prob

// 	double eps = 1e-4;
// 	double kmax = 20;
// 	double F, dFdProb, FdivM; // F = dL/dProb
// 	double n, m;
// 	double p = this->prob;
// 	// Update of prob with Newton Method
// 	for (int k=0; k<kmax; k++)
// 	{
// 		F = dFdProb = 0.0;
// 		for(int t=0; t<this->obs_total.size(); t++)
// 		{
// 			if (this->obs_total[t] >= this->min_obs)
// 			{
// 				m = (double)this->obs_test[t];
// 				n = (double)this->obs_total[t];
// 				F += weights(rows[0],t) * (m/p + (m-n)/(1-p));
// 				dFdProb += weights(rows[0],t) * (-m/p/p + (m-n)/(1-p)/(1-p));
// 			}
// 		}
// 		FdivM = F/dFdProb;
// 		if (FdivM < p)
// 		{
// 			p = p-FdivM;
// 		}
// 		else if (FdivM >= p)
// 		{
// 			p = p/2.0;
// 		}
// 		if(fabs(F)<eps)
// 		{
// 			break;
// 		}
// 	}
// 	this->prob = p;
}

void BinomialTest::update_constrained(const Rcpp::NumericMatrix & weights, const int * rows, double r)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double eps = 1e-4;
	double kmax = 20;
	double F, dFdProb, FdivM; // F = dL/dProb
	double n, m;
	double p = this->prob;
	double pnew;

	// Update of prob with Newton Method
// 	time = clock();
	for (int k=0; k<kmax; k++)
	{
		F = dFdProb = 0.0;
		for(int t=0; t<this->obs_total.size(); t++)
		{
			if (this->obs_total[t] >= this->min_obs)
			{
				m = (double)this->obs_test[t];
				n = (double)this->obs_total[t];
				F += weights(rows[0],t) * (m/p + (m-n)/(1-p));
				dFdProb += weights(rows[0],t) * (-m/p/p + (m-n)/(1-p)/(1-p));
				F += weights(rows[1],t) * (m/(p+r) + (m-n)/(2-p-r));
				dFdProb += weights(rows[1],t) * (-m/(p+r)/(p+r) + (m-n)/(2-p-r)/(2-p-r));
			}
		}
		FdivM = F/dFdProb;
		pnew = p - FdivM;
		if ((pnew >= 0) and (pnew <= 1))
		{
			p = pnew;
		}
		else if (pnew < 0)
		{
			p = p/2.0;
		}
		else if (pnew > 1)
		{
			p = p + (1-p)/2.0;
		}
		if(fabs(F)<eps)
		{
			break;
		}
	}
	this->prob = p;

// 	dtime = clock() - time;

}

double BinomialTest::getLogDensityAt(int test, int total)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logdens;
	if (total < this->min_obs)
	{
		logdens = log(1.0/this->min_obs);
	}
	else
	{
		logdens = R::dbinom(test, total, this->prob, true);
	}
	if (std::isnan(logdens))
	{
		throw nan_detected;
	}
	
	return(logdens);
}

// Getter and Setter ------------------------------------------
DensityName BinomialTest::get_name()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(BINOMIAL_TEST);
}

double BinomialTest::get_prob()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(this->prob);
}

void BinomialTest::set_prob(double prob)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->prob = prob;
}


// ============================================================
// Binomial test context dependent
// ============================================================

// Constructor and Destructor ---------------------------------
BinomialTestContext::BinomialTestContext() { }

BinomialTestContext::BinomialTestContext(const Rcpp::IntegerVector & obs_total, const Rcpp::IntegerVector & obs_test, const Rcpp::IntegerVector & context, Rcpp::NumericVector prob, int min_obs, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs_total = obs_total;
	this->obs_test = obs_test;
	this->context = context;
	this->prob = prob;
	this->min_obs = min_obs;
}

BinomialTestContext::~BinomialTestContext()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
}

// Methods ----------------------------------------------------
void BinomialTestContext::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logdens_min = log(1.0/this->min_obs);
	double prob_context;
	for (int t=0; t<this->obs_total.size(); t++)
	{
		if (this->obs_total[t] < this->min_obs)
		{
			logdens[t] = logdens_min;
		}
		else
		{
			prob_context = this->prob[this->context[t]];
			logdens[t] = R::dbinom(this->obs_test[t], this->obs_total[t], prob_context, true);
		}
		if (std::isnan(logdens[t]))
		{
			throw nan_detected;
		}
	}
} 

void BinomialTestContext::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double prob_context;
	double dens_min = 1.0/this->min_obs;
	for (int t=0; t<this->obs_total.size(); t++)
	{
		if (this->obs_total[t] < this->min_obs)
		{
			dens[t] = dens_min;
		}
		else
		{
			prob_context = this->prob[this->context[t]];
			dens[t] = R::dbinom(this->obs_test[t], this->obs_total[t], prob_context, false);
		}
		if (std::isnan(dens[t]))
		{
			if (verbosity >= 4) Rprintf("obs_test[t=%d] = %d, obs_total[t] = %d, prob_context = %g\n", t, obs_test[t], obs_total[t], prob_context);
			throw nan_detected;
		}
	}
} 

void BinomialTestContext::update(const Rcpp::NumericMatrix & weights, const int * rows)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	// Update prob (p)
	double numerator, denominator;
// 	clock_t time, dtime;
// 	time = clock();

	for (int c=0; c<this->prob.size(); c++)
	{
		numerator = 0.0;
		denominator = 0.0;
		for (int t=0; t<this->obs_total.size(); t++)
		{
			if (this->context[t] == c)
			{
				if (this->obs_total[t] >= this->min_obs)
				{
					numerator += weights(rows[0],t) * (this->obs_test[t]);
					denominator += weights(rows[0],t) * (this->obs_total[t]);
				}
			}
		}
		this->prob[c] = numerator/denominator; // Update prob
		if (this->prob[c] > 1)
		{
			if (verbosity >= 4) Rprintf("prob[c=%d] = %g\n", c, prob[c]);
		  throw nan_detected;
		}
	}
	
// 	dtime = clock() - time;

}

void BinomialTestContext::update_constrained_context(const Rcpp::NumericMatrix & weights, const int * rows, Rcpp::NumericVector rs)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double eps = 1e-4;
	double kmax = 20;
	double F, dFdProb, FdivM; // F = dL/dProb
	double n, m;
	double p, r, pnew;

	// Update of prob with Newton Method
// 	time = clock();
	for (int c=0; c<this->prob.size(); c++)
	{
		p = this->prob[c];
		r = rs[c];
		for (int k=0; k<kmax; k++)
		{
			F = 0.0;
			dFdProb = 0.0;
			for(int t=0; t<this->obs_total.size(); t++)
			{
				if (this->context[t] == c)
				{
					if (this->obs_total[t] >= this->min_obs)
					{
						m = (double)this->obs_test[t];
						n = (double)this->obs_total[t];
						F += weights(rows[0],t) * (m/p + (m-n)/(1-p)) + weights(rows[1],t) * (m/(p+r) + (m-n)/(2-p-r));
						dFdProb += weights(rows[0],t) * (-m/p/p + (m-n)/(1-p)/(1-p)) + weights(rows[1],t) * (-m/(p+r)/(p+r) + (m-n)/(2-p-r)/(2-p-r));
					}
				}
			}
			FdivM = F/dFdProb;
			pnew = p - FdivM;
			if ((pnew >= 0) and (pnew <= 1))
			{
				p = pnew;
			}
			else if (pnew < 0)
			{
				p = p/2.0;
			}
			else if (pnew > 1)
			{
			  p = p + (1-p)/2.0;
			}
			if(fabs(F)<eps)
			{
				break;
			}
		}
		this->prob[c] = p;
		if (this->prob[c] > 1)
		{
			if (verbosity >= 4) Rprintf("prob[c=%d] = %g\n", c, prob[c]);
		  throw nan_detected;
		}
	}

// 	dtime = clock() - time;

}

// Getter and Setter ------------------------------------------
DensityName BinomialTestContext::get_name()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(BINOMIAL_TEST_CONTEXT);
}

Rcpp::NumericVector BinomialTestContext::get_probs()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(Rcpp::clone(this->prob));
}

void BinomialTestContext::set_probs(Rcpp::NumericVector prob)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	for (int c=0; c<this->prob.size(); c++)
	{
		this->prob[c] = prob[c];
	}
}


// ============================================================
// Negative Binomial density
// ============================================================

// Constructor and Destructor ---------------------------------
NegativeBinomial::NegativeBinomial() { }

NegativeBinomial::NegativeBinomial(const Rcpp::IntegerVector & obs, double size, double prob, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs = obs;
	this->size = size;
	this->prob = prob;
	this->lxfactorials = NULL;
	// Precompute the lxfactorials that are used in computing the densities
	this->max_obs = Rcpp::max(obs);
	this->lxfactorials = Rcpp::NumericVector(max_obs+1);
	this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of Calloc
	this->lxfactorials[1] = 0.0;
	for (int j=2; j<=max_obs; j++)
	{
		this->lxfactorials[j] = this->lxfactorials[j-1] + log(j);
	}

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
	}
	// Get the index for each observation
	this->uobsind_per_t = Rcpp::IntegerVector(this->obs.size());
	for (int t=0; t<this->obs.size(); t++)
	{
		this->uobsind_per_t[t] = uobsind_per_obs[this->obs[t]];
	}

}

NegativeBinomial::NegativeBinomial(const Rcpp::IntegerVector & obs, const Rcpp::IntegerVector & obs_unique, const Rcpp::IntegerVector & uobsind_per_t, double size, double prob, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs = obs;
	this->obs_unique = obs_unique;
	this->uobsind_per_t = uobsind_per_t;
	this->size = size;
	this->prob = prob;
	this->lxfactorials = NULL;
	// Precompute the lxfactorials that are used in computing the densities
	this->max_obs = Rcpp::max(obs);
	this->lxfactorials = Rcpp::NumericVector(max_obs+1);
	this->lxfactorials[0] = 0.0;	// Not necessary, already 0 because of Calloc
	this->lxfactorials[1] = 0.0;
	for (int j=2; j<=max_obs; j++)
	{
		this->lxfactorials[j] = this->lxfactorials[j-1] + log(j);
	}
}

NegativeBinomial::~NegativeBinomial()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
}

// Methods ----------------------------------------------------
void NegativeBinomial::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);

	// Limit case
	if (this->size == 0 or this->prob == 1)
	{
		for (int t=0; t<this->obs.size(); t++)
		{
			if (this->obs[t] == 0)
			{
				logdens[t] = 0;
			}
			else
			{
				logdens[t] = -INFINITY;
			}
		}
	}
	else
	{

		// Normal case
		double logp = log(this->prob);
		double log1minusp = log(1-this->prob);
		double lGammaR,lGammaRplusX,lxfactorial;
		double obs_j;
		lGammaR = lgamma(this->size);
		// Select strategy for computing gammas, redundant since obs_unique.size() always < obs.size()
		if (this->obs_unique.size() <= this->obs.size())
		{
			std::vector<double> logdens_per_uobs(this->obs_unique.size());
			for (int j=0; j<this->obs_unique.size(); j++)
			{
				obs_j = this->obs_unique[j];
				logdens_per_uobs[j] = lgamma(this->size + obs_j) - lGammaR - this->lxfactorials[obs_j] + this->size * logp + obs_j * log1minusp;
			}
			for (int t=0; t<this->obs.size(); t++)
			{
				logdens[t] = logdens_per_uobs[this->uobsind_per_t[t]];
				if (std::isnan(logdens[t]))
				{
					throw nan_detected;
				}
			}
		}
		else
		{
			for (int t=0; t<this->obs.size(); t++)
			{
				lGammaRplusX = lgamma(this->size + this->obs[t]);
				lxfactorial = this->lxfactorials[(int) this->obs[t]];
				logdens[t] = lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp;
				if (std::isnan(logdens[t]))
				{
					throw nan_detected;
				}
			}
		}

	}

} 

void NegativeBinomial::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);

	// Limit case
	if (this->size == 0 or this->prob == 1)
	{
		for (int t=0; t<this->obs.size(); t++)
		{
			if (this->obs[t] == 0)
			{
				dens[t] = 1;
			}
			else
			{
				dens[t] = 0;
			}
		}
	}
	else
	{
		
		// Normal case
		double logp = log(this->prob);
		double log1minusp = log(1-this->prob);
		double lGammaR,lGammaRplusX,lxfactorial;
		double obs_j;
		lGammaR = lgamma(this->size);
		// Select strategy for computing gammas, redundant since obs_unique.size() always < obs.size()
		if (this->obs_unique.size() <= this->obs.size())
		{
	// clock_t clocktime = clock(), dtime;
			std::vector<double> dens_per_uobs(this->obs_unique.size());
			for (int j=0; j<this->obs_unique.size(); j++)
			{
				obs_j = this->obs_unique[j];
	// 			dens_per_uobs[j] = R::dnbinom(obs_j, this->size, this->prob, 0); // TOO SLOW!!
				dens_per_uobs[j] = exp( lgamma(this->size + obs_j) - lGammaR - this->lxfactorials[obs_j] + this->size * logp + obs_j * log1minusp );
				if (std::isnan(dens_per_uobs[j]))
				{
					if (verbosity>=4) Rprintf("    size = %g, prob = %g, logp = %g, log1minusp = %g\n", size, prob, logp, log1minusp);
					if (verbosity>=4) Rprintf("    lGammaR = %g, lgamma(size + obs=%d) = %g\n", lGammaR, obs_j, lgamma(size + obs_j));
					throw nan_detected;
				}
			}
	// dtime = clock() - clocktime;
	// Rprintf("dtime = %Lg\n", (long double)dtime);
			for (int t=0; t<this->obs.size(); t++)
			{
				dens[t] = dens_per_uobs[this->uobsind_per_t[t]];
			}
		}
		else
		{
			for (int t=0; t<this->obs.size(); t++)
			{
				lGammaRplusX = lgamma(this->size + this->obs[t]);
				lxfactorial = this->lxfactorials[(int) this->obs[t]];
				dens[t] = exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
				if (std::isnan(dens[t]))
				{
					throw nan_detected;
				}
			}
		}

	}

} 

void NegativeBinomial::calc_CDFs(Rcpp::NumericMatrix::Row & CDF)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR;
	lGammaR=lgamma(this->size);
	std::vector<double> precomputed_CDF(this->max_obs+1);
	double dens;

	// Calculate for j=0
	precomputed_CDF[0] = exp( lgamma(this->size) - lGammaR - this->lxfactorials[0] + this->size * logp );
	// Calculate for j>0
	for (int j=1; j<=this->max_obs; j++)
	{
		dens = exp( lgamma(this->size + j) - lGammaR - this->lxfactorials[j] + this->size * logp + j * log1minusp );
		if (std::isnan(dens))
		{
			throw nan_detected;
		}
		precomputed_CDF[j] = precomputed_CDF[j-1] + dens;
		if (precomputed_CDF[j] >= 1)
		{
			precomputed_CDF[j] = precomputed_CDF[j-1]; 
		}
	}
	for (int t=0; t<this->obs.size(); t++)
	{
		CDF[t] = precomputed_CDF[(int)obs[t]];
		if (std::isnan(CDF[t]))
		{
			throw nan_detected;
		}
	}
}

void NegativeBinomial::calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR;
	lGammaR=lgamma(this->size);
	std::vector<double> precomputed_logCDF(this->max_obs+1);
	double logdens;

	// Calculate for j=0
	precomputed_logCDF[0] = lgamma(this->size) - lGammaR - this->lxfactorials[0] + this->size * logp;
	// Calculate for j>0
	for (int j=1; j<=this->max_obs; j++)
	{
		logdens = lgamma(this->size + j) - lGammaR - this->lxfactorials[j] + this->size * logp + j * log1minusp;
		if (std::isnan(logdens))
		{
			throw nan_detected;
		}
		precomputed_logCDF[j] = log( exp(precomputed_logCDF[j-1]) + exp(logdens) );
		if (precomputed_logCDF[j] >= 0)
		{
			precomputed_logCDF[j] = precomputed_logCDF[j-1]; 
		}
	}
	for (int t=0; t<this->obs.size(); t++)
	{
		logCDF[t] = precomputed_logCDF[(int)obs[t]];
		if (std::isnan(logCDF[t]))
		{
			throw nan_detected;
		}
	}
}

void NegativeBinomial::update(const Rcpp::NumericMatrix & weights, const int * rows)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double eps = 1e-4;
	double kmax = 20;
	double numerator, denominator, size0, DigammaSize, TrigammaSize;
	double F, dFdSize, FdivM;
	double logp = log(this->prob);
	// Update prob (p)
	numerator=denominator=0.0;
// 	clock_t time, dtime;
// 	time = clock();
	for (int t=0; t<this->obs.size(); t++)
	{
		numerator += weights(rows[0],t) * this->size;
		denominator += weights(rows[0],t) * (this->size + this->obs[t]);
	}
	this->prob = numerator/denominator; // Update this->prob
// 	logp = log(this->prob); // Update of size is done with new prob
	
// 	dtime = clock() - time;

	// Update of size with Newton Method
	size0 = this->size;
// 	time = clock();
	// Select strategy for computing digammas
	if (this->max_obs <= this->obs.size())
	{
		std::vector<double> DigammaSizePlusX(this->max_obs+1);
		std::vector<double> TrigammaSizePlusX(this->max_obs+1);
		for (int k=0; k<kmax; k++)
		{
			F=dFdSize=0.0;
			DigammaSize = R::digamma(size0); // boost::math::digamma<>(size0);
			TrigammaSize = R::trigamma(size0); // boost::math::digamma<>(size0);
			// Precompute the digammas by iterating over all possible values of the observation vector
			for (int j=0; j<=this->max_obs; j++)
			{
				DigammaSizePlusX[j] = R::digamma(size0+j);
				TrigammaSizePlusX[j] = R::trigamma(size0+j);
			}
			for(int t=0; t<this->obs.size(); t++)
			{
				if(this->obs[t]==0)
				{
					F += weights(rows[0],t) * logp;
					//dFdSize+=0;
				}
				if(this->obs[t]!=0)
				{
					F += weights(rows[0],t) * (logp - DigammaSize + DigammaSizePlusX[(int)obs[t]]);
					dFdSize += weights(rows[0],t) * (-TrigammaSize + TrigammaSizePlusX[(int)obs[t]]);
				}
			}
			FdivM = F/dFdSize;
// Rprintf("k = %d, F = %g, dFdSize = %g, FdivM = %g, size0 = %g\n", k, F, dFdSize, FdivM, size0);
			if (FdivM < size0)
			{
				size0 = size0-FdivM;
			}
			else if (FdivM >= size0)
			{
				size0 = size0/2.0;
			}
			if(fabs(F)<eps)
			{
				break;
			}
		}
	}
	else
	{
		double DigammaSizePlusX, TrigammaSizePlusX;
		for (int k=0; k<kmax; k++)
		{
			F = dFdSize = 0.0;
			DigammaSize = R::digamma(size0); // boost::math::R::digamma<>(size0);
			TrigammaSize = R::trigamma(size0); // boost::math::R::digamma<>(size0);
			for(int t=0; t<this->obs.size(); t++)
			{
				DigammaSizePlusX = R::digamma(size0+this->obs[t]); //boost::math::digamma<>(size0+this->obs[ti]);
				TrigammaSizePlusX = R::trigamma(size0+this->obs[t]);
				if(this->obs[t]==0)
				{
					F += weights(rows[0],t) * logp;
					//dFdSize+=0;
				}
				if(this->obs[t]!=0)
				{
					F += weights(rows[0],t) * (logp - DigammaSize + DigammaSizePlusX);
					dFdSize += weights(rows[0],t) * (-TrigammaSize + TrigammaSizePlusX);
				}
			}
			FdivM = F/dFdSize;
			if (FdivM < size0)
			{
				size0 = size0-FdivM;
			}
			else if (FdivM >= size0)
			{
				size0 = size0/2.0;
			}
			if(fabs(F)<eps)
			{
				break;
			}
		}
	}
	this->size = size0;

// 	dtime = clock() - time;

}

double NegativeBinomial::getLogDensityAt(int x)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	double logdens;
	// Calculate variance
	double mean = 0, variance = 0;
	for(int t=0; t<this->obs.size(); t++)
	{
		mean += obs[t];
	}
	mean = mean / this->obs.size();
	for(int t=0; t<this->obs.size(); t++)
	{
		variance += pow(obs[t] - mean, 2);
	}
	variance = variance / this->obs.size();
	// Calculate logdensity
	lGammaR = lgamma(this->size);
	lGammaRplusX = lgamma(this->size + x);
	lxfactorial = this->lxfactorials[x];
	logdens = lGammaRplusX - lGammaR - lxfactorial + this->size * logp + x * log1minusp;
	if (std::isnan(logdens))
	{
		throw nan_detected;
	}
	
	return(logdens);
}

// Getter and Setter ------------------------------------------
double NegativeBinomial::get_mean()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return this->size*(1-this->prob)/this->prob;
}

double NegativeBinomial::get_variance()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return this->size*(1-this->prob)/this->prob/this->prob;
}

DensityName NegativeBinomial::get_name()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(NEGATIVE_BINOMIAL);
}

double NegativeBinomial::get_size()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(this->size);
}

double NegativeBinomial::get_prob()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(this->prob);
}


// ============================================================
// Zero Inflation density
// ============================================================

// Constructor and Destructor ---------------------------------
ZeroInflation::ZeroInflation(const Rcpp::IntegerVector & obs, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs = obs;
}

ZeroInflation::~ZeroInflation()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
}

// Methods ----------------------------------------------------
void ZeroInflation::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	for (int t=0; t<this->obs.size(); t++)
	{
		if(obs[t]==0)
		{
			logdens[t] = 0.0;
		};
		if(obs[t]>0)
		{
			logdens[t] = -INFINITY;
		}
	}
}

void ZeroInflation::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	for (int t=0; t<this->obs.size(); t++)
	{
		if(obs[t]==0)
		{
			dens[t] = 1.0;
		}
		if(obs[t]>0)
		{
			dens[t] = 0.0;
		}
	}
}

void ZeroInflation::calc_CDFs(Rcpp::NumericMatrix::Row & CDF)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	for (int t=0; t<this->obs.size(); t++)
	{
		CDF[t] = 1.0;
	}
}

void ZeroInflation::calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	for (int t=0; t<this->obs.size(); t++)
	{
		logCDF[t] = 0.0;
	}
}

void ZeroInflation::update(const Rcpp::NumericMatrix &, const int *)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
}

double ZeroInflation::getLogDensityAt(int x)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double logdens;
	// Calculate logdensity
	if (x == 0)
	{
		logdens = 0;
	}
	else
	{
		logdens = -INFINITY;
	}
	
	return(logdens);
}

// Getter and Setter ------------------------------------------
double ZeroInflation::get_mean()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return 0;
}

double ZeroInflation::get_variance()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return 0;
}

DensityName ZeroInflation::get_name()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(ZERO_INFLATION);
}


// ============================================================
// Beta density
// ============================================================

// Constructor and Destructor ---------------------------------
Beta::Beta() { }

Beta_symmetric::Beta_symmetric() { }

Beta_mirror::Beta_mirror() { }

Beta::Beta(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & logObs, const Rcpp::NumericVector & log1mObs, double a, double b, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs = obs;
	this->logObs = logObs;
	this->log1mObs = log1mObs;
	this->a = a;
	this->b = b;
}

Beta_symmetric::Beta_symmetric(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & logObs, const Rcpp::NumericVector & log1mObs, double a, double b, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs = obs;
	this->logObs = logObs;
	this->log1mObs = log1mObs;
	this->a = a;
	this->b = b;
}

Beta_mirror::Beta_mirror(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & logObs, const Rcpp::NumericVector & log1mObs, double a, double b, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs = obs;
	this->logObs = logObs;
	this->log1mObs = log1mObs;
	this->a = a;
	this->b = b;
}

Beta::~Beta() { }

Beta_symmetric::~Beta_symmetric() { }

Beta_mirror::~Beta_mirror() { }

// Methods ----------------------------------------------------
void Beta::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
} 

void Beta::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double cutoff = 1e10;
	for (int t=0; t<this->obs.size(); t++)
	{
		dens[t] = R::dbeta(this->obs[t], this->a, this->b, 0);
		if (dens[t] > cutoff)
		{
			dens[t] = cutoff;
		}
	}
} 

void Beta_mirror::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double cutoff = 1e10;
	for (int t=0; t<this->obs.size(); t++)
	{
		dens[t] = R::dbeta(this->obs[t], this->a, this->b, 0);
		if (dens[t] > cutoff)
		{
			dens[t] = cutoff;
		}
	}
} 

void Beta_symmetric::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double cutoff = 1e10;
	for (int t=0; t<this->obs.size(); t++)
	{
		dens[t] = R::dbeta(this->obs[t], this->a, this->b, 0);
		if (dens[t] > cutoff)
		{
			dens[t] = cutoff;
		}
	}
} 

void Beta::calc_CDFs(Rcpp::NumericMatrix::Row & CDF)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
}

void Beta::calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
}

void Beta::update(const Rcpp::NumericMatrix & weights, const int * rows)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
}

void Beta_mirror::update(const Rcpp::NumericMatrix & weights, const int * rows)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	// Updates with Newton-Raphson
	double eps = 1e-4;
	double kmax = 20;
	double DiC, TriC;
	double F, dFdx, FdivM;
	double a0 = this->get_a();
	double b0 = this->get_b();

	// A
	for (int k=0; k<kmax; k++)
	{
		F = dFdx = 0.0;
		DiC = - R::digamma(a0) + R::digamma(a0+b0);
		TriC = - R::trigamma(a0) + R::trigamma(a0+b0);
		for(int t=0; t<this->obs.size(); t++)
		{
			F += weights(rows[0],t) * ( DiC + this->logObs[t] );
			F += weights(rows[1],t) * ( DiC + this->log1mObs[t] );
			dFdx += ( weights(rows[0],t) + weights(rows[1],t) ) * TriC;
		}
		FdivM = F/dFdx;
		if (FdivM < a0)
		{
			a0 = a0-FdivM;
		}
		else if (FdivM >= a0)
		{
			a0 = a0/2.0;
		}
		if(fabs(F)<eps)
		{
// Rprintf("k(a0) = %d\n", k);
			break;
		}
	}
	// Artificially restrict to values <= 1 to avoid Inf when obs=0 or 1
// 	Rprintf("a0 = %g\n", a0);
	if (a0 > 1)
	{
		a0 = 1;
	}

	// B
	for (int k=0; k<kmax; k++)
	{
		F = dFdx = 0.0;
		DiC = - R::digamma(b0) + R::digamma(a0+b0);
		TriC = - R::trigamma(b0) + R::trigamma(a0+b0);
		for(int t=0; t<this->obs.size(); t++)
		{
			F += weights(rows[0],t) * ( DiC + this->log1mObs[t] );
			F += weights(rows[1],t) * ( DiC + this->logObs[t] );
			dFdx += ( weights(rows[0],t) + weights(rows[1],t) ) * TriC;
		}
		FdivM = F/dFdx;
		if (FdivM < b0)
		{
			b0 = b0-FdivM;
		}
		else if (FdivM >= b0)
		{
			b0 = b0/2.0;
		}
		if(fabs(F)<eps)
		{
// Rprintf("k(b0) = %d\n", k);
			break;
		}
	}
// 	Rprintf("b0 = %g\n", b0);
	// Artificially restrict to values >= 1 to avoid Inf when obs=0 or 1
	if (b0 < 1)
	{
		b0 = 1;
	}

	this->a = a0;
	this->b = b0;

}

void Beta_symmetric::update(const Rcpp::NumericMatrix & weights, const int * rows)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	// Updates with Newton-Raphson
	double eps = 1e-4;
	double kmax = 20;
	double DiC, TriC;
	double F, dFdx, FdivM;
	double logObs, log1mObs;
	double a0 = this->get_a();

	for (int k=0; k<kmax; k++)
	{
		F = dFdx = 0.0;
		DiC = - 2*R::digamma(a0) + 2*R::digamma(a0+a0);
		TriC = - 2*R::trigamma(a0) + 2*R::trigamma(a0+a0);
		for(int t=0; t<this->obs.size(); t++)
		{
			F += weights(rows[0],t) * ( DiC + this->logObs[t] + this->log1mObs[t] );
			dFdx += weights(rows[0],t) * TriC;
		}
		FdivM = F/dFdx;
		if (FdivM < a0)
		{
			a0 = a0-FdivM;
		}
		else if (FdivM >= a0)
		{
			a0 = a0/2.0;
		}
		if(fabs(F)<eps)
		{
			break;
		}
	}
	// Artificially restrict to values >= 1
	if (a0 < 1)
	{
		a0 = 1;
	}

	this->a = a0;
	this->b = a0;

}

double Beta::getLogDensityAt(double x)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(R::dbeta(x, this->a, this->b, 1));
}

// Getter and Setter ------------------------------------------
double Beta::get_mean()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return this->a / ( this->a + this->b );
}

double Beta::get_variance()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double a = this->a;
	double b = this->b;
	return a*b / ( (a+b)*(a+b) * (a+b+1) );
}

DensityName Beta::get_name()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(BETA);
}

DensityName Beta_symmetric::get_name()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(BETA_SYMMETRIC);
}

DensityName Beta_mirror::get_name()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(BETA_MIRROR);
}

double Beta::get_a()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(this->a);
}

double Beta_mirror::get_a()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(this->a);
}

double Beta_symmetric::get_a()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(this->a);
}

double Beta::get_b()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(this->b);
}

double Beta_mirror::get_b()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(this->b);
}

double Beta_symmetric::get_b()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(this->b);
}

void Beta::set_a(double a)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->a = a;
}

void Beta_mirror::set_a(double a)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->a = a;
}

void Beta_symmetric::set_a(double a)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->a = a;
}

void Beta::set_b(double b)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->b = b;
}

void Beta_mirror::set_b(double b)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->b = b;
}
	
void Beta_symmetric::set_b(double b)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->b = b;
}

// ============================================================
// Multivariate Copula Approximation
// ============================================================

// Constructor and Destructor ---------------------------------
MVCopulaApproximation::MVCopulaApproximation(const Rcpp::IntegerMatrix & obs, const Rcpp::IntegerVector & statedef, const Rcpp::List & emissionParamsList, const Rcpp::NumericMatrix & cor_matrix_inv, const double & cor_matrix_det, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs = obs;
	this->cor_matrix_inv = cor_matrix_inv;
	this->cor_matrix_det = cor_matrix_det;
	// Create marginal distributions
	int ndistr = emissionParamsList.size();
	int stateindex;
	for (int imod=0; imod<ndistr; imod++)
	{
		stateindex = statedef[imod] - 1;
		Rcpp::IntegerMatrix::Column iobs = this->obs(Rcpp::_, imod);
		Rcpp::DataFrame emissionParams = Rcpp::as<Rcpp::DataFrame>(emissionParamsList[imod]);
		Rcpp::CharacterVector emissionTypes = emissionParams["type"];
		Rcpp::NumericVector sizes = Rcpp::as<Rcpp::NumericVector>(emissionParams["size"]);
		Rcpp::NumericVector probs = Rcpp::as<Rcpp::NumericVector>(emissionParams["prob"]);
		std::string dtype = Rcpp::as<std::string>(emissionTypes[stateindex]);
		if (dtype.compare("delta") == 0)
		{
			// Zero Inflation
			ZeroInflation * d = new ZeroInflation(iobs, this->verbosity);
			this->marginals.push_back(d);
		}
		else if (dtype.compare("dnbinom") == 0)
		{
			// Negative Binomial
			NegativeBinomial * d = new NegativeBinomial(iobs, sizes[stateindex], probs[stateindex], this->verbosity);
			this->marginals.push_back(d);
		}
	}
}

MVCopulaApproximation::~MVCopulaApproximation()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	for (int imod=0; imod<this->marginals.size(); imod++)
	{
		delete this->marginals[imod];
	}
}

// Methods ----------------------------------------------------
void MVCopulaApproximation::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
// Rprintf("new state\n");
	// Calculate logdensities for marginals
	Rcpp::NumericMatrix marginals_logdensities = Rcpp::NumericMatrix(this->obs.ncol(), this->obs.nrow());
	Rcpp::NumericMatrix marginals_CDFs = Rcpp::NumericMatrix(this->obs.ncol(), this->obs.nrow());
	for (int imod=0; imod<this->obs.ncol(); imod++)
	{
		Rcpp::NumericMatrix::Row marginals_logdensities_row = marginals_logdensities(imod, Rcpp::_);
		this->marginals[imod]->calc_logdensities(marginals_logdensities_row);
		Rcpp::NumericMatrix::Row marginals_CDFs_row = marginals_CDFs(imod, Rcpp::_);
		this->marginals[imod]->calc_CDFs(marginals_CDFs_row);
	}
	// Calculate multivariate Copula approximation
	double sum, uniform, exponent, exponentTemp;
	Rcpp::NumericVector z = Rcpp::NumericVector(this->obs.ncol());
	for (int t=0; t<this->obs.nrow(); t++)
	{
		sum = 0.0;
		for (int imod=0; imod<this->obs.ncol(); imod++)
		{
			sum += marginals_logdensities(imod,t);
			uniform = marginals_CDFs(imod,t);
			z[imod] = R::qnorm(uniform, 0, 1, 1, 0);
			if (std::isnan(z[imod]))
			{
				throw nan_detected;
			}
// if (t==0)
// {
// 	Rprintf("\nmarginal_logdensities[imod=%d][%d] = %g\n", imod, t, marginals_logdensities[imod][t]);
// 	Rprintf("sum = %g\n", sum);
// 	Rprintf("uniform = %g\n", uniform);
// 	Rprintf("z[imod=%d] = %g\n", imod, z[imod]);
// }
		}
		exponent = 0.0;
		for (int imod=0; imod<this->obs.ncol(); imod++)
		{
			exponentTemp = 0.0;
			for(int jmod=0; jmod<obs.ncol(); jmod++)
			{
				if (std::isinf(z[jmod]))
				{
					exponentTemp = INFINITY;
					break;
				}
				if (imod==jmod)
				{
					exponentTemp += z[jmod] * (this->cor_matrix_inv(jmod, imod) - 1);
				}
				else
				{
					exponentTemp += z[jmod] * this->cor_matrix_inv(jmod, imod);
				}
				if (std::isnan(exponentTemp))
				{
					throw nan_detected;
				}
			}
			if (std::isinf(exponentTemp))
			{
				exponent = INFINITY;
				break;
			}
			else
			{
				exponent += exponentTemp * z[imod];
			}
			if (std::isnan(exponent))
			{
				throw nan_detected;
			}
		}
		logdens[t] = -0.5 * log(this->cor_matrix_det) - 0.5 * exponent + sum;
		if (std::isnan(logdens[t]))
		{
			throw nan_detected;
		}		
// if (t==0)
// {
// 	Rprintf("\nlogdens[%d] = %g\n", t, logdens[t]);
// 	Rprintf("-0.5*exponent = %g\n", -0.5*exponent);
// 	Rprintf("sum = %g\n", sum);
// 	Rprintf("cor_matrix_det = %g\n", this->cor_matrix_det);
// 	Rprintf("-0.5*log(cor_matrix_det) = %g\n", -0.5*log(this->cor_matrix_det));
// }
	}

	// Clean up
}

void MVCopulaApproximation::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->calc_logdensities(dens);

	for (int t=0; t<this->obs.nrow(); t++)
	{
		dens[t] = exp( dens[t] );
	}
}

// Getter and Setter ------------------------------------------
DensityName MVCopulaApproximation::get_name()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(OTHER);
}

	
// ============================================================
// Multivariate Product of Bernoullis
// ============================================================

// Constructor and Destructor ---------------------------------
BernoulliProduct::BernoulliProduct(const Rcpp::NumericMatrix & obs, Rcpp::LogicalVector & binary_states, int verbosity)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	this->verbosity = verbosity;
	this->obs = obs;
	this->binary_states = binary_states;
}

BernoulliProduct::~BernoulliProduct()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
}

// Methods ----------------------------------------------------
void BernoulliProduct::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	double d, mult;
	Rcpp::NumericMatrix tempPost = Rcpp::NumericMatrix(this->obs.nrow(), this->obs.ncol());

	for (int t=0; t<this->obs.nrow(); t++)
	{
		d = 1.0;
		for (int imod=0; imod<this->obs.ncol(); imod++)
		{
			//if state[iN] is such that modification imod is unmodified, multiProb[t][imod] is the univariate posterior of being unmodified. 
			//if state[iN] is such that modification imod is modified, multiProb[t][imod] is the univariate posterior of being modified
			if (binary_states[imod])
			{
				mult = 1-this->obs(t,imod);
			}
			else
			{
				mult = this->obs(t,imod);
			}
			if(mult>=1) mult=0.9999999999999;
			if(mult<=0) mult=0.0000000000001;
			d=d*mult;
		}
		logdens[t] = log(d);
	}
}

// Getter and Setter ------------------------------------------
DensityName BernoulliProduct::get_name()
{
	if (verbosity>=2) Rprintf("    %s\n", __PRETTY_FUNCTION__);
	return(OTHER);
}

