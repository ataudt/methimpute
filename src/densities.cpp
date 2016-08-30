#include "densities.h"

// ============================================================
// Zero-inflated Negative Binomial
// ============================================================

// Constructor and Destructor ---------------------------------
ZiNB::ZiNB(const Rcpp::IntegerVector & obs, double size, double prob, double w)
{
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
}

ZiNB::~ZiNB()
{
}

// Methods ----------------------------------------------------
void ZiNB::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->obs.size())
	{
		std::vector<double> logdens_per_read(this->max_obs+1);
		logdens_per_read[0] = log( this->w + (1-this->w) * exp( lgamma(this->size + 0) - lGammaR - lxfactorials[0] + this->size * logp + 0 * log1minusp ) );
		for (int j=1; j<=this->max_obs; j++)
		{
			logdens_per_read[j] = log(1-this->w) + lgamma(this->size + j) - lGammaR - lxfactorials[j] + this->size * logp + j * log1minusp;
		}
		for (int t=0; t<this->obs.size(); t++)
		{
			logdens[t] = logdens_per_read[(int) this->obs[t]];
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
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR = lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->obs.size())
	{
		std::vector<double> dens_per_read(this->max_obs+1);
		dens_per_read[0] = this->w + (1-this->w) * exp( lgamma(this->size + 0) - lGammaR - lxfactorials[0] + this->size * logp + 0 * log1minusp );
		for (int j=1; j<=this->max_obs; j++)
		{
			dens_per_read[j] = (1-this->w) * exp( lgamma(this->size + j) - lGammaR - lxfactorials[j] + this->size * logp + j * log1minusp );
		}
		for (int t=0; t<this->obs.size(); t++)
		{
			dens[t] = dens_per_read[(int) this->obs[t]];
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

void ZiNB::copy(Density* other)
{
	ZiNB* o = (ZiNB*)other;
	this->prob = o->prob;
	this->size = o->size;
	this->obs = o->obs;
	this->w = o->w;
}

double ZiNB::getLogDensityAt(int x)
{
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
	return (1-this->w)*this->size*(1-this->prob)/this->prob;
}

double ZiNB::get_variance()
{
	return (1-this->w)*this->size*(1-this->prob)/this->prob/this->prob; //TODO: Is this correct?
}

DensityName ZiNB::get_name()
{
	return(ZERO_INFLATED_NEGATIVE_BINOMIAL);
}

double ZiNB::get_size()
{
  return(this->size);
}

double ZiNB::get_prob()
{
  return(this->prob);
}

double ZiNB::get_w()
{
	return(this->w);
}


// ============================================================
// Negative Binomial density
// ============================================================

// Constructor and Destructor ---------------------------------
NegativeBinomial::NegativeBinomial() { }

NegativeBinomial::NegativeBinomial(const Rcpp::IntegerVector & obs, double size, double prob)
{
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
}

NegativeBinomial::~NegativeBinomial()
{
}

// Methods ----------------------------------------------------
void NegativeBinomial::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR = lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->obs.size())
	{
		std::vector<double> logdens_per_read(this->max_obs+1);
		for (int j=0; j<=this->max_obs; j++)
		{
			logdens_per_read[j] = lgamma(this->size + j) - lGammaR - lxfactorials[j] + this->size * logp + j * log1minusp;
		}
		for (int t=0; t<this->obs.size(); t++)
		{
			logdens[t] = logdens_per_read[(int) this->obs[t]];
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

void NegativeBinomial::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	double logp = log(this->prob);
	double log1minusp = log(1-this->prob);
	double lGammaR,lGammaRplusX,lxfactorial;
	lGammaR=lgamma(this->size);
	// Select strategy for computing gammas
	if (this->max_obs <= this->obs.size())
	{
		std::vector<double> dens_per_read(this->max_obs+1);
		for (int j=0; j<=this->max_obs; j++)
		{
			dens_per_read[j] = exp( lgamma(this->size + j) - lGammaR - lxfactorials[j] + this->size * logp + j * log1minusp );
		}
		for (int t=0; t<this->obs.size(); t++)
		{
			dens[t] = dens_per_read[(int) this->obs[t]];
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
			dens[t] = exp( lGammaRplusX - lGammaR - lxfactorial + this->size * logp + this->obs[t] * log1minusp );
			if (std::isnan(dens[t]))
			{
				throw nan_detected;
			}
		}
	}
} 

void NegativeBinomial::calc_CDFs(Rcpp::NumericMatrix::Row & CDF)
{
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

void NegativeBinomial::copy(Density* other)
{
	NegativeBinomial* o = (NegativeBinomial*)other;
	this->prob = o->prob;
	this->size = o->size;
	this->obs = o->obs;
}

double NegativeBinomial::getLogDensityAt(int x)
{
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
	return this->size*(1-this->prob)/this->prob;
}

double NegativeBinomial::get_variance()
{
	return this->size*(1-this->prob)/this->prob/this->prob;
}

DensityName NegativeBinomial::get_name()
{
	return(NEGATIVE_BINOMIAL);
}

double NegativeBinomial::get_size()
{
	return(this->size);
}

double NegativeBinomial::get_prob()
{
	return(this->prob);
}


// ============================================================
// Zero Inflation density
// ============================================================

// Constructor and Destructor ---------------------------------
ZeroInflation::ZeroInflation(const Rcpp::IntegerVector & obs)
{
	this->obs = obs;
}

ZeroInflation::~ZeroInflation()
{
}

// Methods ----------------------------------------------------
void ZeroInflation::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
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

void ZeroInflation::copy(Density*) {}

void ZeroInflation::update(const Rcpp::NumericVector &)
{
}

double ZeroInflation::getLogDensityAt(int x)
{
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
	return 0;
}

double ZeroInflation::get_variance()
{
	return 0;
}

DensityName ZeroInflation::get_name()
{
	return(ZERO_INFLATION);
}


// ============================================================
// Beta density
// ============================================================

// Constructor and Destructor ---------------------------------
Beta::Beta() { }

Beta_symmetric::Beta_symmetric() { }

Beta_mirror::Beta_mirror() { }

Beta::Beta(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & logObs, const Rcpp::NumericVector & log1mObs, const Rcpp::IntegerVector & intervals, double a, double b, int numintervals)
{
	this->obs = obs;
	this->logObs = logObs;
	this->log1mObs = log1mObs;
	this->intervals = intervals;
	this->a = a;
	this->b = b;
	this->numintervals=numintervals;
}

Beta_symmetric::Beta_symmetric(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & logObs, const Rcpp::NumericVector & log1mObs, const Rcpp::IntegerVector & intervals, double a, double b, int numintervals)
{
	this->obs = obs;
	this->logObs = logObs;
	this->log1mObs = log1mObs;
	this->intervals = intervals;
	this->a = a;
	this->b = b;
	this->numintervals=numintervals;
}

Beta_mirror::Beta_mirror(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & logObs, const Rcpp::NumericVector & log1mObs, const Rcpp::IntegerVector & intervals, double a, double b, int numintervals)
{
	this->obs = obs;
	this->logObs = logObs;
	this->log1mObs = log1mObs;
	this->intervals = intervals;
	this->a = a;
	this->b = b;
	this->numintervals=numintervals;
}

Beta::~Beta() { }

Beta_symmetric::~Beta_symmetric() { }

Beta_mirror::~Beta_mirror() { }

// Methods ----------------------------------------------------
void Beta::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
} 

void Beta::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	double binwidth = 1.0 / this->numintervals;
	std::vector<double> probs(this->numintervals);
	for (int i1=0; i1<probs.size(); i1++)
	{
		probs[i1] = R::pbeta(binwidth*(i1+1), this->a, this->b, 1, 0) - R::pbeta(binwidth*i1, this->a, this->b, 1, 0);
	}
	for (int t=0; t<this->intervals.size(); t++)
	{
		dens[t] = probs[this->intervals[t]];
	}
} 

void Beta_mirror::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	double binwidth = 1.0 / this->numintervals;
	std::vector<double> probs(this->numintervals);
	for (int i1=0; i1<probs.size(); i1++)
	{
		probs[i1] = R::pbeta(binwidth*(i1+1), this->a, this->b, 1, 0) - R::pbeta(binwidth*i1, this->a, this->b, 1, 0);
	}
	for (int t=0; t<this->intervals.size(); t++)
	{
		dens[t] = probs[this->intervals[t]];
	}
} 

void Beta_symmetric::calc_densities(Rcpp::NumericMatrix::Row & dens)
{
	double binwidth = 1.0 / this->numintervals;
	std::vector<double> probs(this->numintervals);
	for (int i1=0; i1<probs.size(); i1++)
	{
		probs[i1] = R::pbeta(binwidth*(i1+1), this->a, this->b, 1, 0) - R::pbeta(binwidth*i1, this->a, this->b, 1, 0);
	}
	for (int t=0; t<this->intervals.size(); t++)
	{
		dens[t] = probs[this->intervals[t]];
	}
} 

void Beta::calc_CDFs(Rcpp::NumericMatrix::Row & CDF)
{
}

void Beta::calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF)
{
}

void Beta::update(const Rcpp::NumericMatrix & weights, const int * rows)
{
}

void Beta_mirror::update(const Rcpp::NumericMatrix & weights, const int * rows)
{
	// Updates with Newton-Raphson
	double eps = 1e-4;
	double kmax = 20;
	double a0, b0;
	double DiC, TriC;
	double F, dFdx, FdivM;

	// A
	a0 = this->get_a();
	b0 = this->get_b();
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
			break;
		}
	}
	this->a = a0;

	// B
	a0 = this->get_a();
	b0 = this->get_b();
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
			break;
		}
	}
	this->b = b0;

}

void Beta_symmetric::update(const Rcpp::NumericMatrix & weights, const int * rows)
{
	// Updates with Newton-Raphson
	double eps = 1e-4;
	double kmax = 20;
	double a0;
	double DiC, TriC;
	double F, dFdx, FdivM;
	double logObs, log1mObs;

	a0 = this->get_a();
	for (int k=0; k<kmax; k++)
	{
		F = dFdx = 0.0;
		DiC = - 2*R::digamma(a0) + R::digamma(a0+a0);
		TriC = - 2*R::trigamma(a0) + R::trigamma(a0+a0);
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
	return(R::dbeta(x, this->a, this->b, 1));
}

// Getter and Setter ------------------------------------------
double Beta::get_mean()
{
	return this->a / ( this->a + this->b );
}

double Beta::get_variance()
{
	double a = this->a;
	double b = this->b;
	return a*b / ( (a+b)*(a+b) * (a+b+1) );
}

DensityName Beta::get_name()
{
	return(BETA);
}

DensityName Beta_symmetric::get_name()
{
	return(BETA_SYMMETRIC);
}

DensityName Beta_mirror::get_name()
{
	return(BETA_MIRROR);
}

double Beta::get_a()
{
	return(this->a);
}

double Beta_mirror::get_a()
{
	return(this->a);
}

double Beta_symmetric::get_a()
{
	return(this->a);
}

double Beta::get_b()
{
	return(this->b);
}

double Beta_mirror::get_b()
{
	return(this->b);
}

double Beta_symmetric::get_b()
{
	return(this->b);
}

void Beta::set_a(double a)
{
	this->a = a;
}

void Beta_mirror::set_a(double a)
{
	this->a = a;
}

void Beta_symmetric::set_a(double a)
{
	this->a = a;
}

void Beta::set_b(double b)
{
	this->b = b;
}

void Beta_mirror::set_b(double b)
{
	this->b = b;
}
	
void Beta_symmetric::set_b(double b)
{
	this->b = b;
}

// ============================================================
// Multivariate Copula Approximation
// ============================================================

// Constructor and Destructor ---------------------------------
MVCopulaApproximation::MVCopulaApproximation(const Rcpp::IntegerMatrix & obs, const Rcpp::IntegerVector & statedef, const Rcpp::List & emissionParamsList, const Rcpp::NumericMatrix & cor_matrix_inv, const double & cor_matrix_det)
{
	this->obs = obs;
	this->cor_matrix_inv = cor_matrix_inv;
	this->cor_matrix_det = cor_matrix_det;
	// Create marginal distributions
	int ndist = emissionParamsList.size();
	for (int imod=0; imod<ndist; imod++)
	{
		Rcpp::IntegerMatrix::Column iobs = this->obs(Rcpp::_, imod);
		Rcpp::DataFrame emissionParams = Rcpp::as<Rcpp::DataFrame>(emissionParamsList[imod]);
		Rcpp::CharacterVector emissionTypes = emissionParams["type"];
		Rcpp::NumericVector sizes = Rcpp::as<Rcpp::NumericVector>(emissionParams["size"]);
		Rcpp::NumericVector probs = Rcpp::as<Rcpp::NumericVector>(emissionParams["prob"]);
		NegativeBinomial * d0 = new NegativeBinomial(iobs, sizes[statedef[imod]-1], probs[statedef[imod]-1]);
		this->marginals.push_back(d0);
	}
}

MVCopulaApproximation::~MVCopulaApproximation()
{
	for (int imod=0; imod<this->obs.ncol(); imod++)
	{
		delete this->marginals[imod];
	}
}

// Methods ----------------------------------------------------
void MVCopulaApproximation::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
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
// 		logdens[t] = log(1/sqrt(this->cor_matrix_det)) - 0.5 * exponent + sum;
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
	this->calc_logdensities(dens);

	for (int t=0; t<this->obs.nrow(); t++)
	{
		dens[t] = exp( dens[t] );
	}
}

// Getter and Setter ------------------------------------------
DensityName MVCopulaApproximation::get_name()
{
	return(OTHER);
}

	
// ============================================================
// Multivariate Product of Bernoullis
// ============================================================

// Constructor and Destructor ---------------------------------
BernoulliProduct::BernoulliProduct(const Rcpp::NumericMatrix & obs, Rcpp::LogicalVector & binary_states)
{
	this->obs = obs;
	this->binary_states = binary_states;
}

BernoulliProduct::~BernoulliProduct()
{
}

// Methods ----------------------------------------------------
void BernoulliProduct::calc_logdensities(Rcpp::NumericMatrix::Row & logdens)
{
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
	return(OTHER);
}

