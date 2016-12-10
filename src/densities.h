#include <Rcpp.h>
#include <Rmath.h> // digamma() and qnorm()
#include <vector> // storing density functions in MVCopula

#ifndef DENSITIES_H
#define DENSITIES_H

enum DensityName {ZERO_INFLATION, BINOMIAL_TEST, BINOMIAL_TEST_CONTEXT, NEGATIVE_BINOMIAL, ZERO_INFLATED_NEGATIVE_BINOMIAL, BETA, BETA_MIRROR, BETA_SYMMETRIC, OTHER};

/* custom error handling class */
static class exception_nan: public std::exception
{
  virtual const char* what() const throw()
  {
    return "nan detected";
  }
} nan_detected; // this line creates an object of this class

class Density {
	public:
		// Constructor and Destructor
		virtual ~Density() {};
		// Methods
		virtual void calc_logdensities(Rcpp::NumericMatrix::Row &) {};
		virtual void calc_densities(Rcpp::NumericMatrix::Row &) {};
		virtual void calc_logCDFs(Rcpp::NumericMatrix::Row &) {};
		virtual void calc_CDFs(Rcpp::NumericMatrix::Row &) {};
		virtual void update(const Rcpp::NumericMatrix &, const int * rows) {}; 
		virtual void update_constrained(const Rcpp::NumericMatrix &, const int * rows, double r) {}; 
		virtual void update_constrained_context(const Rcpp::NumericMatrix &, const int * rows, Rcpp::NumericVector rs) {}; 
		virtual double getLogDensityAt(int) { return(0); };
		// Getter and Setter
		virtual DensityName get_name() { return(OTHER); };
		virtual double get_mean() { return(0); };
		virtual double get_variance() { return(0); };
		virtual double get_size() { return(0); };
		virtual double get_prob() { return(0); };
		virtual Rcpp::NumericVector get_probs() { return(0); };
		virtual void set_prob(double prob) {};
		virtual void set_probs(Rcpp::NumericVector prob) {};
		virtual double get_a() { return(0); };
		virtual double get_b() { return(0); };
		virtual void set_a(double a) {};
		virtual void set_b(double b) {};

};  


class ZiNB : public Density {
	public:
		// Constructor and Destructor
		ZiNB();
		ZiNB(const Rcpp::IntegerVector & obs, double size, double prob, double w, int verbosity);
		ZiNB(const Rcpp::IntegerVector & obs, const Rcpp::IntegerVector & obs_unique, const Rcpp::IntegerVector & uobsind_per_obs, double size, double prob, double w, int verbosity);
		~ZiNB();
	
		// Methods
		void calc_logdensities(Rcpp::NumericMatrix::Row & logdens);
		void calc_densities(Rcpp::NumericMatrix::Row & logdens);
		void calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF);
		void calc_CDFs(Rcpp::NumericMatrix::Row & CDF);
		double getLogDensityAt(int x);
	
		// Getter and Setter
		double get_mean();
		double get_variance();
		DensityName get_name();
		double get_size();
		double get_prob();
		double get_w();

	private:
		// Member variables
		int verbosity; ///< verbosity parameter for debugging
		double size; ///< parameter of the distribution
		double prob; ///< parameter of the distribution
		double w; ///< parameter of the distribution
		Rcpp::IntegerVector obs; ///< vector [NDATA] of observations
		Rcpp::IntegerVector obs_unique; ///< vector [?] of unique observations
		Rcpp::IntegerVector uobsind_per_t; ///< index of unique observations for each element in obs
		Rcpp::NumericVector weight; ///< temporary storage for weights in update()
		int max_obs; ///< maximum value in obs
		Rcpp::NumericVector lxfactorials; ///< vector [max_obs] of precomputed factorials (x!)
};

class BinomialTest : public Density {
	public:
		// Constructor and Destructor
		BinomialTest();
		BinomialTest(const Rcpp::IntegerVector & obs_total, const Rcpp::IntegerVector & obs_test, double prob, int min_obs, int verbosity);
		~BinomialTest();

		// Methods
		void calc_logdensities(Rcpp::NumericMatrix::Row & logdens);
		void calc_densities(Rcpp::NumericMatrix::Row & dens);
		void update(const Rcpp::NumericMatrix & weights, const int * rows);
		void update_constrained(const Rcpp::NumericMatrix & weights, const int * rows, double r);
		double getLogDensityAt(int test, int total);

		// Getter and Setter
		DensityName get_name();
		double get_prob();
		void set_prob(double prob);

	private:
		// Member variables
		int verbosity; ///< verbosity parameter for debugging
		double prob; ///< parameter of the distribution
		Rcpp::IntegerVector obs_total; ///< vector [NDATA] of observations
		Rcpp::IntegerVector obs_test; ///< vector [NDATA] of observations
		int min_obs; ///< positions with fewer than min_obs reads will be assigned 0.5
// 		int max_obs; ///< maximum value in obs
// 		Rcpp::NumericVector lxfactorials; ///< vector [max_obs] of precomputed factorials (x!)
};


class BinomialTestContext : public Density {
	public:
		// Constructor and Destructor
		BinomialTestContext();
		BinomialTestContext(const Rcpp::IntegerVector & obs_total, const Rcpp::IntegerVector & obs_test, const Rcpp::IntegerVector & context, Rcpp::NumericVector prob, int min_obs, int verbosity);
		~BinomialTestContext();

		// Methods
		void calc_logdensities(Rcpp::NumericMatrix::Row & logdens);
		void calc_densities(Rcpp::NumericMatrix::Row & dens);
		void update(const Rcpp::NumericMatrix & weights, const int * rows);
		void update_constrained_context(const Rcpp::NumericMatrix & weights, const int * rows, Rcpp::NumericVector rs);

		// Getter and Setter
		DensityName get_name();
		Rcpp::NumericVector get_probs();
		void set_probs(Rcpp::NumericVector prob);

	private:
		// Member variables
		int verbosity; ///< verbosity parameter for debugging
		Rcpp::NumericVector prob; ///< vector with context-dependent parameters of the distribution
		Rcpp::IntegerVector obs_total; ///< vector [NDATA] of observations
		Rcpp::IntegerVector obs_test; ///< vector [NDATA] of observations
		Rcpp::IntegerVector context; ///< vector [NDATA] of contexts
		int min_obs; ///< positions with fewer than min_obs reads will be assigned 0.5
// 		int max_obs; ///< maximum value in obs
// 		Rcpp::NumericVector lxfactorials; ///< vector [max_obs] of precomputed factorials (x!)
};


class NegativeBinomial : public Density {
	public:
		// Constructor and Destructor
		NegativeBinomial();
		NegativeBinomial(const Rcpp::IntegerVector & obs, double size, double prob, int verbosity);
		NegativeBinomial(const Rcpp::IntegerVector & obs, const Rcpp::IntegerVector & obs_unique, const Rcpp::IntegerVector & uobsind_per_obs, double size, double prob, int verbosity);
		~NegativeBinomial();

		// Methods
		void calc_logdensities(Rcpp::NumericMatrix::Row & logdens);
		void calc_densities(Rcpp::NumericMatrix::Row & dens);
		void calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF);
		void calc_CDFs(Rcpp::NumericMatrix::Row & CDF);
		void update(const Rcpp::NumericMatrix & weights, const int * rows);
		double getLogDensityAt(int x);

		// Getter and Setter
		double get_mean();
		double get_variance();
		DensityName get_name();
		double get_size();
		double get_prob();

	private:
		// Member variables
		int verbosity; ///< verbosity parameter for debugging
		double size; ///< parameter of the distribution
		double prob; ///< parameter of the distribution
		Rcpp::IntegerVector obs; ///< vector [NDATA] of observations
		Rcpp::IntegerVector obs_unique; ///< vector [?] of unique observations
		Rcpp::IntegerVector uobsind_per_t; ///< index of unique observations for each element in obs
		int max_obs; ///< maximum value in obs
		Rcpp::NumericVector lxfactorials; ///< vector [max_obs] of precomputed factorials (x!)
};


class ZeroInflation : public Density {
	public:
		// Constructor and Destructor
		ZeroInflation();
		ZeroInflation(const Rcpp::IntegerVector & obs, int verbosity);
		~ZeroInflation();

		// Methods
		void calc_logdensities(Rcpp::NumericMatrix::Row & logdens);
		void calc_densities(Rcpp::NumericMatrix::Row & dens);
		void calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF);
		void calc_CDFs(Rcpp::NumericMatrix::Row & CDF);
		void update(const Rcpp::NumericMatrix & weights, const int * rows);
		double getLogDensityAt(int x);

		// Getters and Setters
		double get_mean();
		double get_variance();
		DensityName get_name();

	private:
		// Member variables
		int verbosity; ///< verbosity parameter for debugging
		Rcpp::IntegerVector obs; ///< vector [NDATA] of observations
};


class Beta : public Density {
	public:
		// Constructor and Destructor
		Beta();
		Beta(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & logObs, const Rcpp::NumericVector & log1mObs, double a, double b, int verbosity);
		~Beta();

		// Methods
		void calc_logdensities(Rcpp::NumericMatrix::Row & logdens);
		void calc_densities(Rcpp::NumericMatrix::Row & dens);
		void calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF);
		void calc_CDFs(Rcpp::NumericMatrix::Row & CDF);
		void update(const Rcpp::NumericMatrix & weights, const int * rows);
		double getLogDensityAt(double x);

		// Getter and Setter
		double get_mean();
		double get_variance();
		DensityName get_name();
		double get_a();
		double get_b();
		void set_a(double a);
		void set_b(double b);

	private:
		// Member variables
		int verbosity; ///< verbosity parameter for debugging
		double a; ///< parameter of the distribution
		double b; ///< parameter of the distribution
		Rcpp::NumericVector obs; ///< vector [NDATA] of observations
		Rcpp::NumericVector logObs; ///< vector [NDATA] of log(observations)
		Rcpp::NumericVector log1mObs; ///< vector [NDATA] of log(1-observations)
};


class Beta_mirror : public Beta {
	public:
		// Constructor and Destructor
		Beta_mirror();
		Beta_mirror(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & logObs, const Rcpp::NumericVector & log1mObs, double a, double b, int verbosity);
		~Beta_mirror();

		// Methods
// 		void calc_logdensities(Rcpp::NumericMatrix::Row & logdens);
		void calc_densities(Rcpp::NumericMatrix::Row & dens);
// 		void calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF);
// 		void calc_CDFs(Rcpp::NumericMatrix::Row & CDF);
		void update(const Rcpp::NumericMatrix & weights, const int * rows);
// 		double getLogDensityAt(double x);

		// Getter and Setter
// 		double get_mean();
// 		double get_variance();
		DensityName get_name();
		double get_a();
		double get_b();
		void set_a(double a);
		void set_b(double b);

	private:
		// Member variables
		int verbosity; ///< verbosity parameter for debugging
		double a; ///< parameter of the distribution
		double b; ///< parameter of the distribution
		Rcpp::NumericVector obs; ///< vector [NDATA] of observations
		Rcpp::NumericVector logObs; ///< vector [NDATA] of log(observations)
		Rcpp::NumericVector log1mObs; ///< vector [NDATA] of log(1-observations)
};


class Beta_symmetric : public Beta {
	public:
		// Constructor and Destructor
		Beta_symmetric();
		Beta_symmetric(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & logObs, const Rcpp::NumericVector & log1mObs, double a, double b, int verbosity);
		~Beta_symmetric();

		// Methods
// 		void calc_logdensities(Rcpp::NumericMatrix::Row & logdens);
		void calc_densities(Rcpp::NumericMatrix::Row & dens);
// 		void calc_logCDFs(Rcpp::NumericMatrix::Row & logCDF);
// 		void calc_CDFs(Rcpp::NumericMatrix::Row & CDF);
		void update(const Rcpp::NumericMatrix & weights, const int * rows);
// 		double getLogDensityAt(double x);

		// Getter and Setter
// 		double get_mean();
// 		double get_variance();
		DensityName get_name();
		double get_a();
		double get_b();
		void set_a(double a);
		void set_b(double b);

	private:
		// Member variables
		int verbosity; ///< verbosity parameter for debugging
		double a; ///< parameter of the distribution
		double b; ///< parameter of the distribution
		Rcpp::NumericVector obs; ///< vector [NDATA] of observations
		Rcpp::NumericVector logObs; ///< vector [NDATA] of log(observations)
		Rcpp::NumericVector log1mObs; ///< vector [NDATA] of log(1-observations)
};


class MVCopulaApproximation : public Density {
	public:
		// Constructor and Destructor
		MVCopulaApproximation(const Rcpp::IntegerMatrix & obs, const Rcpp::IntegerVector & statedef, const Rcpp::List & emissionParamsList, const Rcpp::NumericMatrix & cor_matrix_inv, const double & cor_matrix_det, int verbosity);
		~MVCopulaApproximation();
	
		// Methods
		void calc_logdensities(Rcpp::NumericMatrix::Row & logdens);
		void calc_densities(Rcpp::NumericMatrix::Row & dens);

		// Getters and Setters
		DensityName get_name();

	private:
		// Member variables
		int verbosity; ///< verbosity parameter for debugging
		Rcpp::IntegerMatrix obs; ///< matrix [Nmod x NDATA] of observations
		std::vector<Density*> marginals; ///< vector [Nmod] of marginal distributions
		Rcpp::NumericMatrix cor_matrix_inv; ///< vector with elements of the inverse of the correlation matrix
		double cor_matrix_det; ///< determinant of the correlation matrix
};


class BernoulliProduct : public Density {
	public:
		// Constructor and Destructor
		BernoulliProduct(const Rcpp::NumericMatrix & obs, Rcpp::LogicalVector & binary_states, int verbosity);
		~BernoulliProduct();
		// Methods
		void calc_logdensities(Rcpp::NumericMatrix::Row & logdens);
		// Getters and Setters
		DensityName get_name();

	private:
		// Member variables
		int verbosity; ///< verbosity parameter for debugging
		Rcpp::NumericMatrix obs;
		Rcpp::LogicalVector binary_states;
};


#endif
