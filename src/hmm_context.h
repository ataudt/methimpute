#ifndef HMM_CONTEXT_H
#define HMM_CONTEXT_H

#include "densities.h"
#include <Rcpp.h>
#include <R.h> // R_CheckUserInterrupt()
#include <vector> // storing density functions
#include <time.h> // time(), difftime()
#include <string> // strcmp

#ifdef _OPENMP
#include <omp.h>
#endif

enum whichvariate {UNIVARIATE, MULTIVARIATE};

class HMM_context  {

	public:
		// Constructor and Destructor
		HMM_context();
		// Binomial test context transition
		HMM_context(const Rcpp::IntegerVector & obs_total, const Rcpp::IntegerVector & obs_meth, const Rcpp::IntegerVector & context, const Rcpp::IntegerVector & transitionContext, const Rcpp::NumericVector & distances, Rcpp::NumericVector startProbs_initial, Rcpp::List transProbs_initial, Rcpp::NumericVector transDist, Rcpp::List emissionParams_initial, int min_obs, int verbosity);
		~HMM_context();

		// Member variables

		// Methods
		Rcpp::List forward_backward(double eps, double maxiter, double maxtime);
		Rcpp::List baumWelch(double eps, double maxiter, double maxtime);
		Rcpp::NumericVector calc_weights();
		void calc_weights(Rcpp::NumericVector & weights);

		// Getters and Setters
		int get_NSTATES();
		int get_NDATA();
		void get_posteriors(Rcpp::NumericMatrix & post);
		double get_posterior(int iN, int t);
		double get_density(int iN, int t);
		double get_startProbs(int i);
		double get_loglik();

	private:
		// Member variables
		int verbosity; ///< verbosity level
		int NDATA; ///< length of observed sequence
		int NSTATES; ///< number of states
		int NMOD; ///< number of modifications / marks
		Rcpp::NumericVector logObs; ///< vector [NDATA] of log(observations)
		Rcpp::NumericVector log1mObs; ///< vector [NDATA] of log(1-observations)
		Rcpp::IntegerVector obs_unique; ///< vector [?] of unique observations
		Rcpp::IntegerVector uobsind_per_t; ///< vector [NDATA] of indices of unique observations for each element in obs
		Rcpp::List transProbsList; ///< List with matrices [NSTATES x NSTATES] of transition probabilities
		Rcpp::NumericVector transDist; ///< constant for decay of transition probabilities
		Rcpp::NumericVector transExp; ///< vector [NDATA] with exponential factors for decay of transition probabilities
		Rcpp::IntegerVector transitionContext; ///< vector [NDATA] with transition contexts
		Rcpp::NumericVector startProbs; ///< initial probabilities [NSTATES]
		double loglik; ///< loglikelihood
		Rcpp::NumericVector distances; ///< vector [NDATA] of distances between observations
		Rcpp::NumericVector scalefactoralpha; ///< vector[NDATA] of scaling factors
		Rcpp::NumericMatrix scalealpha; ///< matrix [NDATA x NSTATES] of forward probabilities
		Rcpp::NumericMatrix scalebeta; ///<  matrix [NDATA x NSTATES] of backward probabilities
		Rcpp::NumericMatrix densities; ///< matrix [NSTATES x NDATA] of density values
// 		Rcpp::NumericMatrix tdensities; ///< matrix [NDATA x NSTATES] of density values, for use in multivariate !increases speed, but on cost of RAM usage and that seems to be limiting
		Rcpp::NumericVector sumgamma; ///< vector[NSTATES] of sum of posteriors (gamma values)
		Rcpp::NumericMatrix sumxi; ///< matrix[NSTATES x NSTATES] of xi values
		Rcpp::NumericMatrix gamma; ///< matrix[NSTATES x NDATA] of posteriors
		double dloglik; ///< difference in loglikelihood from one iteration to the next
		time_t baumWelchStartTime_sec; ///< start time of the Baum-Welch in sec
		int baumWelchTime_real; ///< elapsed time from start of the 0th iteration
		int sumdiff_state_last; ///< sum of the difference in the state 1 assignments from one iteration to the next (only univariate case)
		whichvariate xvariate; ///< enum which stores if UNIVARIATE or MULTIVARIATE
		Rcpp::DataFrame emissionParams; ///< DataFrame with emission parameters for the distributions
		Rcpp::List emissionParamsList; ///< List with DataFrames with emission parameters for multivariate
		std::vector<Density*> emissionDensities; ///< density functions for each state
		
		// Methods
		void forward();
		void backward();
		void calc_sumgamma();
		void calc_sumxi();
		void calc_loglikelihood();
		void calc_densities();
		void update_transProbs();
		void update_startProbs();
		void update_transDist();
		void print_uni_iteration(int iteration);
		void print_multi_iteration(int iteration);
		void print_multi_params();
};

#endif
