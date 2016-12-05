#ifndef SCALEHMM_H
#define SCALEHMM_H

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

class ScaleHMM  {

	public:
		// Constructor and Destructor
		ScaleHMM();
		ScaleHMM(const Rcpp::NumericVector & obs, const Rcpp::NumericVector & distances, Rcpp::NumericVector startProbs_initial, Rcpp::NumericMatrix transProbs_initial, double transDist, Rcpp::DataFrame emissionParams_initial, int verbosity);
		ScaleHMM(const Rcpp::IntegerVector & obs, const Rcpp::NumericVector & distances, Rcpp::NumericVector startProbs_initial, Rcpp::NumericMatrix transProbs_initial, double transDist, Rcpp::DataFrame emissionParams_initial, int verbosity);
		ScaleHMM(const Rcpp::IntegerVector & obs_total, const Rcpp::IntegerVector & obs_meth, const Rcpp::IntegerVector & obs_unmeth, const Rcpp::NumericVector & distances, Rcpp::NumericVector startProbs_initial, Rcpp::NumericMatrix transProbs_initial, double transDist, Rcpp::DataFrame emissionParams_initial, int verbosity);
		ScaleHMM(const Rcpp::IntegerMatrix & multi_obs, const Rcpp::NumericVector & distances, Rcpp::NumericVector startProbs_initial, Rcpp::NumericMatrix transProbs_initial, double transDist, Rcpp::List emissionParamsList, int verbosity, const Rcpp::List & cor_mat_inv, const Rcpp::NumericVector & determinant, const Rcpp::DataFrame & statedef);
		~ScaleHMM();

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
		double get_transProbs(int i, int j);
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
		Rcpp::NumericMatrix transProbs; ///< matrix [NSTATES x NSTATES] of transition probabilities
		Rcpp::NumericVector transExp; ///< vector [NDATA] with exponential factors for decay of transition probabilities
// 		double transDist; ///< characteristic decaying constant for the transition probabilities
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
		void print_uni_iteration(int iteration);
		void print_multi_iteration(int iteration);
		void print_multi_params();
};

#endif
