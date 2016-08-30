#include <Rcpp.h>
#include "scalehmm.h"

using namespace Rcpp;

static ScaleHMM * hmm; // declare as static outside the function because we only need one and this enables memory-cleanup on R_CheckUserInterrupt()

// [[Rcpp::export]]
List fitHMM(const IntegerVector & counts, const NumericVector & distances, const List & params, const int & algorithm) {

    // access variables by name
    const NumericVector startProbs_initial = as<NumericVector>(params["startProbs_initial"]);
    const NumericMatrix transProbs_initial = as<NumericMatrix>(params["transProbs_initial"]);
		const double transDist = as<double>(params["transDist"]);
		const DataFrame emissionParams_initial = as<DataFrame>(params["emissionParams_initial"]);
    const double eps = as<double>(params["eps"]);
    const double maxtime = as<double>(params["maxtime"]);
    const double maxiter = as<double>(params["maxiter"]);
		const int verbosity = as<int>(params["verbosity"]);

		// Initialize the HMM
		hmm = new ScaleHMM(counts, distances, startProbs_initial, transProbs_initial, transDist, emissionParams_initial, verbosity);

		// Estimate parameters
		List results;
		std::string error = "";
		if (algorithm == 1)
		{
			try {
					results = hmm->baumWelch(eps, maxiter, maxtime);
			} catch (std::exception & e) {
					if (verbosity>=1) Rprintf("HMM: Error in Baum-Welch: %s\n", e.what());
					error = e.what();
			}
		}
		else if (algorithm == 2)
		{
			try {
					results = hmm->viterbi(eps, maxiter, maxtime);
			} catch (std::exception & e) {
					if (verbosity>=1) Rprintf("HMM: Error in Viterbi: %s\n", e.what());
					error = e.what();
			}
		}
		results.push_back(error, "error");

		// Delete the HMM
		delete hmm;
		hmm = NULL; // assign NULL to defuse the additional delete in on.exit() call

		// Return
		return results;

}

// [[Rcpp::export]]
List fitHMMratio(const NumericVector & ratio, const NumericVector & distances, const List & params, const int & algorithm) {

    // access variables by name
    const NumericVector startProbs_initial = as<NumericVector>(params["startProbs_initial"]);
    const NumericMatrix transProbs_initial = as<NumericMatrix>(params["transProbs_initial"]);
		const double transDist = as<double>(params["transDist"]);
		const DataFrame emissionParams_initial = as<DataFrame>(params["emissionParams_initial"]);
    const double eps = as<double>(params["eps"]);
    const double maxtime = as<double>(params["maxtime"]);
    const double maxiter = as<double>(params["maxiter"]);
		const int verbosity = as<int>(params["verbosity"]);

		// Initialize the HMM
		hmm = new ScaleHMM(ratio, distances, startProbs_initial, transProbs_initial, transDist, emissionParams_initial, verbosity);

		// Estimate parameters
		List results;
		std::string error = "";
		if (algorithm == 1)
		{
			try {
					results = hmm->baumWelch(eps, maxiter, maxtime);
			} catch (std::exception & e) {
					if (verbosity>=1) Rprintf("HMM: Error in Baum-Welch: %s\n", e.what());
					error = e.what();
			}
		}
		else if (algorithm == 2)
		{
			try {
					results = hmm->viterbi(eps, maxiter, maxtime);
			} catch (std::exception & e) {
					if (verbosity>=1) Rprintf("HMM: Error in Viterbi: %s\n", e.what());
					error = e.what();
			}
		}
		results.push_back(error, "error");

		// Delete the HMM
		delete hmm;
		hmm = NULL; // assign NULL to defuse the additional delete in on.exit() call

		// Return
		return results;

}


// [[Rcpp::export]]
List fitMultiHMM(const IntegerMatrix & counts, const NumericVector & distances, const List & params) {

    // access variables by name
    const NumericVector startProbs_initial = as<NumericVector>(params["startProbs_initial"]);
    const NumericMatrix transProbs_initial = as<NumericMatrix>(params["transProbs_initial"]);
		const double transDist = as<double>(params["transDist"]);
		const List emissionParamsList = as<List>(params["emissionParamsList"]);
    const double eps = as<double>(params["eps"]);
    const double maxtime = as<double>(params["maxtime"]);
    const double maxiter = as<double>(params["maxiter"]);
		const int verbosity = as<int>(params["verbosity"]);
		const NumericVector determinant = as<NumericVector>(params["determinant"]);
		const List cor_mat_inv = as<List>(params["correlationMatrixInverse"]);
		const DataFrame statedef = as<DataFrame>(params["statedef"]);

		// Initialize the HMM
		hmm = new ScaleHMM(counts, distances, startProbs_initial, transProbs_initial, transDist, emissionParamsList, verbosity, cor_mat_inv, determinant, statedef);

		// Estimate parameters
		List results;
		std::string error = "";
		try {
				results = hmm->baumWelch(eps, maxiter, maxtime);
		} catch (std::exception & e) {
				if (verbosity>=1) Rprintf("HMM: Error in Baum-Welch: %s\n", e.what());
				error = e.what();
		}
		results.push_back(error, "error");

		// Delete the HMM
		delete hmm;
		hmm = NULL; // assign NULL to defuse the additional delete in on.exit() call

		// Return
		return results;

}

// =======================================================
// This function make a cleanup if anything was left over
// =======================================================
// [[Rcpp::export]]
void cleanup()
{
	delete hmm;
}

