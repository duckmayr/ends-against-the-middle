// [[Rcpp::depends(RcppDist)]]
#include <Rcpp.h>
#include <4beta.h>

// Function for log probability of response given respondent & item parameters
// [[Rcpp::export]]
double logprob(const int choice, const double th, const double a,
               const double d, const Rcpp::NumericVector& t){
    int K = t.size();
    double result = 0, numerator = 0, denominator = 0, tSum = 0;
    for ( int k = 0; k < K; ++k ) {
        tSum += t[k];
        numerator = std::exp(a * (k * (th - d) - tSum));
        numerator += std::exp(a * ((2*K - 1 - k) * (th - d) - tSum));
        denominator += numerator;
        if ( k == choice ) {
            result = numerator;
        }
    }
    return std::log(result) - std::log(denominator);
}

/* This gives the log likelihood of data given a set of posterior samples.
 * Notably, this is set up specifically for dichotomous data;
 * some tweaks would need to be made for polytomous data.
 */
// [[Rcpp::export]]
double ggum_log_likelihood(const Rcpp::IntegerMatrix& data,
                           const Rcpp::NumericMatrix& samples) {
    double res = 0.0;             // Log likelihood
    double prob = 0.0;            // Prob. of response i,j for sample s
    double tmp = 0.0;             // Mean prob. of response i,j across samples
    int n = data.nrow();          // Number of respondents
    int m = data.ncol();          // Number of items
    int S = samples.nrow();       // Number of posterior samples
    double M = std::log(1.0 / S); // Weight for prob. from each sample
    Rcpp::NumericVector tau(2);   // NumericVector for tau, initially w/ all 0s
    double progtot = (double)n*m; // Variables to help track progress
    int progcurrent = 0;
    for ( int j = 0; j < m; ++j ) {    
        for ( int i = 0; i < n; ++i ) {
            progcurrent += 1;
            Rprintf("\r%4.2f", progcurrent / progtot);
            Rcpp::checkUserInterrupt();
            if ( !Rcpp::IntegerVector::is_na(data(i, j)) ) {
                tmp = 0.0;
                for ( int s = 0; s < S; ++s ) {
                    tau[1] = samples(s, n + 2*m + j);
                    prob = logprob(data(i, j), samples(s, i), samples(s, n + j),
                                   samples(s, n + m + j), tau);
                    tmp += std::exp(M + prob);
                }
            } else {
                tmp = 1.0;
            }
            res += std::log(tmp);
        }
    }
    // Return result
    return res;
}

// [[Rcpp::export]]
double cjr_log_likelihood(const Rcpp::IntegerMatrix& data,
                          const Rcpp::NumericMatrix& samples) {
    double res = 0.0;             // Log likelihood
    double prob = 0.0;            // Prob. of response i,j for sample s
    double tmp = 0.0;             // Mean prob. of response i,j across samples
    int n = data.nrow();          // Number of respondents
    int m = data.ncol();          // Number of items
    int S = samples.nrow();       // Number of posterior samples
    double M = std::log(1.0 / S); // Weight for prob. from each sample
    double progtot = (double)n*m; // Variables to help track progress
    int progcurrent = 0;
    for ( int j = 0; j < m; ++j ) {    
        for ( int i = 0; i < n; ++i ) {
            progcurrent += 1;
            Rprintf("\r%4.2f", progcurrent / progtot);
            Rcpp::checkUserInterrupt();
            if ( !Rcpp::IntegerVector::is_na(data(i, j)) ) {
                tmp = 0.0;
                for ( int s = 0; s < S; ++s ) {
                    // MCMCirt1d stores parameters like theta1, ..., thetaN,
                    // alpha1, beta1, ..., alphaM, betaM
                    double mu = samples(s, i) * samples(s,n + 2 * j + 1) -
                        samples(s, n + 2 * j);
                    prob = R::pnorm(mu, 0.0, 1.0, data(i, j), 1);
                    tmp += std::exp(M + prob);
                }
            } else {
                tmp = 1.0;
            }
            res += std::log(tmp);
        }
    }
    Rprintf("\r1.00\n");
    // Return result
    return res;
}

// [[Rcpp::export]]
double mq_log_likelihood(const Rcpp::IntegerMatrix& data,
                         const Rcpp::NumericMatrix& samples,
                         const Rcpp::IntegerVector& time_map) {
    double res = 0.0;             // Log likelihood
    double prob = 0.0;            // Prob. of response i,j for sample s
    double tmp = 0.0;             // Mean prob. of response i,j across samples
    int n = data.nrow();          // Number of respondents
    int m = data.ncol();          // Number of items
    int S = samples.nrow();       // Number of posterior samples
    int T = Rcpp::max(time_map);  // Number of time periods
    double M = std::log(1.0 / S); // Weight for prob. from each sample
    double progtot = (double)n*m; // Variables to help track progress
    int progcurrent = 0;
    for ( int j = 0; j < m; ++j ) {
        int t = time_map[j] - 1;
        for ( int i = 0; i < n; ++i ) {
            progcurrent += 1;
            Rprintf("\r%4.2f", progcurrent / progtot);
            Rcpp::checkUserInterrupt();
            if ( !Rcpp::IntegerVector::is_na(data(i, j)) ) {
                tmp = 0.0;
                for ( int s = 0; s < S; ++s ) {
                    // MCMCirt1d stores parameters like
                    // theta11, ..., theta1T, ..., thetaN1, ..., thetaNT,
                    // alpha1, ..., alphaM, beta1, ..., betaM
                    double mu = samples(s, i * T + t) *
                        samples(s, n * T + m + j) -
                        samples(s, n * T + j);
                    prob = R::pnorm(mu, 0.0, 1.0, data(i, j), 1);
                    tmp += std::exp(M + prob);
                }
            } else {
                tmp = 1.0;
            }
            res += std::log(tmp);
        }
    }
    Rprintf("\r1.00\n");
    // Return result
    return res;
}
