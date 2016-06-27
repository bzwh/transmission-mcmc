#ifndef CHAIN_HPP_INCLUDED
#define CHAIN_HPP_INCLUDED

#include <vector>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
#include "Model.hpp"

class Chain  {
public:
  Chain(int,int,int,int,int,int);
  void run();
  void gen(int);
  int acceptreject();
  void scalefactor(int,double);
  void update_mucov(int);

  Model model;              /// Black box, give it par vector and get back logprior and loglik
  int n_accd;               /// Counter for acceptance
  int n_burn;               /// Length of burn-in period
  int n_samp;               /// Number of post-burn-in samples to take
  int n_thin;               /// Thinning
  gsl_rng* r;               /// PRNG
  int npar;                 /// Number of parameters being fit
  int ready;
  double arate;             /// Acceptance rate
  double sd;                /// Scaling factor
  Eigen::VectorXd par_old;  /// Last accepted set of parameters
  Eigen::VectorXd par_new;  /// Current proposed set of parameters
  Eigen::VectorXd mu_old;   /// Mean parameter values at t-1
  Eigen::VectorXd mu_new;   /// Current means
  Eigen::MatrixXd cov;      /// Variance-covariance matrix of all parameters (empirical)
  Eigen::MatrixXd cov_prop; /// Adapted variance-covariance matrix (rescaled by sd)
  double lhood_old;         /// Old log-likelihood
  double prior_old;         /// Old log-prior
  double lhood_new;         /// Current log-likelihood
  double prior_new;         /// Current log-prior
  int ncov;
  std::ofstream out_cov;                  /// Dump cov matrix
  std::ofstream logdmp;
private:

};

#endif // CHAIN_HPP_INCLUDED

