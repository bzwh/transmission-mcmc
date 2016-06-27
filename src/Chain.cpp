#include <iostream>
#include <cmath>
#include <ctime>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>
#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "Chain.hpp"
#include "Model.hpp"
#include "Mvngen.hpp"

#define COV_DUMP 1

using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;
using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::isfinite;
using std::max;
using std::min;
using std::to_string;
using std::cin;


Chain::Chain(int nb,int ns,int nt,int cno,int mf,int df)  {
  // RNG
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,(cno+1)*time(0));

  // Chain stuff

  lhood_old = std::numeric_limits<double>::infinity();
  prior_old = std::numeric_limits<double>::infinity();
  lhood_new = std::numeric_limits<double>::infinity();
  prior_new = std::numeric_limits<double>::infinity();
  n_accd = 0;
  n_burn = nb;
  n_samp = ns;
  n_thin = nt;
  ncov = 0;
  ready = 0;
  // Model setup - Set paths+filenames, load data & priors, initialise first sample
  model.setup(cno,mf,df);
  if (COV_DUMP)  {
    out_cov.open(model.opath+"cov_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt");
  }
  par_old = model.init_samp(r,prior_old,lhood_old,ready);

  npar = model.npar;
  sd = 2.38*2.38/double(npar);
  par_new = VectorXd::Zero(npar);
  mu_old = par_old;
  mu_new = VectorXd::Zero(npar);
  cov = 1.0e-3*MatrixXd::Identity(npar,npar);
  cov_prop = 1.0e-3*MatrixXd::Identity(npar,npar);

  logdmp.open(model.opath+"logdmp_"+to_string((long long)df)+"_"+ to_string((long long)cno)+".txt");
}


/** \brief
 * \return void
 */
void Chain::run()  {
  int cno = omp_get_thread_num();
  stringstream logstream;
  double arate = 0.0;
  double arate2 = 0.0;
  for (int i=0;i<n_burn+n_samp;++i)  {
    logstream.str(string());
    /* SAMPLE */
    gen(i);
    model.parse(par_new);
    /* TEST */
    int accpd = acceptreject();   // gets model to calc prior and likelihood
    /* UPDATE */
    if (accpd)  {
      ++n_accd;
      update_mucov(i);
      /* update chain */
      par_old = par_new;
      prior_old = prior_new;
      lhood_old = lhood_new;
    }
    arate = n_accd/double(i+1.0);
    arate2 = (i<n_burn) ? 0.0 : ((i-n_burn)*arate2+accpd)/double(i-n_burn+1.0);
    scalefactor(i,arate);
    /* OUTPUT */
    if (i%n_thin==0)  {
      //model.write(par_old,prior_old,lhood_old);
      //if (i<n_burn) model.writeburn(par_old);
      (i<n_burn) ? model.writeburn(par_old)                   // Still warming up
                 : model.write(par_old,prior_old,lhood_old);  // Real thing
    }
    if (i==n_burn)  { // finished adapting
      out_cov << cov/sd << "\n\n" << flush;
      model.out_brn.close();
    }
    logstream << "Chain: " << cno << "-" << to_string((long long) (i+1)) << "\t"
              << sd << "\t" << arate*100.0 << "%\t " << arate2*100 << "%\t "
              << accpd << " " << prior_new << "\t" << lhood_new << "\t";
    if (i%1000==0)  {
      logdmp << logstream.str() << "\n";
    }
  }

  cout << logstream.str() << endl;
  cout << model.opath << endl;
  logdmp.close();
  model.closefiles();
  gsl_rng_free(r);
}



/** \brief Generate new set of parameters.
 *  Applying scale factor to empirical cov matrix
 *  Generate new proposal using Eigen and GSL
 *  Extracting infection times, latent periods and model parameters handled by Model
 * \return void
 */
void Chain::gen(int n)  {  // FIXME separate the non-infecteds too
  if ((n<=n_burn)&&(n<=2*npar||n_accd==0))  { // FIXME conditions....
    // still burning, not many accepted
    cov_prop = 1.0e-3*MatrixXd::Identity(npar,npar);
  }
  else  {
    // Now using empirical cov matrix
    cov_prop = sd*(cov+0.001*MatrixXd::Identity(npar,npar));
  }
  // Generation - gsl providing ugaussian samples, eigen making multinormal
  par_new = par_old+mgen(cov_prop,r);   // multivariate normal proposal
}


/** \brief Test acceptance of current sample - calc prior and lhood
 * \return 1 for accepted, 0 for rejected
 */
int Chain::acceptreject()  {
  // Acceptance flag
  int to_acc = 0;
  double lp = 0.0;
  // Calc prior
  prior_new = model.prior_calc();
  lp += prior_new;
  // Don't bother calculating likelihood if prior already rejects
  if (isfinite(lp))  {
    lhood_new = model.lhood_calc();
    lp += lhood_new;
  }
  else  {
    return(0);
  }

  if (isfinite(lp))  {
    double rr = gsl_rng_uniform(r);
    to_acc = (rr<exp(lp-(lhood_old+prior_old))) ? 1 : 0;
  }
  return(to_acc);
}


/** \brief Iteratively update mean and variance-covariance matrix for next sample. GREEDY
 * \param n int number of iterations
 * \return void
 */
void Chain::update_mucov(int n)  {
  if (n_accd==1)  { // First accepted set, init mean and cov for iterative updates later
  //if (n==1)  { // First accepted set, init mean and cov for iterative updates later
    mu_new = (par_old+par_new)*0.5;
    cov = ((par_old-mu_new)*((par_old+mu_new).transpose())
             +(par_new-mu_new)*((par_new+mu_new).transpose()));
  }
  else if ((n<=n_burn)&&(n_accd>1))  {
    mu_old = mu_new;
    // Update mean
    mu_new = (n_accd*mu_old + par_new)/double(n_accd+1.0);
    // Update covariance
    cov = ((n_accd-1.0)*cov
        + ((n_accd*mu_old*mu_old.transpose()-(n_accd+1.0)*mu_new*mu_new.transpose())
        + par_new*par_new.transpose()))
        / double(n_accd);
  }
  else {
    // finished burning
  }
}


/** \brief Adjust scale factor (external control to call every 100th step during burn-in)
 *  Increase/decrease to keep accpetance rate above 0.2/below 0.4
 * \param i int step number
 * \param arate double current acceptance ratio
 * \return void
 */
void Chain::scalefactor(int i,double arate)  {
  if((i<n_burn)&&(i%100==0))  {
    if (arate<0.2)  {
      sd = max(0.0000002,sd*0.5);      // asfv-pigs ~ 0.02
    }
    else if (arate>0.4)  {
      sd = min(100.0,sd*2.0);
    }
  }
  else  {
    // Finished burning-in
  }
}
