#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>
#include "Chain.hpp"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main(int argc, char* argv[])  {
  // If running from terminal - specify input file name ranges (posterior checks)
  for (int i=0;i<argc;++i)  { cout << i << " " << argv[i] << endl;}
  const int df1 = (argc==3) ? atoi(argv[1]) : 0;
  const int df2 = (argc==3) ? atoi(argv[2]) : 1;

  const int mflag = 2;        // 2:sep inoc/contact. 1:single parameter set. 0:SIR model
  const int nchain = 4;       // num parallel chains
  const int n_samp = 1e8;       // samples
  const int n_burn = 1e8;       // burn in
  const int n_thin = 1e4;       // thinning
  const int lflag = 1;        // log to screen
  for (int dflag=df1;dflag<df2;++dflag)  {        // id of experiment to fit
    #pragma omp parallel num_threads(nchain)
    {
      int cno = omp_get_thread_num();
      stringstream logstream;

      /* Init chain */
      Chain chain(n_burn,n_samp,n_thin,cno,mflag,dflag);
      if(!chain.ready)  {
        logstream << "Could not init: " << cno << "\n";
        cout << logstream.str() << flush;
      }
      else  {
        #pragma omp single
        {
          for (int rms=0;rms<chain.model.n_r;++rms)  {
            logstream << "Room" << rms << ": P" << chain.model.bflags[rms] << " C" << chain.model.cflags[rms] << "\n";
            cout << logstream.str() << flush;
            logstream.str(string());
          }
        }
        logstream << "  Run: " << dflag << ". Chain: " <<  cno << ".\n";
        cout <<  logstream.str() << flush;
        chain.run();
      }
    }
  }
  return(0);
}








