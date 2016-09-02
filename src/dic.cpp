#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <Eigen/dense>

#include "Model.hpp"

//#include "mex.h"

using namespace std;
using namespace Eigen;

double dic(int,int,int);

int main(void)  {
  int mflag = 1;
  ofstream loglik("./outputs/llik.txt");
  double llk = dic(0,mflag,-1);
  loglik << llk << endl;
  loglik.close();
  cout << llk << endl;
  return(0);
}



double dic(int cno,int mf,int df)  {
  Model model;
  model.setup(cno,mf,df);

  vector<double> v;
  ifstream pfile("./input/dic.txt");
  double x = 0.0;
  if (pfile.is_open())  {
    while(true)  {
      pfile >> x;
      if (pfile.eof()) break;
      v.push_back(x);
    }
  }
  else {
    cout << "FILE!" << endl;
    exit(-1);
  }

  cout << v.size() << endl;

  if (v.size()!=model.npar)  {
    cout << "which experiment?!\t" << v.size() << " " << model.npar << endl;
    exit(-1);
  }

  Map<VectorXd> pars(v.data(),model.npar);
  model.parse(pars);
  double logpri = 0.0;
  double loglik = 0.0;
  logpri = model.prior_calc();
  //if (MexIsFinite(logpri))  {
  if (isfinite(logpri))  {
    loglik = model.lhood_calc();
  }

  //cout << logpri << " " << loglik << endl;

  return(loglik);
}




/*
mxArray* getMexArray(const std::vector<double>& v){
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(),mxREAL);
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    return mx;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])  {
  if (nrhs!=1)  {
    mexErrMsgTxt("Three inputs required");
  }
  int cno = 0;
  int dflag = -1;
  int mflag;
  mflag = floor(mxGetScalar(prhs[1]));

  vector<double> llik(1,0.0);
  llik[0] = dic(cno,mflag,dflag);
  plhs[0] = getMexArray(llik);
}
*/






