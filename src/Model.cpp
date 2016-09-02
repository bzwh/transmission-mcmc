#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "Model.hpp"

#define TEST_FLAG 1
#define BETAHACK 1

const double EPSILON = std::numeric_limits<double>::epsilon();

using namespace Eigen;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::isfinite;
using std::max;
using std::min;
using std::to_string;
using std::cin;

#ifdef _WIN32
using std::isinf;
#endif // _WIN32


Model::Model()  {
  bflag = 1;  // max number of pens over all rooms == number of beta parameters
  n_a = 0;
  n_c = 0;
  n_r = 0;
  n_l = 0;
  n_ii = 0;
  n_ti = 0;
  orsel_delay = 1.0;
}


void Model::setup(int cno,int mf,int df)  {
  mflag = mf*2;
  setfns(cno,df);
  load_data();
  load_priors();
  // TODO Probably shouldn't initialise storage - move to setup. Keep purely for chain.
  // Nuisance parameter storage
  tInf = VectorXd::Zero(n_a);     // Infection times
  tinfS = VectorXd::Zero(n_ti);
  lat_P = VectorXd::Zero(n_a);    // Latent period duration
  lat_pS = VectorXd::Zero(n_ii);
  inf_p = VectorXd::Zero(n_ii);   // Infectious period duration
  // Model parameters storage
  plat = VectorXd::Zero(mflag);   // Gamma latent
  pinf = VectorXd::Zero(2);       // Gamma infectious
  beta = VectorXd::Zero(bflag);   // Transmission

  parsamp = VectorXd::Zero(npar);
}


void Model::setfns(int cno,int df)  {
  if (df<10)  {
    const int species_flag = 3;
    switch (species_flag)  {
      case 0:
        dname = "./input/fmd-sheep/ChallengeData_sheep.txt";
        pname = "./input/fmd-sheep/priors_fmdsheep.txt";
        opath = "./outputs/fmd_sheep/";
        break;

      case 1:
        //dname = "./input/fmd-pigs/ChallengeDataPigs-exp1.txt";
        //dname = "./input/fmd-pigs/exp2-challenge_only.txt";
        dname = "./input/fmd-pigs/combined.txt";
        //dname = "./input/fmd-pigs/ChallengeDataPigs-combined.txt";
        pname = "./input/fmd-pigs/priors_fmdv_pigs.txt";
        opath = "./outputs/fmd_pigs/";
        break;

      case 2:
        dname = "./input/asf-ppc/ASFVChallengeData.txt";
        pname = "./input/asf-ppc/priors_asf.txt";
        opath = "./outputs/asf_pigs/";
        orsel_delay = 0.0;
        break;

      case 3:
        dname = "./input/vacc-pigs/ChallengeDataPigsVacc.txt";
        pname = "./input/fmd-pigs/priors_fmdv_pigs.txt";
        opath = "./outputs/vacc_pigs/";
        break;

      case 4:
        dname = "./input/eble-pigs/eble06-pigs.txt";
        pname = "./input/eble-pigs/priors_fmdv_pigs.txt";
        opath = "./outputs/eble_pigs/";
        break;

      default:
        cout << "Wtf are you running?" << endl;
        exit(-1);
        break;
    }
  }
  else  {
    if(1)  {
      dname = "./input/synth-fmd/synth_"+to_string(df)+".txt";
      pname = "./input/synth-fmd/priors.txt";
      opath = "./outputs/synth-fmd/";
    }
    else  {
      dname = "./input/synth-asf/synth_"+to_string(df)+".txt";
      pname = "./input/synth-asf/priors.txt";
      opath = "./outputs/synth-asf/";
      orsel_delay = 0.0;
    }
  }

  if (df!=-1)  {
    out_par.open(opath+"par_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt");  // parameters
    out_lat.open(opath+"latp_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt");  // latent periods
    out_inf.open(opath+"tinf_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt"); // infection times
    //out_ipd.open(opath+"infp_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt"); // infectious periods
    out_lhd.open(opath+"lhd_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt");  // likelihoods
    out_brn.open(opath+"brn_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt");  // likelihoods
    if (!out_par.is_open() || !out_lat.is_open() || !out_inf.is_open())  {
      cout << "Output file error. Check paths" << endl;
      exit(-1);
    }
  }

}

/** \brief Just reading experiment data
 * \return void
 */
void Model::load_data()  {
  ifstream indat(dname.data());
  if (indat.is_open())  {
    indat >> n_a >> n_c >> n_r >> exend; // #Animals, #Contacts, #Rooms, end of experiment
    rooms.resize(n_a,0);  // room number for each animal
    nroom1.resize(n_r,0); // numbers of animals in pen 1 of each room
    nroom2.resize(n_r,0); // numbers of animals in pen 2 of each room
    bflags.resize(n_r,1);
    cflags.resize(n_r,0);
    pflags.resize(n_r,1); // Default to assume there's no inocs, correct in loop
    tMove.resize(n_r,100.0);  // time animals moved from each room - orselpigs
    iTyp.resize(n_a,0);
    // catch tNeg/tPos/tEnd == -1  for non-infections? or just use tNeg
    tNeg.resize(n_a,0);
    tPos.resize(n_a,0);
    tEnd.resize(n_a,0);
    cull.resize(n_a,0);
    for (int i=0;i<n_a;++i)  {
      indat >> rooms[i] >> iTyp[i] >> tNeg[i] >> tPos[i] >> tEnd[i] >> cull[i];
      if (iTyp[i]==3)  {
        bflags[rooms[i]] = 2; // This is a between-pen animal, so need 2 parameters
        bflag = 2;
      }
      else if (iTyp[i]==4)  {
        cflags[rooms[i]] = 1; // This room contains a C2 pig!
      }
    }
  }
  else  {
    cout << "No data" << endl;
    exit(-1);
  }
  rLst.resize(n_r); //

  // Parsing animals. Get room counts, ids. Sort out pos/neg infections
  id_ij.resize(n_a,-1);
  for (int i=0;i<n_a;++i)  {
    int rnumi = rooms[i];
    rLst[rnumi].push_back(i);  //so rLst[room] = {animal ids}
    if (tPos[i]>=0.0)  {  // tPos==-1 for iTyp==0
      id_allinf.push_back(i);
      id_ij[i] = id_allinf.size()-1;
      if (iTyp[i]>1) {
        id_contact.push_back(i);
      }
    }
    if (bflags[rnumi]==2)  {    // Two pens
      (iTyp[i]==3) ? ++nroom2[rnumi] : ++nroom1[rnumi];
    }
    else  {                     // One pen - coopting nroom2 for tmove
      ++nroom1[rnumi];
      if (iTyp[i]==4) ++nroom2[rnumi];
    }
  }
  n_ii = id_allinf.size();  // total number of infected inocs + contacts
  n_l = (mflag) ? n_ii : 0; // number of latent periods to be inferred
  n_ti = id_contact.size(); // number of infection times (ie successful contact transmissions)

  // Extract tMove for Orsel's pigs in relevant rooms. Otw set to 0
  for (int room=0;room<n_r;++room)  {
    if (cflags[room])  {
      for (unsigned int aid=0;aid<rLst[room].size();++aid)  {
        int animalid = rLst[room][aid];
        if ((iTyp[animalid]==2)&&(tNeg[animalid]>=0))  {
          tMove[room] = min(tMove[room],double(tPos[animalid]));
        }
      }
    }
    else  {
      tMove[room] = orsel_delay;
    }
  }
  /* Check for inocs in room and adjust the PIGHACK*/
  for (int aid=0;aid<n_a;++aid)  {
    if (iTyp[aid]<2)  {
      pflags[rooms[aid]] = 0;
    }
  }
  ncount();
  // tI, E, lat (k mu), inf (k mu), beta (B W?)
  npar = n_ti+n_l+mflag+2+bflag;

}


/** \brief Fetch total number of animals in each room [0:n_r-1] on each day [0:exend+1]
 * Stored for each pen separately
 * \return void
 */
void Model::ncount()  {
  deltaN1.resize(n_r,vector<int>(1,0));
  deltaN2.resize(n_r,vector<int>(1,0));
  ntot_pens0 = MatrixXd(n_r,exend+1);
  ntot_pens1 = MatrixXd(n_r,exend+1);
 for(int rm=0;rm<n_r;++rm)  {
    ntot_pens0.row(rm).fill(nroom1[rm]);
    ntot_pens1.row(rm).fill(nroom2[rm]);
    if(cflags[rm]==1)  {
      ntot_pens0.block(rm,0,1,tMove[rm]).array() -= 5;
    }
  }
  for (int t=1;t<exend+1;++t)  { // FIXME FIXME FIXME
    for (int i=0;i<n_a;++i)  {
      if ((cull[i])&&(tEnd[i]<exend))  {  // Culled (not just censored)
        if (t>=tEnd[i])  {
          (iTyp[i]>2) ? ntot_pens1(rooms[i],t) -= 1.0
                      : ntot_pens0(rooms[i],t) -= 1.0;
        }
      }
    }
  }
  deltaN.resize(n_r,vector<int>(1,0));
  for (int rm=0;rm<n_r;++rm)  {
    for (int t=1;t<exend+1;++t)  {
      if (ntot_pens0(rm,t)<ntot_pens0(rm,t-1))  {
        deltaN1[rm].push_back(t);
        deltaN[rm].push_back(t);
      }
      if (ntot_pens1(rm,t)<ntot_pens1(rm,t-1))  {
        deltaN2[rm].push_back(t);
        deltaN[rm].push_back(t);
      }
    }
    deltaN1[rm].push_back(exend);
    deltaN2[rm].push_back(exend);
    deltaN[rm].push_back(exend);
  }

  // FIXME Orsel's pigs not tested on cull date. tEnd is not tCull like it is with Guinat's
  if (TEST_FLAG==2)  {
    cout << endl << ntot_pens0 << endl;
    for (int rm=0;rm<n_r;++rm)  {
      for (int dt=0;dt<deltaN1[rm].size();++dt)  {
        cout << deltaN1[rm][dt] << " ";
      }cout << endl;
    }

    cout << endl << ntot_pens1 << endl;
    for (int rm=0;rm<n_r;++rm)  {
      for (int dt=0;dt<deltaN2[rm].size();++dt)  {
        cout << deltaN2[rm][dt] << " ";
      }cout << endl;
    }

    cout << endl;
    for (int rm=0;rm<n_r;++rm)  {
      for (int dt=0;dt<deltaN[rm].size();++dt)  {
        cout << deltaN[rm][dt] << " ";
      }cout << endl;
    }
  }
}


/** \brief Dump out data structures built from raw data.
 * \return void
 */
void Model::test_data()  {
  for (unsigned int rit=0;rit<rLst.size();++rit)  {
    cout << "Room " << rit << ": ";
    for (auto ait=rLst[rit].begin();ait!=rLst[rit].end();++ait)  {
      cout << *ait << " ";
    } cout << endl;
  }
  cout << "id_allinf: ";
  for (auto iit=id_allinf.begin();iit!=id_allinf.end();++iit)  {
    cout << *iit << " ";
  } cout << endl;
  cout << "id_ij: ";
  for (auto iit=id_ij.begin();iit!=id_ij.end();++iit)  {
    cout << *iit << " ";
  } cout << endl;
  cout << "id_contact: ";
  for (auto iit=id_contact.begin();iit!=id_contact.end();++iit)  {
    cout << *iit << " ";
  } cout << endl;
}


/** \brief Read in priors.
 * Hyperparameters for gamma hyperpriors on gamma latent and infectious periods
 * Values are shapes and scales for distribution of the shape and mean latent/infectious periods.
 * Gamma transsmission parameter.
 * \return void
 */
void Model::load_priors()  {
  prior_latk.resize(2,0.0);
  prior_latm.resize(2,0.0);
  prior_infk.resize(2,0.0);
  prior_infm.resize(2,0.0);
  prior_beta.resize(2,0.0);
  ifstream indat(pname.data());
  if (indat.is_open())  { // Just assuming the file conforms...
    indat >> prior_latk[0] >> prior_latk[1] >> prior_latm[0] >> prior_latm[1]
          >> prior_infk[0] >> prior_infk[1] >> prior_infm[0] >> prior_infm[1]
          >> prior_beta[0] >> prior_beta[1];
  }
  else  {
    cout << "No priors" << endl;
    exit(-1);
  }
}


/** \brief Keep generating sample parameter set until finite likelihood and prior
 * \return int: 0 - everything ok. (-1) failed to init
 */
VectorXd Model::init_samp(gsl_rng* r,double& logprr,double& loglik,int& ready)  {
  int init_counter = 0;
  int give_up = 10000; // Really shouldn't struggle this bad
  // TODO hardcoded... wide enough for all data sets?
  vector<double> shp_l = {0.0,5.0};
  vector<double> mue_l = {0.0,2.0};
  vector<double> shp_i = {0.0,10.0};
  vector<double> mue_i = {0.0,10.0};
  vector<double> bt = {0.0,10.0};


  loglik = std::numeric_limits<double>::infinity();
  logprr = std::numeric_limits<double>::infinity();
  while (isinf(loglik))  {
    logprr = std::numeric_limits<double>::infinity();
    while (isinf(logprr))  { // separate this to gen
      ++init_counter;
      if(init_counter>give_up)  {
        out_par << "did not get finite prior: " << logprr << "\t" << perr_flag << "\t" << loglik << "\t" << lerr_flag << endl;
        MatrixXd mtmp(n_a,2); mtmp << tInf,lat_P;
        cout << mtmp << endl;
        ready = 0;
        return(VectorXd::Zero(npar));
      }
      /* Initial infection times  */
      tInf.fill(0.0);
      tinfS.fill(0.0);
      int it_tinf = 0;
      for (int i=0;i<n_a;++i)  {
        // Didn't get infected(infectious)
        if (tNeg[i]<0.0)  {
          tInf[i] = -1.0;
        }
        else  {
          // Contact infection: start guessing ~ [5,7]?? More dispersed? Dep on data!
          switch (iTyp[i])  {
            case 0:
              tInf[i] = 0.0;
            break;

            case 1:
              // Inoculated at t=0;
              tInf[i] = 0.0;
            break;

            case 2:
              // Contact transmission (pen1)
              tInf[i] = gsl_ran_flat(r,orsel_delay,tPos[i]);  // [SG]
              tinfS[it_tinf++] = tInf[i];
            break;

            case 3:
              // Contact transmission (pen2)
              tInf[i] = gsl_ran_flat(r,orsel_delay,tPos[i]);  // [SG]
              //tInf[i] = orsel_delay+gsl_rng_uniform(r)*(5.0);
              tinfS[it_tinf++] = tInf[i];
            break;

            case 4:
              // Contact transmission (challenge2)
              tInf[i] = gsl_ran_flat(r,tMove[rooms[i]],tPos[i]);  // [SG]
              //tInf[i] = tMove[rooms[i]]+gsl_rng_uniform(r)*(tPos[i]-tMove[rooms[i]])*0.5; // FIXME - C2 pigs?
              tinfS[it_tinf++] = tInf[i];
            break;
          }
        }
      }

      /* Sample initial latent period durations */
      lat_P.fill(0.0);
      lat_pS.fill(0.0);
      int it_lat = 0;
      for (int i=0;i<n_a;++i)  {
        if (id_ij[i]==-1)  {  // Either not infected or a 0Inoc (no latent period to infer)
          lat_P[i] = -1.0;
        }
        else  {
          switch(iTyp[i])  {
            case 0:
              cout << "SHOULD NOT HAPPEN" << endl;
              exit(-1);
              break;
            case 1:
              //lat_P[i] = gsl_rng_uniform(r)*(tPos[i]);
              lat_P[i] = gsl_ran_flat(r,tNeg[i],tPos[i]);
              lat_pS[it_lat++] = lat_P[i];
              break;
            default:
              lat_P[i] = (tInf[i]>tNeg[i]) ? gsl_ran_flat(r,tInf[i],tPos[i])-tInf[i]
                                           : gsl_ran_flat(r,tNeg[i]-tInf[i],tPos[i]-tInf[i]); // FIXME INIT ISSUE
              lat_pS[it_lat++] = lat_P[i];
              break;
          }
        }
      }

      /* Sample initial infectious period durations */
      inf_p.fill(0.0);
      for (int j=0;j<n_ii;++j)  {
        int i = id_allinf[j];
        if (tNeg[i]<0)  { // ignore when no infection?
          cout << "THIS REALLY SHOULDNT HAPPEN - negative dpi" << endl;
          exit(-1);
        }
        inf_p[j] = tEnd[i]-tInf[i]-lat_pS[j];// potentially censored, definitely binned
      }

      /* Sample initial latent period parameters  */
      plat.fill(0.0);
      switch (mflag)  {
        case 4:
          plat[0] = gsl_ran_flat(r,shp_l[0],shp_l[1]);  // Contact shape
          plat[1] = gsl_ran_flat(r,mue_l[0],mue_l[1]);  // Contact mean
          plat[2] = gsl_ran_flat(r,shp_l[0],shp_l[1]);  // Inoc shape
          plat[3] = gsl_ran_flat(r,mue_l[0],mue_l[1]);  // Inoc mean
          break;
        case 2:
          plat[0] = gsl_ran_flat(r,shp_l[0],shp_l[1]);  // Combined shape
          plat[1] = gsl_ran_flat(r,mue_l[0],mue_l[1]);  // Combined mean
          break;
        case 0:
          break;
        default:
          break;
      }

      /* Sample initial infectious period parameters  */
      pinf.fill(0.0);
      pinf[0] = gsl_ran_flat(r,shp_i[0],shp_i[1]);
      pinf[1] = gsl_ran_flat(r,mue_i[0],mue_i[1]);

      /* Sample initial transmission parameters */
      beta.fill(0.0);
      for (int bi=0;bi<bflag;++bi)  {
        beta[bi] = gsl_ran_flat(r,bt[0],bt[1]);
      }

      if (mflag)  {
        parsamp << tinfS , lat_pS , plat , pinf , beta;
      }
      else  {
        parsamp << tinfS , pinf , beta;
      }

      //cout << "Calc prior...";// << flush;
      logprr = prior_calc();
      //cout << prior_old << endl;
    }
    //cout << "Calc lhood...";// << flush;
    //cout << init_counter << endl;
    loglik = lhood_calc();
    //cout << lhood_old << endl;
    //std::cin.get();


  }
  if (TEST_FLAG)  {
    stringstream lstream;
    //lstream <<"INIT " << omp_get_thread_num() << "- samp gen'd in: "<<init_counter<<"\n";
    cout << lstream.str() << flush;
  }
  ready = 1;
  return(parsamp);
}


/** \brief Split big vector of paramters into infection times, latent/infectious periods, etc
 * \param par_vec VectorXd
 * \return void
 */
void Model::parse(VectorXd par_vec)  {
  if (par_vec.size()!=npar)  {
    cout << "how the fuck have you got the wrong number of parameters?" << endl;
    exit(-1);
  }
  parsamp = par_vec;
  /* Extract infection times */
  tinfS = parsamp.head(n_ti);            // infection times of contacts only
  tInf.fill(0.0);                       // Inocs assumed to be infected at t=0
  // id_contact tells us which i \in [0,n_a) to populate
  for (int j=0;j<n_ti;++j)  {
    int i = id_contact[j];
    tInf[i] = tinfS[j];
  }
  /* Extract latent periods */
  if (mflag)  {
    lat_pS = parsamp.segment(n_ti,n_l);     // unobserved latent periods
    lat_P.fill(-1.0);
    for (int j=0;j<n_ii;++j)  {
      int i = id_allinf[j];
      lat_P[i] = lat_pS[j];
    }
  }
  /* Extract infectious periods */
  for (int j=0;j<n_ii;++j)  {
    int i = id_allinf[j];
    inf_p[j] = tEnd[i]-tInf[i]-lat_P[i];
  }
  /* Extract Parameters */
  plat = parsamp.segment(n_ti+n_l,mflag);    // shape and mean latent period inocs and contacts
  pinf = parsamp.segment(n_ti+n_l+mflag,2);  // shape and mean infectious period
  beta = parsamp.tail(bflag);
}

/** \brief Wanna have a gander at the burn-in.
 * Get a better look at convergence from dispersed starting points.
 * \return void
 */
void Model::writeburn(VectorXd pars)  {
  // all the model parameters
  out_brn << pars.transpose() << "\n";
}


/** \brief Calculate prior.
 * Gamma hyperpriors on gamma distributed shapes and means of latent and infectious periods
 * Exponential prior on transmission parameter - not particularly informative, stop going too large
 * Read from file defined in constructor
 * \return double
 */
double Model::prior_calc()  {
  double ptmp = 0.0;
  lerr_flag = 0.0;
  perr_flag = 0.0;
  /* latent period parameters combined or contacts first  */
  if (mflag)  {
    ptmp += log(gsl_ran_gamma_pdf(plat[0],prior_latk[0],prior_latk[1]))
         +  log(gsl_ran_gamma_pdf(plat[1],prior_latm[0],prior_latm[1]));
    if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag = 1;if(TEST_FLAG==2)cout<<"1P "<<flush;return(ptmp); }

    /* latent period parameters for separated inoculated (same priors)  */
    if (mflag==4)  {
      ptmp += log(gsl_ran_gamma_pdf(plat[2],prior_latk[0],prior_latk[1]))
           +  log(gsl_ran_gamma_pdf(plat[3],prior_latm[0],prior_latm[1]));
    }
    if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag = 2;if(TEST_FLAG==2)cout<<"2P "<<flush;return(ptmp); }
  }

  /* 3 infectious period parameters */
  ptmp += log(gsl_ran_gamma_pdf(pinf[0],prior_infk[0],prior_infk[1]))
       +  log(gsl_ran_gamma_pdf(pinf[1],prior_infm[0],prior_infm[1]));
  if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag = 3;if(TEST_FLAG==2)cout<<"3P "<<flush;return(ptmp); }

  /* 4 transmission parameters  */
  for (int bi=0;bi<bflag;++bi)  {
    ptmp += log(gsl_ran_gamma_pdf(beta[bi],prior_beta[0],prior_beta[1]));
  }
  if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag = 4;if(TEST_FLAG==2)cout<<"4P "<<flush;return(ptmp); }

  /* 5 Constrain infection time between 0/orsel_delay/tMove and first positive  */
  for (auto jit=id_contact.begin();jit!=id_contact.end();++jit)  {
    // directly grab indices of contacts rather than looping through iTyp
    int i = *jit;
    if (iTyp[i]==4)  {
      // Orsel C2 pig. Transmission only occurs [tMove,tPos]
      if (mflag)  {
        ptmp += log(gsl_ran_flat_pdf(tInf[i],tMove[rooms[i]],tPos[i]));
      }
      else  {
        ptmp += log(gsl_ran_flat_pdf(tInf[i],tNeg[i],tPos[i]));
      }
    }
    else  {
      // Normal contact (poss within or between pen). Transmission occurs [intro,tPos]
      if (mflag)  {
        ptmp += log(gsl_ran_flat_pdf(tInf[i],orsel_delay,tPos[i]));
      }
      else  { // Latent periods assumed zero... so stronger constraints. Issue for inocs though?
        ptmp += log(gsl_ran_flat_pdf(tInf[i],tNeg[i],tPos[i]));
      }
    }
    if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag=5+i/100.0;if(TEST_FLAG==2)cout<<"5P:"<<i<<" ";return(ptmp); }
  }

  /* 6 Constrain latent period between tNeg-tI and tPos-tI: E gets from tI to tN but before tP  */
  for (int j=0;j<n_l;++j)  {
    int i = id_allinf[j]; // this is animal id for [0:nA]
    if (BETAHACK)  {
      double ee = (tPos[i]-tNeg[i]);
      ptmp += log(gsl_ran_beta_pdf((lat_P[i]+tInf[i]-tNeg[i])/ee,1.01,1.01));
    }
    else  {
      ptmp += log(gsl_ran_flat_pdf(lat_P[i],tNeg[i]-tInf[i],tPos[i]-tInf[i]));
    }
    if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag = 6+j/100.0;if(TEST_FLAG==2)cout<<"\t\t\t\t\t6P-"<<i<<"\t"<<tInf[i] << " " <<lat_pS[j]<<" "<<flush;return(ptmp); }
    // Wee bit of early return... can skip one loop but not so efficient if decent acceptance rates?
    if (isinf(ptmp))  {
      if (TEST_FLAG==2)  {cout << endl;}
      return(ptmp);
    }
  }

  if ((TEST_FLAG==2)&&(isinf(ptmp)))  { if(TEST_FLAG==2)cout<<endl<<tInf<<"!"<<endl<<lat_pS<<endl<<endl; }
  return(ptmp);
}
// TODO test and early return on inf prior


/** \brief Count number of infectious animals in each pen at time contact animal i is infected.
 * Only call this when contact animal i was actually infected
 * \param i int The id of animal being infected
 * \return void
 */
void Model::ni_calc(int i,std::vector<double>& ninf)  {
  // Infected animal i, want number of still infectious in room when i got infected
  int room = rooms[i];
  for (auto jit=rLst[room].begin();jit!=rLst[room].end();++jit)  {
    int j = *jit;   // animal id for each animal in this room
    if (i==j)  {
      continue; // obv don't want to count self. not that it'd actually do anything...
    }
    else  {
      // Check each other animal in this room
      if (tNeg[j]<0)  {
        // j was never infected/infectious, so no contribution
        continue;
      }
      else  {
        // j was infected/infectious at some point
        // Check j was infectious && not yet dead
        if (iTyp[j]==0)  {
          // Assumed infectious source
          if ((tInf[i]<=tMove[room]))  {
            ninf[0] += 1.0;
          }
        }
        else if ((tInf[i]>(tInf[j]+lat_P[j])) && (tInf[i]<tEnd[j]))  { // TODO tEnd[j]+1? t_FN?
          (iTyp[j]==3) ? ninf[1]+=1.0 : ninf[0]+=1.0; // Check which pen...
        }
      }
    }
  }
}


/** \brief Integrating total number of infectious up to time contact animal i is infected
 * For infected i, sums up to t=tInf. Non-infected sums to t=tEnd
 * \param i int
 * \return vector<double>
 */
void Model::ns_calc(int i,std::vector<double>& ninfsum)  {
  // If infected - want to stop at tInf. If suscpetible - stop at end of experiment/cull(?)
  double iStop = (tNeg[i]==-1) ? tEnd[i] : tInf[i];
  // Animal i, want total infectious time of animals in room up to istop
  int room = rooms[i];

  // Count up total infectious time each j has with our animal i
  for (auto jit=rLst[room].begin();jit!=rLst[room].end();++jit)  {
    int j = *jit; // animal id for each in this room
    if ((tNeg[j]<0)||(i==j))  {
      // j was never infected, no contribution. and i can't infect itself...
      continue;
    }
    else  {
      switch (iTyp[j])  {
        case 0: // Inoc with no data - assuming to be infectious for full window [0,tMove]
          if (iTyp[i]==2)  { // C1 pig
            (iStop>tMove[room]) ? ninfsum[0]+=tMove[room]-orsel_delay //exposed to inocs for full time
                                : ninfsum[0]+=iStop-orsel_delay;      //infected before move
          }
          // else an inoc or C2 - neither affected by source
          break;

        case 1: // Normal inoc
          if (iStop>(lat_P[j]+inf_p[id_ij[j]]))  {
            // j recovered before i infected.
            ninfsum[0] += inf_p[id_ij[j]];
          }
          else  {
            ninfsum[0] += max(iStop-lat_P[j],0.0);
          }
          break;

        case 2: // Challenge from pen 1
          //if (iStop>(lat_P[j]+inf_p[id_ij[j]]))  {  // FIXME should have the tInf[j] too?
          if (iStop>tEnd[j])  {
            ninfsum[0] += inf_p[id_ij[j]];
          }
          else  {
            ninfsum[0] += max(iStop-(tInf[j]+lat_P[j]),0.0);
          }
          break;

        case 3: // Challenge from pen 2
          if (iStop>tEnd[j])  {
            ninfsum[1] += inf_p[id_ij[j]];
          }
          else  {
            ninfsum[1] += max(iStop-(tInf[j]+lat_P[j]),0.0);
          }
          break;

        case 4: // C2 challenge
          if (iStop>tEnd[j])  {
            ninfsum[1] += inf_p[id_ij[j]];
          }
          else  {
            ninfsum[0] += max(iStop-(tInf[j]+lat_P[j]),0.0);
          }
          break;
      }
    }
  }
}


/** \brief Calculate probability of infection for each i \in n_a at t=tInf[i]
 *  Loop through all animals. iTyp=={0,1} immediately makes prob_inf==1
 *  Infected contacts call ni_calc and ns_calc
 * \param prob_inf std::vector<double>&
 * \return void
 */
void Model::probinf_calc(std::vector<double>& prob_inf)  {
  // For each animal in this room - calculate probability of infection
  for (int i=0;i<n_a;++i)  {
    if (iTyp[i]>1)  {
      vector<double> n_pen = {double(nroom1[rooms[i]]),double(nroom2[rooms[i]])};
      vector<double> ninf_sum(bflags[rooms[i]],0.0);
      ns_calc(i,ninf_sum);

      vector<double> ninf(bflags[rooms[i]],0.0);
      if (id_ij[i]!=-1)  {
        ni_calc(i,ninf);
      }

      if (bflags[rooms[i]]==1)  {  // No second pen.
        // Already guaranteed to be a contact
        if (tNeg[i]<0)  {
          prob_inf[i] =  exp(-(beta[0]*ninf_sum[0]/n_pen[0]));
        }
        else  {
          prob_inf[i] = (beta[0]*ninf[0]/n_pen[0])
                      *  exp(-(beta[0]*ninf_sum[0]/n_pen[0]));
        }
      }

      else  {// multiple pens (experiment - not individual rooms...)
        switch (iTyp[i])  {
          case 2: // Within-pen contact
            if (tNeg[i]<0)  { // No transmission
              prob_inf[i] =  exp( -(beta[0]*ninf_sum[0]/n_pen[0] + beta[1]*ninf_sum[1]/(n_pen[0]+n_pen[1])) );
            }
            else  {           // Transmission
              prob_inf[i] = (beta[0]*ninf[0]/n_pen[0] + beta[1]*ninf[1]/(n_pen[0]+n_pen[1]))
                          *  exp(-(beta[0]*ninf_sum[0]/n_pen[0]+beta[1]*ninf_sum[1]/(n_pen[0]+n_pen[1])));
            }
            break;

          case 3:   // Between-pen transmission
            if (tNeg[i]<0)  { // Not infected contact
              prob_inf[i] =  exp(-(beta[0]*ninf_sum[1]/n_pen[1]+beta[1]*ninf_sum[0]/(n_pen[0]+n_pen[1])));
            }
            else  { // Infected
              prob_inf[i] = (beta[0]*ninf[1]/n_pen[1]+beta[1]*ninf[0]/(n_pen[0]+n_pen[1]))
                          *  exp(-(beta[0]*ninf_sum[1]/n_pen[1]+beta[1]*ninf_sum[0]/(n_pen[0]+n_pen[1])));
            }
            break;

          default:
            cout << "probinf_calc - have 2 pens but not iTyp 2 or 3" << endl;
            exit(-1);
            break;
        }
      }
    }
    else  {
      // Not contact!
      prob_inf[i] = 1.0;
    }
  }
}



/** \brief Calculate log-likelihood.
 * First calc probs of infection (move outside?)
 * Then contributions from gamma distributed latent and infectious periods.
 * \return double
 */
double Model::lhood_calc()  {
  lerr_flag = 0.0;
  double loglik = 0.0;
  // Get probabilities of infection for all contacts
  vector<double> prob_inf(n_a,0.0);
  //probinf_calc(prob_inf);
  //vector<double> prob_inf2(n_a,0.0);
  probinf_calc2(prob_inf);


  //Contributions from infection times for contact transmissionsn
  for (int i=0;i<n_a;++i)  {  // prob_inf already handles if transmission occurs or not
    if (pflags[rooms[i]])  // FIXME EXP2 HACK - ITYP>2
      loglik += (iTyp[i]>2) ? log(prob_inf[i]) : 0.0;
    else
      loglik += (iTyp[i]>1) ? log(prob_inf[i]) : 0.0;
    if ((TEST_FLAG)&&(isinf(loglik)))  { lerr_flag=1+i/100.0;if(TEST_FLAG==2)cout<<"1L ";return(loglik); }
  }

  // Contributions from latent periods
  for (auto iit=id_allinf.begin();iit!=id_allinf.end();++iit)  {
    int i = *iit;
    // Contribution from latent periods (probably prettier way to do this but fuck it)
    switch (mflag)  {
      case 0:
        // No latent periods inferred... so no contributions
        break;
      case 2:
        // All inoc and contacts together
        loglik += log(gsl_ran_gamma_pdf(lat_P[i],plat[0],plat[1]/plat[0]));
        break;
      case 4:
        // Separated contributions from inocs and contacts
        switch(iTyp[i])  {
          case 0:
            // Censored inoc pig - no contribution
            break;
          case 1:
            // Normal inoc
            loglik += log(gsl_ran_gamma_pdf(lat_P[i],plat[2],plat[3]/plat[2]));
            break;
          default:
            // Contacts (case 2,3,4)
            loglik += log(gsl_ran_gamma_pdf(lat_P[i],plat[0],plat[1]/plat[0]));
            break;
        }
        break;
    }
    if ((TEST_FLAG)&&(isinf(loglik)))  {
      lerr_flag=2+i/100.0;
      if (lat_P[i]<0.0) lerr_flag+=0.001;
      if(TEST_FLAG==2)cout<<"2L ";return(loglik);}
  }

  // Contributions from infectious periods - potentially censored & definitely binned
  /* Extract infectious periods */
  for (int j=0;j<n_ii;++j)  {
    int i = id_allinf[j];
    switch (cull[i])  {
      case 0:
        loglik += log( gsl_cdf_gamma_P(inf_p[j]+1.0,pinf[0],pinf[1]/pinf[0])
                      -gsl_cdf_gamma_P(inf_p[j],pinf[0],pinf[1]/pinf[0]));
        break;
      case 1:
        loglik += log(gsl_cdf_gamma_Q(inf_p[j],pinf[0],pinf[1]/pinf[0]));
        break;
      default:
        cout << "LHOOD - cull status not zero or one WTF?" << endl;
        exit(-1);
        break;
    }
    if ((TEST_FLAG)&&(isinf(loglik)))  { lerr_flag=3+i/100.0;if(TEST_FLAG==2)cout<<"3L " << endl;return(loglik); }
  }
  // Finished
  return(loglik);
}


/** \brief Writing current chain state to relevant files.
 *  Infection times and latent period durations in separate files.
 *  Model parameters all together in own file
 * \return void
 */
void Model::write(VectorXd pars,double logprr,double loglik)  {
  // infection times for the (3) contacts
  out_inf << pars.head(n_ti).transpose() << "\n";
  // latent periods for all infecteds (inocs and contacts)
  out_lat << pars.segment(n_ti,n_l).transpose() << "\n";
  // all the model parameters
  out_par << pars.tail(mflag+2+bflag).transpose() << "\n";
  // dump prior and likelihoods
  // FIXME perr_flag does not correspond to the dumped prior and likelihood yet
  out_lhd << logprr<< "\t" << loglik<< "\t" << perr_flag << "\t" << lerr_flag << "\n";
  // infectious periods, well, sort of. Not inferred, just tEnd-latp-tInf
  //out_ipd << inf_p.transpose() << "\n";
}


/** \brief Close ofstreams.
 * \return void
 */
void Model::closefiles()  {
  out_par.close();
  out_lat.close();
  out_inf.close();
  //out_ipd.close();
  out_lhd.close();
  out_brn.close();
}




/** \brief Integrating I(t)/N(t) for t=[0,tInf ? tEnd]
 * For infected i, integrates up to t=tInf. Non-infected sums to t=tEnd
 * Discretised by wheverever N changes. But only caught in each room...
 * \param i int
 * \return vector<double>
 */
void Model::ns_calc2(int i,std::vector<double>& ninfsum)  {
  double iS = (tNeg[i]<0) ? tEnd[i] : tInf[i];
  int room = rooms[i];
  for (int asd=0;asd<deltaN[room].size()-1;++asd)  {
    double tau1 = deltaN[room][asd];
    double tau2 = deltaN[room][asd+1];
    double iStop = min(tau2,iS);//(tNeg[i]==-1) ? tEnd[i] : tInf[i];
    for (auto jit=rLst[room].begin();jit!=rLst[room].end();++jit)  {
      int j = *jit; // animal id for each in this room
      if ((tNeg[j]<0)||(i==j))  {
        // j was never infected, no contribution. and i can't infect itself...
        continue;
      }
      else  {
        double ijStop = min(iStop,double(tEnd[j]));
        double iStart = max(tau1,tInf[j]+lat_P[j]);

        if (bflags[room]==1)  {
          ninfsum[0] += max(ijStop-iStart,0.0)/double(ntot_pens0(room,floor(tau1)));
        }
        else  {
          if (iTyp[i]==3)  {
            /* Second pen being infected by... */
            if (iTyp[j]==3)  {  // Same pen
              ninfsum[0] += max(ijStop-iStart,0.0)/double(ntot_pens1(room,floor(tau1)));
            }
            else  {             // Other pen
              ninfsum[1] += max(ijStop-iStart,0.0)/double(ntot_pens0(room,floor(tau1))+ntot_pens1(room,floor(tau1)));
            }
          }
          else  {
            /* First pen (inocs and C1) being infected by... */
            if (iTyp[j]==3)  {  // Other pen
              ninfsum[1] += max(ijStop-iStart,0.0)/double(ntot_pens0(room,floor(tau1)));
            }
            else  {             // Same pen
              ninfsum[0] += max(ijStop-iStart,0.0)/double(ntot_pens0(room,floor(tau1)));
            }
          }
        }
      }
    }
    if (tau2>iS)
      break;
  }
}

void Model::probinf_calc2(std::vector<double>& prob_inf)  {
  // For each animal in this room - calculate probability of infection
  for (int i=0;i<n_a;++i)  {
    if (iTyp[i]>1)  {
      vector<double> n_pen = {double(nroom1[rooms[i]]),double(nroom2[rooms[i]])};
      vector<double> ninf_sum(bflags[rooms[i]],0.0);
      ns_calc2(i,ninf_sum);
      vector<double> ninf(bflags[rooms[i]],0.0);

      if (id_ij[i]!=-1)  {
        ni_calc(i,ninf);
      }
      if (bflags[rooms[i]]==1)  {  // No second pen.
        // Already guaranteed to be a contact
        if (tNeg[i]<0)  {
          prob_inf[i] =  exp(-(beta[0]*ninf_sum[0]));
        }
        else  {
          prob_inf[i] = (beta[0]*ninf[0]/ntot_pens0(rooms[i],floor(tInf[i])))
                      *  exp(-(beta[0]*ninf_sum[0]));
        }
      }
      else  {
        switch(iTyp[i])  {
          case 2:
            prob_inf[i] = (beta[0]*ninf[0]/ntot_pens0(rooms[i],floor(tInf[i]))
                          +beta[1]*ninf[1]/(ntot_pens0(rooms[i],floor(tInf[i]))+ntot_pens1(rooms[i],floor(tInf[i]))))
                        *  exp(-(beta[0]*ninf_sum[0]+beta[1]*ninf_sum[1]));
            break;
          case 3:
            prob_inf[i] = (beta[0]*ninf[1]/ntot_pens1(rooms[i],floor(tInf[i]))
                          +beta[1]*ninf[0]/(ntot_pens0(rooms[i],floor(tInf[i]))+ntot_pens1(rooms[i],floor(tInf[i]))))
                        *  exp(-(beta[0]*ninf_sum[0]+beta[1]*ninf_sum[1]));
            break;
          default:
            cout << "asdfg" << endl;
            exit(-1);
            break;
        }
      }
    }
    else  {
      // Not contact!
      prob_inf[i] = 1.0;
    }
  }
}
