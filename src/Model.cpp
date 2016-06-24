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

// Set to day when inocs introduced to susc
#define ORSELHACK 1.0
#define TEST_FLAG 1

const double EPSILON = std::numeric_limits<double>::epsilon();

using namespace Eigen;
using std::ifstream;
using std::ofstream;
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
}

void Model::setup(int cno,int mf,int df)  {
  mflag = mf*2;
  setfns(cno,df);
  load_data();
  load_priors();
}


void Model::setfns(int cno,int df)  {
  const int species_flag = 2;
  switch (species_flag)  {
    case 0:
      dname = "./input/fmdv-sheep/ChallengeData_sheep.txt";
      pname = "./input/fmdv-sheep/priors_fmdsheep.txt";
      opath = "./outputs/fmd_sheep1/";
      break;

    case 1:
      dname = "./input/ChallengeDataPigs-exp1.txt";
      pname = "./input/priors_fmdv_pigs.txt";
      opath = "./outputs/fmd_pigs/";
      break;

    case 2:
      dname = "./input/asfv-ppc/ASFVChallengeData.txt";
      pname = "./input/asfv-ppc/priors_asf.txt";
      opath = "./outputs/asf_pigs/";
      break;

    default:
      cout << "Wtf are you running?" << endl;
      exit(-1);
      break;
  }

  out_par.open(opath+"par_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt");  // parameters
  out_lat.open(opath+"latp_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt");  // latent periods
  out_inf.open(opath+"tinf_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt"); // infection times
  out_ipd.open(opath+"infp_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt"); // infectious periods
  out_lhd.open(opath+"lhd_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt");  // likelihoods
  out_brn.open(opath+"brn_"+to_string((long long)df) + "_" + to_string((long long)cno)+".txt");  // likelihoods
  if (!out_par.is_open() || !out_lat.is_open() || !out_inf.is_open() || !out_ipd.is_open())  {
    cout << "Output file error. Check paths" << endl;
    exit(-1);
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
    (iTyp[i]==3) ? ++nroom2[rnumi] : ++nroom1[rnumi];
  }
  n_ii = id_allinf.size();  // total number of infected inocs + contacts
  n_l = (mflag) ? n_ii : 0; // number of latent periods to be inferred
  n_ti = id_contact.size(); // number of infection times (ie successful contact transmissions)

  // Extract tMove for Orsel's pigs in relevant rooms. Otw set to 0
  for (int room=0;room<n_r;++room)  {
    if (cflags[room])  {
      for (unsigned int aid=0;aid<rLst[room].size();++aid)  {
        int animalid = rLst[room][aid];
        if (iTyp[animalid]==2)  {
          tMove[room] = min(tMove[room],double(tPos[animalid]));
        }
      }
    }
    else  {
      tMove[room] = ORSELHACK;
    }
  }
  // tI, E, lat (k mu), inf (k mu), beta (B W?)
  npar = n_ti+n_l+mflag+2+bflag;

  /* number of infection times, latent periods, also latent periods
  cout << "Animals: " << n_a << "\nContacts: " << n_c << "\nRooms: " << n_r << endl;
  cout << n_ti << " " << n_ii << " " << id_allinf.size() << endl;*/
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
  int give_up = 1000000; // Really shouldn't struggle this bad
  // FIXME hardcoded... wide enough for all data sets?
  // Uniform bounds on initial infection times
  double timin = ORSELHACK;
  /*vector<double> shp_l = {0.0,25.0};
  vector<double> mue_l = {0.0,10.0};
  vector<double> shp_i = {0.0,25.0};
  vector<double> mue_i = {0.0,15.0};
  vector<double> bt = {0.0,25.0};      // uniform bounds on transmission parameter parameters
  */
  vector<double> shp_l = {0.0,5.0};
  vector<double> mue_l = {0.0,10.0};
  vector<double> shp_i = {0.0,25.0};
  vector<double> mue_i = {0.0,15.0};
  vector<double> bt = {0.0,5.0};
  /*vector<double> shp_l = {0.0,2.5};
  vector<double> mue_l = {0.0,2.5};
  vector<double> shp_i = {0.0,5.0};
  vector<double> mue_i = {0.0,10.0};
  vector<double> bt = {0.0,10.0};      // uniform bounds on transmission parameter parameters. shouldnt be the same...
*/
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
  loglik = std::numeric_limits<double>::infinity();
  logprr = std::numeric_limits<double>::infinity();
  while (isinf(loglik))  {
    logprr = std::numeric_limits<double>::infinity();
    while (isinf(logprr))  { // separate this to gen
      // FIXME Rewrite initial sample to be gen'd directly in to par_old and parse()
      ++init_counter;
      // Initial infection times
      // FIXME Shouldn't initialise storage - move to setup. Keep purely for chain.
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
              tInf[i] = timin+gsl_rng_uniform(r)*(tPos[i]-timin);
              tinfS[it_tinf++] = tInf[i];
            break;

            case 3:
              // Contact transmission (pen2)
              tInf[i] = timin+gsl_rng_uniform(r)*(tPos[i]-timin);
              tinfS[it_tinf++] = tInf[i];
            break;

            case 4:
              // Contact transmission (challenge2)
              tInf[i] = tMove[rooms[i]]+gsl_rng_uniform(r)*(tPos[i]-tMove[rooms[i]]); // FIXME - C2 pigs?
              tinfS[it_tinf++] = tInf[i];
            break;
          }
        }
      }

      // Sample initial latent period durations
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
              lat_P[i] = tNeg[i]+gsl_rng_uniform(r)*(tPos[i]-tNeg[i]);
              lat_pS[it_lat++] = lat_P[i];
              break;
            case 2:
              lat_P[i] = (tInf[i]>tNeg[i]) ? gsl_rng_uniform(r)*(tPos[i]-tInf[i])
                                           : tNeg[i]-tInf[i] + gsl_rng_uniform(r)*(tPos[i]-tNeg[i]);
              lat_pS[it_lat++] = lat_P[i];
              break;
            case 3:
              lat_P[i] = (tInf[i]>tNeg[i]) ? gsl_rng_uniform(r)*(tPos[i]-tInf[i])
                                           : tNeg[i]-tInf[i] + gsl_rng_uniform(r)*(tPos[i]-tNeg[i]);
              lat_pS[it_lat++] = lat_P[i];
              break;
            case 4:
              lat_P[i] = (tInf[i]>tNeg[i]) ? gsl_rng_uniform(r)*(tPos[i]-tInf[i])
                                           : tNeg[i]-tInf[i] + gsl_rng_uniform(r)*(tPos[i]-tNeg[i]);
              lat_pS[it_lat++] = lat_P[i];
              break;
            default :
              break;
          }
        }
      }

      // Sample initial infectious period durations
      inf_p.fill(0.0);
      for (int j=0;j<n_ii;++j)  {
        int i = id_allinf[j];
        if (tNeg[i]<0)  { // ignore when no infection?
          cout << "THIS REALLY SHOULDNT HAPPEN - negative dpi" << endl;
          exit(-1);
        }
        inf_p[j] = tEnd[i]-tInf[i]-lat_pS[j];// potentially censored, definitely binned
      }

      // Sample initial latent period parameters
      plat.fill(0.0);
      switch (mflag)  {
        case 4:
          plat[0] = shp_l[0]+gsl_rng_uniform(r)*(shp_l[1]-shp_l[0]);
          plat[1] = mue_l[0]+gsl_rng_uniform(r)*(mue_l[1]-mue_l[0]);
          plat[2] = shp_l[0]+gsl_rng_uniform(r)*(shp_l[1]-shp_l[0]);
          plat[3] = mue_l[0]+gsl_rng_uniform(r)*(mue_l[1]-mue_l[0]);
          break;
        case 2:
          plat[0] = shp_l[0]+gsl_rng_uniform(r)*(shp_l[1]-shp_l[0]);
          plat[1] = mue_l[0]+gsl_rng_uniform(r)*(mue_l[1]-mue_l[0]);
          break;
        case 0:
          break;
        default:
          break;
      }

      // Sample initial infectious period parameters
      pinf.fill(0.0);
      pinf[0] = gsl_rng_uniform(r)*(shp_i[1]-shp_i[0])+shp_i[0];
      pinf[1] = gsl_rng_uniform(r)*(mue_i[1]-mue_i[0])+mue_i[0];

      // Sample initial transmission parameters
      beta.fill(0.0);
      for (int bi=0;bi<bflag;++bi)  {
        beta[bi] = bt[0] + gsl_rng_uniform(r)*(bt[1]-bt[0]);
      }

/*      cout << tinfS.size() << " " << lat_pS.size() << " "
           << plat.size()  << " " << pinf.size()  << "\n"
           << par_old.size() << endl;*/
      if (mflag)  {
        parsamp << tinfS , lat_pS , plat , pinf , beta;
      }
      else  {
        parsamp << tinfS , pinf , beta;
      }

      //cout << "Calc prior...";// << flush;
      logprr = prior_calc();
      //cout << prior_old << endl;
      if(init_counter>give_up)  {
        out_par << "did not get finite prior: "+opath << endl;
        ready = 0;
        return(VectorXd::Zero(npar));
      }
    }
    //cout << "Calc lhood...";// << flush;
    loglik = lhood_calc();
    //cout << lhood_old << endl;
    //std::cin.get();
  }
  if (TEST_FLAG)  { cout <<"INIT " << omp_get_thread_num() << "- samp gen'd in: "<<init_counter<<"\n"; }
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
  // Extract infection times
  tinfS = parsamp.head(n_ti);            // infection times of contacts only
  tInf.fill(0.0);                       // Inocs assumed to be infected at t=0
  // id_contact tells us which i \in [0,n_a) to populate
  for (int j=0;j<n_ti;++j)  {
    int i = id_contact[j];
    tInf[i] = tinfS[j];
  }

  // Extract latent periods
  if (mflag)  {
    lat_pS = parsamp.segment(n_ti,n_l);     // unobserved latent periods
    lat_P.fill(-1.0);
    for (int j=0;j<n_ii;++j)  {
      int i = id_allinf[j];//FIXME id_allinf or id_ij?
      lat_P[i] = lat_pS[j];
    }
  }

  // Extract infectious periods
  for (int j=0;j<n_ii;++j)  {
    int i = id_allinf[j];
    inf_p[j] = tEnd[i]-tInf[i]-lat_pS[j];
  }

  // Extract Parameters
  plat = parsamp.segment(n_ti+n_l,mflag);    // shape and mean latent period inocs and contacts
  pinf = parsamp.segment(n_ti+n_l+mflag,2);  // shape and mean infectious period
  beta = parsamp.tail(bflag);
  //out_acc << prior_old << " " << lhood_old << " " << par_old.transpose() << endl;
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
  perr_flag=0;
  // latent period parameters combined or contacts first
  if (mflag)  {
    ptmp += log(gsl_ran_gamma_pdf(plat[0],prior_latk[0],prior_latk[1]))
         +  log(gsl_ran_gamma_pdf(plat[1],prior_latm[0],prior_latm[1]));
    if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag = 1;if(TEST_FLAG==2)cout<<"1P "<<flush;return(ptmp); }

    // latent period parameters for separated inoculated (same priors)
    if (mflag==4)  {
      ptmp += log(gsl_ran_gamma_pdf(plat[2],prior_latk[0],prior_latk[1]))
           +  log(gsl_ran_gamma_pdf(plat[3],prior_latm[0],prior_latm[1]));
    }
    if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag = 2;if(TEST_FLAG==2)cout<<"2P "<<flush;return(ptmp); }
  }

  // infectious period parameters
  ptmp += log(gsl_ran_gamma_pdf(pinf[0],prior_infk[0],prior_infk[1]))
       +  log(gsl_ran_gamma_pdf(pinf[1],prior_infm[0],prior_infm[1]));
  if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag = 3;if(TEST_FLAG==2)cout<<"3P "<<flush;return(ptmp); }

  // transmission parameters
  for (int bi=0;bi<bflag;++bi)  {
    //ptmp += log(gsl_ran_flat_pdf(beta[bi],prior_beta[0],prior_beta[1]));    // Testing flat prior...
    ptmp += log(gsl_ran_gamma_pdf(beta[bi],prior_beta[0],prior_beta[1]));
  }
  if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag = 4;if(TEST_FLAG==2)cout<<"4P "<<flush;return(ptmp); }

  // Constrain latent period between tNeg-tI and tPos-tI: E gets from tI to tN but before tP
  for (int j=0;j<n_l;++j)  {
    int i = id_allinf[j];
    ptmp += log(gsl_ran_flat_pdf(lat_P[i],max(0.0,tNeg[i]-tInf[i]),tPos[i]-tInf[i]));
    if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag = 5;if(TEST_FLAG==2)cout<<"\t\t\t\t\t5P-"<<i<<"\t"<<tInf[i] << " " <<lat_pS[j]<<" "<<flush;return(ptmp); }
    // Wee bit of early return... can skip one loop but not so efficient if decent acceptance rates?
    if (isinf(ptmp))  {
      if (TEST_FLAG==2)  {cout << endl;}
      return(ptmp);
    }
  }

  // Constrain infection time between 0/ORSELHACK/tMove and first positive
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
        ptmp += log(gsl_ran_flat_pdf(tInf[i],ORSELHACK,tPos[i]));
      }
      else  { // Latent periods assumed zero... so stronger constraints. Issue for inocs though?
        ptmp += log(gsl_ran_flat_pdf(tInf[i],tNeg[i],tPos[i]));
      }
    }
    if ((TEST_FLAG)&&(isinf(ptmp)))  { perr_flag=6;if(TEST_FLAG==2)cout<<"6P:"<<i<<" ";return(ptmp); }
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
  if (tNeg[i]<0)  { // Should never trigger...  never infected?
    cout << "WTF ni_calc" << endl;
    exit(-1);
  }
  // Animal i, want number of infectious left in room when i got infected
  int room = rooms[i];
  for (auto jit=rLst[room].begin();jit!=rLst[room].end();++jit)  {
    int j = *jit;   // animal id for each animal in this room
    if (i==j)  {
      continue; // obv don't want to count self. not that it'd actually do anything...
    }
    else  {
      // Check each other animal in this room
      if (tNeg[j]<0)  {
        // j was never infected/infectious, so no contirbution
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
        //else if ((tInf[i]>(tInf[j]+lat_pS[id_ij[j]])) && (tInf[i]<tEnd[j]))  {
        else if ((tInf[i]>(tInf[j]+lat_P[j])) && (tInf[i]<tEnd[j]))  {
          // iTyp is 2 for pen1 contact and 3 for pen2 contact
          (iTyp[j]==3) ? ninf[1]+=1.0 : ninf[0]+=1.0;
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
    int j = *jit; // animal id for each in same room as i
    if ((tNeg[j]<0)||(i==j))  {
      // j was never infected, no contribution
      continue;
    }
    else  {
        // Censored Inoc is assumed infectious from t==0
      switch (iTyp[j])  {
        case 0: // Inoc with no data - assuming to be infectious for full window [0,tMove]
          if (iTyp[i]==2)  { // C1 pig
            (iStop>tMove[room]) ? ninfsum[0]+=tMove[room]-ORSELHACK //exposed to inocs for full time
                                : ninfsum[0]+=iStop-ORSELHACK;      //infected before move
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
          if (iStop>(lat_P[j]+inf_p[id_ij[j]]))  {
            ninfsum[0] += inf_p[id_ij[j]];
          }
          else  {
            ninfsum[0] += max(iStop-(tInf[j]+lat_P[j]),0.0);
          }
          break;

        case 3: // Challenge from pen 2
          if (iStop>(lat_P[j]+inf_p[id_ij[j]]))  {
            ninfsum[1] += inf_p[id_ij[j]];
          }
          else  {
            ninfsum[1] += max(iStop-(tInf[j]+lat_P[j]),0.0);
          }
          break;

        case 4: // C2 challenge
          ninfsum[0] += max(iStop-(tInf[j]+lat_P[j]),0.0);
          break;
      }
/*        double istart = max(tInf[j]+lat_pS[id_ij[j]],double(tMove[room]));
        if (iTyp[j]==3)  {
          ninfsum[1] += max(istop-istart, 0.0);
        }
        else  {
          ninfsum[0] += max(istop-istart, 0.0);
        }*/
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
  for (int rm=rooms.front();rm<=rooms.back();++rm)  { // know rooms are incremental. otw iterator.
    int i_beg = rLst[rm].front();   // first animal index in this room
    int i_end = rLst[rm].back();    // last animal index
    vector<double> n_pen = {double(nroom1[rm]),double(nroom2[rm])};

    // For each animal in this room - calculate probability of infection
    for (int i=i_beg;i<=i_end;++i)  {
      if (iTyp[i]<2)  { // So iTyp 0 or 1 - inoculated
        prob_inf[i] = 1.0;
      }
      // Contact animals not nec infected, calc probs
      else  {
        // numbers of infectious in each pen when i is infected
        vector<double> ninf(bflags[rooms[i]],0.0);
        if (id_ij[i]!=-1)  { // Definitely infected - (tPos>=0)
          ni_calc(i,ninf);
        }
        // Cumulative infectious time for all animals in pen until i is infected
        // Need this whether or not i got infected
        vector<double> ninf_sum(bflags[rooms[i]],0.0);
        ns_calc(i,ninf_sum);

        // Now calculate probability of infection at sampled tInf
        if (bflags[rooms[i]]==1)  {  // No second pen
          switch (iTyp[i])  {
            case 2:
              // Normal contact
              if (tNeg[i]<0)  { // No transmission
                prob_inf[i] =  exp(-(beta[0]*ninf_sum[0]/n_pen[0]));
              }
              else  {           // Transmission
                prob_inf[i] = (beta[0]*ninf[0]/n_pen[0])
                            *  exp(-(beta[0]*ninf_sum[0]/n_pen[0]));
              }
              break;

            case 3:
              cout << "WTF lhood - bflags==1 but have iTyp==3?!" << endl;
              exit(-1);
              break;

            case 4:
              // Within-pen contact
              if (tNeg[i]<0)  { // No transmission
                prob_inf[i] =  exp(-(beta[0]*ninf_sum[0]/n_pen[0]));
              }
              else  {           // Transmission
                prob_inf[i] = (beta[0]*ninf[0]/n_pen[0])
                            *  exp(-(beta[0]*ninf_sum[0]/n_pen[0]));
              }
              break;

            default:
              cout << "WTF" << endl;
              exit(-1);
              break;
          }
          if ((TEST_FLAG)&&(prob_inf[i]==0.0))  {
            if(TEST_FLAG==2)cout << i << " " << ninf[0] << " " << ninf_sum[0] << endl;
          }
        }
        else  {
          // multiple pens (experiment - not individual rooms...)
          switch (iTyp[i])  {
            case 2:
              // Within-pen contact
              if (tNeg[i]<0)  { // No transmission
                prob_inf[i] =  exp( -(beta[0]*ninf_sum[0]/n_pen[0] + beta[1]*ninf_sum[1]/(n_pen[0]+n_pen[1])) );
              }
              else  {           // Transmission
                prob_inf[i] = (beta[0]*ninf[0]/n_pen[0] + beta[1]*ninf[1]/(n_pen[0]+n_pen[1]))
                            *  exp(-(beta[0]*ninf_sum[0]/n_pen[0]+beta[1]*ninf_sum[1]/(n_pen[0]+n_pen[1])));
              }
              break;

            case 3:
              // Between-pen transmission
              if (tNeg[i]<0)  { // Not infected contact
                prob_inf[i] =  exp(-(beta[0]*ninf_sum[1]/n_pen[1]+beta[1]*ninf_sum[0]/(n_pen[0]+n_pen[1])));
              }
              else  { // Infected
                prob_inf[i] = (beta[0]*ninf[1]/n_pen[1]+beta[1]*ninf[0]/(n_pen[0]+n_pen[1]))
                            *  exp(-(beta[0]*ninf_sum[1]/n_pen[1]+beta[1]*ninf_sum[0]/(n_pen[0]+n_pen[1])));
              }
              //cout << " " << betaw*ninf[1]/n_pen[1] << "\t" << betab*ninf[0]/(n_pen[0]+n_pen[1]) << endl;
              break;

            case 4:
              // Within-pen contact
              cout << "WTF LHOOD - bflags==2 && iTyp==4 ie between pen transmission on Orsel C2 pig?" << endl;
              exit(-1);
              break;

            default:
              cout << "WTF" << endl;
              exit(-1);
              break;
          }
        }
      }
    }
  }
}



/** \brief Calculate log-likelihood.
 * First calc probs of infection (move outside?)
 * Then contributions from gamma distributed latent and infectious periods.
 * \return double
 */
double Model::lhood_calc()  {
  double loglik = 0.0;
  // Get probabilities of infection for all contacts
  vector<double> prob_inf(n_a,0.0);
  probinf_calc(prob_inf);
  //Contributions from infection times for contact transmissionsn
  for (int i=0;i<n_a;++i)  {  // prob_inf already handles if transmission occurs or not
    loglik += (iTyp[i]>1) ? log(prob_inf[i]) : 0.0;
  }
  if ((TEST_FLAG)&&(isinf(loglik)))  { if(TEST_FLAG==2)cout<<"1L "; }

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
  }
  if ((TEST_FLAG)&&(isinf(loglik)))  { if(TEST_FLAG==2)cout<<"2L "; }

  // Contributions from infectious periods - potentially censored & definitely binned
  for (auto iit=id_allinf.begin();iit!=id_allinf.end();++iit)  {
    int i = *iit;
    switch (cull[i])  {
      case 0:
        // Binned - tEnd is last positive. infectiousness bounded < tEnd+1 (first negative)
        loglik += log(gsl_cdf_gamma_P(inf_p[id_ij[i]]+1.0,pinf[0],pinf[1]/pinf[0])
                     -gsl_cdf_gamma_P(inf_p[id_ij[i]],pinf[0],pinf[1]/pinf[0]));
        break;
      case 1:
        loglik += log(1.0-gsl_cdf_gamma_P(inf_p[id_ij[i]],pinf[0],pinf[1]/pinf[0]));
        break;
      default:
        // Should not happen!!
        cout << "LHOOD - cull status not zero or one WTF?" << endl;
        cin.get();
        break;
    }
  }
  if ((TEST_FLAG)&&(isinf(loglik)))  { if(TEST_FLAG==2)cout<<"3L " << endl; }

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
  out_lhd << logprr<< " " << perr_flag << " " << loglik<< "\n";
  // infectious periods, well, sort of. Not inferred, just tEnd-latp-tInf
  out_ipd << inf_p.transpose() << "\n";
}


/** \brief Close ofstreams.
 * \return void
 */
void Model::closefiles()  {
  out_par.close();
  out_lat.close();
  out_inf.close();
  out_ipd.close();
  out_lhd.close();
  out_brn.close();
}
