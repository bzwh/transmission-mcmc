#ifndef MODEL_HPP_INCLUDED
#define MODEL_HPP_INCLUDED

#include <fstream>
#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
class Model  {
public:
  Model();
  void setup(int,int,int);
  void setfns(int,int);
  void load_data();
  void ncount();
  Eigen::VectorXd init_samp(gsl_rng*,double&,double&,int&);
  void load_priors();
  void parse(Eigen::VectorXd);
  void ni_calc(int,std::vector<double>&);
  void ns_calc(int,std::vector<double>&);
  void ns_calc2(int,std::vector<double>&);
  void probinf_calc(std::vector<double>&);
  void probinf_calc2(std::vector<double>&);
  double lhood_calc();
  double prior_calc();

  void write(Eigen::VectorXd,double,double);
  void writeburn(Eigen::VectorXd);
  void closefiles();
  void test_data();

  // File names and paths for input and outputs
  std::string dname;                      /// Data filename
  std::string pname;                      /// Priors filename
  std::string opath;                      /// Path to write outputs
  std::ofstream out_par;                  /// Dump whole parameter vector
  std::ofstream out_lat;                  /// Dump latent periods
  std::ofstream out_inf;                  /// Dump infection times
  std::ofstream out_ipd;                  /// Dump infectious periods
  std::ofstream out_lhd;                  /// Dump likelihoods
  std::ofstream out_brn;

  // Data stuffs
  int bflag;                              /// data includes multi-pen setup
  int mflag;                              /// a toggle for separated and combined latent periods
  int n_a;                                /// total number of animals
  int n_c;                                /// number of contacts
  int n_r;                                /// number of rooms
  int n_l;                                /// number latent periods to infer (ignore eg iTyp==0)
  int n_ii;                               /// number of infecteds - inoc+successful contacts
  int n_ti;                               /// number infection times (n_l - #inoc ?) to infer
  int exend;                              /// length of experiment
  double orsel_delay;                     /// Intro delay for first challenge. Orsel default =1.0
  std::vector<int> bflags;                /// number of pens in each room
  std::vector<int> cflags;                /// indicates C2 challenge (orsel pigs)
  std::vector<int> pflags;
  std::vector<int> id_allinf;   /// all infecteds. for lat * infectious periods. size n_l - rescaled to n_a
  std::vector<int> id_ij;       /// Converts i \in [0,n_a] to j \in [0,n_ii] for access to lat_p/inf_p
  std::vector<int> id_contact;  /// +ve contacts for infection times. x[k \in n_i] = {i \in [0,n_a]}
  std::vector< std::vector<int> > rLst;   /// rLst[room id] = {animal IDs}
  std::vector<int> iTyp;                  /// 1-inoc 2-contact 3-between pen contact
  std::vector<int> tNeg;                  /// time of last negative
  std::vector<int> tPos;                  /// time of first positive
  std::vector<int> tEnd;                  /// time removed/recovered(?)
  std::vector<int> cull;                  /// culled or not - ie censored infectious period
  std::vector<double> tMove;              /// time animals moved (for each room) - 0 if not moved
  std::vector<int> rooms;                 /// rooms[i] room id for animal i
  std::vector<int> nroom1;                /// number of animals in each room's main pen
  std::vector<int> nroom2;                /// number of animals in each room's 2nd pen(?)
  std::vector< std::vector<int> > deltaN;
  Eigen::MatrixXd ntot_pens0;
  std::vector< std::vector<int> > deltaN1;
  Eigen::MatrixXd ntot_pens1;
  std::vector< std::vector<int> > deltaN2;
  int npar;                               /// Number of parameters being fit

  // Priors
  std::vector<double> prior_latk; /// shape and scale of gamma latent shape
  std::vector<double> prior_latm; /// shape and scale of gamma latent mean
  std::vector<double> prior_infk; /// shape and scale of gamma infectious shape
  std::vector<double> prior_infm; /// shape and scale of gamma infectious mean
  std::vector<double> prior_beta; /// shape and scale of gamma transmission parameter
  double perr_flag;
  double lerr_flag;

  // Time infected
  Eigen::VectorXd tInf;     /// Time infected for all animals. 0 for inocs.
  Eigen::VectorXd tinfS;    /// (n_i) purely for convenience parsing gen'd pars. Immediately populates tInf, never used directly.
  // Latent period duration
  Eigen::VectorXd lat_P;    /// Latent period duration for all animals. -1 when not inferred/irrelevant
  Eigen::VectorXd lat_pS;   /// (n_l) unobserved latent periods. tEnd-t_inf-lat_p = infectious_p
  Eigen::VectorXd inf_p;    /// (n_l) infectious periods - possibly censored, def binned
  // Shape and mean parameters of gamma distributed latent/infectious periods
  Eigen::VectorXd plat;     /// Shape and mean parameters for latent period.
  Eigen::VectorXd pinf;     /// Shape and mean parameters for infectious period.
  // Transmission parameters - (w)ithin pen and (b)etween pen
  Eigen::VectorXd beta;     /// Transmisssion parameters - within and between pens

  Eigen::VectorXd parsamp;
private:
};


#endif // MODEL_HPP_INCLUDED
