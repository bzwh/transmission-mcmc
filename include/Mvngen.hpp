#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace Eigen;




/** \brief Given variance-covariance matrix, sample from zero mean multinormal distribution.
 * Decomposes the matrix first using dcmp
 * \param sigma MatrixXd
 * \param r gsl_rng*
 * \return VectorXd
 */
VectorXd mgen(const MatrixXd&, gsl_rng*);


/** \brief Decompose given matrix. Tries Choleksy, if fails does eigendecomp instead
 * \param sigma const MatrixXd&
 * \return MatrixXd
 */
VectorXd mgen_dcmp(const MatrixXd&, gsl_rng*);


MatrixXd dcmp(const MatrixXd&);
