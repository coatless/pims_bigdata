#include <RcppEigen.h>
#include <omp.h>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/Sparse>
using namespace Rcpp;
using namespace Eigen;

// Type definition to save space
typedef Eigen::Triplet<double> T;

// Dependency / flag sets
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]


// [[Rcpp::export]]
double nll_eigen(double x,
                  const Eigen::MatrixXd& d, 
                  const Eigen::VectorXd& z,
                  double n){

  // Create the correlation matrix
  Eigen::MatrixXd corr_matrix = (-1/x * d).array().exp();
  
  // -----
  // Chol decomposition
  
  // Compute the Cholesky decomposition of A
  Eigen::LLT< Eigen::MatrixXd > llt(corr_matrix); 
  
  // Retrieve U matrix in the decomposition call it Q
  Eigen::MatrixXd Q = llt.matrixU(); 

  // Obtain the log determinant
  double logdet = 2.0 * Q.diagonal().array().log().sum();
  
  // Eigen library rox
  Eigen::VectorXd bs = Q.transpose().triangularView<Lower>().solve(z);

  // Inner product x^t %*% x
  double distval = bs.dot(bs);

  // -----  
  // Calculate and return y
  return (n * log(2.0 * M_PI) + n * log(distval / n) + logdet + n)/2.0;
} 

//' @title First Tapering Implementation
//' @description The function implements the first tapering method described
//' in the Covariance Tapering paper.
//' @param x A \code{double}
//' @param n A \code{double} that indicates the matrix row/col.
//' @param good_dists A \code{VectorXd} that contains the distances which meet tapering requirements
//' @param taps A \code{VectorXd} results of applying the taper function
//' @param ia A \code{VectorXd} containing the ordered row indices
//' @param ja A \code{VectorXd} containing the ordered column indices
//' @param z A \code{VectorXd} element position in matrix 
//' @param rescol An \code{int} that contains how many non-zero elements should be in a column. 
//' @return A \code{double}
//' @details Currently the EIGEN sparse matrix cholesky decomposition differs from that of SPAM. 
//' As a result, the forward solve step finds a different vector. This in turn returns a higher distval.
//' @examples
//' data(anom1962)
//' d = rdist_earth1(loc)
//' setup.eigen = make_tapersetup_eigen(d,taprange = 50)
//' val = nll_1taper(20, setup.eigen$n,
//'                  setup.eigen$good.dists,
//'                  setup.eigen$taps,
//'                  setup.eigen$ia,
//'                  setup.eigen$ja,
//'                  z,
//'                  setup.eigen$rescol)
// [[Rcpp::export]]
double nll_1taper(double x, double n,
                  const Eigen::VectorXd& good_dists, 
                  const Eigen::VectorXd& taps,
                  const Eigen::VectorXd& ia,
                  const Eigen::VectorXd& ja,
                  const Eigen::VectorXd& z,
                  unsigned int rescol){
  
  // -----
  // Creation and filling of a sparse matrix

  // Define sparse matrix with ColMajor orientation
  Eigen::SparseMatrix<double> corr_matrix_Taper(n,n);
  
  // Begin filling sparse matrix
  std::vector<T> tripletList;
  
  // Guess how many entries are non-empty per column
  tripletList.reserve(rescol);
    
  // Fill matrix with pairs
  for(unsigned int i = 0; i < good_dists.size(); i++){
    tripletList.push_back( T( ia(i), ja(i),  // i,j
                             exp(-good_dists(i)/x)*taps(i) // v
                             ) );
  }
  
  corr_matrix_Taper.setFromTriplets(tripletList.begin(), tripletList.end());
  // End fill statement
    
  // Release memory
  tripletList.clear();
  

  // -----
  // Chol decomposition
  
  // compute the Cholesky decomposition of A
  Eigen::SimplicialLLT< Eigen::SparseMatrix<double> > llt_sparse(corr_matrix_Taper); 
  
  // retrieve U matrix in the decomposition call it Q
  Eigen::SparseMatrix<double> Q = llt_sparse.matrixU(); 
    
  // -----
  // Taking a determinant 
  
  // Setup the SparseLU solver
  Eigen::SparseLU< Eigen::SparseMatrix<double> > slu;
  
  // Compute the derivative
  slu.compute(Q);
  
  // Retrieve the log modulus
  double logdet = 2.0 * slu.logAbsDeterminant(); 
  
  // -----  
  // Forward solve  
  VectorXd bs = Q.triangularView<Lower>().solve(z);
  
  // Inner product x^t %*% x
  double distval = bs.dot(bs);
    
  // -----  
  // Calculate and return y
  return (n * log(2.0 * M_PI) + n * log(distval / n) + logdet + n)/2.0;
} 

//' @title Second Tapering Implementation
//' @description The function implements the second tapering method described
//' in the Covariance Tapering paper.
//' @param x A \code{double}
//' @param n A \code{double} that indicates the matrix row/col.
//' @param good_dists A \code{VectorXd} that contains the distances which meet tapering requirements
//' @param taps A \code{VectorXd} results of applying the taper function
//' @param ia A \code{VectorXd} containing the ordered row indices
//' @param ja A \code{VectorXd} containing the ordered column indices
//' @param z A \code{VectorXd} element position in matrix 
//' @param rescol An \code{int} that contains how many non-zero elements should be in a column. 
//' @return A \code{double}
//' @details Currently the EIGEN sparse matrix cholesky decomposition differs from that of SPAM. 
//' As a result, the backsolve and forwardsolve step finds a different vector. 
//' This in turn returns a higher distval.
//' @examples
//' data(anom1962)
//' d = rdist_earth1(loc)
//' setup.eigen = make_tapersetup_eigen(d,taprange = 50)
//' val = nll_2taper(20, setup.eigen$n,
//'                  setup.eigen$good.dists,
//'                  setup.eigen$taps,
//'                  setup.eigen$ia,
//'                  setup.eigen$ja,
//'                  z,
//'                  setup.eigen$rescol)
// [[Rcpp::export]]
double nll_2taper(double x, double n,
                  const Eigen::VectorXd& good_dists, 
                  const Eigen::VectorXd& taps,
                  const Eigen::VectorXd& ia,
                  const Eigen::VectorXd& ja,
                  const Eigen::VectorXd& z,
                  unsigned int rescol){

  // -----
  // Creation and filling of a sparse matrix
  
  // Define sparse matrix with ColMajor orientation
  Eigen::SparseMatrix<double> corr_matrix_Taper(n,n);
  
  // Begin filling sparse matrix
  std::vector<T> tripletList;
  
  // Guess how many non-empty elements exist per column
  // Better guesses => faster timing
  tripletList.reserve(rescol);
    
  // Fill matrix with entires
  for(unsigned int i = 0; i < good_dists.size(); i++){
    tripletList.push_back( 
                           T( 
                              ia(i), ja(i),  // i, j matrix positions
                              exp(-good_dists(i)/x)*taps(i) // value
                             ) 
                         );
  }
  
  corr_matrix_Taper.setFromTriplets(tripletList.begin(), tripletList.end());
  // End matrix fill statement
  
  // Release memory
  tripletList.clear();
  
  // -----
  // Chol decomposition
  
  // Compute the Cholesky decomposition of A
  Eigen::SimplicialLLT< Eigen::SparseMatrix<double> > llt_sparse(corr_matrix_Taper); 
  
  // Retrieve U matrix in the decomposition call it Q
  Eigen::SparseMatrix<double> Q = llt_sparse.matrixU(); 
   
  // -----
  // Taking a determinant 
  
  // Setup the SparseLU solver
  SparseLU< SparseMatrix<double> > slu;
  
  // Compute the derivative
  slu.compute(Q);
  
  // Retrieve the log modulus
  double logdet = 2.0 * slu.logAbsDeterminant(); 
  
    
  // -----
  
  // ----------------------------------------------------
  //
  // This section is new when compared to nll1_taper
  //
  // ----------------------------------------------------

  // Backward solve (upper triangular) on a Forward Solve (lower triangular)
  
  // Create an identity matrix
  Eigen::MatrixXd I;
  I.setIdentity(n, n);

  // Access triangular views and solve.
  Eigen::MatrixXd inv = Q.triangularView<Eigen::Upper>().solve(
                            Q.triangularView<Eigen::Lower>().solve(I)
                         );

  // -----
  // Refill matrix with updated entries based on inv
  
  // Guess how many non-empty elements exist per column
  // Better guesses => faster timing
  tripletList.reserve(rescol);
    
  // Fill matrix with pairs
  for(unsigned int k = 0; k < ia.size(); k++){
    unsigned int i = ia(k), j = ja(k); // temp storage
    
    tripletList.push_back( T( i, j,  // i,j old positions
                             inv(i,j)*taps(k) // new v
                             ) );
  }
  
  corr_matrix_Taper.setFromTriplets(tripletList.begin(), tripletList.end());
  // End fill statement

  // -----  
  // Distance calculation z^T * Corr * z
  double distval = z.transpose() * corr_matrix_Taper * z;

  // -----
  // Calculate and return y
  return (n * log(2.0 * M_PI) + n * log(distval / n) + logdet + n)/2.0;
}

//' @describeIn nll_1taper
//' @param cores An \code{int} that indicates the number of cores to divide the task over.
// [[Rcpp::export]]
Eigen::VectorXd nll_1taper_parallel(
                  Eigen::VectorXd x, double n,
                  const Eigen::VectorXd& good_dists, 
                  const Eigen::VectorXd& taps,
                  const Eigen::VectorXd& ia,
                  const Eigen::VectorXd& ja,
                  const Eigen::VectorXd& z,
                  unsigned int rescol,
                  unsigned int cores){
  unsigned int i;
  
  #pragma omp parallel for num_threads(cores) private(i)
  for(i = 0; i < x.size(); i++){
    x(i) = nll_1taper(x(i), n, good_dists, taps, ia, ja, z, rescol);
  }
  
  return x;
} 

//' @describeIn nll_2taper
//' @param cores An \code{int} that indicates the number of cores to divide the task over.
// [[Rcpp::export]]
Eigen::VectorXd nll_2taper_parallel(
                  Eigen::VectorXd x, double n,
                  const Eigen::VectorXd& good_dists, 
                  const Eigen::VectorXd& taps,
                  const Eigen::VectorXd& ia,
                  const Eigen::VectorXd& ja,
                  const Eigen::VectorXd& z,
                  unsigned int rescol,
                  unsigned int cores){
  
  // Declare variables
  unsigned int i;
  
  // Request parallelization for the loop
  // Parallelize using N cores
  // Keep the looping variable (i) private
  #pragma omp parallel for num_threads(cores) private(i)
  for(i = 0; i < x.size(); i++){
    x(i) = nll_2taper(x(i), n, good_dists, taps, ia, ja, z, rescol);
  }
  
  return x;
}
