#ifndef BASICOPERATOR_H_
#define BASICOPERATOR_H_

#include <complex>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Core>

#include "macros.h"
#include "ranlxs.h"
#include "input_parms.h"
#include "tmLQCD.h"
#include "typedefs.h"

// from twisted mass code
#ifdef __cplusplus                                                              
extern "C"                                                                      
{                                                                               
#endif /* __cplusplus */     

#include "ranlxd.h"
#include "su3.h"
#include "su3adj.h"
#include "io/spinor.h"
#include "start.h"

#ifdef __cplusplus                                                              
}                                                                               
#endif  


#include "mpi.h"
#include "omp.h"

#undef complex

extern "C" void zgemm_(char* TRANSA, char* TRANSB, const int* M,
                       const int* N, const int* K, std::complex<double>* alpha, std::complex<double>* A,
                       const int* LDA, std::complex<double>* B, const int* LDB, std::complex<double>* beta,
                       std::complex<double>* C, const int* LDC); 

namespace LapH {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class distillery{
  // ------------------------ METHODS ------------------------------------------
public:
  // standard constructor - doesn't do much -> Work is done in this->initialise
  distillery(){
    V = NULL;
    random_vector = NULL;
    memory_allocated = false;
  };
  // standard deconstructor
  ~distillery(){ 
    if(memory_allocated == true)
      this->clean();
  };

  // ---------------------------------------------------------------------------
  // ------------------------- INTERFACE ---------------------------------------
  // Initialisation:
  // allocation of memory, reading eigenvectors, generating randomvector
  // input: in_param -> all necessary input parameters (see struct definition)
  void initialise(const input_parameter& in_param);
  // cleaning everything up and freeing all memory
  void clean();
  // setting everything to zero again, and reding new eigenvectors
  // only difference to initialise is that no new memory is allocated
  // input: in_param -> all necessary input parameters (see struct definition)
  void reset_all(const input_parameter& in_param);
  // setting the perambulators to zero again, e.g. for a new rnd vector run
  // and generating a new random vector
  // input: seed -> new seed to generate new random vector
  //        random_vector_number -> just the id of random vector in progress
  void reset_perambulator_and_randomvector(const size_t rnd_id);
  // creates all sources for inversions 
  // output:  sources -> complex vector to store the sources 
  void create_source(std::complex<double>** source);
  // multiplies a propagator with V^dagger and stores it in one column of
  // of the perambulator. The column number is computed with dil_t, dil_e,
  // and dil_d
  // input: propagator -> the propagator contains the inverted Dirac
  //                      operator
  void add_to_perambulator(const std::complex<double>* const * const propagator);
  // writes the perambulator to disk
  // the indices from slow to fast: t -> Dirac -> eigenvectors :: This ordering
  // is important because it speeds up the multiplication with V^dagger which
  // is blockdiagonal in Dirac space
  // input: rnd_id -> id of the current randomvector
  void write_perambulator_to_disk(const size_t rnd_id);
  // ----------------------- END INTERFACE -------------------------------------
  // ---------------------------------------------------------------------------

private:
  DISALLOW_COPY_AND_ASSIGN(distillery); // no copy nor assignment
  void read_eigenvectors();
  size_t create_eigenvector_chunks(std::vector<ev_chunk>& ev_chunks);
  void set_random_vector(const size_t rnd_id);
  void set_sink_random_vector(const size_t rnd_id, const size_t sink_id, 
                              const size_t rnd_id_si, Eigen::VectorXcd& out);
  void read_random_vector();
  void copy_to_V(const std::complex<double>* const eigen_vec, const int t,
                 const int nev);
  int write_random_vector_to_disk(const size_t rnd_id);
  // -------------------------- DATA -------------------------------------------
  Eigen::MatrixXcd* V;                // memory eigensystem
  std::vector<std::vector<
    std::vector<Eigen::MatrixXcd> > > 
                        perambulator; // memory perambulator
  Eigen::VectorXcd* random_vector;    // memory random vector
  std::vector<
      std::vector<Eigen::VectorXcd> >
                    random_vector_si; // memory random vector
  input_parameter param;              // all necessary input parameters
  tmLQCD_mpi_params* tmLQCD_params;   // tmLQCD MPI parameters
  // some flags
  bool memory_allocated;              // is all memory allocated?

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

}

#endif /* BASICOPERATOR_H_ */
