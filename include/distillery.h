#ifndef BASICOPERATOR_H_
#define BASICOPERATOR_H_

#include <array>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include "boost/multi_array.hpp"

#include "macros.h"
#include "input_parms.h"
#include "random_vector.h"

namespace LapH {

// typedefs for boost multidimensional array - no need fo nasty new and delete
typedef boost::multi_array<Eigen::MatrixXcd, 2> array_type;
typedef std::complex<double> cmplx;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class distillery{
  // ------------------------ METHODS ------------------------------------------
public:
  // standard constructor - Construction is done in initialise
  distillery(const input_parameter& in_param) : 
                          V(), rnd_vec_source(), rnd_vec_sink(), 
                          perambulator(), param(in_param) {
    initialise();
  };
  // standard deconstructor
  ~distillery() {};

  // ---------------------------------------------------------------------------
  // ------------------------- INTERFACE ---------------------------------------
  // setting everything to zero again, and reding new eigenvectors
  // only difference to initialise is that no new memory is allocated
  // input: in_param -> all necessary input parameters (see struct definition)
  void reset_all(const input_parameter& in_param);
  // setting the perambulators to zero again, e.g. for a new rnd vector run
  // and generating a new random vector
  // input: rnd_id -> id of the new randomvector
  void reset_perambulator_and_randomvector(const size_t rnd_id);
  // creates the source for inversions -- param.dilution_type_so[2] sources are
  // created simultaneously. This happens because the sources in Dirac space
  // are very similar and it is faster to create these sources at once
  // input:   dil_t -> counts the number of inversions in time
  //          dil_e -> counts the number of inversions in eigenvector space
  // output:  source -> complex vector to store the sources 
  void create_source(const size_t dil_t, const size_t dil_e, 
                     std::complex<double>** source);
  // multiplies a propagator with V^dagger and stores it in one column of
  // of the perambulator. The column number is computed with dil_t, dil_e,
  // and dil_d
  // input:   dil_t -> counts the number of inversions in time
  //          dil_e -> counts the number of inversions in eigenvector space
  //          dil_d -> counts the number of inversions in Dirac space
  //          propagator -> the propagator contains the inverted Dirac
  //                        operator
  void add_to_perambulator(const size_t dil_t, 
                           const size_t dil_e, 
                           const size_t dil_d, 
                           const std::complex<double>* const propagator);
  // writes the perambulator to disk
  // the indices from slow to fast: t -> Dirac -> eigenvectors :: This ordering
  // is important because it speeds up the multiplication with V^dagger which
  // is blockdiagonal in Dirac space
  void write_perambulator_to_disk();
  // ----------------------- END INTERFACE -------------------------------------
  // ---------------------------------------------------------------------------

private:
  DISALLOW_COPY_ASSIGN_MOVE(distillery); // no copy, assignment nor move
  // Initialisation:
  // allocation of memory, reading eigenvectors, generating randomvector
  void initialise();
  // reading eigenvectors
  void read_eigenvectors();

  // -------------------------- DATA -------------------------------------------
  std::vector<Eigen::MatrixXcd> V; // memory eigensystem
  random_vector rnd_vec_source;    // memory random vector
  std::vector<random_vector> rnd_vec_sink;      // memory random vector
  // memory perambulator - this is a two-dimensional array!!
  array_type perambulator;   
 
  input_parameter param;           // all necessary input parameters

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

}

#endif /* BASICOPERATOR_H_ */
