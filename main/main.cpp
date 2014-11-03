#include "distillery.h"
#include "input_parms.h"

extern "C" {
  #include "pre_wrapper.h"
}

int main(int argc, char *argv[]){

  pre_wrapper(argc, argv);

  // input parameter -> see struct definition for input parameter description
  // TODO: write class that can read input file!
  LapH::input_parameter param;
  param.parse_input_file(argc, argv);

  Eigen::setNbThreads(4);
  Eigen::initParallel();

  // memory for sources
  // see discription of distillation::create_source for explanation for 
  // two-dimensional source vector!
  std::complex<double>** source = 
    new std::complex<double>*[param.dilution_size_so[2]];
  for(size_t dil_d = 0; dil_d < param.dilution_size_so[2]; ++dil_d){
    source[dil_d] = 
              new std::complex<double>[3*4*param.Lt*param.Ls*param.Ls*param.Ls];
  }
  // memory for propagator 
  std::complex<double>* propagator = 
              new std::complex<double>[3*4*param.Lt*param.Ls*param.Ls*param.Ls];

  // initialising distillation class
  LapH::distillery dis(param);

  // loop over configurations
  for(size_t config = param.config; config <= param.config; config += 8){
    tmLQCD_read_gauge(config);
    std::cout << "processing config: " << config << "\n" << std::endl;
    // resetting everything
    if(config != param.config) {
      param.config = config;
      dis.reset_all(param);
    }
    // loop over random vectors
    for(size_t rnd_id = 0; rnd_id < param.nb_rnd; ++rnd_id) {
      // loop over all inversions
      for(size_t dil_t = 0; dil_t < param.dilution_size_so[0]; ++dil_t){
        for(size_t dil_e = 0; dil_e < param.dilution_size_so[1]; ++dil_e){         
          // creating sources --------------------------------------------------
          dis.create_source(dil_t, dil_e, source); // --------------------------
          for(size_t dil_d = 0; dil_d < param.dilution_size_so[2]; ++dil_d){

            std::cout << dil_t << "\t" << dil_e << "\t" << dil_d << std::endl;
   
      	    // I assume here for the moment that the first operator defined 
            // in the input file is the one to be used.
      	    // I also assume that this is _not_ a doublet type of operator
      	    invert((double *) propagator, (double *) (source[dil_d]));

            // add the propagatro to the perambulator --------------------------
            dis.add_to_perambulator(dil_t, dil_e, dil_d, propagator); // -------
            
          }
        }
      } // end of loop over inversions
      // writing perambulator to disk ------------------------------------------
      dis.write_perambulator_to_disk(); // -------------------------------------
      // creating new random vector --------------------------------------------
      if(rnd_id < param.nb_rnd-1)
        dis.reset_perambulator_and_randomvector(rnd_id+1);
    } // end of loop over random vectors
  } // end of loop over configurations

	return 0;
}

