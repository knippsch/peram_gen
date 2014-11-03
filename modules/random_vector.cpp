#include "random_vector.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::random_vector::set(const int seed, const size_t length) {

  // resizing the random vector
  vec.resize(length);

  // initialisation of the rando vector to create Z2 random vector
  rlxs_init(0, seed);
  std::vector<float> rnd(2*length);
  ranlxs(&(rnd[0]), 2*length);

  // generating a Z_2 source
  for(size_t i = 0; i < length; ++i){
    const double sqrt2 = 0.5*sqrt(2.0);
    double re, im;
    if (rnd[2*i] < 0.5)
      re = -sqrt2;
    else
      re = sqrt2;
    if (rnd[2*i + 1] < 0.5)
      im = -sqrt2;
    else
      im = sqrt2;
    vec[i] = cmplx(re, im);
  }

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::random_vector::set(const int seed, const size_t length,
                              const std::string& filename) {
  set(seed, length);
  write_random_vector(filename);

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::random_vector::write_random_vector(
                                       const std::string& filename) const {
  // writing random vector to file
  FILE* fp = NULL;
  if((fp = fopen(filename.c_str(), "wb")) == NULL){
    std::cout << "failed to open file to write random vector: " 
              << filename << "\n" << std::endl;
    exit(0);
  }   
  int check_read_in = fwrite(&(vec[0]), sizeof(cmplx), vec.size(), fp);
  if(check_read_in !=  (int) vec.size())
    std::cout << "It seems that not all data are written to: "
              << filename.c_str() << "\n" << std::endl;
  fclose(fp);

} 
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::random_vector::read_random_vector(const std::string& filename) {

//  char infile[400];
//  FILE *fp = NULL;
//  const int Lt = param.Lt;
//  const int verbose = param.verbose;
//  const int number_of_rnd_vec = 1;
//  const int number_of_eigen_vec = param.nb_ev;
//  const int rnd_vec_length = number_of_eigen_vec * 4;
//  // memory for reading random vectors
//  cmplx* rnd_vec_read =
//      new cmplx[Lt*rnd_vec_length];
//  std::string filename = "/data/LapHs/peram_gen/new/test/randomvector";
//
//  if(verbose) printf("reading random vectors from files:\n");
//  else printf("\treading random vectors:");
//  int check_read_in = 0;
//
//  for(int rnd_vec_i = 0; rnd_vec_i < number_of_rnd_vec; ++rnd_vec_i){
//    // read random vector
//    sprintf(infile, "%s.%03d.u.Tf.%04d", filename.c_str(), 0, 1000);
//    if(verbose) printf("\tread file: %s\n", infile);
//    if((fp = fopen(infile, "rb")) == NULL){
//      std::cout << "failed to open file: " << infile << "\n" << std::endl;
//      exit(0);
//    }   
//    check_read_in += fread(rnd_vec_read, sizeof(cmplx),
//        Lt*rnd_vec_length, fp);
//    // copy into matrix structure
//    for(int t = 0; t < Lt; t++){
//      for(int row_i = 0; row_i < rnd_vec_length; ++row_i){
//        random_vector[t](row_i) = rnd_vec_read[t * rnd_vec_length + row_i];
//      }   
//    }
//    fclose(fp);
//  }   
//  delete[] rnd_vec_read;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
 
