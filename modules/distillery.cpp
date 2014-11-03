#include "./distillery.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::initialise() {

  param.print_options();
  const size_t Ls = param.Ls;
  const size_t Lt = param.Lt;
  const size_t dim_row = Ls * Ls * Ls * 3;
  const size_t cnfg = param.config;

  // Initializing memory for eigenvectors, perambulator and random vector
  V.resize(Lt);
  perambulator.resize(boost::extents[param.nb_sinks][Lt]);
  for(size_t t = 0; t < Lt; ++t){
    V[t] = Eigen::MatrixXcd::Zero(dim_row, param.nb_ev);
    for(size_t sinks = 0; sinks < param.nb_sinks; sinks++)
      perambulator[sinks][t] = Eigen::MatrixXcd::Zero(4 * param.nb_ev,
                                                    param.dilution_size_so[0] * 
                                                    param.dilution_size_so[1] * 
                                                    param.dilution_size_so[2] );
  }
  // reading eigenvectors from disk
  read_eigenvectors(); 
  // generating random vector for source
  rnd_vec_source.set((param.seed[0])^cnfg, Lt*4*param.nb_ev, 
                      param.create_out_filename("randomvector", 0));
  // generating random vector for sink
  size_t nb_st_si = param.nb_stoch_sinks;
  rnd_vec_sink.resize(nb_st_si);
  for(size_t sinks = 0; sinks < param.nb_stoch_sinks; sinks++){
    const size_t i = (param.nb_sinks - nb_st_si) + sinks;
    // length of the random vector
    const size_t length = 12*Lt*param.dilution_sink[i].dilution_size_si[1];
    rnd_vec_sink[sinks].set((param.dilution_sink[i].seed)^cnfg, length, 
                            param.create_out_filename("randomvector", 0));
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::reset_all(const LapH::input_parameter& in_param){

  param = in_param;
  const size_t Ls = param.Ls;
  const size_t Lt = param.Lt;
  const size_t dim_row = Ls * Ls * Ls * 3;

  for(size_t t = 0; t < Lt; ++t){
    V[t] = Eigen::MatrixXcd::Zero(dim_row, param.nb_ev);
    for(size_t sinks = 0; sinks < param.nb_sinks; sinks++)
      perambulator[sinks][t] = Eigen::MatrixXcd::Zero(4 * param.nb_ev,
                                                    param.dilution_size_so[0] * 
                                                    param.dilution_size_so[1] * 
                                                    param.dilution_size_so[2] );
  }
  read_eigenvectors();  
  rnd_vec_source.set((param.seed[0])^param.config, Lt*4*param.nb_ev, 
                      param.create_out_filename("randomvector", 0));
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::reset_perambulator_and_randomvector(const size_t rnd_id){

  // setting perambulator and random vector to zero 
  for(size_t t = 0; t < param.Lt; ++t){
    for(size_t sinks = 0; sinks < param.nb_sinks; sinks++)
      perambulator[sinks][t] = Eigen::MatrixXcd::Zero(4 * param.nb_ev,
                                                    param.dilution_size_so[0] * 
                                                    param.dilution_size_so[1] * 
                                                    param.dilution_size_so[2] );
  }
  // create new random vector
  rnd_vec_source.set((param.seed[rnd_id])^param.config, param.Lt*4*param.nb_ev, 
                      param.create_out_filename("randomvector", 0));

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Creates the dilution lookup table. 
// input:  dilution_size  -> size of dilution, e.g. block size or interlace size
//         dilution_entry -> index of the dilution
//         type           -> just the type of the dilution
// output: lookup         -> the final lookup table for dilution
static void create_dilution_lookup(const size_t nb_of_nonzero_entries, 
                                   const size_t nb_inversions, 
                                   const size_t dilution_entry, 
                                   const std::string& type, size_t lookup[]){
  if(type == "B" || type == "N")
    for(size_t i = 0; i < nb_of_nonzero_entries; ++i)
      lookup[i] = nb_of_nonzero_entries * dilution_entry + i;
  else if(type == "I" || type == "F")
    for(size_t i = 0; i < nb_of_nonzero_entries; ++i)
      lookup[i] = dilution_entry + i * nb_inversions;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::create_source(const size_t dil_t, const size_t dil_e,
                                     cmplx** source) {
  
  if(dil_t >= param.dilution_size_so[0] ||
     dil_e >= param.dilution_size_so[1] ){
        std::cout << "dilution is out of bounds in \"create_source\"" <<
        std::endl;
        exit(0);
  }

  const size_t Lt = param.Lt;
  const size_t Ls = param.Ls;
  const size_t number_of_eigen_vec = param.nb_ev;
  const size_t Vs = Ls*Ls*Ls;
  const size_t dim_row = Vs*3;

  // setting source to zero 
  for(size_t dil_d = 0; dil_d < param.dilution_size_so[2]; ++dil_d)
    for(size_t i = 0; i < 12*Vs*Lt; ++i)
      source[dil_d][i] = {0.0, 0.0};

  // indices of timeslices with non-zero entries
  size_t nb_of_nonzero_t = Lt/param.dilution_size_so[0];  
  // TODO: think about smart pointer here!
  std::vector<size_t> t_index(nb_of_nonzero_t);
  create_dilution_lookup(nb_of_nonzero_t, param.dilution_size_so[0], 
                         dil_t, param.dilution_type_so[0], &(t_index[0]));

  // indices of eigenvectors to be combined
  size_t nb_of_ev_combined = number_of_eigen_vec/param.dilution_size_so[1];
  std::vector<size_t> ev_index(nb_of_ev_combined);
  create_dilution_lookup(nb_of_ev_combined, param.dilution_size_so[1], 
                         dil_e, param.dilution_type_so[1], &(ev_index[0]));

  // indices of Dirac components to be combined 
  size_t nb_of_dirac_combined = 4/param.dilution_size_so[2];
  boost::multi_array<size_t, 2> 
       d_index(boost::extents[param.dilution_size_so[2]][nb_of_dirac_combined]);
  for(size_t dil_d = 0; dil_d < param.dilution_size_so[2]; ++dil_d)
    create_dilution_lookup(nb_of_dirac_combined, param.dilution_size_so[2],
                        dil_d, param.dilution_type_so[2], &(d_index[dil_d][0]));

  // creating the source 
  // running over nonzero timeslices
  for(size_t t = 0; t < nb_of_nonzero_t; ++t){ 
    // intermidiate memory
    Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(dim_row, 4);
    const size_t time = t_index[t];         // helper index
    // building source on one timeslice
    for(size_t ev = 0; ev < nb_of_ev_combined; ++ev){
      const size_t ev_h = ev_index[ev] * 4; // helper index 
      for(size_t d = 0; d < 4; ++d){
        S.col(d) += rnd_vec_source[time*4*number_of_eigen_vec + ev_h + d] * 
                    (V[time]).col(ev_index[ev]); 
      }
    }
    // copy the created source into output array
    size_t t_h = time*Vs; // helper index
    for(size_t x = 0; x < Vs; ++x){
      const size_t x_h  = (t_h + x)*12;  // helper index
      const size_t x_h2 = x*3;           // helper index
      for(size_t d2 = 0; d2 < param.dilution_size_so[2]; ++d2){
        for(size_t d3 = 0; d3 < nb_of_dirac_combined; ++d3){
          const size_t d_h = x_h + d_index[d2][d3]*3; // helper index
          for(size_t c = 0; c < 3; ++c){
            source[d2][d_h + c] = S(x_h2 + c, d_index[d2][d3]);
          }
        }
      }    
    }
  } // end of loop over nonzero timeslices

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::add_to_perambulator(const size_t dil_t, 
                          const size_t dil_e, const size_t dil_d,
                          const cmplx* const propagator) { 
 
  if(dil_t >= param.dilution_size_so[0] ||
     dil_e >= param.dilution_size_so[1] ||
     dil_d >= param.dilution_size_so[2] ){
        std::cout << "dilution is out of bounds in \"add_to_perambulator\"" <<
        std::endl;
        exit(0);
  }
 
  const size_t Lt = param.Lt;
  const size_t Ls = param.Ls;
  const size_t vec_length = Ls*Ls*Ls*3;
  const size_t nb_ev = param.nb_ev;
  // computing the column index of perambulator
  const size_t col = dil_t*param.dilution_size_so[1]*param.dilution_size_so[2] +
                     dil_e*param.dilution_size_so[2] + 
                     dil_d;

  // intermediate memory for multiplication with V^dagger
  Eigen::MatrixXcd vec = Eigen::MatrixXcd::Zero(vec_length, 4);

  // running over sink time index
  for(size_t t = 0; t < Lt; ++t){ 
    // running over all indices on one timeslice 
    // -> resorting propagator and copyying it into vec
    const size_t t_h = 12*Ls*Ls*Ls*t; // helper index
    for(size_t x = 0; x < Ls*Ls*Ls; ++x){ // spatial
      const size_t x_h = t_h + 12*x;  // helper index
      for(size_t d = 0; d < 4; ++d){      // Dirac
        const size_t d_h = x_h + 3*d; // helper index
        for(size_t c = 0; c < 3; ++c){    // colour
           vec(3*x + c, d) = propagator[d_h + c];
        }
      }
    }

    // running over different sinks - smeared and stochastic
    // TODO: Implementation of stochastic sinks!
    for(size_t sinks = 0; sinks < param.nb_sinks; sinks++){
      // CHANGED for easier comparison with existing code!!!!!!!!
      Eigen::MatrixXcd vec1 = Eigen::MatrixXcd::Zero(nb_ev, 4);
      for(size_t d = 0; d < 4; ++d)    
        vec1 = V[t].adjoint()*vec;
      for(size_t c = 0; c < nb_ev; ++c)  
        for(size_t d = 0; d < 4; ++d)    
          (perambulator[sinks][t])(4*c+d, col) = vec1(c, d);
    }
  } // end of loop over sink time index

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::write_perambulator_to_disk() {

  for(size_t sinks = 0; sinks < param.nb_sinks; sinks++){
    FILE *fp = NULL;
    const size_t Lt = param.Lt;
    const size_t verbose = param.verbose;
    const size_t number_of_eigen_vec = param.nb_ev;
    const size_t number_of_inversions = param.dilution_size_so[0] *
                                        param.dilution_size_so[1] *
                                        param.dilution_size_so[2] ;
    const size_t size_perambulator_entry = number_of_inversions * Lt * 4 * 
                                           number_of_eigen_vec;
  
    // memory for reading perambulators
    std::vector<cmplx> perambulator_write(size_perambulator_entry);
    // copy perambulator into writing array
    for(size_t t = 0; t < Lt; ++t){
      size_t t_h = t*4*number_of_eigen_vec; // helper index
      for(size_t row_i = 0; row_i < 4 * number_of_eigen_vec; ++row_i){
        size_t row_i_h = (row_i + t_h) * number_of_inversions; // helper index
        for(size_t col_i = 0; col_i < number_of_inversions; ++col_i){
          perambulator_write[row_i_h + col_i] = 
                perambulator[sinks][t](row_i, col_i);
        }
      }
    }
    // create perambulator file name
    std::string outfile = param.create_out_filename("perambulator", sinks);
    // writing data
    if((fp = fopen(outfile.c_str(), "wb")) == NULL){
      std::cout << "failed to open file: " << outfile << "\n" << std::endl;
      exit(0);
    }
    if(verbose) printf("\tread file: %s\n", outfile.c_str());
    int check_read_in = 
            fwrite(&(perambulator_write[0]), sizeof(cmplx),
                   size_perambulator_entry, fp);
    if(check_read_in != (int) size_perambulator_entry)
      std::cout << "It seems that not all data are written to: "
                << outfile.c_str() << "\n" << std::endl;
    fclose(fp);
  }
}
////////////////////////////////////////////////////////////////////////////////
/////////////////////// private methods ////////////////////////////////////////
void LapH::distillery::read_eigenvectors(){

  const size_t Ls = param.Ls;
  const size_t Lt = param.Lt;
  const size_t dim_row = Ls * Ls * Ls * 3;
  const size_t verbose = param.verbose;
  const size_t number_of_eigen_vec = param.nb_ev;

  //buffer for read in
  cmplx* eigen_vec = new cmplx[dim_row];

  if(verbose) printf("reading eigen vectors from files:\n");
  else printf("\treading eigenvectors\n");
  fflush(stdout);

  for(size_t t = 0; t < Lt; ++t){
    //setting up file
    char name[200];
    sprintf(name, "%s/eigenvectors.%04d.%03d", 
                  param.inpath_ev.c_str(), (int) param.config, (int) t);
    if(verbose) std::cout << "Reading file: " << name << std::endl;
    std::ifstream infile(name, std::ifstream::binary);
    if (infile) {
      for (size_t nev = 0; nev < number_of_eigen_vec; ++nev) {
        infile.read( (char*) eigen_vec, 2*dim_row*sizeof(double));
        for(size_t nrow = 0; nrow < dim_row; ++nrow)
          (V[t])(nrow, nev) = eigen_vec[nrow];
      }
    }
    else {
      std::cout << "eigenvector file does not exist!!!\n" << std::endl;
      exit(0);
    }
    infile.close();

    // small test of trace and sum over the eigen vector matrix!
    if(verbose){
      std::cout << "trace of V^d*V on t = " << t << ":\t"
          << (V[t].adjoint() * V[t]).trace() << std::endl;
      std::cout << "sum over all entries of V^d*V on t = " << t << ":\t"
          << (V[t].adjoint() * V[t]).sum() << std::endl;
    }
  }
  delete[] eigen_vec;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------












