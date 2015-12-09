#include "./distillery.h"


// -----------------------------------------------------------------------------
// swapping endianess ----------------------------------------------------------
template <typename T> T swap_endian(T u) {

  union {
    T u;
    unsigned char u8[sizeof(T)];
  } source, dest;

  source.u = u;

  for (size_t k = 0; k < sizeof(T); k++)
    dest.u8[k] = source.u8[sizeof(T) - k - 1];

  return dest.u;
}
// -----------------------------------------------------------------------------
inline std::complex<double> swap_complex(std::complex<double> val){
  return std::complex<double>(swap_endian<double>(std::real(val)),
                              swap_endian<double>(std::imag(val)));
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::initialise(const LapH::input_parameter& in_param) {

  param = in_param;
  const size_t Ls = param.Ls;
  const size_t Lt = param.Lt;
  const size_t nb_of_sinks = param.nb_of_sinks;

  // initialisation of mpi tmLQCD parameters
  tmLQCD_params = new tmLQCD_mpi_params;
  tmLQCD_get_mpi_params(tmLQCD_params);

  const size_t T = Lt/tmLQCD_params->nproc_t;
  const size_t dim_row = Ls * Ls * Ls * 3 / (tmLQCD_params->nproc_x * 
                                             tmLQCD_params->nproc_y * 
                                             tmLQCD_params->nproc_z);

  // Initializing memory for eigenvectors, perambulator and random vector
  V = new Eigen::MatrixXcd[T];
  for(size_t t = 0; t < T; ++t)
    V[t] = Eigen::MatrixXcd::Zero(dim_row, param.nb_ev);
  
  perambulator.resize(nb_of_sinks);
  for(size_t nbs = 0; nbs < nb_of_sinks; nbs++){
    perambulator[nbs].resize(T);
    // TODO: usually the sink is only smeared so the size is 4 time the number
    //       of eigen vectors. If the stochastik sink is applied only a block
    //       scheme in spatial coordinates is implied, the rest is fully
    //       dilutet!
    size_t sink_size = 4*param.nb_ev;
    if (!param.dilution_type_si[nbs][1].compare("B"))
      sink_size = 12*param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1];
    for(size_t t = 0; t < T; ++t){
      perambulator[nbs][t] = Eigen::MatrixXcd::Zero(sink_size,
                                                    param.dilution_size_so[0] * 
                                                    param.dilution_size_so[1] * 
                                                    param.dilution_size_so[2] );
    }
  }
  random_vector = new Eigen::VectorXcd[Lt];
  for(size_t t = 0; t < Lt; ++t)
    random_vector[t] = Eigen::VectorXcd::Zero(4 * param.nb_ev);


  // reading eigenvectors from disk
  read_eigenvectors(); 
  // generating random vector for source
  set_random_vector(param.rnd_id[0]);

  // generating random vectors for sink
  random_vector_si.resize(param.nb_of_sinks);
  const size_t volume = 12 * Ls * Ls * Ls * Lt / (tmLQCD_params->nproc_t *
                               tmLQCD_params->nproc_x * tmLQCD_params->nproc_y * 
                               tmLQCD_params->nproc_z);
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){
    if (param.dilution_type_si[nbs][1].compare("B"))
      continue;
    random_vector_si[nbs] = Eigen::VectorXcd::Zero(volume);
    set_sink_random_vector(param.rnd_id[0], nbs, random_vector_si[nbs]);
  }

  // is everything allocated?
  memory_allocated = true;
  
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::clean() {

  const size_t T = param.Lt/tmLQCD_params->nproc_t;

  // matrices are resized to size 0
  for(size_t t = 0; t < T; ++t){
    (V[t]).resize(0, 0);
    (random_vector[t]).resize(0);
  }
  for(size_t t = 0; t < param.Lt; ++t){
    (random_vector[t]).resize(0);
  }
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){
    for(size_t t = 0; t < T; ++t){
      (perambulator[nbs][t]).resize(0, 0);
    }
    perambulator[nbs].resize(0);
    (random_vector_si[nbs]).resize(0);
  }
  perambulator.resize(0);
  random_vector_si.resize(0);

  delete[] V;
  delete[] random_vector;
  V = NULL;
  random_vector = NULL;

  memory_allocated = false;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::reset_all(const LapH::input_parameter& in_param){

  param = in_param;
  const size_t Ls = param.Ls;
  const size_t Lt = param.Lt;

  const size_t T = Lt/tmLQCD_params->nproc_t;
  const size_t dim_row = Ls * Ls * Ls * 3 / (tmLQCD_params->nproc_x * 
                                             tmLQCD_params->nproc_y * 
                                             tmLQCD_params->nproc_z);
  const size_t volume = 12 * Ls * Ls * Ls * Lt / (tmLQCD_params->nproc_t *
                               tmLQCD_params->nproc_x * tmLQCD_params->nproc_y * 
                               tmLQCD_params->nproc_z);
  // eigen vectors
  for(size_t t = 0; t < T; ++t){
    V[t] = Eigen::MatrixXcd::Zero(dim_row, param.nb_ev);
  }
  read_eigenvectors();  
  // perambulator
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){
    // TODO: usually the sink is only smeared so the size is 4 time the number
    //       of eigen vectors. If the stochastik sink is applied only a block
    //       scheme in spatial coordinates is implied, the rest is fully
    //       dilutet!
    size_t sink_size = 4*param.nb_ev;
    if (!param.dilution_type_si[nbs][1].compare("B"))
      sink_size = 12*param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1];
    for(size_t t = 0; t < T; ++t)
      perambulator[nbs][t] = Eigen::MatrixXcd::Zero(sink_size,
                                                    param.dilution_size_so[0] * 
                                                    param.dilution_size_so[1] * 
                                                    param.dilution_size_so[2] );
  }
  // source random vector
  for(size_t t = 0; t < Lt; ++t)
    random_vector[t] = Eigen::VectorXcd::Zero(4 * param.nb_ev);
  set_random_vector(param.rnd_id[0]);
  // sink random vectors
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){
    if (param.dilution_type_si[nbs][1].compare("B"))
      continue;
    random_vector_si[nbs] = Eigen::VectorXcd::Zero(volume);
    set_sink_random_vector(param.rnd_id[0], nbs, random_vector_si[nbs]);
  }

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::reset_perambulator_and_randomvector(const size_t rnd_id){

  const size_t Ls = param.Ls;
  const size_t Lt = param.Lt;
  const size_t T = param.Lt/tmLQCD_params->nproc_t;
  const size_t volume = 12 * Ls * Ls * Ls * Lt / (tmLQCD_params->nproc_t *
                               tmLQCD_params->nproc_x * tmLQCD_params->nproc_y * 
                               tmLQCD_params->nproc_z);

  // setting perambulator and random vector to zero 
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){
    // TODO: usually the sink is only smeared so the size is 4 time the number
    //       of eigen vectors. If the stochastik sink is applied only a block
    //       scheme in spatial coordinates is implied, the rest is fully
    //       dilutet!
    size_t sink_size = 4*param.nb_ev;
    if (!param.dilution_type_si[nbs][1].compare("B"))
      sink_size = 12*param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1];
    for(size_t t = 0; t < T; ++t)
      perambulator[nbs][t] = Eigen::MatrixXcd::Zero(sink_size,
                                                    param.dilution_size_so[0] * 
                                                    param.dilution_size_so[1] * 
                                                    param.dilution_size_so[2] );
  }
  for(size_t t = 0; t < param.Lt; ++t)
    random_vector[t] = Eigen::VectorXcd::Zero(4 * param.nb_ev);

  // resetting random seed and generating new random vector
  set_random_vector(rnd_id); // source
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){ // sinks
    if (param.dilution_type_si[nbs][1].compare("B"))
      continue;
    random_vector_si[nbs] = Eigen::VectorXcd::Zero(volume);
    set_sink_random_vector(param.rnd_id[rnd_id], nbs, random_vector_si[nbs]);
  }

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
                                   const std::string type, 
                                   size_t* const lookup){
  if(type == "B" || type == "N")
    for(size_t i = 0; i < nb_of_nonzero_entries; ++i)
      lookup[i] = nb_of_nonzero_entries * dilution_entry + i;
  else if(type == "I" || type == "F")
    for(size_t i = 0; i < nb_of_nonzero_entries; ++i)
      lookup[i] = dilution_entry + i * nb_inversions;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::create_source(std::complex<double>** source) {

  MPI_Barrier(MPI_COMM_WORLD); 
  double time1 = MPI_Wtime();

  const size_t Lt = param.Lt;
  const size_t Ls = param.Ls;
  const size_t T = Lt/tmLQCD_params->nproc_t;
  const size_t X = Ls/tmLQCD_params->nproc_x;
  const size_t Y = Ls/tmLQCD_params->nproc_y;
  const size_t Z = Ls/tmLQCD_params->nproc_z;
  const size_t Vs = X * Y * Z;
  const size_t dim_row = Vs * 3;
  const size_t length = 4*T*dim_row;
  const size_t number_of_eigen_vec = param.nb_ev;
  const size_t nb_of_inversions = param.dilution_size_so[0] * 
                                  param.dilution_size_so[1] *
                                  param.dilution_size_so[2];
  // setting source to zero 
  for(size_t inversion = 0; inversion < nb_of_inversions; inversion++)
    for(size_t i = 0; i < length; i++)
      source[inversion][i] = std::complex<double>(0.0, 0.0);

  // indices of timeslices with non-zero entries
  size_t nb_of_nonzero_t = Lt/param.dilution_size_so[0];  
  // TODO: think about smart pointer here!
  size_t* t_index = new size_t[nb_of_nonzero_t];

  // indices of eigenvectors to be combined
  size_t nb_of_ev_combined = number_of_eigen_vec/param.dilution_size_so[1];
  size_t* ev_index = new size_t[nb_of_ev_combined];

  // indices of Dirac components to be combined 
  // TODO: This is needed only once - could be part of class
  size_t nb_of_dirac_combined = 4/param.dilution_size_so[2];
  size_t** d_index = new size_t*[param.dilution_size_so[2]];
  for(size_t i = 0; i < param.dilution_size_so[2]; ++i)
    d_index[i] = new size_t[nb_of_dirac_combined];

  for(size_t dil_t = 0; dil_t < param.dilution_size_so[0]; ++dil_t){
    // lookup in time
    create_dilution_lookup(nb_of_nonzero_t, param.dilution_size_so[0], 
                             dil_t, param.dilution_type_so[0], t_index);
    for(size_t dil_e = 0; dil_e < param.dilution_size_so[1]; ++dil_e){    

      // lookup in eigenspace
      create_dilution_lookup(nb_of_ev_combined, param.dilution_size_so[1], 
                             dil_e, param.dilution_type_so[1], ev_index);
      // lookup in Diracspace
      for(size_t dil_d = 0; dil_d < param.dilution_size_so[2]; ++dil_d)
        create_dilution_lookup(nb_of_dirac_combined, param.dilution_size_so[2],
                             dil_d, param.dilution_type_so[2], d_index[dil_d]); 

      size_t inversion_number = (dil_t * param.dilution_size_so[1] + dil_e) 
                                                  * param.dilution_size_so[2];
      // creating the source ---------------------------------------------------
      for(size_t t = 0; t < nb_of_nonzero_t; ++t){ 

        size_t time = t_index[t];                         
        size_t t_start = T*tmLQCD_params->proc_coords[0]; // for if statement
        if((t_start <= time) && (time < (t_start + T))) {

          // intermediate memory for source
          Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(dim_row, 4);
          size_t time_id = time % T; 

          // building source on one timeslice
          for(size_t ev = 0; ev < nb_of_ev_combined; ++ev){
            size_t ev_h = ev_index[ev] * 4; // helper index 
            for(size_t d = 0; d < 4; ++d){
              S.col(d) += random_vector[time](ev_h+d) *
                              (V[time_id]).col(ev_index[ev]); 
            }
          }

          // copy the created source into output array
          size_t t_h = time_id*Vs; // helper index
          for(size_t x = 0; x < Vs; ++x){
            size_t x_h  = (t_h + x)*12;    // helper index
            size_t x_h2 = x*3;            // helper index
            for(size_t d2 = 0; d2 < param.dilution_size_so[2]; ++d2){
              for(size_t d3 = 0; d3 < nb_of_dirac_combined; ++d3){
                size_t d_h = x_h + d_index[d2][d3]*3; // helper index
                for(size_t c = 0; c < 3; ++c){
                  source[inversion_number + d2][d_h + c] = 
                                            S(x_h2 + c, d_index[d2][d3]);
          }}}}
        } // end of if-statement
      }   // end of loop over nonzero timeslices
  }}    // end of loops over eigenvector and Diracdilution
  // freeing memory
  delete[] t_index;
  delete[] ev_index;
  for(size_t i = 0; i < param.dilution_size_so[2]; ++i)
    delete[] d_index[i];
  delete[] d_index;
  t_index = NULL;
  ev_index = NULL;
  d_index = NULL;

  MPI_Barrier(MPI_COMM_WORLD);
  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  time1 = MPI_Wtime() - time1;
  if(myid == 0)
    std::cout << "\tTime for source creation: " << time1 << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::add_to_perambulator(
                         const std::complex<double>* const * const propagator) { 
  
  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Barrier(MPI_COMM_WORLD); 
  double time1 = MPI_Wtime();
  Eigen::setNbThreads(1);
  if(myid == 0){
    std::cout << "\tBuilding Perambulator" << std::endl;
    std::cout << "\tEigen will use " <<  Eigen::nbThreads() 
              << " threads!" << std::endl;
  }

  const size_t Lt = param.Lt;
  const size_t Ls = param.Ls;
  const size_t T = Lt/tmLQCD_params->nproc_t;
  const size_t X = Ls/tmLQCD_params->nproc_x;
  const size_t Y = Ls/tmLQCD_params->nproc_y;
  const size_t Z = Ls/tmLQCD_params->nproc_z;
  const size_t Vs = X * Y * Z;
  const size_t dim_row = Vs * 3;
  const size_t nb_ev = param.nb_ev;
  const size_t nb_of_inversions = param.dilution_size_so[0] * 
                                  param.dilution_size_so[1] *
                                  param.dilution_size_so[2];
  const size_t px = tmLQCD_params->proc_coords[1];
  const size_t py = tmLQCD_params->proc_coords[2];
  const size_t pz = tmLQCD_params->proc_coords[3];


  // running over all sinks
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){
    // checking if smeared or stochastic sink must be computed
    if (!param.dilution_type_si[nbs][1].compare("F")){ // smeared sink
      // running over sink time index
      #pragma omp parallel for
      for(size_t t = 0; t < T; ++t){ 
        Eigen::MatrixXcd vec = Eigen::MatrixXcd::Zero(dim_row, 
                                                            4*nb_of_inversions);
        // running over inversions
        for(size_t col = 0; col < nb_of_inversions; col++){
          // running over all indices on one timeslice 
          // -> resorting propagator and copying it into vec
          size_t t_h = 12*Vs*t; // helper index
          for(size_t x = 0; x < Vs; ++x){ // spatial
            size_t x_h = t_h + 12*x;  // helper index
            for(size_t d = 0; d < 4; ++d){      // Dirac
              size_t d_h = x_h + 3*d; // helper index
              for(size_t c = 0; c < 3; ++c){    // colour
                 vec(3*x + c, col*4 + d) = propagator[col][d_h + c];
          }}}
        }
    
        // multiplication with V^dagger and storing result in perambulator
        Eigen::MatrixXcd vec1 = Eigen::MatrixXcd::Zero(nb_ev, 
                                                            4*nb_of_inversions);
        vec1 = ((V[t]).adjoint())*vec;
        for(size_t col = 0; col < nb_of_inversions; col++)
          for(size_t ev = 0; ev < nb_ev; ++ev)  
            for(size_t d = 0; d < 4; ++d)    
              (perambulator[nbs][t])(4*ev+d, col) = vec1(ev, 4*col+d);
    
      } // end of loop over sink time index
    }
    else { // stochastic sink
      for(size_t t = 0; t < T; ++t){ 
        // resorting propagator to make the dilution easier
        Eigen::MatrixXcd M = 
               Eigen::MatrixXcd::Zero(Vs, 12*nb_of_inversions);
        // running over inversions
        for(size_t col = 0; col < nb_of_inversions; col++){
          // running over all indices on one timeslice 
          // -> resorting propagator and copying it into vec
          size_t t_h = 12*Vs*t; // helper index
          for(size_t x = 0; x < Vs; ++x){ // spatial
            size_t x_h = t_h + 12*x;  // helper index
            for(size_t d = 0; d < 4; ++d){      // Dirac
              size_t d_h = x_h + 3*d; // helper index
              for(size_t c = 0; c < 3; ++c){    // colour
                M(x, 12*col + 3*d + c) = propagator[col][d_h + c];
          }}}
        }

        // building the diluted random vector matrix
        size_t dil_size = param.dilution_size_si[nbs][1];
        Eigen::MatrixXcd R = Eigen::MatrixXcd::Zero(Vs, 
                                                 12*dil_size*dil_size*dil_size);
        for(size_t x = 0; x < X; x++){
          for(size_t y = 0; y < Y; y++){
            for(size_t z = 0; z < Z; z++){
              size_t xx = X*px + x;
              size_t x_quadrant = 0;
              size_t x_h = Ls/dil_size;
              for(size_t x1 = 0; x1 < dil_size; x1++){
                if(xx >= x1*x_h && xx < (x1+1)*x_h){
                  x_quadrant = x1;
                  break;
                }
              }
              size_t yy = Y*py + y;
              size_t y_quadrant = 0;
              size_t y_h = Ls/dil_size;
              for(size_t y1 = 0; y1 < dil_size; y1++){
                if(yy >= y1*y_h && yy < (y1+1)*y_h){
                  y_quadrant = y1;
                  break;
                }
              }
              size_t zz = Z*pz + z;
              size_t z_quadrant = 0;
              size_t z_h = Ls/dil_size;
              for(size_t z1 = 0; z1 < dil_size; z1++){
                if(zz >= z1*z_h && zz < (z1+1)*z_h){
                  z_quadrant = z1;
                  break;
                }
              }
              size_t col = (x_quadrant*dil_size + y_quadrant) * dil_size + 
                                                                     z_quadrant;
              size_t row = x*Y*Z + y*Z + z;
              R.block(row, 12*col, 1, 12) = 
                 (random_vector_si[nbs].segment(12*Vs*t + row, 12)).transpose();
            }
          }
        }
        Eigen::MatrixXcd res = R.adjoint()*M;
        // resorting res into correct propagator order
        for(size_t col = 0; col < nb_of_inversions; col++)
          for(size_t ev = 0; ev < dil_size*dil_size*dil_size; ++ev)  
            for(size_t d = 0; d < 4; ++d)    
              for(size_t c = 0; c < 3; ++c)    
                (perambulator[nbs][t])(12*ev+3*d+c, col) = 
                                                          res(ev, 12*col+3*d+c);

      } // end of loop over sink time index
    } // end of if else desision of sink type
  } // end of loop over sinks

  MPI_Barrier(MPI_COMM_WORLD);
  time1 = MPI_Wtime() - time1;
  if(myid == 0)
    std::cout << "\tTime for perambulator creation: " 
              << time1 << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::write_perambulator_to_disk(const size_t rnd_id) {

  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Barrier(MPI_COMM_WORLD); 
  double time1 = MPI_Wtime();

  char outfile[400];
  FILE *fp = NULL;
  const size_t Lt = param.Lt;
  const size_t T = Lt/tmLQCD_params->nproc_t;
  const size_t verbose = param.verbose;
  const size_t number_of_inversions = param.dilution_size_so[0] *
                                      param.dilution_size_so[1] *
                                      param.dilution_size_so[2] ;


  // loop over all sinks
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){
    const size_t number_of_eigen_vec = param.nb_ev;
    size_t sink_size = 4*number_of_eigen_vec;
    if (!param.dilution_type_si[nbs][1].compare("B"))
      sink_size = 12*param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1];

    const size_t size_perambulator_entry = number_of_inversions*Lt*sink_size;

    // memory for reading perambulators and setting it to zero
    std::vector<std::complex<double> > perambulator_write(
                       size_perambulator_entry, std::complex<double>(0.0, 0.0));

    // copy perambulator into writing array
    for(size_t t = 0; t < T; ++t){
      size_t t_global = T*tmLQCD_params->proc_coords[0] + t;
      size_t t_h = t_global*sink_size; // helper index
      for(size_t row_i = 0; row_i < sink_size; ++row_i){
        size_t row_i_h = (row_i + t_h) * number_of_inversions; // helper index
        for(size_t col_i = 0; col_i < number_of_inversions; ++col_i){
          perambulator_write[row_i_h + col_i] = 
                                             perambulator[nbs][t](row_i, col_i);
        }
      }
    }
    // Summing all the perambulator data on process 0
    if(myid == 0)
      MPI_Reduce(MPI_IN_PLACE, &(perambulator_write[0]), 
                 2*size_perambulator_entry, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    else 
      MPI_Reduce(&(perambulator_write[0]), 0,
                 2*size_perambulator_entry, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
    // process 0 writes the data to disk
    if(myid == 0){
      // changing endianess
      for (size_t i = 0; i < size_perambulator_entry; i++)
        //perambulator_write[i] = swap_complex(perambulator_write[i]);
        perambulator_write[i] = perambulator_write[i];
      // data path
      std::string filename = param.outpath + "/perambulator";
      if(verbose) printf("writing perambulators from files:\n");
      else printf("\twriting perambulator\n");
      // create perambulator file name
      sprintf(outfile,
          "%s.rndvecnb%02d.%s.Tso%s%04d.Vso%s%04d.Dso%s%01d.Tsi%s%04d.Ssi%s%04d" \
          ".Dsi%s%01d.Csi%s%01d.smeared%01d.%05d",
          filename.c_str(), (int) param.rnd_id[rnd_id], param.quarktype.c_str(), 
          param.dilution_type_so[0].c_str(), (int) param.dilution_size_so[0], 
          param.dilution_type_so[1].c_str(), (int) param.dilution_size_so[1], 
          param.dilution_type_so[2].c_str(), (int) param.dilution_size_so[2],
          param.dilution_type_si[nbs][0].c_str(), (int) param.dilution_size_si[nbs][0], 
          param.dilution_type_si[nbs][1].c_str(), (int) param.dilution_size_si[nbs][1], 
          param.dilution_type_si[nbs][2].c_str(), (int) param.dilution_size_si[nbs][2], 
          param.dilution_type_si[nbs][3].c_str(), (int) param.dilution_size_si[nbs][3], 
          1, (int) param.config);
     
      // writing data
      if((fp = fopen(outfile, "wb")) == NULL){
        std::cout << "failed to open file: " << outfile << "\n" << std::endl;
        exit(0);
      }
      if(verbose) printf("\tread file: %s\n", outfile);
      fwrite(&(perambulator_write[0]), sizeof(std::complex<double>),
          size_perambulator_entry, fp);
      fclose(fp);
     
    } // if for process 0 ends here
  } // for loop over sinks ends here
  MPI_Barrier(MPI_COMM_WORLD);
  time1 = MPI_Wtime() - time1;
  if(myid == 0)
    std::cout << "\tTime for perambulator writing: " 
              << time1 << std::endl;
}
////////////////////////////////////////////////////////////////////////////////
/////////////////////// private methods ////////////////////////////////////////
void LapH::distillery::set_random_vector(const size_t rnd_id) {

  double sqrt2 = 0.5*sqrt(2.0);
  double re, im;
  rlxs_init(0, (param.seed[rnd_id])^param.config);
  int rnd_length = 2*param.Lt*param.nb_ev*4;
  float* rnd = new float[rnd_length];
  ranlxs(rnd, rnd_length);

  // generating a Z_2 source
  for(size_t t = 0; t < param.Lt; ++t ){ 
    random_vector[t].setRandom(4 * param.nb_ev);
    for(size_t i = 0; i < 4 * param.nb_ev; ++i){
      if (rnd[8 * param.nb_ev * t + 2*i] < 0.5)
        re = -sqrt2;
      else 
        re = sqrt2;
      if (rnd[8 * param.nb_ev * t + 2*i + 1] < 0.5)
        im = -sqrt2;
      else 
        im = sqrt2;
      random_vector[t](i) = std::complex<double>(re, im);
    }
  }

  // writing random vector to disc
  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if(myid == 0)
    write_random_vector_to_disk(rnd_id);

  delete[] rnd;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::set_sink_random_vector(const size_t rnd_id, 
                                  const size_t sink_id, Eigen::VectorXcd& out) {

  rlxd_init(2, ((param.seed[rnd_id])^param.config)^(param.seed_si[sink_id]));
  // Reinterpret cast is NOT good style but necessary here to connect the data
  // here with the structs of tmLQCD!
  random_spinor_field_lexic(reinterpret_cast<spinor*>(&out[0]), 1, RN_Z2);

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::read_random_vector() {

  char infile[400];
  FILE *fp = NULL;
  const int Lt = param.Lt;
  const int verbose = param.verbose;
  const int number_of_rnd_vec = 1;
  const int number_of_eigen_vec = param.nb_ev;
  const int rnd_vec_length = number_of_eigen_vec * 4;
  // memory for reading random vectors
  std::complex<double>* rnd_vec_read =
      new std::complex<double>[Lt*rnd_vec_length];
  std::string filename = "/data/LapHs/peram_gen/new/test/randomvector";

  if(verbose) printf("reading random vectors from files:\n");
  else printf("\treading random vectors:");
  int check_read_in = 0;

  for(int rnd_vec_i = 0; rnd_vec_i < number_of_rnd_vec; ++rnd_vec_i){
    // read random vector
    sprintf(infile, "%s.%03d.u.Tf.%04d", filename.c_str(), 0, 1000);
    if(verbose) printf("\tread file: %s\n", infile);
    if((fp = fopen(infile, "rb")) == NULL){
      std::cout << "failed to open file: " << infile << "\n" << std::endl;
      exit(0);
    }   
    check_read_in += fread(rnd_vec_read, sizeof(std::complex<double>),
        Lt*rnd_vec_length, fp);
    // copy into matrix structure
    for(int t = 0; t < Lt; t++){
      for(int row_i = 0; row_i < rnd_vec_length; ++row_i){
        random_vector[t](row_i) = rnd_vec_read[t * rnd_vec_length + row_i];
      }   
    }
    fclose(fp);
  }   
  delete[] rnd_vec_read;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::write_random_vector_to_disk(size_t rnd_id){

  char outfile[400];
  FILE *fp = NULL;
  const size_t Lt = param.Lt;
  const size_t verbose = param.verbose;
  const size_t number_of_eigen_vec = param.nb_ev;
  const size_t rnd_vec_length = Lt * number_of_eigen_vec * 4;

  // copy from eigen structure into intermediate memory
  std::complex<double>* rnd_vec_write =
      new std::complex<double>[rnd_vec_length];
  for(size_t t = 0; t < Lt; ++t)
    for(size_t row_i = 0; row_i < 4 * number_of_eigen_vec; ++row_i)
      //rnd_vec_write[row_i + t * rnd_vec_length/Lt] = 
      //                                    swap_complex(random_vector[t](row_i));
      rnd_vec_write[row_i + t * rnd_vec_length/Lt] = 
                                          random_vector[t](row_i);

  // creating name and path of outfile
  std::string filename = param.outpath + "/randomvector";
  if(verbose) printf("writing random vector to files:\n");
  else printf("\twriting random vector\n");
  int check_read_in = 0;
  // TODO: Must be changed to write stochastic sink random vectors
  sprintf(outfile, "%s.rndvecnb%02d.%s.nbev%04d.%04d", filename.c_str(), 
                    (int) param.rnd_id[rnd_id], param.quarktype.c_str(), 
                    (int) param.nb_ev, (int) param.config);
  if(verbose) printf("\twrite to file: %s\n", outfile);
  if((fp = fopen(outfile, "wb")) == NULL){
    std::cout << "failed to open file to write random vector: " 
              << outfile << "\n" << std::endl;
    exit(0);
  }   
  check_read_in += fwrite(rnd_vec_write, sizeof(std::complex<double>),
                          rnd_vec_length, fp);
  if(check_read_in != (int) rnd_vec_length){
    std::cout << "failed to write random vector: "
              << outfile << "\n" << std::endl;
    exit(0);
  }   
  fclose(fp);
  // delete intermediate memory
  delete[] rnd_vec_write;

} 
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::copy_to_V(const std::complex<double>* const eigen_vec, 
                                 const int t, const int nev){
  const size_t Ls = param.Ls;
  const size_t X = Ls/tmLQCD_params->nproc_x;
  const size_t Y = Ls/tmLQCD_params->nproc_y;
  const size_t Z = Ls/tmLQCD_params->nproc_z;
  const size_t nproc_y = tmLQCD_params->nproc_y;
  const size_t nproc_z = tmLQCD_params->nproc_z;
  const size_t px = tmLQCD_params->proc_coords[1];
  const size_t py = tmLQCD_params->proc_coords[2];
  const size_t pz = tmLQCD_params->proc_coords[3];

  size_t j = 0;    
  for(size_t x = 0; x < X; x++){
    for(size_t y = 0; y < Y; y++){
      for(size_t z = 0; z < Z; z++){
        for(size_t c = 0; c < 3; c++){

          size_t i = 3*( ((X*px + x)*Y*nproc_y + (Y*py + y))*Z*nproc_z  
                                                             +  Z*pz + z) + c;
          // byte swap to change endianess
          //(V[t])(j, nev) = swap_complex(eigen_vec[i]);
          (V[t])(j, nev) = eigen_vec[i];
          j++;

        }
      }
    }
  }

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::read_eigenvectors(){

  MPI_Barrier(MPI_COMM_WORLD); 
  double time1 = MPI_Wtime();

  const size_t Ls = param.Ls;
  const size_t Lt = param.Lt;
  const size_t verbose = param.verbose;
  const size_t number_of_eigen_vec = param.nb_ev;

  const size_t T = Lt/tmLQCD_params->nproc_t;
  const size_t dim_row = Ls * Ls * Ls * 3;

  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  //buffer for read in
  std::complex<double>* eigen_vec = new std::complex<double>[dim_row];

  if(myid == 0){
    if(verbose) printf("reading eigen vectors from files:\n");
    else printf("\treading eigenvectors\n");
    fflush(stdout);
  }

  // running over all timeslices on this process
  for(size_t t = 0; t < T; t++){

    const int real_t = T*tmLQCD_params->proc_coords[0] + t;
 
    //setting up file
    char name[200];
    sprintf(name, "%s/eigenvectors.%04d.%03d", 
                  param.inpath_ev.c_str(), (int) param.config, real_t);
    if(verbose) std::cout << "Reading file: " << name << std::endl;
    std::ifstream infile(name, std::ifstream::binary);
    if (infile) {
      for (size_t nev = 0; nev < number_of_eigen_vec; ++nev) {
        // reading the full vector
        infile.read( (char*) eigen_vec, 2*dim_row*sizeof(double));
        // copying the correct components into V
        copy_to_V(eigen_vec, t, nev);
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

  MPI_Barrier(MPI_COMM_WORLD);
  time1 = MPI_Wtime() - time1;
  if(myid == 0)
    std::cout << "\tTime for eigenvector reading: " << time1 << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------












