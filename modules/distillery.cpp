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
    perambulator[nbs].resize(param.nb_of_sink_rnd_vec[nbs]);
    for(size_t nbr = 0; nbr < param.nb_of_sink_rnd_vec[nbs]; nbr++){
      perambulator[nbs][nbr].resize(T);
      // TODO: usually the sink is only smeared so the size is 4 time the number
      //       of eigen vectors. If the stochastik sink is applied only a block
      //       scheme in spatial coordinates is implied, the rest is fully
      //       dilutet!
      size_t sink_size = 4*param.nb_ev;
      if (param.dilution_type_si[nbs][1].compare("F"))
        sink_size = 12*param.dilution_size_si[nbs][1] *
                       param.dilution_size_si[nbs][1] *
                       param.dilution_size_si[nbs][1];
      for(size_t t = 0; t < T; ++t){
        perambulator[nbs][nbr][t] = Eigen::MatrixXcd::Zero(sink_size,
                                                    param.dilution_size_so[0] * 
                                                    param.dilution_size_so[1] * 
                                                    param.dilution_size_so[2] );
      }
    }
  }
  random_vector = new Eigen::VectorXcd[Lt];
  for(size_t t = 0; t < Lt; ++t)
    random_vector[t] = Eigen::VectorXcd::Zero(4 * param.nb_ev);


  // reading eigenvectors from disk
  read_eigenvectors(); 
  // generating random vector for source
  set_random_vector(0);

  // generating random vectors for sink
  random_vector_si.resize(param.nb_of_sinks);
  const size_t volume = 12 * Ls * Ls * Ls * Lt / (tmLQCD_params->nproc_t *
                               tmLQCD_params->nproc_x * tmLQCD_params->nproc_y *
                               tmLQCD_params->nproc_z);
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){
    if (!param.dilution_type_si[nbs][1].compare("F"))
      continue;
    random_vector_si[nbs].resize(param.nb_of_sink_rnd_vec[nbs]);
    for(size_t nbr = 0; nbr < param.nb_of_sink_rnd_vec[nbs]; nbr++){
      random_vector_si[nbs][nbr] = Eigen::VectorXcd::Zero(volume);
      set_sink_random_vector(0, nbs, nbr, random_vector_si[nbs][nbr]);
    }
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
    for(size_t nbr = 0; nbr < param.nb_of_sink_rnd_vec[nbs]; nbr++){
      for(size_t t = 0; t < T; ++t){
        (perambulator[nbs][nbr][t]).resize(0, 0);
      }
      (perambulator[nbs][nbr]).resize(0);
      if (!param.dilution_type_si[nbs][1].compare("F"))
        continue;
      (random_vector_si[nbs][nbr]).resize(0);
    }
    (perambulator[nbs]).resize(0);
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
    size_t sink_size = 4*param.nb_ev;
    if (param.dilution_type_si[nbs][1].compare("F"))
      sink_size = 12*param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1];
    for(size_t nbr = 0; nbr < param.nb_of_sink_rnd_vec[nbs]; nbr++){
      // TODO: usually the sink is only smeared so the size is 4 time the number
      //       of eigen vectors. If the stochastik sink is applied only a block
      //       scheme in spatial coordinates is implied, the rest is fully
      //       dilutet!
      for(size_t t = 0; t < T; ++t)
        perambulator[nbs][nbr][t] = Eigen::MatrixXcd::Zero(sink_size,
                                                    param.dilution_size_so[0] * 
                                                    param.dilution_size_so[1] * 
                                                    param.dilution_size_so[2] );
    }
  }
  // source random vector
  for(size_t t = 0; t < Lt; ++t)
    random_vector[t] = Eigen::VectorXcd::Zero(4 * param.nb_ev);
  set_random_vector(0);
  // sink random vectors
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){
    if (!param.dilution_type_si[nbs][1].compare("F"))
      continue;
    for(size_t nbr = 0; nbr < param.nb_of_sink_rnd_vec[nbs]; nbr++){
      random_vector_si[nbs][nbr] = Eigen::VectorXcd::Zero(volume);
      set_sink_random_vector(0, nbs, nbr, random_vector_si[nbs][nbr]);
    }
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
    if (param.dilution_type_si[nbs][1].compare("F"))
      sink_size = 12*param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1] *
                     param.dilution_size_si[nbs][1];
    for(size_t nbr = 0; nbr < param.nb_of_sink_rnd_vec[nbs]; nbr++){
      for(size_t t = 0; t < T; ++t)
        perambulator[nbs][nbr][t] = Eigen::MatrixXcd::Zero(sink_size,
                                                    param.dilution_size_so[0] * 
                                                    param.dilution_size_so[1] * 
                                                    param.dilution_size_so[2] );
    }
  }
  for(size_t t = 0; t < param.Lt; ++t)
    random_vector[t] = Eigen::VectorXcd::Zero(4 * param.nb_ev);

  // resetting random seed and generating new random vector
  set_random_vector(rnd_id); // source
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){ // sinks
    if (!param.dilution_type_si[nbs][1].compare("F"))
      continue;
    for(size_t nbr = 0; nbr < param.nb_of_sink_rnd_vec[nbs]; nbr++){
      random_vector_si[nbs][nbr] = Eigen::VectorXcd::Zero(volume);
      set_sink_random_vector(rnd_id, nbs, nbr, 
                             random_vector_si[nbs][nbr]);
    }
  }

}
// -----------------------------------------------------------------------------
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
static size_t create_sink_dilution_index(const size_t xx, const size_t Ls, 
                                         const size_t size, 
                                         const std::string type){
  if(type == "B"){
    size_t x_h = Ls/size;
    for(size_t x1 = 0; x1 < size; x1++)
      if(xx >= x1*x_h && xx < (x1+1)*x_h)
        return x1;
  }
  else if(type == "I"){
    for(size_t x1 = 0; x1 < size; x1++){
      if(xx%size == x1){
        return x1;
      }
    }
  }
  return 0;
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
      for(size_t t = 0; t < T; ++t){ 
        std::vector<std::complex<double> > vec(dim_row*4*nb_of_inversions);
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
                 vec.at(4*nb_of_inversions*(3*x + c) + col*4 + d) = 
                                                       propagator[col][d_h + c];
          }}}
        }
    
        // multiplication with V^dagger and storing result in perambulator
        std::vector<std::complex<double> > vec1(nb_ev*4*nb_of_inversions);
        std::complex<double> zero(0.0, 0.0);
        std::complex<double> one(1.0, 0.0);
        std::string type = "C";
        const int M = nb_ev;
        const int N = 4*nb_of_inversions;
        const int K = dim_row;
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, 
                                           Eigen::RowMajor> Vt = V[t].adjoint();
        zgemm_(&type[0], &type[0], &M, &N, &K, &one, Vt.data(), 
               &K, &vec[0], &N, &zero, &vec1[0], &M);

        for(size_t col = 0; col < nb_of_inversions; col++)
          for(size_t ev = 0; ev < nb_ev; ++ev)  
            for(size_t d = 0; d < 4; ++d)    
              (perambulator[nbs][0][t])(4*ev+d, col) = 
                                       std::conj(vec1.at(ev + nb_ev*(4*col+d)));
    
      } // end of loop over sink time index
    }
    else { // stochastic sink
      for(size_t nbr = 0; nbr < param.nb_of_sink_rnd_vec[nbs]; nbr++){
        for(size_t t = 0; t < T; ++t){ 
          // resorting propagator to make the dilution easier
          std::vector<std::complex<double> > vec(Vs*12*nb_of_inversions);
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
                  vec.at(nb_of_inversions*(12*x + 3*d + c) + col) = 
                                                       propagator[col][d_h + c];
            }}}
          }
          // building the diluted random vector matrix
          size_t dil_size = param.dilution_size_si[nbs][1];
          std::string dil_type = param.dilution_type_si[nbs][1];
          Eigen::MatrixXcd R = Eigen::MatrixXcd::Zero(12*Vs, 
                                                 12*dil_size*dil_size*dil_size);
          for(size_t x = 0; x < X; x++){
            for(size_t y = 0; y < Y; y++){
              for(size_t z = 0; z < Z; z++){
  
                size_t x_glb = X*px + x;
                size_t x_pos = create_sink_dilution_index(
                                                 x_glb, Ls, dil_size, dil_type);
                size_t y_glb = Y*py + y;
                size_t y_pos = create_sink_dilution_index(
                                                 y_glb, Ls, dil_size, dil_type);
                size_t z_glb = Z*pz + z;
                size_t z_pos = create_sink_dilution_index(
                                                 z_glb, Ls, dil_size, dil_type);
                size_t col = (x_pos*dil_size + y_pos) * dil_size + z_pos;
                size_t row = 12*(x*Y*Z + y*Z + z);
                R.block(row, 12*col, 12, 12) = 
                    (random_vector_si[nbs][nbr].segment(12*Vs*t + row, 12)).
                                                                   asDiagonal();
              }
            }
          }
          // multiplication with R^dagger and storing result in perambulator
          std::vector<std::complex<double> > 
                           vec1(12*dil_size*dil_size*dil_size*nb_of_inversions);
          std::complex<double> zero(0.0, 0.0);
          std::complex<double> one(1.0, 0.0);
          std::string type = "C";
          const int M = 12*dil_size*dil_size*dil_size;
          const int N = nb_of_inversions;
          const int K = 12*Vs;
          Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, 
                                              Eigen::RowMajor> Rt = R.adjoint();
          zgemm_(&type[0], &type[0], &M, &N, &K, &one, Rt.data(), &K, &vec[0], 
                                                       &N, &zero, &vec1[0], &M);
          // resorting res into correct propagator order
          for(size_t col = 0; col < nb_of_inversions; col++)
            for(size_t ev = 0; ev < 12*dil_size*dil_size*dil_size; ++ev)  
              (perambulator[nbs][nbr][t])(ev, col) = 
                   std::conj(vec1.at(ev + dil_size*dil_size*dil_size * 12*col));
        } // end of loop over sink time index
      } // end of loop over random sinks
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

  int myid = 0, nb_procs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
  MPI_Barrier(MPI_COMM_WORLD); 
  double time1 = MPI_Wtime();

  char outfile[400];
  const size_t Lt = param.Lt;
  const size_t T = Lt/tmLQCD_params->nproc_t;

  // create communicator for processes sharing the same time slices ------------
  MPI_Comm ts_comm;
  MPI_Comm_split(MPI_COMM_WORLD, tmLQCD_params->proc_coords[0], myid, &ts_comm);
  int myid_t = 0;
  MPI_Comm_rank(ts_comm, &myid_t);
  // create communicator for writing the data in parallel
  MPI_Comm write_comm;
  MPI_Comm_split(MPI_COMM_WORLD, myid_t, myid, &write_comm);
  int myid_write = 0;
  MPI_Comm_rank(write_comm, &myid_write);


  // loops over all sinksa -----------------------------------------------------
  for(size_t nbs = 0; nbs < param.nb_of_sinks; nbs++){
    for(size_t nbr = 0; nbr < param.nb_of_sink_rnd_vec[nbs]; nbr++){

      // memory for writing - The resizing is necesseray because only the 
      // relevant processes do have the correct memory
      std::vector<std::complex<double> > write_peram;
      if(myid_t == 0)
        write_peram.resize(perambulator[nbs][nbr][0].size()*T);

      // Summing all perambulator data -----------------------------------------
      for(size_t t = 0; t < T; ++t){
        // convert Eigen matrix to raw data pointer
        std::complex<double>* peram_send = &(perambulator[nbs][nbr][t](0,0));
        if(myid_t == 0)
          MPI_Reduce(MPI_IN_PLACE, peram_send, 
                     2*perambulator[nbs][nbr][t].size(), MPI_DOUBLE, MPI_SUM, 
                     0, ts_comm);
        else 
          MPI_Reduce(peram_send, 0,
                     2*perambulator[nbs][nbr][t].size(), MPI_DOUBLE, MPI_SUM, 
                     0, ts_comm);
        MPI_Barrier(MPI_COMM_WORLD); 
        // copy to writing array and swap endianness if neccessary -------------
        if(myid_t == 0){
          if(param.endianness == "little")
            for(int row = 0; row < perambulator[nbs][nbr][t].rows(); row++){
              int row_off = perambulator[nbs][nbr][t].size()*t + 
                            row * perambulator[nbs][nbr][t].cols();
              for(int col = 0; col < perambulator[nbs][nbr][t].cols(); col++) 
                write_peram[row_off + col] = 
                                (perambulator[nbs][nbr][t](row, col));
            }
          else
            for(int row = 0; row < perambulator[nbs][nbr][t].rows(); row++){
              int row_off = perambulator[nbs][nbr][t].size()*t + 
                            row * perambulator[nbs][nbr][t].cols();
              for(int col = 0; col < perambulator[nbs][nbr][t].cols(); col++) 
                write_peram[row_off + col] = 
                              swap_complex(perambulator[nbs][nbr][t](row, col));
            }
        }
      }
      // creating perambulator file name ---------------------------------------
      std::string filename = param.outpath + "/perambulator";
      sprintf(outfile,
        "%s.rndvecnb%02d.%s.Tso%s%04d.Vso%s%04d.Dso%s%01d.Tsi%s%04d." \
        "Ssi%s%04d.Dsi%s%01d.Csi%s%01d.smeared%01d.%05d",
        filename.c_str(), (int) param.rnd_id[rnd_id], param.quarktype.c_str(),
        param.dilution_type_so[0].c_str(), (int) param.dilution_size_so[0], 
        param.dilution_type_so[1].c_str(), (int) param.dilution_size_so[1], 
        param.dilution_type_so[2].c_str(), (int) param.dilution_size_so[2],
        param.dilution_type_si[nbs][0].c_str(), 
        (int) param.dilution_size_si[nbs][0], 
        param.dilution_type_si[nbs][1].c_str(), 
        (int) param.dilution_size_si[nbs][1], 
        param.dilution_type_si[nbs][2].c_str(), 
        (int) param.dilution_size_si[nbs][2], 
        param.dilution_type_si[nbs][3].c_str(), 
        (int) param.dilution_size_si[nbs][3], 
        (int) nbr, (int) param.config);
      if(myid==0) printf("\twriting perambulator to file: %s\n", outfile);

      // writing data in parallel, but only with processes wich have the data --
      if(myid_t == 0){
        MPI_File fh;
        MPI_Status status;
        MPI_File_open(write_comm, outfile, 
                      MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
        int my_offset = myid_write*write_peram.size()*2*sizeof(double);
        MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
        MPI_File_write(fh, &(write_peram[0]), 2*write_peram.size(), 
                       MPI_DOUBLE, &status);
        MPI_File_close(&fh);
        // check that all data have been written
        int size_of_written_data = 0;
        MPI_Get_count(&status, MPI_DOUBLE, &size_of_written_data);
        if(size_of_written_data != int(2*write_peram.size()))
          std::cout << "\n\n\t\tERROR\n\n" 
                    << "\t\tNot all perambulator data were written for " 
                    << outfile << "\n\n" << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);

    } // for loop over random vectors ends here
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

  // writing random vector to disc ---------------------------------------------
  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int kill_program = 0; // flag to check if data were correctly written
  if(myid == 0)
    kill_program = write_random_vector_to_disk(rnd_id);
  MPI_Bcast(&kill_program, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(kill_program){ // kill program in case of write failure of rnd vecs
    MPI_Finalize();
    exit(0);
  }

  delete[] rnd;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::set_sink_random_vector(const size_t rnd_id, 
          const size_t sink_id, const size_t rnd_id_si, Eigen::VectorXcd& out) {

  rlxd_init(2, 
       ((param.seed[rnd_id])^param.config)^(param.seed_si[sink_id][rnd_id_si]));
  // Reinterpret cast is NOT good style but necessary here to connect the data
  // here with the structs of tmLQCD!
  random_spinor_field_lexic(reinterpret_cast<spinor*>(&out[0]), 1, RN_Z2);

  
  char outfile[400];
  FILE *fp = NULL;

  // writing data to disk
  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int kill_program = 0; // flag to check if data were correctly written
  if(myid == 0){  
    std::string filename = "./randomvector_sink";
    sprintf(outfile, "%s.%zu_%zu", filename.c_str(), rnd_id, sink_id);
  
    if((fp = fopen(outfile, "wb")) == NULL){
      std::cout << "\n\nfailed to open file to write random vector: " 
                << outfile << "\n" << std::endl;
      kill_program = 1;
    }   
    int check_write_size = fwrite(&out[0], sizeof(std::complex<double>), 
                                  out.size(), fp);
    if(check_write_size != out.size()){
      std::cout << "\n\ncould not write all data for sink rnd vector!\n: " 
                << outfile << "\n" << std::endl;
      kill_program = 1;
    }   
    fclose(fp);
  }
  MPI_Bcast(&kill_program, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(kill_program){ // kill program in case of write failure of rnd vecs
    MPI_Finalize();
    exit(0);
  }

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
    if(check_read_in != Lt*rnd_vec_length){
      std::cout << "Failed to read random vectors - program aborted!" 
                << std::endl;
      MPI_Finalize();
      exit(0);
    }
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
int LapH::distillery::write_random_vector_to_disk(size_t rnd_id){

  char outfile[400];
  FILE *fp = NULL;
  const size_t Lt = param.Lt;
  const size_t verbose = param.verbose;
  const size_t number_of_eigen_vec = param.nb_ev;
  const size_t rnd_vec_length = Lt * number_of_eigen_vec * 4;

  // copy from eigen structure into intermediate memory
  std::complex<double>* rnd_vec_write =
      new std::complex<double>[rnd_vec_length];
  if(param.endianness == "little")
    for(size_t t = 0; t < Lt; ++t)
      for(size_t row_i = 0; row_i < 4 * number_of_eigen_vec; ++row_i)
        rnd_vec_write[row_i + t * rnd_vec_length/Lt] = 
                                            random_vector[t](row_i);
  else
    for(size_t t = 0; t < Lt; ++t)
      for(size_t row_i = 0; row_i < 4 * number_of_eigen_vec; ++row_i)
        rnd_vec_write[row_i + t * rnd_vec_length/Lt] = 
                                            swap_complex(random_vector[t](row_i));
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
  int kill_program = 0;
  if(check_read_in != (int) rnd_vec_length){
    kill_program = 1;
    std::cout << "\n\nfailed to write random vector: "
              << outfile << "\n" << std::endl;
  }   
  fclose(fp);
  // delete intermediate memory
  delete[] rnd_vec_write;
  return kill_program;

} 
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::read_eigenvectors(){

  MPI_Barrier(MPI_COMM_WORLD); 
  double time1 = MPI_Wtime(), time2 = MPI_Wtime(), time3 = MPI_Wtime();

  const int Ls = param.Ls;
  const int Lt = param.Lt;
  const int verbose = param.verbose;
  const int number_of_eigen_vec = param.nb_ev;

  const int T = Lt/tmLQCD_params->nproc_t;
  const int X = Ls/tmLQCD_params->nproc_x;
  const int Y = Ls/tmLQCD_params->nproc_y;
  const int Z = Ls/tmLQCD_params->nproc_z;

  const int nproc_x = tmLQCD_params->nproc_x;
  const int nproc_y = tmLQCD_params->nproc_y;
  const int nproc_z = tmLQCD_params->nproc_z;

  MPI_File fh;
  int file_open_error;

  int myid = 0, nb_ranks = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_ranks);

  if(myid == 0){
    if(verbose) printf("reading eigen vectors from files:\n");
    else printf("\treading eigenvectors\n");
    fflush(stdout);
  }

  // create communicator to be able to read diff. timeslices on diff. procs
  MPI_Comm ts_comm;
  MPI_Comm_split(MPI_COMM_WORLD, tmLQCD_params->proc_coords[0], 
                 myid, &ts_comm );
  int myid_t = 0, nb_procs_t = 0;
  MPI_Comm_rank(ts_comm, &myid_t);
  MPI_Comm_size(ts_comm, &nb_procs_t);

  // variables for checking trace and sum of v^daggerv
  std::complex<double> trace_s(.0,.0), trace_r(.0,.0), 
                       sum_r(.0,.0), sum_s(.0,.0);

  // running over all timeslices on this process -------------------------------
  for(int t = 0; t < T; t++){

    time2 = MPI_Wtime();// TODO: Just for testing
    // setting up filename
    const int real_t = T*tmLQCD_params->proc_coords[0] + t;
    char name[200];
    sprintf(name, "%s/eigenvectors.%04d.%03d", 
                  param.inpath_ev.c_str(), (int) param.config, real_t);
    if(myid_t == 0) std::cout << "Reading file: " << name << std::endl;
    // open file and check if it worked
    file_open_error = MPI_File_open(ts_comm, name, MPI_MODE_RDONLY, 
                                    MPI_INFO_NULL, &fh);
    if (file_open_error != MPI_SUCCESS) {
  
      char error_string[BUFSIZ];
      int length_of_error_string, error_class;
  
      MPI_Error_class(file_open_error, &error_class);
      MPI_Error_string(error_class, error_string, &length_of_error_string);
      printf("%3d: %s\n", myid, error_string);
  
      MPI_Error_string(file_open_error, error_string, &length_of_error_string);
      printf("%3d: %s\n", myid, error_string);
  
      MPI_Abort(ts_comm, file_open_error);
    }
    
    // reading and distributing data -------------------------------------------
    if(myid_t == 0) std::cout << "starting to read eigenvectors!" << std::endl;
      time3 = MPI_Wtime();// TODO: Just for testing
MPI_Barrier(ts_comm);
// creating new MPI datatype ---------------------------------------------------
MPI_Datatype newtype;

////////////////////////////////////////////////////////////////////////////////

// Old Try

int array_of_sizes[4] = {number_of_eigen_vec, Ls, Ls, 2*3*Ls};
int array_of_subsizes[4] = {1, nproc_x, nproc_y, nproc_z};

int distribs[4] = {MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK};
int dargs[4] = {MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG};

//for(int id = 0; id < nb_ranks; id++){
//  if(id == myid){
//    std::cout << "I am process: " << myid << std::endl;
//    std::cout << "t: "   << T*tmLQCD_params->proc_coords[0] + t << std::endl;
//    std::cout << distribs[0] << " " << distribs[1] << " "
//              << distribs[2] << " "<< distribs[3] << " "
//              << dargs[0] << " " << dargs[1] << " " 
//              << dargs[2] << " " << dargs[3] << " " << std::endl;
//    std::cout << "\n\n" << std::endl;
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//}

MPI_Type_create_darray(nb_procs_t, myid_t, 4, array_of_sizes, distribs, dargs, array_of_subsizes, MPI_ORDER_C, MPI_DOUBLE, &newtype);
MPI_Type_commit(&newtype);


////////////////////////////////////////////////////////////////////////////////

// New Try












char datarep[] = "native";
MPI_File_set_view(fh, 0, MPI_DOUBLE, newtype, datarep, MPI_INFO_NULL);
MPI_Barrier(ts_comm);


// reading data ---------------------------------------------------------------
std::vector<std::complex<double> > eigen_vec(number_of_eigen_vec*3*Z*Y*Z, std::complex<double>(0.0,0.0));
MPI_Status status;
MPI_File_read_all(fh, eigen_vec.data(), 2*3*Z*Y*X*number_of_eigen_vec, MPI_DOUBLE, &status);
MPI_File_close(&fh);


// check if there was a problem while reading data
if(status.MPI_ERROR != MPI_SUCCESS)
  std::cout << "\n\n\nError while reading on process: " << myid 
            << " " << status.MPI_SOURCE << " " << status.MPI_TAG
            << " " << status.MPI_ERROR << "\n\n" << std::endl;
// TODO: just a rest of how many data were written
int size_of_written_data = 0;
MPI_Get_count(&status, MPI_DOUBLE, &size_of_written_data);
std::cout << "process " << myid << " read " << size_of_written_data 
          << " bytes! It should have read " << number_of_eigen_vec*2*3*Z*Y*X 
          << std::endl;

copy_to_V(eigen_vec, t);

for(int id = 0; id < nb_ranks; id++){
  if(id == myid && real_t == 0){
    std::cout << "\n\n\t myid: " << myid << std::endl;
    std::cout << "T: "   << T*tmLQCD_params->proc_coords[0] 
              << "\tX: " << X*tmLQCD_params->proc_coords[1]
              << "\tY: " << Y*tmLQCD_params->proc_coords[2] 
              << "\tZ: " << Z*tmLQCD_params->proc_coords[3] << std::endl;
    for(int row = 0; row < V[t].rows(); row++)
      for(int col = 0; col < V[t].cols(); col++)
        std::cout << V[t](row, col) << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

    // TODO: This is merely a test
    MPI_Barrier(ts_comm);
    if(myid == 0)
      std::cout << "\t\tTime for reading eigenvectors on timeslice " 
                << real_t << " : " << MPI_Wtime() - time3 << std::endl;
    // compute the trace and sum of v^daggerv --------------
    trace_s += (V[t].adjoint() * V[t]).trace();
    sum_s += (V[t].adjoint() * V[t]).sum();
    MPI_Barrier(ts_comm);
    // TODO: This is merely a test
    if(myid == 0)
      std::cout << "\tTime for EV reading on timeslice " << real_t 
                << " : " << MPI_Wtime() - time2 << std::endl;
  }

  // checking trace and sum of v^dagger v ---------------------------------------
  MPI_Allreduce(&trace_s, &trace_r, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&sum_s, &sum_r, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if((fabs(trace_r.real() - Lt*number_of_eigen_vec) > 10e-6) ||
     (fabs(sum_r.real()   - Lt*number_of_eigen_vec) > 10e-6) ||
     (fabs(trace_r.imag()) > 10e-6) || (fabs(sum_r.imag()) > 10e-6) ){
    if(myid == 0)
      std::cout << "\n\nTrace of sum of V^daggerV is not correct! "
                << "- abort program\n Sum = " << sum_r << " Trace = " 
                << trace_r << "\n" << std::endl;
    MPI_Finalize();
    exit(212);
  }
  // taking the time of reading eigenvectors -----------------------------------
  MPI_Barrier(MPI_COMM_WORLD);
  time1 = MPI_Wtime() - time1;
  if(myid == 0)
    std::cout << "\tTime for eigenvector reading: " << time1 << std::endl;

// TODO: Test to check if data are correctly read
if(myid == 0) std::cout << "\n\n" << std::endl;
for(int id = 0; id < nb_ranks; id++){
  if(id == myid){
    std::cout << "I am process: " << myid << std::endl;
    std::cout << "T: "   << T*tmLQCD_params->proc_coords[0] 
              << "\tX: " << X*tmLQCD_params->proc_coords[1]
              << "\tY: " << Y*tmLQCD_params->proc_coords[2] 
              << "\tZ: " << Z*tmLQCD_params->proc_coords[3] << std::endl;
    for(int t = 0; t < T; t++)
      std::cout << "t: " << T*tmLQCD_params->proc_coords[0] + t
                << "\tSum: "   << (V[t].adjoint() * V[t]).sum() 
                << "\tTrace: " << (V[t].adjoint() * V[t]).trace() << std::endl;
    std::cout << "\n\n" << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

  MPI_Finalize();
  exit(0);
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::distillery::copy_to_V(
                            const std::vector<std::complex<double> >& eigen_vec, 
                            const int t){
  const size_t Ls = param.Ls;
  const size_t number_of_eigen_vec = param.nb_ev;
  const size_t X = Ls/tmLQCD_params->nproc_x;
  const size_t Y = Ls/tmLQCD_params->nproc_y;
  const size_t Z = Ls/tmLQCD_params->nproc_z;
  const size_t dim_row = X*Y*Z;

  for(size_t nev = 0; nev < number_of_eigen_vec; nev++){
    size_t j = 0;    
    for(size_t x = 0; x < X; x++){
      for(size_t y = 0; y < Y; y++){
        for(size_t z = 0; z < Z; z++){
          for(size_t c = 0; c < 3; c++){
  
            size_t i = 3*(dim_row*nev + x*Y*Z + y*Z + z) + c;
            if(param.endianness == "little")
              (V[t])(j, nev) = eigen_vec.at(i);
            else
              (V[t])(j, nev) = swap_complex(eigen_vec[i]);
            j++;
  
          }
        }
      }
    }
  }
}
//// -----------------------------------------------------------------------------
//// -----------------------------------------------------------------------------
//void LapH::distillery::read_eigenvectors(){
//
//  MPI_Barrier(MPI_COMM_WORLD);
//  double time1 = MPI_Wtime();
//  double time2;
//
//  const size_t Ls = param.Ls;
//  const size_t Lt = param.Lt;
//  const size_t verbose = param.verbose;
//  const size_t number_of_eigen_vec = param.nb_ev;
//
//  const size_t T = Lt/tmLQCD_params->nproc_t;
//  const size_t X = Ls/tmLQCD_params->nproc_x;
//  const size_t Y = Ls/tmLQCD_params->nproc_y;
//  const size_t Z = Ls/tmLQCD_params->nproc_z;
//  const size_t nproc_x = tmLQCD_params->nproc_x;
//  const size_t nproc_y = tmLQCD_params->nproc_y;
//  const size_t nproc_z = tmLQCD_params->nproc_z;
//  const size_t px = tmLQCD_params->proc_coords[1];
//  const size_t py = tmLQCD_params->proc_coords[2];
//  const size_t pz = tmLQCD_params->proc_coords[3];
//
//  const size_t dim_row = Ls * Ls * Ls * 3;
//
//  MPI_Offset my_offset;
//  MPI_File fh;
//  int file_open_error;
//  MPI_Status* status = new MPI_Status[nproc_x*nproc_y*nproc_z];
//
//  int myid = 0;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//
//  // determine how many MPI reads it will take to read the eigenvectors
//  std::vector<ev_chunk> ev_chunks;
//  const size_t max_buflen = create_eigenvector_chunks(ev_chunks);
//
//  //buffer for read in
//  std::vector<std::complex<double> > eigen_vec(max_buflen);
//
//  if(myid == 0){
//    if(verbose) printf("reading eigen vectors from files:\n");
//    else printf("\treading eigenvectors\n");
//    fflush(stdout);
//  }
//
//  // create communicator for this timeslice amongst all processes which 
//  // also contain it
//  MPI_Comm ts_comm;
//  MPI_Comm_split(MPI_COMM_WORLD, tmLQCD_params->proc_coords[0], 
//                 myid, &ts_comm );
//
//  // variables for checking trace and sum of v^daggerv
//  std::complex<double> trace_s(.0,.0), trace_r(.0,.0), sum_r(.0,.0), sum_s(.0,.0);
//  // running over all timeslices on this process
//  for(size_t t = 0; t < T; t++){
//    MPI_Barrier(MPI_COMM_WORLD);
//    // setting up filename
//    const int real_t = T*tmLQCD_params->proc_coords[0] + t;
//
//    char name[200];
//    snprintf(name, 200, "%s/eigenvectors.%04d.%03d",
//                  param.inpath_ev.c_str(), (int) param.config, real_t);
//    if(verbose) std::cout << "Reading file: " << name << std::endl;
//    // open file and check if it worked
//    MPI_Barrier(ts_comm);
//    file_open_error = MPI_File_open(ts_comm, name,
//                          MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
//    if (file_open_error != MPI_SUCCESS) {
//      char error_string[BUFSIZ];
//      int length_of_error_string, error_class;
//
//      MPI_Error_class(file_open_error, &error_class);
//      MPI_Error_string(error_class, error_string, &length_of_error_string);
//      printf("%3d: %s\n", myid, error_string);
//
//      MPI_Error_string(file_open_error, error_string, &length_of_error_string);
//      printf("%3d: %s\n", myid, error_string);
//
//      MPI_Abort(MPI_COMM_WORLD, file_open_error);
//    }
//    
//    // XLC has no support for range loops
//    for(size_t i_chunk = 0; i_chunk < ev_chunks.size(); ++i_chunk){
//      if(verbose){
//        std::cout << "read_eigenvectors: Process: " << myid << " reading eigenvectors " << ev_chunks[i_chunk].offset 
//                  << " to " << ev_chunks[i_chunk].offset+ev_chunks[i_chunk].stride-1 << " on timeslice " << real_t << std::endl;
//      }
//      MPI_Barrier(ts_comm);
//      time2 = MPI_Wtime();
//
//      // reset all statuses to not terminate in case MPI_File_read_all does not properly set the status struct
//      // (which seems to be the case...)
//      for(int i_status = 0; i_status < nproc_x*nproc_y*nproc_z; ++i_status  ) status[i_status].MPI_ERROR = MPI_SUCCESS;
//      
//      // the offset argument is generally of type size_t and does not overflow,
//      // the count argument however potentially overflows twice, once because it's an int
//      // and a second time when it gets converted into the number of bytes (also an int...)
//      MPI_File_read_at_all(fh, 2*dim_row*ev_chunks[i_chunk].offset, &eigen_vec[0], 2*dim_row*ev_chunks[i_chunk].stride, MPI_DOUBLE, status);
//      
//      bool do_abort = false;
//      for(int i_status = 0; i_status < nproc_x*nproc_y*nproc_z; ++i_status  ) {
//        if( status[i_status].MPI_ERROR != MPI_SUCCESS ){
//          std::cout << "read_eigenvectors: File read error " << status[i_status].MPI_ERROR <<  " on process " << myid << " status element " << i_status
//                    << " local timeslice: " << t << " global timeslice: " << real_t << std::endl;
//          do_abort = true;
//        }
//      }
//      // delay abort up to here to collect more errors
//      if(do_abort) MPI_Abort( MPI_COMM_WORLD, 234 );
//      
//      MPI_Barrier(ts_comm);
//      if(verbose) std::cout << "read_eigenvectors: Time for reading " << 2*dim_row*ev_chunks[i_chunk].stride*sizeof(double)/(1024*1024)
//                            << " MB: " << MPI_Wtime()-time2 << " seconds. Process: " << myid << std::endl;
//      
//      time2=MPI_Wtime();
//      for (size_t nev = ev_chunks[i_chunk].offset; nev < (ev_chunks[i_chunk].offset+ev_chunks[i_chunk].stride); ++nev) {
//        // copying the correct components into V
//        copy_to_V(&eigen_vec[(nev-ev_chunks[i_chunk].offset)*dim_row], t, nev);
//      }
//      MPI_Barrier(ts_comm);
//      if(verbose) std::cout << "read_eigenvectors: Time for copy_to_V of " << ev_chunks[i_chunk].stride << " eigenvectors: "
//                            << MPI_Wtime()-time2 << " seconds" << std::endl;
//    }
//    MPI_Barrier(ts_comm);
//    MPI_File_close(&fh); // this works fine on all processes
//    std::cout << "Closed filehandle: " << myid << std::endl;  
//
//    // computing trace and sum for check! --------------------------------------
//    MPI_Barrier(ts_comm);
//    time2 = MPI_Wtime();
//    if(verbose) 
//    std::cout << "read_eigenvectors: computing local V^d V" << std::endl;
//    trace_s += (V[t].adjoint() * V[t]).trace();
//    sum_s += (V[t].adjoint() * V[t]).sum();
//    if(verbose) 
//      std::cout << "read_eigenvectors: time for v^dagger v: " 
//                << MPI_Wtime()-time2 << " seconds on Process: " << myid 
//                << std::endl;
//    MPI_Barrier(ts_comm);    // it seems like this barrier is not reached by 70+-30 out of 512 processes...
//     /* by removing the trace and sum, execution can be made to continue
//      * this seems to indicate that something about V[t] is not right
//      * which makes the computation very slow on some processes
//      * */
//    std::cout << "Freed ts_comm communicator: " << myid << std::endl;  
//  }
//
//  // checking trace nd sum of v^dagger v ---------------------------------------
//  std::cout << "Communicating trace and sum " << myid << std::endl;
//  MPI_Barrier(MPI_COMM_WORLD);
//  MPI_Allreduce(&trace_s, &trace_r, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  MPI_Allreduce(&sum_s, &sum_r, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  if((fabs(trace_r.real() - Lt*number_of_eigen_vec) > 10e-6) ||
//     (fabs(sum_r.real()   - Lt*number_of_eigen_vec) > 10e-6) ||
//     (fabs(trace_r.imag()) > 10e-6) || (fabs(sum_r.imag()) > 10e-6) ){
//    if(myid == 0)
//      std::cout << "\n\nTrace of sum of V^daggerV is not correct! "
//                << "- abort program\n Sum = " << sum_r << " Trace = " 
//                << trace_r << "\n" << std::endl;
//    MPI_Finalize();
//    exit(212);
//  }
//  
//  // taking the time -----------------------------------------------------------
//  MPI_Barrier(MPI_COMM_WORLD);
//  time1 = MPI_Wtime() - time1;
//  if(myid == 0) std::cout << "\tTotal time for eigenvector reading: " 
//                          << time1 << std::endl;
//  free(status);
//
//  // TODO: just a test
//  MPI_Finalize();
//  exit(0);
//}
//// -----------------------------------------------------------------------------
//// -----------------------------------------------------------------------------
//size_t LapH::distillery::create_eigenvector_chunks(std::vector<ev_chunk>& ev_chunks){
//  const size_t Ls = param.Ls;
//  const size_t dim_row = Ls*Ls*Ls*3;
//  const size_t number_of_eigen_vec = param.nb_ev;
//
//  // number of times that MPI_File_read_all has to be called (at least) to avoid int overflow due to ridiculous MPI interface
//  // the int overflow happens internally when the "count" argument is converted from the MPI datatype to bytes
//  // one therefore has to take into account the number of bytes rather tha the number of doubles...
//  size_t nchunks       = (int)ceil( ((double) 2*number_of_eigen_vec*dim_row*sizeof(double) )/INT_MAX );
//  size_t nev_in_chunk  = number_of_eigen_vec / nchunks;
//
//  size_t nev_remaining = number_of_eigen_vec;
//  size_t nev_processed = 0;
//  // if the integer division above had a remainder, the last chunk will process
//  // fewer eigenvectors, we thus construct a vector of tuples to represent this
//  for(size_t i = 0; i < nchunks; ++i){
//    size_t nev_this_chunk = nev_remaining > nev_in_chunk ? nev_in_chunk : nev_remaining;
//    ev_chunk this_chunk = { nev_processed, nev_this_chunk };
//    ev_chunks.push_back( this_chunk );
//    nev_processed += nev_this_chunk;
//    nev_remaining -= nev_this_chunk;
//  }
//  if( nev_remaining != 0 ){
//    std::cout << "Problem with eigenvector reading logic. Aborting!" << std::cout;
//    MPI_Abort(MPI_COMM_WORLD, 111);
//  }
//  // we return the maximum buffer size required for reading these chunks (in units of complex double)
//  return(nev_in_chunk*dim_row);
//}
//// -----------------------------------------------------------------------------
//// -----------------------------------------------------------------------------


//void LapH::distillery::read_eigenvectors(){
//
//  MPI_Barrier(MPI_COMM_WORLD); 
//  double time1 = MPI_Wtime();
//
//  const size_t Ls = param.Ls;
//  const size_t Lt = param.Lt;
//  const size_t verbose = param.verbose;
//  const size_t number_of_eigen_vec = param.nb_ev;
//
//  const size_t T = Lt/tmLQCD_params->nproc_t;
//  const size_t dim_row = Ls * Ls * Ls * 3;
//
//  int myid = 0;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//
//  //buffer for read in
//  std::complex<double>* eigen_vec = new std::complex<double>[dim_row];
//
//  if(myid == 0){
//    if(verbose) printf("reading eigen vectors from files:\n");
//    else printf("\treading eigenvectors\n");
//    fflush(stdout);
//  }
//  // variables for checking trace and sum of v^daggerv
//  std::complex<double> trace_s(.0,.0), trace_r(.0,.0), sum_r(.0,.0), sum_s(.0,.0);
//  // running over all timeslices on this process
//  for(size_t t = 0; t < T; t++){
//
//    const int real_t = T*tmLQCD_params->proc_coords[0] + t;
// 
//    //setting up file
//    char name[200];
//    sprintf(name, "%s/eigenvectors.%04d.%03d", 
//                  param.inpath_ev.c_str(), (int) param.config, real_t);
//    if(verbose) std::cout << "Reading file: " << name << std::endl;
//    std::ifstream infile(name, std::ifstream::binary);
//    if (infile) {
//      for (size_t nev = 0; nev < number_of_eigen_vec; ++nev) {
//        // reading the full vector
//        infile.read((char*) eigen_vec, 2*dim_row*sizeof(double));
//
//        //if(infile.gcount() != int(2*dim_row*sizeof(double))){
//        if(!infile){
//          std::cout << "\n\nreading eigenvectors failed at t = " << real_t 
//                    << " and ev = " << nev << std::endl;
//          MPI_Finalize();
//          exit(0);
//        }
//        // copying the correct components into V
//        copy_to_V(eigen_vec, t, nev);
//      }
//    }
//    else {
//      std::cout << "eigenvector file does not exist!!!\n" << std::endl;
//      exit(0);
//    }
//    infile.close();
//    // computing trace and sum for check!
//    trace_s += (V[t].adjoint() * V[t]).trace();
//    sum_s += (V[t].adjoint() * V[t]).sum();
//  }
//  delete[] eigen_vec;
//
//  // checking the trace and the sum of v^daggerv -------------------------------
//  MPI_Allreduce(&trace_s, &trace_r, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  MPI_Allreduce(&sum_s, &sum_r, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  if((fabs(trace_r.real() - Lt*number_of_eigen_vec) > 10e-6) ||
//     (fabs(sum_r.real()   - Lt*number_of_eigen_vec) > 10e-6) ||
//     (fabs(trace_r.imag()) > 10e-6) || (fabs(sum_r.imag()) > 10e-6) ){
//    if(myid == 0)
//      std::cout << "\n\nTrace of sum of V^daggerV is not correct! "
//                << "- abort program\n Sum = " << sum_r << " Trace = " 
//                << trace_r << "\n" << std::endl;
//    MPI_Finalize();
//    exit(0);
//  }
//
//  MPI_Barrier(MPI_COMM_WORLD);
//  time1 = MPI_Wtime() - time1;
//  if(myid == 0)
//    std::cout << "\tTime for eigenvector reading: " << time1 << std::endl;
//}
//// -----------------------------------------------------------------------------
//// -----------------------------------------------------------------------------












