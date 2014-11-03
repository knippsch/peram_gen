#include "input_parms.h"
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// swap function - must be in namespace because it is a friend function
namespace LapH {
  void swap(input_parameter& first, input_parameter& second){

  first.dilution_size_so.swap(second.dilution_size_so); 
  first.dilution_type_so.swap(second.dilution_type_so);
  first.dilution_sink.swap(second.dilution_sink);

  std::swap(first.config,         second.config);
  std::swap(first.nb_config,      second.nb_config);
  std::swap(first.delta_config,   second.delta_config);
  std::swap(first.Ls,             second.Ls);
  std::swap(first.Lt,             second.Lt);
  std::swap(first.nb_ev,          second.nb_ev);
  std::swap(first.nb_sinks,       second.nb_sinks);
  std::swap(first.nb_stoch_sinks, second.nb_stoch_sinks);
  std::swap(first.verbose,        second.verbose);

  std::swap(first.nb_rnd, second.nb_rnd);
  first.rnd_id.swap(second.rnd_id);
  first.seed.swap(second.seed);

  std::swap(first.quarktype,         second.quarktype);
  std::swap(first.outpath,           second.outpath);
  std::swap(first.inpath_ev,         second.inpath_ev);
  std::swap(first.peram_file_name,   second.peram_file_name);
  std::swap(first.rnd_vec_file_name, second.rnd_vec_file_name);

  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// deep copy constructor
LapH::input_parameter::input_parameter(const input_parameter& other){

  dilution_size_so = other.dilution_size_so;
  dilution_type_so = other.dilution_type_so;

  dilution_sink = other.dilution_sink;

  config = other.config;
  nb_config = other.nb_config;
  delta_config = other.delta_config;
  Ls = other.Ls;
  Lt = other.Lt;
  nb_ev = other.nb_ev;
  nb_sinks = other.nb_sinks;
  nb_stoch_sinks = other.nb_stoch_sinks;
  verbose = other.verbose;

  nb_rnd = other.nb_rnd;
  rnd_id = other.rnd_id;
  seed = other.seed;
  
  quarktype = other.quarktype;
  outpath = other.outpath;
  inpath_ev = other.inpath_ev;
  peram_file_name = other.peram_file_name;
  rnd_vec_file_name = other.rnd_vec_file_name;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
static void exit_check(){
  std::cerr << "\n\tcheck size of dilution space\n" << std::endl;
  exit(0);
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// checking that the dilution parameters are correctly set in infile
// input:  type       -> gives the type of the dilution scheme
//         max_size   -> maximal size of space
// output: nb_dil_vec -> sets (in case of full and no dilution) and checks the 
//                       number of dilution vectors
static void check_dilution_input(const std::string type, const size_t max_size, 
                                 const size_t& nb_dil_vec){
  // check for type
  if( (type.compare("F") != 0) && 
      (type.compare("I") != 0) && 
      (type.compare("B") != 0) && 
      (type.compare("N") != 0) ) {
        std::cerr << "Inversion type has to be one of \"F\", \"I\"," \
                     " \"B\" or \"N\"." << std::endl;
        std::cerr << "Aborting..." << std::endl;
        exit(0);
  }
  // check and set number of inversions in corresponding space
  if((max_size % nb_dil_vec) != 0)
    exit_check();

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::input_parameter::check_input_parameters() const {

  if(config > 100000) {
    std::cout << "Please check whether configuration number is correct: " 
              << config << std::endl;
    exit(0);
  }
  if(Ls > 500) {
    std::cout << "Please check whether Ls is correct: " <<  Ls << std::endl;
    exit(0);
  }
  if(Lt > 500) {
    std::cout << "Please check whether Lt is correct: " <<  Lt << std::endl;
    exit(0);
  }
  if(nb_ev > 5000) {
    std::cout << "Please check whether number of eigenvectors is correct: " 
              <<  nb_ev << std::endl;
    exit(0);
  }
  if(nb_rnd > 1000) {
    std::cout << "Please check whether starting number of randomvector" \
                 " is correct: " <<  nb_rnd << std::endl;
    exit(0);
  }
  if((quarktype.compare("u") != 0) && 
     (quarktype.compare("d") != 0) && 
     (quarktype.compare("s") != 0) && 
     (quarktype.compare("c") != 0) ) {
       std::cerr << "Quarktype has to be one of \"u\", \"d\", \"s\" or \"c\"." 
                 << std::endl;
       std::cerr << "Aborting..." << std::endl;
       exit(0);
  }
  check_dilution_input(dilution_type_so[0], Lt, dilution_size_so[0]); 
  check_dilution_input(dilution_type_so[1], nb_ev, dilution_size_so[1]); 
  check_dilution_input(dilution_type_so[2], 4, dilution_size_so[2]);
  for(size_t sinks = 0; sinks < nb_sinks; sinks++){
    check_dilution_input(dilution_sink[sinks].dilution_type_si[0], Lt, 
                         dilution_sink[sinks].dilution_size_si[0]); 
    check_dilution_input(dilution_sink[sinks].dilution_type_si[1], 
                         Ls*Ls*Ls, 
                         dilution_sink[sinks].dilution_size_si[1]); 
    check_dilution_input(dilution_sink[sinks].dilution_type_si[2], 4, 
                         dilution_sink[sinks].dilution_size_si[2]); 
    check_dilution_input(dilution_sink[sinks].dilution_type_si[3], 3, 
                         dilution_sink[sinks].dilution_size_si[3]); 
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// TODO: Discription
std::string LapH::input_parameter::create_out_filename(const std::string& type, 
                                               const int sink_id) const {
  char name[500];
  std::string filename;
  if(type.compare("randomvector") == 0){       // for randomvector -------------
    sprintf(name, "%s/randomvector.rndvecnb%02d.%s.nbev%04d.%04d", 
            outpath.c_str(), (int) nb_rnd, quarktype.c_str(), (int) nb_ev, 
            (int) config);
    filename = name;
    printf("Writing randomvector to file: %s\n", filename.c_str());
  }
  else if (type.compare("perambulator") == 0){ // for perambulator -------------
    sprintf(name,
        "%s/perambulator.rndvecnb%02d.%s.Tso%s%04d.Vso%s%04d.Dso%s%01d." \
        "Tsi%s%04d.Ssi%s%04d.Dsi%s%01d.Csi%s%01d.smeared%01d.%05d",
        outpath.c_str(), (int) nb_rnd, quarktype.c_str(), 
        dilution_type_so[0].c_str(), (int) dilution_size_so[0], 
        dilution_type_so[1].c_str(), (int) dilution_size_so[1], 
        dilution_type_so[2].c_str(), (int) dilution_size_so[2],
        dilution_sink[sink_id].dilution_type_si[0].c_str(), 
        (int) dilution_sink[sink_id].dilution_size_si[0], 
        dilution_sink[sink_id].dilution_type_si[1].c_str(), 
        (int) dilution_sink[sink_id].dilution_size_si[1], 
        dilution_sink[sink_id].dilution_type_si[2].c_str(), 
        (int) dilution_sink[sink_id].dilution_size_si[2], 
        dilution_sink[sink_id].dilution_type_si[3].c_str(), 
        (int) dilution_sink[sink_id].dilution_size_si[3], 
        1, (int) config);
    filename = name;
    printf("Writing perambulator to file: %s\n", filename.c_str());
  }
  else{
    std::cout << "\n\n\tCould not find type of which a filename should be" \
                 "genereated!\n\n" << std::endl;
    exit(0);
  }
  // TODO: check if files already exist and rename just in case!!! 
  return filename;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::input_parameter::parse_input_file(int argc, char *argv[]) {

  int opt = -1;
  int reader = 0;
  char infilename[200];
  char readin[256];
  FILE* infile = NULL;

  // search for command line option and put filename in "infilename"
  for(int i = 0; i < argc; ++i) {
    if(std::strcmp(argv[i], "-LapHsin") == 0) {
      opt = i+1;
      break;
    }
  }
  if(opt < 0) {
    std::cout << "No input file specified, trying infile.in" << std::endl;
    sprintf(infilename, "infile.in");
  } else {
    sprintf(infilename, "%s", argv[opt]);
    std::cout << "Trying input file " << infilename << std::endl;
  }
  // open file for reading
  if ((infile = fopen(infilename, "r")) == NULL ) {
    std::cerr << "Could not open file " << infilename << std::endl;
    std::cerr << "Aborting..." << std::endl;
    exit(-10);
  }

  // scan infile and check arguments ----------------------------------------
  // configs
  reader += fscanf(infile, "config = %zu \n", &config);
  reader += fscanf(infile, "total number of configs = %zu \n", &nb_config);
  reader += fscanf(infile, "distance between configs = %zu \n", &delta_config);
  // spatial extend
  reader += fscanf(infile, "Ls = %zu \n", &Ls);
  // temporal extend
  reader += fscanf(infile, "Lt = %zu \n", &Lt);
  // number of eigenvectors
  reader += fscanf(infile, "nb_ev = %zu \n", &nb_ev);
  // starting number of randomvectors
  // TODO: must be changed to allow for several random vectors
  reader += fscanf(infile, "nb_rnd = %zu \n", &nb_rnd);
  // seed and id for randomvectors
  rnd_id.resize(nb_rnd);
  seed.resize(nb_rnd);
  for (size_t i = 0; i < nb_rnd; i++)
    reader += fscanf(infile, "id %i seed %i\n", &(rnd_id[i]), &(seed[i]));
  // verbosity
  reader += fscanf(infile, "verbose = %zu\n", &( verbose));
  // quarktype
  reader += fscanf(infile, "quarktype = %255s\n", readin);
  quarktype.assign(readin);

  // SOURCE --------------------------------------------------------------------
  // type and number of inversions for the source in time ----------------------
  reader += fscanf(infile, "inversion_source_type_t = %255s\n", readin);
  dilution_type_so[0].assign(readin);
  reader += fscanf(infile, "inversion_source_number_t = %zu\n", 
                   &(dilution_size_so[0]));
  // type and number of inversions for the source in eigenvector space ---------
  reader += fscanf(infile, "inversion_source_type_v = %255s\n", readin);
  dilution_type_so[1].assign(readin);
  reader += fscanf(infile, "inversion_source_number_v = %zu\n", 
                   &(dilution_size_so[1]));
  // type and number of inversions for the soure in Dirac space ----------------
  reader += fscanf(infile, "inversion_source_type_d = %255s\n", readin);
  dilution_type_so[2].assign(readin);
  reader += fscanf(infile, "inversion_source_number_d = %zu\n", 
                   &(dilution_size_so[2]));

  // SINK ----------------------------------------------------------------------
  // number of sinks in total and number of stochastik sinks
  reader += fscanf(infile, "nb_sinks = %zu\n", &nb_sinks); 
  reader += fscanf(infile, "nb_stoch_sinks = %zu\n", &nb_stoch_sinks);
  dilution_sink.resize(nb_sinks); 
  for(size_t sinks = 0; sinks < nb_sinks; sinks++){ // reading all sinks

    reader += fscanf(infile, "%d %d %d\n", &(dilution_sink[sinks].id),
                                           &(dilution_sink[sinks].nb_vectors),
                                           &(dilution_sink[sinks].seed));
    
    // type and number of dilution vectors for the sink in time ----------------
    reader += fscanf(infile, "inversion_sink_type_t = %255s\n", readin);
    dilution_sink[sinks].dilution_type_si[0].assign(readin);
    reader += fscanf(infile, "inversion_sink_number_t = %zu\n", 
                     &(dilution_sink[sinks].dilution_size_si[0]));
    // type and number of dilution vectors for the sink in space ---------------
    reader += fscanf(infile, "inversion_sink_type_s = %255s\n", readin);
    dilution_sink[sinks].dilution_type_si[1].assign(readin);
    reader += fscanf(infile, "inversion_sink_number_s = %zu\n", 
                     &(dilution_sink[sinks].dilution_size_si[1]));
    // type and number of dilution vectors for the sink in Dirac space ---------
    reader += fscanf(infile, "inversion_sink_type_d = %255s\n", readin);
    dilution_sink[sinks].dilution_type_si[2].assign(readin);
    reader += fscanf(infile, "inversion_sink_number_d = %zu\n", 
                     &(dilution_sink[sinks].dilution_size_si[2]));
    // type and number of dilution vectors for the sink in colour space --------
    reader += fscanf(infile, "inversion_sink_type_c = %255s\n", readin);
    dilution_sink[sinks].dilution_type_si[3].assign(readin);
    reader += fscanf(infile, "inversion_sink_number_c = %zu\n", 
                     &(dilution_sink[sinks].dilution_size_si[3]));

  }
  // output path for the perambulators and the randomvectors
  reader += fscanf(infile, "outpath = %255s\n", readin);
  outpath.assign(readin);
  // input path for the eigensystems
  reader += fscanf(infile, "inpath_ev = %255s\n", readin);
  inpath_ev.assign(readin);

  // close input file
  fclose(infile);
  
  // checking input parameters
  check_input_parameters();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::input_parameter::print_options() const {

  std::cout << "config = " <<  config << std::endl;
  std::cout << "Ls = " <<  Ls << ", Lt = " <<  Lt << std::endl;
  std::cout << "nb_ev = " <<  nb_ev << ", nb_rnd = " <<  nb_rnd << std::endl;
  for (size_t i = 0; i < nb_rnd; i++)
    printf("\tid = %d, seed = %d\n", rnd_id[i], seed[i]); 
  std::cout << "verbose = " <<  verbose << std::endl;
  std::cout << "quarktype = " <<  quarktype << std::endl;
  std::cout << "inversion source time: " <<  dilution_type_so[0] << " " 
            <<  dilution_size_so[0] << std::endl;
  std::cout << "inversion source eigenspace: " <<  dilution_type_so[1] << " " 
            <<  dilution_size_so[1] << std::endl;
  std::cout << "inversion source Dirac: " <<  dilution_type_so[2] << " " 
            <<  dilution_size_so[2] << std::endl;

  for(size_t sinks = 0; sinks < nb_sinks; sinks++){ // running over all sinks
    printf("\nid = %d, nb_vectors = %d, seed = %d\n", 
           dilution_sink[sinks].id, dilution_sink[sinks].nb_vectors,
           dilution_sink[sinks].seed);
    std::cout << "inversion sink time: "
              <<  dilution_sink[sinks].dilution_type_si[0] << " " 
              <<  dilution_sink[sinks].dilution_size_si[0] << std::endl;
    std::cout << "inversion sink spatial space: " 
              <<  dilution_sink[sinks].dilution_type_si[1] << " " 
              <<  dilution_sink[sinks].dilution_size_si[1] << std::endl;
    std::cout << "inversion sink Dirac: "
              <<  dilution_sink[sinks].dilution_type_si[2] << " " 
              <<  dilution_sink[sinks].dilution_size_si[2] << std::endl;
    std::cout << "inversion sink color: "
              <<  dilution_sink[sinks].dilution_type_si[3] << " " 
              <<  dilution_sink[sinks].dilution_size_si[3] << std::endl;
  }
  std::cout << "output path: " <<  outpath << std::endl;
  std::cout << "input path ev: " <<  inpath_ev << std::endl;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------



