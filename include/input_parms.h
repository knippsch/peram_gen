#ifndef INPUT_PARMS_H_
#define INPUT_PARMS_H_

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace LapH {

struct sink_parameter{

  // nb_vectors means here either the number of eigenvectors in case of a 
  // smeared sink or the number of dilution vectors in spatial components
  int id, nb_vectors, seed; 
  // SINK dilution sizes and types
  // 0 = time, 1 = x,y,z, 2 = Dirac, 3 = color
  std::array<size_t, 4> dilution_size_si;      
  // NOTE: Please use capital Latters I, B, F, S!! All of these numbers are 
  // NOT the blocksize or interlace sizes, but the number of inversions in 
  // the corresponding space 
  std::array<std::string, 4> dilution_type_si; 

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class input_parameter{
public:
  // ---------------------- METHODS ------------------------------------------
  // standard constructor - initialisation is done by calling parse_input_file
  input_parameter() {};  
  // deep copy constructor
  input_parameter(const input_parameter& other);
  // swap function to avoid code duplication - function is friend to swap
  // private data
  friend void swap(input_parameter& first, input_parameter& second);
  // assignment operator - note that hard copy is much better here than 
  // reference because the copy constructor can do all the work very fast!
  input_parameter& operator=(input_parameter other){
    swap(*this, other);
    return *this;
  }
  // for c++11: move constructor - uses just swap function
  input_parameter(input_parameter&& other) :
          input_parameter() { // initialize via default constructor, C++11 only
    swap(*this, other);
  }
  // standard destructor
  ~input_parameter() {};  

public:
  // --------------------- INTERFACE -----------------------------------------
  void parse_input_file(int argc, char *argv[]);
  void print_options() const;
  // creates the names for outfiles, type = perambulator, randomvector
  // TODO: Name generation should be put in seperate class and it must ne 
  // checked at the beginning if files do exist!!
  std::string create_out_filename(const std::string& type, const int sink_id) const; 
  // -------------------- END INTERFACE --------------------------------------

private:
  // checks the input parameters for consistency
  void check_input_parameters() const;

public: // TODO: should be changed to private and should be initialised with 0
  // -------------------------- DATA -----------------------------------------
  // SOURCE dilution sizes and types, 0 = time, 1 = eigevector space, 2 = Dirac
  std::array<size_t, 3> dilution_size_so; 
  // NOTE: Please use capital Latters I, B, F!! All of these numbers are NOT
  // the blocksize or interlace sizes, but the number of inversions in the 
  // corresponding space 
  std::array<std::string, 3> dilution_type_so; 

  // contains all the sink parameter
  std::vector<sink_parameter> dilution_sink; 

  size_t config;                   // configuration number
  size_t nb_config;                // total number of configs to be processed
  size_t delta_config;             // distance between configurations
  size_t Ls;                       // spatial extend
  size_t Lt;                       // temporal extend
  size_t nb_ev;                    // number of eigenvectors
  size_t nb_sinks;                 // total number of different sinks
  size_t nb_stoch_sinks;           // number of stochastic sinks
  size_t verbose;                  // displaying more informations

  size_t nb_rnd;                   // random vector id of vector in process
  std::vector<int> rnd_id;         // id of random vectors
  std::vector<int> seed;           // seeds for random vector generation

  std::string quarktype;           // quark type: u,d,s,c -> for naming outfiles
  std::string outpath;             // path to write everything
  std::string inpath_ev;           // path to read eigenvectors
  std::string peram_file_name;     // perambulator file name
  std::string rnd_vec_file_name;   // random vektor file name
  // ------------------------- END DATA ---------------------------------------

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // namespace ends here

#endif /* INPUT_PARMS_H_ */

