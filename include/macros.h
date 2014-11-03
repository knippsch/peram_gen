// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_ASSIGN_MOVE(TypeName) \
  TypeName(const TypeName&);\
  TypeName(TypeName&&);\
  void operator=(const TypeName&)           

