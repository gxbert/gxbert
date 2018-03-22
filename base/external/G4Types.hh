#ifndef G4TYPES_HH
#define G4TYPES_HH

#ifdef WIN32
#else
  #define G4DLLEXPORT
  #define G4DLLIMPORT
  #define G4GLOB_DLL
#endif

#include <complex>

// Definitions for Thread Local Storage
//
//#include "tls.hh"
#  define G4ThreadLocalStatic static
#  define G4ThreadLocal 

// Typedefs to decouple from library classes
// Typedefs for numeric types
//
typedef double G4double;
typedef float G4float;
typedef int G4int;
typedef bool G4bool;
typedef long G4long;
typedef std::complex<G4double> G4complex;

// Forward declation of void type argument for usage in direct object
// persistency to define fake default constructors
//
class __void__;

#endif /* G4TYPES_HH */
