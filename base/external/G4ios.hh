#ifndef included_G4ios
#define included_G4ios

#include "G4Types.hh"

#include <iostream>

  extern std::ostream& GXcout;
  extern std::ostream& GXcerr;
#define G4cin std::cin
#define G4endl std::endl

#define G4cout GXcout 
#define G4cerr GXcerr 

#endif
