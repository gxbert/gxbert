// Global Constants and typedefs
//
// History:
// 30.06.95 P.Kent - Created
// 16.02.96 G.Cosmo - Added inclusion of "templates.hh"
// 03.03.96 M.Maire - Added inclusion of "G4PhysicalConstants.hh"
// 08.11.96 G.Cosmo - Added cbrt() definition and G4ApplicationState enum type
// 29.11.96 G.Cosmo - Added typedef of HepBoolean to G4bool
// 22.10.97 M.Maire - Moved PhysicalConstants at the end of the file
// 04.12.97 G.Cosmo,E.Tcherniaev - Migrated to CLHEP
// 26.08.98 J.Allison,E.Tcherniaev - Introduced min/max/sqr/abs functions
// 22.09.98 G.Cosmo - Removed min/max/sqr/abs functions and replaced with
//                    inclusion of CLHEP/config/TemplateFunctions.h for CLHEP-1.3
// 15.12.99 G.Garcia - Included min, max definitions for NT with ISO standard
// 15.06.01 G.Cosmo - Removed cbrt() definition

#ifndef GLOBALS_HH
#define GLOBALS_HH

#include "G4ios.hh"

#ifndef FALSE
  #define FALSE 0
#endif
#ifndef TRUE
  #define TRUE 1
#endif

//#include <algorithm>  // Retrieve definitions of min/max

// Include base types
#include "G4Types.hh"

// Get definition of G4String
//#include "G4String.hh"
typedef std::string G4String;

// Includes some additional definitions: sqr, G4SwapPtr, G4SwapObj.
//#include "templates.hh"
#include <limits>
#ifndef DBL_MIN   /* Min decimal value of a double */
#define DBL_MIN   std::numeric_limits<double>::min()  // 2.2250738585072014e-308
#endif
inline int G4lrint(double ad) {
  return (ad>0) ? static_cast<int>(ad+.5) : static_cast<int>(ad-.5);
}


// Global error function
#include "G4ExceptionSeverity.hh"

typedef std::ostringstream GXExceptionDescription;
#define G4ExceptionDescription GXExceptionDescription

void GXException(const char* originOfException,
                 const char* exceptionCode,
                             G4ExceptionSeverity severity,
                 const char* comments);

void GXException(const char* originOfException,
                 const char* exceptionCode,
                             G4ExceptionSeverity severity,
                 G4String& comments);
void GXException(const char* originOfException,
                 const char* exceptionCode,
                             G4ExceptionSeverity severity,
                 G4ExceptionDescription & description);
void GXException(const char* originOfException,
                 const char* exceptionCode,
                             G4ExceptionSeverity severity,
                 G4ExceptionDescription & description,
                 const char* comments);

#define G4Exception GXException

#endif /* GLOBALS_HH */

