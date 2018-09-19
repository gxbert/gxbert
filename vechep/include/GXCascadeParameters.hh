//
// @File: GXCascadeParameters.h
//
// 20180919  Guilherme Lima -- Created, based on M.Kelsey's G4CascadeParameters
//

#ifndef GXCascadeParameters_hh
#define GXCascadeParameters_hh 1

#include <iostream>
using std::endl;

namespace gxbert {

inline namespace
GXBERT_IMPL_NAMESPACE {

//TK for GXBERT
//class GXCascadeParamMessenger;

class GXCascadeParameters {

private:	// Singleton -- no public constructor

  GXCascadeParameters();
  ~GXCascadeParameters()
  { }

public:

  // Singleton
  static const GXCascadeParameters* Instance()
  {
    if (!fpInstance) {
      fpInstance = new GXCascadeParameters;
      //TK for GXBERT: Exclude to use G4AutoDelete
      //G4AutoDelete::Register(fpInstance);
    }
    return fpInstance;
  }

  void DumpConfig(std::ostream& os) const;

  // Top-level configuration flags
  static int verbose()              { return Instance()->VERBOSE_LEVEL; }
  static bool checkConservation()   { return Instance()->CHECK_ECONS; }
  static bool usePreCompound()      { return Instance()->USE_PRECOMPOUND; }
  static bool doCoalescence()       { return Instance()->DO_COALESCENCE; }
  static bool showHistory()         { return Instance()->SHOW_HISTORY; }
  static bool use3BodyMom()	      { return Instance()->USE_3BODYMOM; }
  static bool usePhaseSpace()       { return Instance()->USE_PHASESPACE; }
  static double piNAbsorption()     { return Instance()->PIN_ABSORPTION; }
  static const std::string& randomFile() { return Instance()->RANDOM_FILE; }

  // Nuclear structure parameters
  static bool useTwoParam()      { return Instance()->TWOPARAM_RADIUS; }
  static double radiusScale()    { return Instance()->RADIUS_SCALE; }	
  static double radiusSmall()    { return Instance()->RADIUS_SMALL; }
  static double radiusAlpha()    { return Instance()->RADIUS_ALPHA; }
  static double radiusTrailing() { return Instance()->RADIUS_TRAILING; }
  static double fermiScale()     { return Instance()->FERMI_SCALE; }
  static double xsecScale()      { return Instance()->XSEC_SCALE; }
  static double gammaQDScale()   { return Instance()->GAMMAQD_SCALE; }

  // Final-state clustering cuts
  static double dpMaxDoublet() { return Instance()->DPMAX_DOUBLET; }
  static double dpMaxTriplet() { return Instance()->DPMAX_TRIPLET; }
  static double dpMaxAlpha()   { return Instance()->DPMAX_ALPHA; }

  static void DumpConfiguration(std::ostream& os) { Instance()->DumpConfig(os); }

private:
  void Initialize();		// Fill parameter values from envvar strings

  // Environment variable values, null pointers mean not set
  const char* G4CASCADE_VERBOSE;
  const char* G4CASCADE_CHECK_ECONS;
  const char* G4CASCADE_USE_PRECOMPOUND;
  const char* G4CASCADE_DO_COALESCENCE;
  const char* G4CASCADE_SHOW_HISTORY;
  const char* G4CASCADE_USE_3BODYMOM;
  const char* G4CASCADE_USE_PHASESPACE;
  const char* G4CASCADE_PIN_ABSORPTION;
  const char* G4CASCADE_RANDOM_FILE;
  const char* G4NUCMODEL_USE_BEST;
  const char* G4NUCMODEL_RAD_2PAR;
  const char* G4NUCMODEL_RAD_SCALE;
  const char* G4NUCMODEL_RAD_SMALL;
  const char* G4NUCMODEL_RAD_ALPHA;
  const char* G4NUCMODEL_RAD_TRAILING;
  const char* G4NUCMODEL_FERMI_SCALE;
  const char* G4NUCMODEL_XSEC_SCALE;
  const char* G4NUCMODEL_GAMMAQD;
  const char* DPMAX_2CLUSTER;
  const char* DPMAX_3CLUSTER;
  const char* DPMAX_4CLUSTER;

  int VERBOSE_LEVEL;		// Top-level configuration flags
  bool CHECK_ECONS;
  bool USE_PRECOMPOUND;
  bool DO_COALESCENCE;
  bool SHOW_HISTORY;
  bool USE_3BODYMOM;
  bool USE_PHASESPACE;
  double PIN_ABSORPTION;
  std::string RANDOM_FILE;

  bool BEST_PAR;		// Nuclear structure parameters
//BEST_PAR has been used in a project on hold.
//Currently setting BEST_PAR or GXNUCMODEL_USE_BEST does not improve physics performance.
//Developer can get more information about this from cascade/test/README

  bool TWOPARAM_RADIUS;
  double RADIUS_SCALE;	
  double RADIUS_SMALL;
  double RADIUS_ALPHA;
  double RADIUS_TRAILING;
  double FERMI_SCALE;
  double XSEC_SCALE;
  double GAMMAQD_SCALE;

  double DPMAX_DOUBLET;	// Final-state clustering cuts
  double DPMAX_TRIPLET;
  double DPMAX_ALPHA;

//TK for GXBERT
//  GXCascadeParamMessenger* messenger;		// For access via UI commands
//  friend class GXCascadeParamMessenger;

  static GXCascadeParameters* fpInstance;
};

GXCascadeParameters* GXCascadeParameters::fpInstance = 0;

// Constructor initializes everything once

#define OLD_RADIUS_UNITS (3.3836/1.2)		// Used with NucModel params

GXCascadeParameters::GXCascadeParameters()
  : G4CASCADE_VERBOSE(getenv("G4CASCADE_VERBOSE")),
    G4CASCADE_CHECK_ECONS(getenv("G4CASCADE_CHECK_ECONS")),
    G4CASCADE_USE_PRECOMPOUND(getenv("G4CASCADE_USE_PRECOMPOUND")),
    G4CASCADE_DO_COALESCENCE(getenv("G4CASCADE_DO_COALESCENCE")),
    G4CASCADE_SHOW_HISTORY(getenv("G4CASCADE_SHOW_HISTORY")),
    G4CASCADE_USE_3BODYMOM(getenv("G4CASCADE_USE_3BODYMOM")),
    G4CASCADE_USE_PHASESPACE(getenv("G4CASCADE_USE_PHASESPACE")),
    G4CASCADE_PIN_ABSORPTION(getenv("G4CASCADE_PIN_ABSORPTION")),
    G4CASCADE_RANDOM_FILE(getenv("G4CASCADE_RANDOM_FILE")),
    G4NUCMODEL_USE_BEST(getenv("G4NUCMODEL_USE_BEST")),
    G4NUCMODEL_RAD_2PAR(getenv("G4NUCMODEL_RAD_2PAR")),
    G4NUCMODEL_RAD_SCALE(getenv("G4NUCMODEL_RAD_SCALE")),
    G4NUCMODEL_RAD_SMALL(getenv("G4NUCMODEL_RAD_SMALL")),
    G4NUCMODEL_RAD_ALPHA(getenv("G4NUCMODEL_RAD_ALPHA")),
    G4NUCMODEL_RAD_TRAILING(getenv("G4NUCMODEL_RAD_TRAILING")),
    G4NUCMODEL_FERMI_SCALE(getenv("G4NUCMODEL_FERMI_SCALE")),
    G4NUCMODEL_XSEC_SCALE(getenv("G4NUCMODEL_XSEC_SCALE")),
    G4NUCMODEL_GAMMAQD(getenv("G4NUCMODEL_GAMMAQD")),
    DPMAX_2CLUSTER(getenv("DPMAX_2CLUSTER")),
    DPMAX_3CLUSTER(getenv("DPMAX_3CLUSTER")),
    //TK for GXBERT
    //DPMAX_4CLUSTER(getenv("DPMAX_4CLUSTER")),
    DPMAX_4CLUSTER(getenv("DPMAX_4CLUSTER"))
    //messenger(0) 
    {
  //TK for GXBERT
  //messenger = new G4CascadeParamMessenger(this);
  Initialize();
}

void GXCascadeParameters::Initialize() {
  VERBOSE_LEVEL = (G4CASCADE_VERBOSE ? atoi(G4CASCADE_VERBOSE) : 0);
  CHECK_ECONS = (0!=G4CASCADE_CHECK_ECONS);
  USE_PRECOMPOUND = (0!=G4CASCADE_USE_PRECOMPOUND &&
		     G4CASCADE_USE_PRECOMPOUND[0]!='0');
  DO_COALESCENCE = (0==G4CASCADE_DO_COALESCENCE ||
		    G4CASCADE_DO_COALESCENCE[0]!='0');
  SHOW_HISTORY = (0!=G4CASCADE_SHOW_HISTORY);
  USE_3BODYMOM = (0!=G4CASCADE_USE_3BODYMOM);
  USE_PHASESPACE = (0!=G4CASCADE_USE_PHASESPACE &&
		    G4CASCADE_USE_PHASESPACE[0]!='0');
  PIN_ABSORPTION = (G4CASCADE_PIN_ABSORPTION ? strtod(G4CASCADE_PIN_ABSORPTION,0)
		    : 0.);
  RANDOM_FILE = (G4CASCADE_RANDOM_FILE ? G4CASCADE_RANDOM_FILE : "");
  BEST_PAR = (0!=G4NUCMODEL_USE_BEST);
  TWOPARAM_RADIUS = (0!=G4NUCMODEL_RAD_2PAR);
  RADIUS_SCALE = (G4NUCMODEL_RAD_SCALE ? strtod(G4NUCMODEL_RAD_SCALE,0)
  		  : (BEST_PAR?1.0:OLD_RADIUS_UNITS));
  //TK for GXBERT
  //if ( G4NUCMODEL_RAD_SCALE == 0 && BEST_PAR == 0 ) HDP.DeveloperGet("BERT_RADIUS_SCALE",RADIUS_SCALE);
  RADIUS_SMALL = ((G4NUCMODEL_RAD_SMALL ? strtod(G4NUCMODEL_RAD_SMALL,0)
		   : (BEST_PAR?1.992:(8.0/OLD_RADIUS_UNITS))) * RADIUS_SCALE);
  RADIUS_ALPHA = (G4NUCMODEL_RAD_ALPHA ? strtod(G4NUCMODEL_RAD_ALPHA,0)
		  : (BEST_PAR?0.84:0.70));
  RADIUS_TRAILING = ((G4NUCMODEL_RAD_TRAILING ? strtod(G4NUCMODEL_RAD_TRAILING,0)
		      : 0.) * RADIUS_SCALE);
  //TK for GXBERT
  //if ( G4NUCMODEL_RAD_TRAILING == 0 ) HDP.DeveloperGet("BERT_RAD_TRAILING",RADIUS_TRAILING),RADIUS_TRAILING*=RADIUS_SCALE;
  FERMI_SCALE = ((G4NUCMODEL_FERMI_SCALE ? strtod(G4NUCMODEL_FERMI_SCALE,0)
		  : (BEST_PAR?0.685:(1.932/OLD_RADIUS_UNITS))) * RADIUS_SCALE);
  //TK for GXBERT
  //if ( G4NUCMODEL_FERMI_SCALE == 0 && BEST_PAR == 0 ) HDP.DeveloperGet("BERT_FERMI_SCALE",FERMI_SCALE),FERMI_SCALE*=RADIUS_SCALE;
  XSEC_SCALE = (G4NUCMODEL_XSEC_SCALE ? strtod(G4NUCMODEL_XSEC_SCALE,0)
  		: (BEST_PAR?0.1:1.0) );
  //TK for GXBERT
  //if ( G4NUCMODEL_XSEC_SCALE == 0 && BEST_PAR == 0 ) HDP.DeveloperGet("BERT_XSEC_SCALE",XSEC_SCALE);
  GAMMAQD_SCALE = (G4NUCMODEL_GAMMAQD?strtod(G4NUCMODEL_GAMMAQD,0):1.);
  DPMAX_DOUBLET = (DPMAX_2CLUSTER ? strtod(DPMAX_2CLUSTER,0) : 0.090);
  DPMAX_TRIPLET = (DPMAX_3CLUSTER ? strtod(DPMAX_3CLUSTER,0) : 0.108);
  DPMAX_ALPHA = (DPMAX_4CLUSTER ? strtod(DPMAX_4CLUSTER,0) : 0.115);
}


// Report any non-default parameters (used by G4CascadeInterface)

void GXCascadeParameters::DumpConfig(std::ostream& os) const
{
  if (G4CASCADE_VERBOSE)
    os << "G4CASCADE_VERBOSE = " << G4CASCADE_VERBOSE << endl;
  if (G4CASCADE_CHECK_ECONS)
    os << "G4CASCADE_CHECK_ECONS = " << G4CASCADE_CHECK_ECONS << endl;
  if (G4CASCADE_USE_PRECOMPOUND)
    os << "G4CASCADE_USE_PRECOMPOUND = " << G4CASCADE_USE_PRECOMPOUND << endl;
  if (G4CASCADE_DO_COALESCENCE)
    os << "G4CASCADE_DO_COALESCENCE = " << G4CASCADE_DO_COALESCENCE << endl;
  if (G4CASCADE_PIN_ABSORPTION)
    os << "G4CASCADE_PIN_ABSORPTION = " << G4CASCADE_PIN_ABSORPTION << endl;
  if (G4CASCADE_SHOW_HISTORY)
    os << "G4CASCADE_SHOW_HISTORY = " << G4CASCADE_SHOW_HISTORY << endl;
  if (G4CASCADE_USE_3BODYMOM)
    os << "G4CASCADE_USE_3BODYMOM = " << G4CASCADE_USE_3BODYMOM << endl;
  if (G4CASCADE_USE_PHASESPACE)
    os << "G4CASCADE_USE_PHASESPACE = " << G4CASCADE_USE_PHASESPACE << endl;
  if (G4CASCADE_RANDOM_FILE)
    os << "G4CASCADE_RANDOM_FILE = " << G4CASCADE_RANDOM_FILE << endl;
  if (G4NUCMODEL_USE_BEST)
    os << "G4NUCMODEL_USE_BEST = " << G4NUCMODEL_USE_BEST << endl;
  if (G4NUCMODEL_RAD_2PAR)
    os << "G4NUCMODEL_RAD_2PAR = " << G4NUCMODEL_RAD_2PAR << endl;
  if (G4NUCMODEL_RAD_SCALE)
    os << "G4NUCMODEL_RAD_SCALE = " << G4NUCMODEL_RAD_SCALE << endl;
  if (G4NUCMODEL_RAD_SMALL)
    os << "G4NUCMODEL_RAD_SMALL = " << G4NUCMODEL_RAD_SMALL << endl;
  if (G4NUCMODEL_RAD_ALPHA)
    os << "G4NUCMODEL_RAD_ALPHA = " << G4NUCMODEL_RAD_ALPHA << endl;
  if (G4NUCMODEL_RAD_TRAILING)
    os << "G4NUCMODEL_RAD_TRAILING = " << G4NUCMODEL_RAD_TRAILING << endl;
  if (G4NUCMODEL_FERMI_SCALE)
    os << "G4NUCMODEL_FERMI_SCALE = " << G4NUCMODEL_FERMI_SCALE << endl;
  if (G4NUCMODEL_XSEC_SCALE)
    os << "G4NUCMODEL_XSEC_SCALE = " << G4NUCMODEL_XSEC_SCALE << endl;
  if (G4NUCMODEL_GAMMAQD)
    os << "G4NUCMODEL_GAMMAQD = " << G4NUCMODEL_GAMMAQD << endl;
  if (DPMAX_2CLUSTER)
    os << "DPMAX_2CLUSTER = " << DPMAX_2CLUSTER << endl;
  if (DPMAX_3CLUSTER)
    os << "DPMAX_3CLUSTER = " << DPMAX_3CLUSTER << endl;
  if (DPMAX_4CLUSTER)
    os << "DPMAX_4CLUSTER = " << DPMAX_4CLUSTER << endl;
}

} // GXBERT_IMPL_NAMESPACE
} // gxbert namespace

#endif // GXCascadeParameters_hh
