//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: GXException.cc 67970 2013-03-13 10:10:06Z gcosmo $
//
// 
// ----------------------------------------------------------------------
// GXException
//
// Global error function prints string to G4cerr (or G4cout in case of
// warning). May abort program according to severity.
// ----------------------------------------------------------------------

#include "globals.hh"
//#include "G4ios.hh"
//#include "G4String.hh"
//#include "GXExceptionSeverity.hh"
//typedef std::ostringstream GXExceptionDescription;
//#include "G4StateManager.hh"

void GXException(const char* originOfException,
                 const char* exceptionCode,
                             G4ExceptionSeverity severity,
                 const char* description)
{
  //G4VExceptionHandler* exceptionHandler
  //  = G4StateManager::GetStateManager()->GetExceptionHandler();
  G4bool toBeAborted = true;
  //if(exceptionHandler)
  //{
  //  toBeAborted = exceptionHandler
  //   ->Notify(originOfException,exceptionCode,severity,description);
  //}
  //else
  //{
    static const G4String es_banner
      = "\n-------- EEEE ------- GXException-START -------- EEEE -------\n";
    static const G4String ee_banner
      = "\n-------- EEEE -------- GXException-END --------- EEEE -------\n";
    static const G4String ws_banner
      = "\n-------- WWWW ------- GXException-START -------- WWWW -------\n";
    static const G4String we_banner
      = "\n-------- WWWW -------- GXException-END --------- WWWW -------\n";
    std::ostringstream message;
    message << "\n*** ExceptionHandler is not defined ***\n"
            << "*** GXException : " << exceptionCode << G4endl
            << "      issued by : " << originOfException << G4endl
            << description << G4endl;
    switch(severity)
    {
     case FatalException:
      G4cerr << es_banner << message.str() << "*** Fatal Exception ***"
             << ee_banner << G4endl;
      break;
     case FatalErrorInArgument:
      G4cerr << es_banner << message.str() << "*** Fatal Error In Argument ***"
             << ee_banner << G4endl;
      break;
     case RunMustBeAborted:
      G4cerr << es_banner << message.str() << "*** Run Must Be Aborted ***"
             << ee_banner << G4endl;
      break;
     case EventMustBeAborted:
      G4cerr << es_banner << message.str() << "*** Event Must Be Aborted ***"
             << ee_banner << G4endl;
      break;
     default:
      G4cout << ws_banner << message.str()
             << "*** This is just a warning message. ***"
             << we_banner << G4endl;
      toBeAborted = false;
      break;
    }
  //}
  if(toBeAborted)
  {
   //if(G4StateManager::GetStateManager()->SetNewState(G4State_Abort))
   //{
     G4cerr << G4endl << "*** GXException: Aborting execution ***" << G4endl;
     abort();
   //}
   //else
   //{
     G4cerr << G4endl << "*** GXException: Abortion suppressed ***"
            << G4endl << "*** No guarantee for further execution ***" << G4endl;
   //}
  }
}

void GXException(const char* originOfException,
                 const char* exceptionCode,
                 G4ExceptionSeverity severity,
                 GXExceptionDescription & description)
{
  G4String des = description.str();
  GXException(originOfException, exceptionCode, severity, des.c_str());
}

void GXException(const char* originOfException,
                 const char* exceptionCode,
                 G4ExceptionSeverity severity,
                 GXExceptionDescription & description,
                 const char* comments)
{
  description << comments << G4endl;
  GXException(originOfException, exceptionCode, severity, description);
}
