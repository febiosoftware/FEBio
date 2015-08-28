
#ifndef KEYGEN_GENERATE_H
#define KEYGEN_GENERATE_H 1

//
// generate.h
//
// Copyright (c) 2010 - University of Utah Software Development Center
//

#include "common.h"

//
// generateKey()
//
//   Produces a 32 character string 'hash' of the input parameters.
//   This string serves as a key to 'validate' the software and encodes
//   the input information for retrieval when needed.
//
//   * 'name' - arbitary string representing the person's name
//   * 'company' - arbitary string representing the person's company's name
//
//   * 'expYear'  - These three integers represent the date that the 
//   * 'expMonth'   key expires.  Year is YYYY, Month is 1-12, and day is 
//   * 'expDay'     1-31.
//
//   * 'productIdentifier' - 
//
//

EXPORT
std::string
generateKey( const std::string & name, 
             const std::string & company,
             const std::string & productIdentifier,
             bool                verbose = false );

EXPORT
std::string
generateKey( const std::string & name, 
             const std::string & company,
             const std::string & productIdentifier,
             int                 expYear,
             int                 expMonth,
             int                 expDay,
             bool                verbose = false );


//
// checkDate()
//
//    Makes sure that month = 1-12, day = 1-31, and year = 2010-.
//    Throws exception (std::string) if this is not the case.
//
EXPORT
void
checkDate( const std::string & month,
           const std::string & day,
           const std::string & year );

#endif // KEYGEN_GENERATE_H
