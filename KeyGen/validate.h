
#ifndef KEYGEN_VALIDATE_H
#define KEYGEN_VALIDATE_H 1

//
// validate.h
//
// Copyright (c) 2010 - University of Utah Software Development Center
//

#include "common.h"

////////////////////////////////////////////////////////////////////

EXPORT 
enum ValidationResult { 
	Passed, 
	Key_Length_Error,
	Failed,
	Registration_Error,
	Product_Id_Error,  
	Expired
};

EXPORT std::string validationResultToString( const ValidationResult type );

////////////////////////////////////////////////////////////////////

EXPORT ValidationResult validateKey( const std::string & key,
                                     const std::string & name,
                                     const std::string & company,
                                     const std::string & productId );

//
//  getExpirationTime()/getTimeOfGeneration()
//
//      Returns NULL if the key does not expire.  If 'verbose', prints out the expiration date.
//
EXPORT struct tm * getTimeOfGeneration ( const std::string & key, bool verbose = false );
EXPORT struct tm * getExpirationTime   ( const std::string & key, bool verbose = false );

EXPORT bool        hasExpirationDate( const std::string & key );

EXPORT std::string getProductIdentifier( const std::string & key );

#endif // !KEYGEN_VALIDATE_H
