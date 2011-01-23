///////////////////////////////////////////////////////////////////////////////
//         ______  ______    ______
//        |  ____||  ____|  /  ____|
//        | |     | |      |  /       _______    _____    _____
//        | |__   | |___   | |       /  ___  \  /  ___|  /  __ \
//        |  __|  |  ___|  | |      |  /   \  ||  /     |  /  \ |
//        | |     | |      | |      | |     | || |      | |   /_/
//        | |     | |____  |  \____ |  \___/  || |      |  \____
//        |_|     |______|  \______| \_______/ |_|       \______|
//
///////////////////////////////////////////////////////////////////////////////
//! The FECore library is a collection of tools that simplify the developement
//! of FE software. 
//-----------------------------------------------------------------------------

#ifndef _FECORE_H_04222008_
#define _FECORE_H_04222008_

//-----------------------------------------------------------------------------
//! The FECore namespace encapsulates all classes that belong to the FECore library
namespace FECore
{
	// retrieve version numbers
	void get_version(int& version, int& subversion);

	// retrieve version number string
	const char* get_version_string();

} // namespace FECore

#endif // _FECORE_H_04222008_
