// fem.h: interface for the FEM class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _FEM_H_07012006_
#define _FEM_H_07012006_

#include "FEBioLib/FEBioModel.h"

//-----------------------------------------------------------------------------
//! The Finite Element Model class. 

//! This class stores solver parameters, geometry data, material data, and 
//! other data that is needed to solve the FE problem.
//! FEBio is designed to solve finite element problems. All the finite element
//! data is collected here in this class. This class also provides
//! routines to initalize, input, output and update the FE data. Although this
//! class provides the main solve routine it does not really solve anything.
//! The actual solving is done by one of the classes derived from the FESolver class.

class FEM : public FEBioModel
{
public:
	//! constructor - sets default variables
	FEM();

	//! destructor
	virtual ~FEM();

	//! read the configuration file
	bool Configure(const char* szfile);

	//! Restart from restart point
	bool Restart(const char* szfile);

	//! check for user-interruption
	virtual void CheckInterruption();

public:
	virtual void PushState();
	virtual void PopState ();

protected:
	void ShallowCopy(FEM& fem);

public:
	bool	m_bInterruptable;	//!< true if this model can be interrupted with ctrl+c
};

#endif // _FEM_H_07012006_
