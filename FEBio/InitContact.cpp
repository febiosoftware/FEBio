#include "stdafx.h"
#include "fem.h"

//-----------------------------------------------------------------------------
//! Initializes contact data

bool FEM::InitContact()
{
	// loop over all contact interfaces
	for (int i=0; i<ContactInterfaces(); ++i)
	{
		// get the contact interface
		FEContactInterface& ci = *m_CI[i];

		// initializes contact interface data
		ci.Init();
	}

	return true;
}
