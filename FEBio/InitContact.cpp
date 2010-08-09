#include "stdafx.h"
#include "fem.h"

//-----------------------------------------------------------------------------
//! Initializes contact data

bool FEM::InitContact()
{
	// set contact flag
	m_bcontact = (ContactInterfaces() > 0);

	// only proceed if there is contact
	if (!m_bcontact) return true;

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
