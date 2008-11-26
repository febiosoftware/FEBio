// FEElementLibrary.h: interface for the FEElementLibrary class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEELEMENTLIBRARY_H__3DB47576_A8D2_48BC_A48A_FD247DD84B43__INCLUDED_)
#define AFX_FEELEMENTLIBRARY_H__3DB47576_A8D2_48BC_A48A_FD247DD84B43__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/vector.h"

class FEElement;
class FEElementTraits;

//-----------------------------------------------------------------------------
//! This class stores the different element traits classes

//! The purpose of this class is to store all the different element traits classes.
//! It are these traits classes that define the different element types. All different
//! traits must be registered before they can be assigned to elements.

class FEElementLibrary  
{
public:
	//! constructor
	FEElementLibrary();

	//! destructor
	virtual ~FEElementLibrary();

	//! Function to register an element traits class
	int RegisterTraits(FEElementTraits* ptrait);

	//! Assign a traits class to an element
	static void SetElementTraits(FEElement& el, int id);

protected:
	static ptr_vector<FEElementTraits>	m_Traits;	//!< pointer to registered element traits
};

extern FEElementLibrary	elem_lib;

#endif // !defined(AFX_FEELEMENTLIBRARY_H__3DB47576_A8D2_48BC_A48A_FD247DD84B43__INCLUDED_)
