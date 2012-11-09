// FEElementLibrary.cpp: implementation of the FEElementLibrary class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElementLibrary.h"
#include "FEElement.h"

std::vector<FEElementTraits*> FEElementLibrary::m_Traits;

FEElementLibrary	elem_lib;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEElementLibrary::FEElementLibrary()
{
	static bool bfirst = true;

	if (bfirst)
	{
		// register element types
		int n;
		n = RegisterTraits(new FEHex8G8  ); assert(n==FE_HEX8G8  );
		n = RegisterTraits(new FEHex8RI  ); assert(n==FE_HEX8RI  );
		n = RegisterTraits(new FEHex8G1  ); assert(n==FE_HEX8G1  );
		n = RegisterTraits(new FETet4G1  ); assert(n==FE_TET4G1  );
		n = RegisterTraits(new FETet4G4  ); assert(n==FE_TET4G4  );
		n = RegisterTraits(new FEPenta6G6); assert(n==FE_PENTA6G6);
		n = RegisterTraits(new FETet10G4 ); assert(n==FE_TET10G4 );
		n = RegisterTraits(new FETet10G8 ); assert(n==FE_TET10G8 );
		n = RegisterTraits(new FEHex20G27); assert(n==FE_HEX20G27);
		n = RegisterTraits(new FEQuad4G4 ); assert(n==FE_QUAD4G4 );
		n = RegisterTraits(new FEQuad4NI ); assert(n==FE_QUAD4NI );
		n = RegisterTraits(new FETri3G1  ); assert(n==FE_TRI3G1  );
		n = RegisterTraits(new FETri3G3  ); assert(n==FE_TRI3G3  );
		n = RegisterTraits(new FETri3NI  ); assert(n==FE_TRI3NI  );
		n = RegisterTraits(new FETri6G3  ); assert(n==FE_TRI6G3  );
		n = RegisterTraits(new FETri6G4  ); assert(n==FE_TRI6G4  );
		n = RegisterTraits(new FETri6G7  ); assert(n==FE_TRI6G7  );
		n = RegisterTraits(new FETri6GL7 ); assert(n==FE_TRI6GL7 );
		n = RegisterTraits(new FETri6NI  ); assert(n==FE_TRI6NI  );
		n = RegisterTraits(new FEQuad8G9 ); assert(n==FE_QUAD8G9 );
		n = RegisterTraits(new FEShellQuadElementTraits  ); assert(n==FE_SHELL_QUAD);
		n = RegisterTraits(new FEShellTriElementTraits   ); assert(n==FE_SHELL_TRI);
		n = RegisterTraits(new FETrussElementTraits      ); assert(n==FE_TRUSS);
		n = RegisterTraits(new FEDiscreteElementTraits   ); assert(n==FE_DISCRETE);
		bfirst = false;
	}
}

FEElementLibrary::~FEElementLibrary()
{
	for (size_t i=0; i<m_Traits.size(); ++i) delete m_Traits[i];
	m_Traits.clear();
}

int FEElementLibrary::RegisterTraits(FEElementTraits* ptrait) 
{ 
	m_Traits.push_back(ptrait); 
	return (m_Traits.size()-1);
}

void FEElementLibrary::SetElementTraits(FEElement& el, int nid)
{
	el.SetTraits( m_Traits[nid] );
}
