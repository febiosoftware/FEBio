// FEElementLibrary.cpp: implementation of the FEElementLibrary class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElementLibrary.h"
#include "FEElement.h"

FEElementLibrary* FEElementLibrary::m_pThis = 0;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEElementLibrary* FEElementLibrary::GetInstance()
{
	if (m_pThis == 0)
	{
		m_pThis = new FEElementLibrary;

		// register element types
		int n;
		n = m_pThis->RegisterTraits(new FEHex8G8   ); assert(n==FE_HEX8G8   );
		n = m_pThis->RegisterTraits(new FEHex8RI   ); assert(n==FE_HEX8RI   );
		n = m_pThis->RegisterTraits(new FEHex8G1   ); assert(n==FE_HEX8G1   );
		n = m_pThis->RegisterTraits(new FETet4G1   ); assert(n==FE_TET4G1   );
		n = m_pThis->RegisterTraits(new FETet4G4   ); assert(n==FE_TET4G4   );
		n = m_pThis->RegisterTraits(new FEPenta6G6 ); assert(n==FE_PENTA6G6 );
		n = m_pThis->RegisterTraits(new FETet10G4  ); assert(n==FE_TET10G4  );
		n = m_pThis->RegisterTraits(new FETet10G8  ); assert(n==FE_TET10G8  );
		n = m_pThis->RegisterTraits(new FETet10GL11); assert(n==FE_TET10GL11);
		n = m_pThis->RegisterTraits(new FETet15G8  ); assert(n==FE_TET10G8  );
		n = m_pThis->RegisterTraits(new FEHex20G27 ); assert(n==FE_HEX20G27 );
		n = m_pThis->RegisterTraits(new FEQuad4G4  ); assert(n==FE_QUAD4G4  );
		n = m_pThis->RegisterTraits(new FEQuad4NI  ); assert(n==FE_QUAD4NI  );
		n = m_pThis->RegisterTraits(new FETri3G1   ); assert(n==FE_TRI3G1   );
		n = m_pThis->RegisterTraits(new FETri3G3   ); assert(n==FE_TRI3G3   );
		n = m_pThis->RegisterTraits(new FETri3NI   ); assert(n==FE_TRI3NI   );
		n = m_pThis->RegisterTraits(new FETri6G3   ); assert(n==FE_TRI6G3   );
		n = m_pThis->RegisterTraits(new FETri6G4   ); assert(n==FE_TRI6G4   );
		n = m_pThis->RegisterTraits(new FETri6G7   ); assert(n==FE_TRI6G7   );
		n = m_pThis->RegisterTraits(new FETri6GL7  ); assert(n==FE_TRI6GL7  );
		n = m_pThis->RegisterTraits(new FETri6NI   ); assert(n==FE_TRI6NI   );
		n = m_pThis->RegisterTraits(new FETri7G7   ); assert(n==FE_TRI7G7   );
		n = m_pThis->RegisterTraits(new FEQuad8G9  ); assert(n==FE_QUAD8G9  );
		n = m_pThis->RegisterTraits(new FEShellQuadElementTraits  ); assert(n==FE_SHELL_QUAD);
		n = m_pThis->RegisterTraits(new FEShellTriElementTraits   ); assert(n==FE_SHELL_TRI);
		n = m_pThis->RegisterTraits(new FETrussElementTraits      ); assert(n==FE_TRUSS);
		n = m_pThis->RegisterTraits(new FEDiscreteElementTraits   ); assert(n==FE_DISCRETE);
	}
	return m_pThis;
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
	el.SetTraits( GetInstance()->m_Traits[nid] );
}

FEElementTraits* FEElementLibrary::GetElementTraits(int ntype)
{ 
	return GetInstance()->m_Traits[ntype]; 
}
