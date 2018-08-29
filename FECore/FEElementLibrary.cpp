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
		n = m_pThis->RegisterTraits(new FEHex8G8    ); assert(n==FE_HEX8G8   );
		n = m_pThis->RegisterTraits(new FEHex8RI    ); assert(n==FE_HEX8RI   );
		n = m_pThis->RegisterTraits(new FEHex8G1    ); assert(n==FE_HEX8G1   );
		n = m_pThis->RegisterTraits(new FETet4G1    ); assert(n==FE_TET4G1   );
		n = m_pThis->RegisterTraits(new FETet4G4    ); assert(n==FE_TET4G4   );
		n = m_pThis->RegisterTraits(new FEPenta6G6  ); assert(n==FE_PENTA6G6 );
        n = m_pThis->RegisterTraits(new FETet10G1   ); assert(n==FE_TET10G1  );
		n = m_pThis->RegisterTraits(new FETet10G4   ); assert(n==FE_TET10G4  );
		n = m_pThis->RegisterTraits(new FETet10G8   ); assert(n==FE_TET10G8  );
		n = m_pThis->RegisterTraits(new FETet10GL11 ); assert(n==FE_TET10GL11);
		n = m_pThis->RegisterTraits(new FETet10G4RI1); assert(n==FE_TET10G4RI1);
		n = m_pThis->RegisterTraits(new FETet10G8RI4); assert(n==FE_TET10G8RI4); 
		n = m_pThis->RegisterTraits(new FETet15G4   ); assert(n==FE_TET15G4  );
		n = m_pThis->RegisterTraits(new FETet15G8   ); assert(n==FE_TET15G8  );
		n = m_pThis->RegisterTraits(new FETet15G11  ); assert(n==FE_TET15G11 );
		n = m_pThis->RegisterTraits(new FETet15G15  ); assert(n==FE_TET15G15 );
		n = m_pThis->RegisterTraits(new FETet15G15RI4); assert(n==FE_TET15G15RI4);
		n = m_pThis->RegisterTraits(new FETet20G15  ); assert(n==FE_TET20G15 );
        n = m_pThis->RegisterTraits(new FEHex20G8   ); assert(n==FE_HEX20G8  );
		n = m_pThis->RegisterTraits(new FEHex20G27  ); assert(n==FE_HEX20G27 );
		n = m_pThis->RegisterTraits(new FEHex27G27  ); assert(n==FE_HEX27G27 );
        n = m_pThis->RegisterTraits(new FEPenta15G8 ); assert(n==FE_PENTA15G8);
        n = m_pThis->RegisterTraits(new FEPenta15G21); assert(n==FE_PENTA15G21);
		n = m_pThis->RegisterTraits(new FEPyra5G8   ); assert(n==FE_PYRA5G8  );
		n = m_pThis->RegisterTraits(new FEQuad4G4   ); assert(n==FE_QUAD4G4  );
		n = m_pThis->RegisterTraits(new FEQuad4NI   ); assert(n==FE_QUAD4NI  );
		n = m_pThis->RegisterTraits(new FETri3G1    ); assert(n==FE_TRI3G1   );
		n = m_pThis->RegisterTraits(new FETri3G3    ); assert(n==FE_TRI3G3   );
		n = m_pThis->RegisterTraits(new FETri3G7    ); assert(n==FE_TRI3G7   );
		n = m_pThis->RegisterTraits(new FETri3NI    ); assert(n==FE_TRI3NI   );
		n = m_pThis->RegisterTraits(new FETri6G3    ); assert(n==FE_TRI6G3   );
		n = m_pThis->RegisterTraits(new FETri6G4    ); assert(n==FE_TRI6G4   );
		n = m_pThis->RegisterTraits(new FETri6G7    ); assert(n==FE_TRI6G7   );
		n = m_pThis->RegisterTraits(new FETri6mG7   ); assert(n==FE_TRI6MG7  );
		n = m_pThis->RegisterTraits(new FETri6GL7   ); assert(n==FE_TRI6GL7  );
		n = m_pThis->RegisterTraits(new FETri6NI    ); assert(n==FE_TRI6NI   );
		n = m_pThis->RegisterTraits(new FETri7G3    ); assert(n==FE_TRI7G3   );
		n = m_pThis->RegisterTraits(new FETri7G4    ); assert(n==FE_TRI7G4   );
		n = m_pThis->RegisterTraits(new FETri7G7    ); assert(n==FE_TRI7G7   );
		n = m_pThis->RegisterTraits(new FETri7GL7   ); assert(n==FE_TRI7GL7  );
		n = m_pThis->RegisterTraits(new FETri10G7   ); assert(n==FE_TRI10G7  );
		n = m_pThis->RegisterTraits(new FETri10G12  ); assert(n==FE_TRI10G12 );
		n = m_pThis->RegisterTraits(new FEQuad8G9   ); assert(n==FE_QUAD8G9  );
        n = m_pThis->RegisterTraits(new FEQuad8NI   ); assert(n==FE_QUAD8NI  );
		n = m_pThis->RegisterTraits(new FEQuad9G9   ); assert(n==FE_QUAD9G9  );
        n = m_pThis->RegisterTraits(new FEQuad9NI   ); assert(n==FE_QUAD9NI  );
        n = m_pThis->RegisterTraits(new FEShellQuad4G8  ); assert(n==FE_SHELL_QUAD4G8 );
        n = m_pThis->RegisterTraits(new FEShellQuad4G12 ); assert(n==FE_SHELL_QUAD4G12);
        n = m_pThis->RegisterTraits(new FEShellQuad8G18 ); assert(n==FE_SHELL_QUAD8G18);
        n = m_pThis->RegisterTraits(new FEShellQuad8G27 ); assert(n==FE_SHELL_QUAD8G27);
        n = m_pThis->RegisterTraits(new FEShellTri3G6   ); assert(n==FE_SHELL_TRI3G6);
        n = m_pThis->RegisterTraits(new FEShellTri3G9   ); assert(n==FE_SHELL_TRI3G9);
        n = m_pThis->RegisterTraits(new FEShellTri6G14  ); assert(n==FE_SHELL_TRI6G14);
        n = m_pThis->RegisterTraits(new FEShellTri6G21  ); assert(n==FE_SHELL_TRI6G21);
        n = m_pThis->RegisterTraits(new FETrussElementTraits      ); assert(n==FE_TRUSS);
		n = m_pThis->RegisterTraits(new FEDiscreteElementTraits   ); assert(n==FE_DISCRETE);
		n = m_pThis->RegisterTraits(new FE2DTri3G1  ); assert(n==FE2D_TRI3G1 );
		n = m_pThis->RegisterTraits(new FE2DTri6G3  ); assert(n==FE2D_TRI6G3 );
		n = m_pThis->RegisterTraits(new FE2DQuad4G4 ); assert(n==FE2D_QUAD4G4);
		n = m_pThis->RegisterTraits(new FE2DQuad8G9 ); assert(n==FE2D_QUAD8G9);
		n = m_pThis->RegisterTraits(new FE2DQuad9G9 ); assert(n==FE2D_QUAD9G9);
		n = m_pThis->RegisterTraits(new FELine2G1   ); assert(n==FE_LINE2G1);
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
	return ((int)m_Traits.size()-1);
}

void FEElementLibrary::SetElementTraits(FEElement& el, int nid)
{
	el.SetTraits( GetInstance()->m_Traits[nid] );
}

FEElementTraits* FEElementLibrary::GetElementTraits(int ntype)
{ 
	return GetInstance()->m_Traits[ntype]; 
}

FE_Element_Shape FEElementLibrary::GetElementShape(int ntype)
{
	FEElementLibrary* p = GetInstance();
	if ((ntype < 0)||(ntype >= p->m_Traits.size())) return FE_ELEM_INVALID_SHAPE;
	return p->m_Traits[ntype]->Shape();
}

//! return the element class of a given element type
FE_Element_Class FEElementLibrary::GetElementClass(int ntype)
{
	FEElementLibrary* p = GetInstance();
	if ((ntype < 0) || (ntype >= p->m_Traits.size())) return FE_ELEM_INVALID_CLASS;
	return p->m_Traits[ntype]->Class();
}

bool FEElementLibrary::IsValid(const FE_Element_Spec& c)
{
	if (c.eclass == FE_ELEM_INVALID_CLASS) return false;
	if (c.eshape == FE_ELEM_INVALID_SHAPE) return false;
	if (c.etype  == FE_ELEM_INVALID_TYPE ) return false;
	if (c.eclass != GetElementClass(c.etype)) return false;
	if (c.eshape != GetElementShape(c.etype)) return false;
	return true;
}
