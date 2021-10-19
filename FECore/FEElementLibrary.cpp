/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEElementLibrary.h"
#include "FEElement.h"
#include "FESolidElementShape.h"
#include "FESurfaceElementShape.h"

FEElementLibrary* FEElementLibrary::m_pThis = 0;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//! initialize library
void FEElementLibrary::Initialize()
{
	// Calling GetInstance will initialize the static pointer
	if (m_pThis == 0) GetInstance();
}

FEElementLibrary* FEElementLibrary::GetInstance()
{
	if (m_pThis == 0)
	{
		m_pThis = new FEElementLibrary;

		int n;
		// register element shapes (must be done before types!)
		n = m_pThis->RegisterShape(new FETet4   ); assert(n == ET_TET4   );
		n = m_pThis->RegisterShape(new FETet5   ); assert(n == ET_TET5   );
		n = m_pThis->RegisterShape(new FETet10  ); assert(n == ET_TET10  );
		n = m_pThis->RegisterShape(new FETet15  ); assert(n == ET_TET15  );
		n = m_pThis->RegisterShape(new FETet20  ); assert(n == ET_TET20  );
		n = m_pThis->RegisterShape(new FEPenta6 ); assert(n == ET_PENTA6 );
		n = m_pThis->RegisterShape(new FEPenta15); assert(n == ET_PENTA15);
		n = m_pThis->RegisterShape(new FEHex8   ); assert(n == ET_HEX8   );
		n = m_pThis->RegisterShape(new FEHex20  ); assert(n == ET_HEX20  );
		n = m_pThis->RegisterShape(new FEHex27  ); assert(n == ET_HEX27  );
		n = m_pThis->RegisterShape(new FEPyra5  ); assert(n == ET_PYRA5  );
        n = m_pThis->RegisterShape(new FEPyra13 ); assert(n == ET_PYRA13 );
		n = m_pThis->RegisterShape(new FEQuad4  ); assert(n == ET_QUAD4  );
		n = m_pThis->RegisterShape(new FEQuad8  ); assert(n == ET_QUAD8  );
		n = m_pThis->RegisterShape(new FEQuad9  ); assert(n == ET_QUAD9  );
		n = m_pThis->RegisterShape(new FETri3   ); assert(n == ET_TRI3   );
		n = m_pThis->RegisterShape(new FETri6   ); assert(n == ET_TRI6   );
		n = m_pThis->RegisterShape(new FETri7   ); assert(n == ET_TRI7   );
		n = m_pThis->RegisterShape(new FETri10  ); assert(n == ET_TRI10  );

		// register element types
		n = m_pThis->RegisterTraits(new FEHex8G8    ); assert(n==FE_HEX8G8   );
		n = m_pThis->RegisterTraits(new FEHex8RI    ); assert(n==FE_HEX8RI   );
		n = m_pThis->RegisterTraits(new FEHex8G1    ); assert(n==FE_HEX8G1   );
		n = m_pThis->RegisterTraits(new FETet4G1    ); assert(n==FE_TET4G1   );
		n = m_pThis->RegisterTraits(new FETet4G4    ); assert(n==FE_TET4G4   );
		n = m_pThis->RegisterTraits(new FETet5G4    ); assert(n==FE_TET5G4   );
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
        n = m_pThis->RegisterTraits(new FEPyra13G8  ); assert(n==FE_PYRA13G8 );
		n = m_pThis->RegisterTraits(new FEQuad4G4   ); assert(n==FE_QUAD4G4  );
		n = m_pThis->RegisterTraits(new FEQuad4NI   ); assert(n==FE_QUAD4NI  );
		n = m_pThis->RegisterTraits(new FETri3G1    ); assert(n==FE_TRI3G1   );
		n = m_pThis->RegisterTraits(new FETri3G3    ); assert(n==FE_TRI3G3   );
		n = m_pThis->RegisterTraits(new FETri3G7    ); assert(n==FE_TRI3G7   );
		n = m_pThis->RegisterTraits(new FETri3NI    ); assert(n==FE_TRI3NI   );
		n = m_pThis->RegisterTraits(new FETri6G3    ); assert(n==FE_TRI6G3   );
		n = m_pThis->RegisterTraits(new FETri6G4    ); assert(n==FE_TRI6G4   );
		n = m_pThis->RegisterTraits(new FETri6G7    ); assert(n==FE_TRI6G7   );
//		n = m_pThis->RegisterTraits(new FETri6mG7   ); assert(n==FE_TRI6MG7  );
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

int FEElementLibrary::RegisterShape(FEElementShape* pshape)
{
	m_Shape.push_back(pshape);
	return ((int)m_Shape.size() - 1);
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

//! return element shape class
FEElementShape* FEElementLibrary::GetElementShapeClass(FE_Element_Shape eshape)
{
	return GetInstance()->m_Shape[eshape];
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

//! get the element spec from the type
FE_Element_Spec FEElementLibrary::GetElementSpecFromType(FE_Element_Type elemType)
{
	FE_Element_Spec espec;
	espec.etype = elemType;
	if (elemType != FE_ELEM_INVALID_TYPE)
	{
		espec.eclass = GetElementClass(elemType);
		espec.eshape = GetElementShape(elemType);
	}
	return espec;
}
