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
		n = RegisterTraits(new FEHexElementTraits       ); assert(n==FE_HEX);
		n = RegisterTraits(new FEHex20ElementTraits     ); assert(n==FE_HEX20);
		n = RegisterTraits(new FERIHexElementTraits     ); assert(n==FE_RIHEX);
		n = RegisterTraits(new FEUDFHexElementTraits    ); assert(n==FE_UDGHEX);
		n = RegisterTraits(new FETetElementTraits       ); assert(n==FE_TET);
		n = RegisterTraits(new FETet10ElementTraits     ); assert(n==FE_TET10);
		n = RegisterTraits(new FEPentaElementTraits     ); assert(n==FE_PENTA);
		n = RegisterTraits(new FEG1TetElementTraits     ); assert(n==FE_TETG1);
		n = RegisterTraits(new FEQuadElementTraits      ); assert(n==FE_QUAD);
		n = RegisterTraits(new FENIQuadElementTraits    ); assert(n==FE_NIQUAD);
		n = RegisterTraits(new FETriElementTraits       ); assert(n==FE_TRI);
		n = RegisterTraits(new FENITriElementTraits     ); assert(n==FE_NITRI);
		n = RegisterTraits(new FEShellQuadElementTraits ); assert(n==FE_SHELL_QUAD);
		n = RegisterTraits(new FEShellTriElementTraits  ); assert(n==FE_SHELL_TRI);
		n = RegisterTraits(new FETrussElementTraits     ); assert(n==FE_TRUSS);
		n = RegisterTraits(new FEDiscreteElementTraits  ); assert(n==FE_DISCRETE);
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
