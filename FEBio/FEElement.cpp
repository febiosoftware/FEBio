// FEElement.cpp: implementation of the FEElement class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElement.h"
#include <math.h>

FESolidElement::~FESolidElement()
{
	if (m_pT)
	{
		// clean up material data
		int nint = GaussPoints();
		assert(m_State.size() == nint);
		for (int i=0; i<nint; ++i) if (m_State[i]) delete m_State[i];
		zero(m_State);
	}
}

FEShellElement::~FEShellElement()
{
	if (m_pT)
	{
		// clean up material data
		int nint = GaussPoints();
		assert(m_State.size() == nint);
		for (int i=0; i<nint; ++i) if (m_State[i]) delete m_State[i];
		zero(m_State);
	}
}

FETrussElement::~FETrussElement()
{
	if (m_pT)
	{
		// clean up material data
		int nint = GaussPoints();
		assert(m_State.size() == nint);
		for (int i=0; i<nint; ++i) if (m_State[i]) delete m_State[i];
		zero(m_State);
	}
}
