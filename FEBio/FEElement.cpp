// FEElement.cpp: implementation of the FEElement class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElement.h"
#include <math.h>

FESolidElement::~FESolidElement()
{
	// clean up material data
	int nint = GaussPoints();
	assert(m_State.size() == nint);
	for (int i=0; i<nint; ++i) if (m_State[i]) delete m_State[i];
	m_State.set(0);
}

FEShellElement::~FEShellElement()
{
	// clean up material data
	int nint = GaussPoints();
	assert(m_State.size() == nint);
	for (int i=0; i<nint; ++i) if (m_State[i]) delete m_State[i];
	m_State.set(0);
}

FETrussElement::~FETrussElement()
{
	// clean up material data
	int nint = GaussPoints();
	assert(m_State.size() == nint);
	for (int i=0; i<nint; ++i) if (m_State[i]) delete m_State[i];
	m_State.set(0);
}
