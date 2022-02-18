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
#include "FEElementList.h"
#include "FEMesh.h"
#include "FEDomain.h"

FEElement& FEElementList::iterator::operator*()
{ 
	return m_pmesh->Domain(m_ndom).ElementRef(m_nel); 
}

FEElement* FEElementList::iterator::operator->()
{
	return &m_pmesh->Domain(m_ndom).ElementRef(m_nel);
}

FECORE_API FEElementList::iterator::operator FEElement* ()
{
	return &m_pmesh->Domain(m_ndom).ElementRef(m_nel);
}

void FEElementList::iterator::operator ++ ()
{
	if (m_pmesh && (m_ndom >= 0) && (m_nel >= 0))
	{
		m_nel++;
		if (m_nel >= m_pmesh->Domain(m_ndom).Elements())
		{
			m_ndom++;
			m_nel = 0;
			if (m_ndom >= m_pmesh->Domains())
			{
				m_ndom = -1;
				m_nel = -1;
			}
		}
	}
}
