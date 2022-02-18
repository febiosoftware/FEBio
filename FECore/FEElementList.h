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



#pragma once
#include "fecore_api.h"

class FEMesh;
class FEElement;

//-----------------------------------------------------------------------------
//! utitlity class for accessing all elements without having to go throug the domains
class FEElementList
{
public:
	class iterator
	{
	public:
		iterator() { m_pmesh = 0; m_ndom = -1; m_nel = -1; }
		iterator(FEMesh* pm) { m_pmesh = pm; m_ndom = 0; m_nel = 0; }

		FECORE_API FEElement& operator*();

		FECORE_API FEElement* operator->();

		FECORE_API operator FEElement* ();

		FECORE_API void operator ++ ();

		bool operator != (const iterator& it)
		{
			return ((m_ndom!=it.m_ndom)||(m_nel != it.m_nel));
		}

	public:
		FEMesh*		m_pmesh;	// pointer to mesh
		int			m_ndom;		// domain index
		int			m_nel;		// element index
	};

public:
	FEElementList(FEMesh& m) : m_mesh(m){}

	iterator begin() { return iterator(&m_mesh); }
	iterator end() { return iterator(); }

private:
	FEMesh&	m_mesh;
};
