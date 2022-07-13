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
#include <vector>
#include "fecore_api.h"

//-----------------------------------------------------------------------------
class FEMesh;
class FESurface;
class FEElement;

//-----------------------------------------------------------------------------
//! This class finds for each element the neighbouring elements
//!
class FECORE_API FEElemElemList
{
public:
	//! constructor
	FEElemElemList(void);

	//! destructor
	~FEElemElemList(void);

	//! create the element-element list
	bool Create(FEMesh* pmesh);

	//! create the element-element list for a surface
	bool Create(const FESurface* psurf);

	//! Find the j-th neighbor element of element n
	FEElement* Neighbor(int n, int j) { return m_pel[ m_ref[n] + j]; }

	//! Find the j-th neighbor element of element n
	int NeighborIndex(int n, int j) { return m_peli[m_ref[n] + j]; }

	//! Return the size of the neighbor vector
	int NeighborSize() { return m_pel.size()/m_ref.size(); }

protected:
	//! Initialization
	void Init();

protected:
	std::vector<int>		m_ref;		//!< start index into pel and peli array
	std::vector<FEElement*>	m_pel;	//!< list of all neighbouring elements (or 0 if no neighbor)
	std::vector<int>		m_peli;	//!< indices of neighbor elems (or -1 if no neighbor)
	FEMesh*	m_pmesh;				//!< pointer to mesh that created this list
};
