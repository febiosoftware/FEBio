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
#include "FELevelStructure.h"

//-----------------------------------------------------------------------------
//! This class implements an algoritm that calculates a permutation of 
//! the node numbering in order to obtain a bandwidth reduced stiffness matrix

//! The algorithm comes from "An algorithm for reducing the bandwidth and 
//! profile of a sparse matrix", by N.E.Gibbs e.a. It applies the algorithm
//! on the node numberings in stead of the actual sparse matrix since that
//! was easier to implement :). In the future I would like to extend it to
//! work with the actual sparse matrix. 

class FECORE_API FENodeReorder
{

public:
	//! default constructor
	FENodeReorder();

	//! destructor
	virtual ~FENodeReorder();

	//! calculates the permutation vector
	void Apply(FEMesh& m, std::vector<int>& P);
};
