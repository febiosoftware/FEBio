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

class FEModel;

//-----------------------------------------------------------------------------
//! This class represents a global system array. It provides functions to assemble
//! local (element) vectors into this array
//! TODO: remove FEModel dependency!
class FECORE_API FEGlobalVector
{
public:
	//! constructor
	FEGlobalVector(FEModel& fem, std::vector<double>& R, std::vector<double>& Fr);

	//! destructor
	virtual ~FEGlobalVector();

	//! Assemble the element vector into this global vector
	virtual void Assemble(std::vector<int>& en, std::vector<int>& elm, std::vector<double>& fe, bool bdom = false);

	//! Assemble into this global vector
	virtual void Assemble(std::vector<int>& lm, std::vector<double>& fe);

	//! assemble a nodel value
	virtual void Assemble(int node, int dof, double f);
    
	//! access operator
	double& operator [] (int i) { return m_R[i]; }

	//! Get the FE model
	FEModel& GetFEModel() { return m_fem; }

	//! get the size of the vector
	int Size() const { return (int) m_R.size(); }

	operator std::vector<double>& () { return m_R; }

protected:
	FEModel&			m_fem;	//!< model
	std::vector<double>&		m_R;	//!< residual
	std::vector<double>&		m_Fr;	//!< nodal reaction forces \todo I want to remove this
};
