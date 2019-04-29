/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include <FECore/FESurfacePairConstraint.h>
#include "febiomech_api.h"

class FEModel;
class FESolver;
class FEGlobalMatrix;

// Macauley bracket
#define MBRACKET(x) ((x)>=0? (x): 0)

// Heavyside function
#define HEAVYSIDE(x) ((x)>=0?1:0)


//-----------------------------------------------------------------------------
//! This is the base class for contact interfaces
class FEBIOMECH_API FEContactInterface : public FESurfacePairConstraint
{
public:
	//! constructor
	FEContactInterface(FEModel* pfem);

	//! destructor
	~FEContactInterface();

	//! serialize data to archive
	void Serialize(DumpStream& ar);

public:
	// The Residual function evaluates the "forces" that contribute to the residual of the system
	virtual void Residual(FEGlobalVector& R, const FETimeInfo& tp) = 0;

	// Evaluates the contriubtion to the stiffness matrix
	virtual void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) = 0;

	// Performs an augmentation step
	virtual bool Augment(int naug, const FETimeInfo& tp) = 0;

protected:
	//! don't call the default constructor
	FEContactInterface() : FESurfacePairConstraint(0){}

	//! auto-penalty calculation
	double AutoPenalty(FESurfaceElement& el, FESurface& s);

public:
	bool	m_blaugon;	//!< augmented lagrangian flag
};
