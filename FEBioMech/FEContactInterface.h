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
#include <FECore/FESurfacePairConstraint.h>
#include "febiomech_api.h"

class FEModel;
class FESolver;
class FEGlobalMatrix;
class FEContactSurface;

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
	void Serialize(DumpStream& ar) override;

public:
	// The LoadVector function evaluates the "forces" that contribute to the residual of the system
	virtual void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) = 0;

	// Evaluates the contriubtion to the stiffness matrix
	virtual void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) = 0;

protected:
	//! don't call the default constructor
	FEContactInterface() : FESurfacePairConstraint(0){}

	//! auto-penalty calculation
	double AutoPenalty(FESurfaceElement& el, FESurface& s);

	// serialize the element pointers
	void SerializeElementPointers(FEContactSurface& ss, FEContactSurface& ms, DumpStream& ar);

    //! cale the penalty factor during Lagrange augmentation
    double GetPenaltyScaleFactor();
    
public:
	int		m_laugon;	//!< contact enforcement method
    double  m_psf;      //!< penalty scale factor during Lagrange augmentation
    double  m_psfmax;   //!< max allowable penalty scale factor during laugon

	DECLARE_FECORE_CLASS();
};
