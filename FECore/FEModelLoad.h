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
#include "FEStepComponent.h"
#include "FEGlobalVector.h"
#include <FECore/FEDofList.h>

//-----------------------------------------------------------------------------
class FELinearSystem;

//-----------------------------------------------------------------------------
//! This class is the base class for all classes that affect the state of the model
//! and contribute directly to the residual and the global stiffness matrix. This
//! includes most boundary loads, body loads, contact, etc.
class FECORE_API FEModelLoad : public FEStepComponent
{
	FECORE_SUPER_CLASS(FELOAD_ID)
	FECORE_BASE_CLASS(FEModelLoad)

public:
	//! constructor
	FEModelLoad(FEModel* pfem);

	const FEDofList& GetDofList() const;
	
	void Serialize(DumpStream& ar) override;

public:
	// all classes derived from this base class must implement
	// the following functions.

	//! evaluate the contribution to the external load vector
	virtual void LoadVector(FEGlobalVector& R);

	//! evaluate the contribution to the global stiffness matrix
	virtual void StiffnessMatrix(FELinearSystem& LS);

protected:
	FEDofList	m_dof;
};
