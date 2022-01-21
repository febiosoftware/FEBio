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
#include "FEModelComponent.h"

//-----------------------------------------------------------------------------
// A Step component is a model component that can be assigned to a step. 
// It adds a mechanism for activating and deactivating the component.
class FECORE_API FEStepComponent : public FEModelComponent
{
public:
	FEStepComponent(FEModel* fem);

	//-----------------------------------------------------------------------------------
	//! This function checks if the component is active in the current step. 
	bool IsActive() const;

	//-----------------------------------------------------------------------------------
	//! Activate the component.
	//! This function is called during the step initialization, right before the step is solved.
	//! This function can be used to initialize any data that could depend on the model state. 
	//! Data allocation and initialization of data that does not depend on the model state should
	//! be done in Init().
	virtual void Activate();

	//-----------------------------------------------------------------------------------
	//! Deactivate the component
	virtual void Deactivate();

public:
	//! serialization
	void Serialize(DumpStream& ar);

private:
	bool		m_bactive;	//!< flag indicating whether the component is active
};
