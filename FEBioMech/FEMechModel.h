/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include <FECore/FEModel.h>
#include "febiomech_api.h"

//---------------------------------------------------------------------------------------
class FERigidSystem;

//---------------------------------------------------------------------------------------
// This class extends the basic FEModel class by adding a rigid body system
class FEBIOMECH_API FEMechModel : public FEModel
{
public:
	FEMechModel();

	// get the rigid system
	FERigidSystem* GetRigidSystem();

	// clear all model data
	void Clear() override;

	// initialize the rigid system
	bool InitRigidSystem() override;

	// model activation
	void Activate() override;

	// reset
	bool Reset() override;

	//! Initialize shells
	void InitShells() override;

	// find a parameter value
	FEParamValue GetParameterValue(const ParamString& param) override;

	//! serialize data for restarts
	void SerializeGeometry(DumpStream& ar) override;

	//! Build the matrix profile for this model
	void BuildMatrixProfile(FEGlobalMatrix& G, bool breset) override;

private:
	FERigidSystem*	m_prs;

	DECLARE_FECORE_CLASS();
};
