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
#include "FECoreBase.h"
#include "FETimeInfo.h"

//-----------------------------------------------------------------------------
//! forward declaration of the FEModel class.
//! All classes inherited from FEModelComponent should take the model as a parameter
//! to the constructor.
class FENodeSet;
class FEMesh;

//-----------------------------------------------------------------------------
//! This class serves as a base class for many of the FECore classes. It defines
//! activation and deactivation functions which are used in multi-step analyses to determine which
//! components are active during an analysis.
//! A model component is basically anything that affects the state of a model.
//! For instance, boundary conditions, loads, contact definitions, etc.
class FECORE_API FEModelComponent : public FECoreBase
{
public:
	//! constructor
	FEModelComponent(FEModel* fem);

	//! destructor
	virtual ~FEModelComponent();

	//-----------------------------------------------------------------------------------
	//! Update the component
	//! This is called whenever the model is updated, i.e. the primary variables were updated.
	virtual void Update();

public: // some convenience functions (to pull data from FEModel without the need to include)
	double CurrentTime() const;
	double CurrentTimeIncrement() const;
	double GetGlobalConstant(const char* sz) const;
	int GetDOFIndex(const char* szvar, int n) const;
	int GetDOFIndex(const char* szdof) const;
	const FETimeInfo& GetTimeInfo() const;
	void AttachLoadController(const char* szparam, int lc);
	void AttachLoadController(void* pd, int lc);

	//! Get the model's mesh
	FEMesh& GetMesh();
};
