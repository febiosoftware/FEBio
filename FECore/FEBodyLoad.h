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
#include "FEModelLoad.h"
#include "FEDomain.h"
#include "FEDomainList.h"

//-----------------------------------------------------------------------------
// forward declaration of FEModel class
class FEModel;
class FELinearSystem;

//-----------------------------------------------------------------------------
//! Base class for body-loads
class FECORE_API FEBodyLoad : public FEModelLoad
{
	FECORE_SUPER_CLASS

public:
	FEBodyLoad(FEModel* pfem);
	virtual ~FEBodyLoad();

	//! initialization
	bool Init() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

public:
	//! return number of domains this load is applied to
	int Domains() const;

	//! return a domain 
	FEDomain* Domain(int i);

	//! add a domain to which to apply this load
	void SetDomainList(FEElementSet* elset);

	//! get the domain list
	FEDomainList& GetDomainList();

public: // This should be overridden by derived classes

	//! Evaluate force vector
	virtual void ForceVector(FEGlobalVector& R);

	//! evaluate stiffness matrix
	virtual void StiffnessMatrix(FELinearSystem& S);

public:
	// NOTE: Work in progress! Working on integrating body loads as a model loads
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

private:
	FEDomainList	m_dom;	//!< list of domains to which to apply the body load
};
