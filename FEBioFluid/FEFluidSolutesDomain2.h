/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include <FECore/FESolidDomain.h>
#include "FEFluidDomain3D.h"
#include "FESolutesDomain.h"
#include "FEFluidSolutesMaterial2.h"
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! domain described by 3D volumetric elements
//!
class FEBIOFLUID_API FEFluidSolutesDomain2 : public FEFluidDomain3D, public FESolutesDomain
{
public:
	enum {
		FLUID_DOMAIN = 0,
		SOLUTES_DOMAIN = 1
	};

public:
	//! constructor
	FEFluidSolutesDomain2(FEModel* pfem);
	~FEFluidSolutesDomain2();

	void SetActiveDomain(int n);

	const FEDofList& GetDOFList() const override;

	void Serialize(DumpStream& ar) override;

public: // overrides from FEDomain

	//! get the material
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pm) override;

	// initialization
	bool Init() override;

	// reset 
	void Reset() override;

	void Activate() override;

	void InitMaterialPoints() override;

public: // overrides of FESolidDomain
	void PreSolveUpdate(const FETimeInfo& tp) override;

	void Update(const FETimeInfo& tp) override;

private:
	int			m_activeDomain;	// 0 = fluid, 1 = solutes
	FEFluidSolutesMaterial2*     m_pMat;
};
