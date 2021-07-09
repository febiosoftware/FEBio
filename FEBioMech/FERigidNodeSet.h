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
#include <FECore/FEBoundaryCondition.h>
#include "febiomech_api.h"

class FENodeSet;

//-----------------------------------------------------------------------------
//! rigid node set
class FEBIOMECH_API FERigidNodeSet : public FEBoundaryCondition
{
public:
	enum SHELL_BC {
		HINGED_SHELL,
		CLAMPED_SHELL
	};

public:
	FERigidNodeSet(FEModel* pfem);
	FERigidNodeSet(const FERigidNodeSet& rs);
	void operator = (const FERigidNodeSet& rs);

	bool Init() override;

	void Serialize(DumpStream& ar) override;

	void Activate() override;
	void Deactivate() override;

	void SetRigidMaterialID(int rid);

	void SetShellBC(SHELL_BC bc);

	// copy data from another class
	void CopyFrom(FEBoundaryCondition* pbc) override;

public: // from FEModelComponent
	void SetNodeSet(FENodeSet* ns) override;

private: // parameters
	int			m_rigidMat;		//!< rigid body's material
	int			m_nshellBC;		//!< flag defining how shells are attached (0=hinged, 1=clamped)

private:
	FENodeSet* m_nodeSet;

	DECLARE_FECORE_CLASS();
};
