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
#include "FEBodyForce.h"
#include <FECore/FESolidElement.h>

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEPointBodyForce : public FEBodyForce
{
public:
	FEPointBodyForce(FEModel* pfem);

	vec3d force(FEMaterialPoint& mp) override;
	double divforce(FEMaterialPoint& mp) override;
	mat3d stiffness(FEMaterialPoint& mp) override;

	void Serialize(DumpStream& ar) override;

	bool Init() override;
	void Update() override;

public:
	double	m_a, m_b;           //!< coefficients of exponential decay of body force
	vec3d	m_rc;               //!< center point of body force
	
	int		m_inode;            //!< node number of center of body force, or -1 if not a node

	bool	m_brigid;           //!< flag if center point is located within a rigid body

	FESolidElement* m_pel;		//!< element in which point m_r0 lies
	double			m_rs[3];	//!< isoparametric coordinates

	DECLARE_FECORE_CLASS();
};
