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
#include <FEBioMech/FEBodyForce.h>

//! This body load emulates the effect of a moving frame
class FEFluidMovingFrameLoad : public FEBodyForce
{
public:
	FEFluidMovingFrameLoad(FEModel* fem);

	void Activate() override;

	void PrepStep() override;

	vec3d force(FEMaterialPoint& pt) override;

	mat3d stiffness(FEMaterialPoint& pt) override;

private:
	vec3d	m_w;	// angular velocity of frame
	vec3d	m_a;	// linear acceleration of frame
	
	quatd	m_q;
	vec3d	m_wp;
	vec3d	m_al;

	DECLARE_FECORE_CLASS();
};
