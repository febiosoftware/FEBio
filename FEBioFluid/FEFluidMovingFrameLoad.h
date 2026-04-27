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
#include <FECore/quatd.h>
#include <FECore/FEFunction1D.h>

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
	vec3d	m_at;	// linear acceleration of moving frame in fixed frame
	vec3d	m_wi;	// angular velocity input
	int		m_frame;	// fixed (=0), or moving (=1) frame for input angular velocity

private:
	vec3d	m_ap, m_a, m_A;
	vec3d	m_wp, m_w, m_Wp, m_W;
	vec3d	m_alt, m_alp, m_al, m_Alt, m_Alp, m_Al;
	quatd	m_qt, m_qp, m_q;

	// output parameters
	vec3d m_wt; // angular velocity in fixed frame (output)
	vec3d m_Wt; // angular velocity in moving frame (output)
	vec3d m_rt; // rotational vector in fixed frame (output)

	DECLARE_FECORE_CLASS();
};
