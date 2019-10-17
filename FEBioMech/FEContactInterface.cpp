/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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



#include "stdafx.h"
#include "FEContactInterface.h"
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"
#include "FECore/FESolver.h"

BEGIN_FECORE_CLASS(FEContactInterface, FESurfacePairConstraint)
	ADD_PARAMETER(m_laugon, "laugon");
END_FECORE_CLASS();

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEContactInterface::FEContactInterface(FEModel* pfem) : FESurfacePairConstraint(pfem)
{
	m_laugon = 0;	// penalty method by default
}

FEContactInterface::~FEContactInterface()
{

}

//-----------------------------------------------------------------------------
//! This function calculates a contact penalty parameter based on the 
//! material and geometrical properties of the slave and master surfaces
//!
double FEContactInterface::AutoPenalty(FESurfaceElement& el, FESurface &s)
{
	// get the mesh
	FEMesh& m = GetFEModel()->GetMesh();

	// get the element this surface element belongs to
	FEElement* pe = el.m_elem[0];
	if (pe == nullptr) return 0.0;

	// extract the elastic material
	FEElasticMaterial* pme = GetFEModel()->GetMaterial(pe->GetMatID())->ExtractProperty<FEElasticMaterial>();
	if (pme == nullptr) return 0.0;

	// get a material point
	FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

	// setup the material point
	pt.m_F = mat3dd(1.0);
	pt.m_J = 1;
	pt.m_s.zero();

	// get the tangent (stiffness) and it inverse (compliance) at this point
	tens4ds S = pme->Tangent(mp);
	tens4ds C = S.inverse();

	// evaluate element surface normal at parametric center
	vec3d t[2];
	s.CoBaseVectors0(el, 0, 0, t);
	vec3d n = t[0] ^ t[1];
	n.unit();
		
	// evaluate normal component of the compliance matrix
	// (equivalent to inverse of Young's modulus along n)
	double eps = 1./(n*(vdotTdotv(n, C, n)*n));

	// get the area of the surface element
	double A = s.FaceArea(el);

	// get the volume of the volume element
	double V = m.ElementVolume(*pe);

	return eps*A/V;
}

//-----------------------------------------------------------------------------
void FEContactInterface::Serialize(DumpStream& ar)
{
	// store base class
	FESurfacePairConstraint::Serialize(ar);

	// save parameters
	ar & m_laugon;
}
