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
#include "stdafx.h"
#include "FEInitialPreStrain.h"
#include <FECore/FEDomain.h>
#include <FECore/FEMesh.h>
#include "FEConstPrestrain.h"


BEGIN_FECORE_CLASS(FEInitialPreStrain, FEInitialCondition)
	ADD_PARAMETER(m_binit , "init" );
	ADD_PARAMETER(m_breset, "reset");
END_FECORE_CLASS();

FEInitialPreStrain::FEInitialPreStrain(FEModel* pfem) : FEInitialCondition(pfem)
{
	m_binit = true;
	m_breset = true;
}

void FEInitialPreStrain::Activate()
{
	FEInitialCondition::Activate();

	FEMesh& mesh = GetMesh();

	// loop over all the domains
	int ND = mesh.Domains();
	for (int i=0; i<ND; ++i)
	{
		FEDomain& dom = mesh.Domain(i);

		// see if this material is the right type
		FEPrestrainMaterial* pmat = dynamic_cast<FEPrestrainMaterial*>(dom.GetMaterial());
		if (pmat)
		{
			// get the prestrain gradient
			FEPrestrainGradient* Fp = pmat->PrestrainGradientProperty();
			if (Fp)
			{
				// loop over the elements
				int NE = dom.Elements();
				for (int i=0; i<NE; ++i)
				{
					FEElement& el = dom.ElementRef(i);

					// loop over the integration points
					const int nint = el.GaussPoints();
					for (int n=0; n<nint; ++n)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(n);
						FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
						FEPrestrainMaterialPoint& pp = *mp.ExtractData<FEPrestrainMaterialPoint>();

						// get the deformation gradient
						mat3d& F = pt->m_F;

						// initialize prestrain material point data
						if (m_binit) Fp->Initialize(F, mp);
						else
						{
							mat3d Fc = F*pp.PrestrainCorrection();
							pp.setPrestrainCorrection(Fc);
						}

						// reset deformation gradient
						F.unit();
					}
				}
			}
		}
	}

	// reset the nodal coordinates
	if (m_breset)
	{
		int dofX = GetDOFIndex("x");
		int dofY = GetDOFIndex("y");
		int dofZ = GetDOFIndex("z");
		int NN = mesh.Nodes();
		for (int i=0; i<NN; ++i)
		{
			FENode& node = mesh.Node(i);
			node.m_r0 = node.m_rt;
			node.set(dofX, 0.0);
			node.set(dofY, 0.0);
			node.set(dofZ, 0.0);
		}
	}
}
