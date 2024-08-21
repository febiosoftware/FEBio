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
#include "FEContactGapCriterion.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include "FEContactInterface.h"
#include "FEContactSurface.h"

BEGIN_FECORE_CLASS(FEContactGapCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_gap, "gap");
END_FECORE_CLASS();

FEContactGapCriterion::FEContactGapCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_gap = 0;
}

FEMeshAdaptorSelection FEContactGapCriterion::GetElementSelection(FEElementSet* elset)
{
	FEMeshAdaptorSelection sel;

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
		if (pci && pci->IsActive())
		{
			FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(pci->GetPrimarySurface());
			if (pcs)
			{
				for (int j = 0; j < pcs->Elements(); ++j)
				{
					FESurfaceElement& face = pcs->Element(j);
					if (face.isActive())
					{
						int nint = face.GaussPoints();
						double max_gap = m_gap;
						for (int k = 0; k < nint; ++k)
						{
							FEMaterialPoint* mp = face.GetMaterialPoint(k);
							FEContactMaterialPoint* cp = dynamic_cast<FEContactMaterialPoint*>(mp);
							if (cp)
							{
								if (cp->m_gap > max_gap) max_gap = cp->m_gap;
							}
						}

						if (max_gap > m_gap)
						{
							FEElement* pe = face.m_elem[0].pe;
							if (pe)
							{
								// TODO: Add check for duplicate element addition 
								sel.push_back({ pe->GetID(), max_gap});
							}
						}
					}
				}
			}
		}
	}

	return sel;
}
