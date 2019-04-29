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
#include "FEDomain.h"
#include "FEMaterial.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FEDomain::FEDomain(int nclass, FEModel* fem) : FEMeshPartition(nclass, fem)
{

}

//-----------------------------------------------------------------------------
void FEDomain::SetMaterial(FEMaterial* pm)
{
	assert(pm);
	if (pm) pm->AddDomain(this);
}

//-----------------------------------------------------------------------------
void FEDomain::SetMatID(int mid)
{
	ForEachElement([=](FEElement& el) { el.SetMatID(mid); });
}

//-----------------------------------------------------------------------------
// This routine allocates the material point data for the element's integration points.
// Currently, this has to be called after the elements have been assigned a type (since this
// determines how many integration point an element gets). 
void FEDomain::CreateMaterialPointData()
{
	FEMaterial* pmat = GetMaterial();
	if (pmat) ForEachElement([=](FEElement& el) {
		for (int k = 0; k<el.GaussPoints(); ++k) el.SetMaterialPointData(pmat->CreateMaterialPointData(), k);
	});
}

//-----------------------------------------------------------------------------
// serialization
void FEDomain::Serialize(DumpStream& ar)
{
	FEMeshPartition::Serialize(ar);

	if (ar.IsShallow())
	{
		int NEL = Elements();
		for (int i = 0; i < NEL; ++i)
		{
			FEElement& el = ElementRef(i);
			el.Serialize(ar);
			int nint = el.GaussPoints();
			for (int j = 0; j < nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
		}
	}
	else
	{
		if (ar.IsSaving())
		{
			FEMaterial* mat = GetMaterial();
			ar << mat;

			int NEL = Elements();
			ar << NEL;
			for (int i = 0; i < NEL; ++i)
			{
				FEElement& el = ElementRef(i);
				el.Serialize(ar);
				int nint = el.GaussPoints();
				for (int j = 0; j < nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
		else
		{
			FEMaterial* pmat = 0;
			ar >> pmat;
			SetMaterial(pmat);

			int NEL = 0;
			ar >> NEL;
			Create(NEL);
			for (int i = 0; i < NEL; ++i)
			{
				FEElement& el = ElementRef(i);
				el.Serialize(ar);
				int nint = el.GaussPoints();
				for (int j = 0; j < nint; ++j)
				{
					el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
					el.GetMaterialPoint(j)->Serialize(ar);
				}
			}
		}
	}
}
