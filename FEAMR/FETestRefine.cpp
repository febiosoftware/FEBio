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
#include "FETestRefine.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEMeshTopo.h>
#include <FECore/FEFixedBC.h>
#include <FECore/FESurface.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FETestRefine, FERefineMesh)
END_FECORE_CLASS();

FETestRefine::FETestRefine(FEModel* fem) : FERefineMesh(fem)
{
}

struct TRI
{
	int n[3];
};

bool FETestRefine::RefineMesh()
{
	FEModel& fem = *GetFEModel();

	if (m_topo == nullptr) return false;
	FEMeshTopo& topo = *m_topo;

	FEMesh& mesh = fem.GetMesh();

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(0));
	int NEL = dom.Elements();
	int N0 = mesh.Nodes();

	// now we recreate the domains
	const int NDOM = mesh.Domains();
	for (int i = 0; i < NDOM; ++i)
	{
		// get the old domain
		FEDomain& oldDom = mesh.Domain(i);
		int NE0 = oldDom.Elements();

		// create a copy of old domain (since we want to retain the old domain)
		FEDomain* newDom = fecore_new<FESolidDomain>(oldDom.GetTypeStr(), &fem);
		newDom->Create(NE0, FEElementLibrary::GetElementSpecFromType(FE_TET4G4));
		for (int j = 0; j < NE0; ++j)
		{
			FEElement& el0 = oldDom.ElementRef(j);
			FEElement& el1 = newDom->ElementRef(j);
			for (int k = 0; k < el0.Nodes(); ++k) el1.m_node[k] = el0.m_node[k];
		}

		// reallocate the old domain
		oldDom.Create(NE0, FEElementLibrary::GetElementSpecFromType(FE_TET4G4));

		// set new element nodes
		int nel = 0;
		for (int j = 0; j < NE0; ++j)
		{
			FEElement& el0 = newDom->ElementRef(j);
			FEElement& el1 = oldDom.ElementRef(nel++);

			el1.m_node[0] = el0.m_node[0];
			el1.m_node[1] = el0.m_node[1];
			el1.m_node[2] = el0.m_node[2];
			el1.m_node[3] = el0.m_node[3];
		}

		// we don't need this anymore
		delete newDom;
	}
	mesh.RebuildLUT();

	// re-init domains
	for (int i = 0; i < NDOM; ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.CreateMaterialPointData();
		dom.Reset();	// NOTE: we need to call this to actually call the Init function on the material points.
		dom.Init();
		dom.Activate();
	}

	// recreate element sets
	for (int i = 0; i < mesh.ElementSets(); ++i)
	{
		FEElementSet& eset = mesh.ElementSet(i);

		// get the domain list
		// NOTE: Don't get the reference, since then the same reference
		// is passed to Create below, which causes problems.
		FEDomainList domList = eset.GetDomainList();
		if (domList.IsEmpty()) { throw std::runtime_error("Error in FEMMGRemesh!"); }

		// recreate the element set from the domain list
		eset.Create(domList);
	}

	// recreate surfaces
	int faceMark = 1;
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);
		int NF0 = surf.Elements();

		vector<int> faceList = topo.FaceIndexList(surf);
		assert(faceList.size() == NF0);

		vector<TRI> tri(NF0);
		for (int j = 0; j < NF0; ++j)
		{
			FESurfaceElement& el = surf.Element(j);
			TRI t;
			t.n[0] = el.m_node[0];
			t.n[1] = el.m_node[1];
			t.n[2] = el.m_node[2];
			tri[j] = t;
		}

		int NF1 = NF0;
		surf.Create(NF1);
		int nf = 0;
		for (int j = 0; j < NF0; ++j)
		{
			TRI& t = tri[j];

			FESurfaceElement& fj = surf.Element(nf++);
			fj.SetType(FE_TRI3G3);
			fj.m_node[0] = t.n[0];
			fj.m_node[1] = t.n[1];
			fj.m_node[2] = t.n[2];
		}

		surf.CreateMaterialPointData();
		surf.Init();

		// also update the facet set if the surface has one
		FEFacetSet* fset = surf.GetFacetSet();
		if (fset)
		{
			fset->Create(surf);
		}

		faceMark++;
	}

	// re-activate the model
	UpdateModel();

	return true;
}
