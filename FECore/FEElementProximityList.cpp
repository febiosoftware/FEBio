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
#include "FEElementProximityList.h"
#include "FEMeshTopo.h"
#include "FEDomain.h"
#include "log.h"
#include <iostream>

FEElementProximityList::FEElementProximityList()
{

}

FEElementProximityList::~FEElementProximityList()
{
	delete m_lut;
}

bool FEElementProximityList::Create(FEMesh& mesh, double R)
{
	m_lut = new FEElementLUT(mesh);

	FEMeshTopo topo;
	cout << "Evaluating mesh topology...";
	if (topo.Create(&mesh) == false)
	{
		cout << "Failed building mesh topo.";
		return false;
	}
	cout << "Evaluating element proximity...";
	int NE = topo.Elements();
	m_EPL.resize(NE);
	double eplmin = 1e9;
	double eplmax = 0;
	double eplavg = 0;
	for (int i = 0; i < NE; ++i) {
		std::vector<FEElement*> epl = topo.ElementProximityList(i, R, false);
		m_EPL[i] = epl;
		double epls = epl.size();
		eplmin = std::min(eplmin, epls);
		eplmax = std::max(eplmax, epls);
		eplavg += epls;
	}
	eplavg /= NE;
	printf("Done.");
	printf("Number of neighboring elements:\n");
	printf("-------------------------------\n");
	printf("Min. = %g, Max. = %g, Avg. = %g\n", eplmin, eplmax, eplavg);
}

bool FEElementProximityList::Create(FEMesh& mesh, FEDomainList& domainList, double R)
{
	m_lut = new FEElementLUT(mesh, domainList);
	int NE = (int)m_lut->Size();
	if (NE == 0) return false;

	FEMeshTopo topo;
	cout << "Evaluating mesh topology...";
	if (topo.Create(&mesh) == false)
	{
		cout << "Failed building mesh topo.";
		return false;
	}
	cout << "Evaluating element proximity...";
	m_EPL.resize(NE);
	double eplmin = 1e9;
	double eplmax = 0;
	double eplavg = 0;
	for (int n = 0; n < domainList.size(); ++n) {
		FEDomain& dom = *domainList[n];
		int ne = dom.Elements();
		for (int i = 0; i < ne; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			int m = m_lut->FindIndex(el.GetID());
			std::vector<FEElement*> epl = topo.ElementProximityList(m, R, false);
			m_EPL[m] = epl;
			double epls = epl.size();
			eplmin = std::min(eplmin, epls);
			eplmax = std::max(eplmax, epls);
			eplavg += epls;
		}
	}
	eplavg /= NE;
	printf("Done.");
	printf("Number of neighboring elements:\n");
	printf("-------------------------------\n");
	printf("Min. = %g, Max. = %g, Avg. = %g\n", eplmin, eplmax, eplavg);
}
