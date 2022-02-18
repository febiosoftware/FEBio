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
#include "FERefineMesh.h"

class FEHexRefine : public FERefineMesh
{
public:
	FEHexRefine(FEModel* fem);

	bool Init() override;

	bool RefineMesh() override;

protected:
	void BuildSplitLists(FEModel& fem);
	void UpdateNewNodes(FEModel& fem);
	void FindHangingNodes(FEModel& fem);
	void BuildNewDomains(FEModel& fem);
	void UpdateNodeSet(FENodeSet& nset);
	bool UpdateSurface(FESurface& surf);

private:
	int		m_elemRefine;		// max nr of elements to refine per step
	double	m_maxValue;

	vector<int>	m_elemList;
	vector<int>	m_edgeList;	// list of edge flags to see whether the edge was split
	vector<int>	m_faceList;	// list of face flags to see whether the face was split
	int			m_N0;
	int			m_NC;
	int			m_NN;

	int	m_splitElems;
	int m_splitFaces;
	int m_splitEdges;
	int	m_hangingNodes;

	FEMeshAdaptorCriterion*	m_criterion;

	DECLARE_FECORE_CLASS();
};
