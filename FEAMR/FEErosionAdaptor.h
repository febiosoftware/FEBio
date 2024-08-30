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
#include <FECore/FEMeshAdaptor.h>
#include <FECore/FEMeshAdaptorCriterion.h>

class FEMeshTopo;

class FEErosionAdaptor : public FEMeshAdaptor
{
	enum SurfaceErodeOption {
		DONT_ERODE,		// don't erode surface elements (default)
		ERODE,			// erode surface elements
		GROW,			// add newly exposed surface elements
		RECONSTRUCT		// reconstruct surfaces
	};

public:
	FEErosionAdaptor(FEModel* fem);

	bool Apply(int iteration) override;

private:
	void RemoveIslands(FEMeshTopo& topo);
	void DeactivateOrphanedNodes();
	void ErodeSurfaces();
	void GrowErodedSurfaces(FEMeshTopo& topo);
	void ReconstructSurfaces(FEMeshTopo& topo);
	void UpdateLinearConstraints();

private:
	int		m_maxIters;			// max iterations per time step
	bool	m_bremoveIslands;	// remove disconnected elements
	int		m_maxelem;			// the max nr of elements to erode per adaptation iteration
	int		m_nsort;			// sort option (0 = none, 1 = smallest to largest, 2 = largest to smallest)
	int		m_erodeSurfaces;	// option to erode surfaces

	FEMeshAdaptorCriterion* m_criterion;

	DECLARE_FECORE_CLASS()
};
