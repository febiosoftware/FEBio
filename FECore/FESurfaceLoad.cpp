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
#include "FESurfaceLoad.h"
#include "FEMesh.h"
#include "DumpStream.h"
#include "FEModel.h"

FESurfaceLoad::FESurfaceLoad(FEModel* pfem) : FEModelLoad(pfem)
{
	m_psurf = 0;
}

FESurfaceLoad::~FESurfaceLoad(void)
{

}

//! Set the surface to apply the load to
void FESurfaceLoad::SetSurface(FESurface* ps)
{
	m_psurf = ps; 
}

bool FESurfaceLoad::Init()
{
	if (m_psurf == 0) return false;
	if (m_psurf->Init() == false) return false;
	return FEModelLoad::Init();
}

void FESurfaceLoad::Serialize(DumpStream& ar)
{
	FEModelLoad::Serialize(ar);
	if (ar.IsShallow()) return;

	ar & m_psurf;

	// the mesh manages surfaces for surface loads
	if (m_psurf && ar.IsLoading())
	{
		FEMesh* pm = m_psurf->GetMesh();
		pm->AddSurface(m_psurf);
	}
}

void FESurfaceLoad::ForceMeshUpdate()
{
	GetFEModel()->SetMeshUpdateFlag(true);
}
