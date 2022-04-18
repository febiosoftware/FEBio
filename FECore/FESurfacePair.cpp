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
#include "FESurfacePair.h"
#include "FEFacetSet.h"
#include "FEMesh.h"
#include "DumpStream.h"

//--------------------------------------------------------
FESurfacePair::FESurfacePair(FEMesh* pm) : m_mesh(pm)
{
	m_surface1 = 0;
	m_surface2 = 0;
}

void FESurfacePair::SetName(const std::string& name)
{
	m_name = name;
}

const std::string& FESurfacePair::GetName() const
{
	return m_name;
}

FEFacetSet* FESurfacePair::GetPrimarySurface()
{
	return m_surface1;
}

void FESurfacePair::SetPrimarySurface(FEFacetSet* pf)
{
	m_surface1 = pf;
}

FEFacetSet* FESurfacePair::GetSecondarySurface()
{
	return m_surface2;
}

void FESurfacePair::SetSecondarySurface(FEFacetSet* pf)
{
	m_surface2 = pf;
}

void FESurfacePair::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	ar& m_surface1;
	ar& m_surface2;
	ar& m_mesh;
}
