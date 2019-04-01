/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FETractionLoad.h"
#include "FECore/FEModel.h"

//=============================================================================
BEGIN_FECORE_CLASS(FETractionLoad, FESurfaceTraction)
	ADD_PARAMETER(m_scale   , "scale");
	ADD_PARAMETER(m_traction, "traction");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FETractionLoad::FETractionLoad(FEModel* pfem) : FESurfaceTraction(pfem)
{
	m_scale = 1.0;
	m_traction = vec3d(0, 0, 0);

	// Since the traction is deformation independent, we need to set the linear flag
	SetLinear(true);
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETractionLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_traction.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
vec3d FETractionLoad::Traction(const FESurfaceMaterialPoint& pt)
{
	vec3d t = m_traction(pt)*m_scale;
	return t;
}
