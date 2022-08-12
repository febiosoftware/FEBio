/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.

 See Copyright-FEBio.txt for details.

 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEPrescribedFluidPressure.h"
#include "FEBioFluid.h"
#include <FECore/FESurface.h>
#include <FECore/FEModel.h>
#include "FEFluid.h"

 //-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEPrescribedFluidPressure, FEPrescribedSurface)
	ADD_PARAMETER(m_p0, "pressure");
END_FECORE_CLASS();


 //! constructor
FEPrescribedFluidPressure::FEPrescribedFluidPressure(FEModel* fem) : FEPrescribedSurface(fem)
{
	m_p0 = 0.0;
}

//! initialize
bool FEPrescribedFluidPressure::Init()
{
    // get the dof index
    m_dofEF = GetDOFIndex(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION), 0);
    if (m_dofEF < 0) return false;
    SetDOFList(m_dofEF);

    if (FEPrescribedSurface::Init() == false) return false;

    // get fluid from first surface element
    // assuming the entire surface bounds the same fluid
    FESurface& surf = *GetSurface();
    FESurfaceElement& el = surf.Element(0);
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;

    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    m_pfluid = pm->ExtractProperty<FEFluidMaterial>();
    if (m_pfluid == nullptr) return false;

    return true;
}

void FEPrescribedFluidPressure::GetNodalValues(int nodelid, std::vector<double>& val)
{
    FESurface& surf = *GetSurface();
    FENode& node = surf.Node(nodelid);

    // calculate the dilatation
    double e = node.get(m_dofEF);
    bool good = m_pfluid->Dilatation(0, m_p0, 0, e);
    assert(good);

    val[0] = e;
}

void FEPrescribedFluidPressure::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement
    assert(false);
}
