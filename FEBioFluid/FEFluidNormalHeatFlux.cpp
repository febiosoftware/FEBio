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
#include "FEFluidNormalHeatFlux.h"
#include <FECore/FEElemElemList.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FEGlobalVector.h>
#include <FECore/log.h>
#include <FECore/LinearSolver.h>
#include "FEBioThermoFluid.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidNormalHeatFlux, FESurfaceLoad)
    ADD_PARAMETER(m_flux    , "flux");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidNormalHeatFlux::FEFluidNormalHeatFlux(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_flux = 0.0;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidNormalHeatFlux::Init()
{
    // get the degrees of freedom
    m_dof.Clear();
    m_dof.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE));
    if (m_dof.IsEmpty()) return false;

    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the prescribed normal velocity
void FEFluidNormalHeatFlux::LoadVector(FEGlobalVector& R)
{
    FESurface& surf = GetSurface();

    // evaluate the integral
    surf.LoadVector(R, m_dof, false, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {
        
        // evaluate pressure at this material point
        double q = m_flux(pt);

        double J = (pt.dxr ^ pt.dxs).norm();

        double H = dof_a.shape;

        val[0] = H*q*J;
    });
}
