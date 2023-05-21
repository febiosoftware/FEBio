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

#include "stdafx.h"
#include "FEInitialFluidVelocity.h"
#include "FEBioFluid.h"
#include <FECore/FEMaterialPoint.h>
#include <FECore/FENode.h>

BEGIN_FECORE_CLASS(FEInitialFluidVelocity, FENodalIC)
    ADD_PARAMETER(m_v0, "value")->setUnits(UNIT_VELOCITY);
END_FECORE_CLASS();

FEInitialFluidVelocity::FEInitialFluidVelocity(FEModel* fem) : FENodalIC(fem)
{
    m_v0 = vec3d(0, 0, 0);
}

// set the initial value
void FEInitialFluidVelocity::SetValue(const vec3d& v0)
{
    m_v0 = v0;
}

// initialization
bool FEInitialFluidVelocity::Init()
{
    FEDofList dofs(GetFEModel());
    if (dofs.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY)) == false) return false;
    SetDOFList(dofs);
    return true;
}

// return the values for node i
void FEInitialFluidVelocity::GetNodalValues(int inode, std::vector<double>& values)
{
    assert(values.size() == 3);
    
    const FENodeSet& nset = *GetNodeSet();
    const FENode& node = *nset.Node(inode);
    
    FEMaterialPoint mp;
    mp.m_r0 = node.m_r0;
    mp.m_index = inode;
    
    vec3d v0 = m_v0(mp);
    
    values[0] = v0.x;
    values[1] = v0.y;
    values[2] = v0.z;
}
