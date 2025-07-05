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
#include "FEInitialDisplacement.h"
#include "FEBioMech.h"
#include <FECore/FEMaterialPoint.h>
#include <FECore/FENode.h>

BEGIN_FECORE_CLASS(FEInitialDisplacement, FENodalIC)
    ADD_PARAMETER(m_u0, "value")->setUnits(UNIT_LENGTH);
END_FECORE_CLASS();

FEInitialDisplacement::FEInitialDisplacement(FEModel* fem) : FENodalIC(fem)
{
    m_u0 = vec3d(0, 0, 0);
}

// set the initial value
void FEInitialDisplacement::SetValue(const vec3d& u0)
{
    m_u0 = u0;
}

// initialization
bool FEInitialDisplacement::Init()
{
    FEDofList dofs(GetFEModel());
    if (dofs.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT)) == false) return false;
    SetDOFList(dofs);
    return true;
}

// return the values for node i
void FEInitialDisplacement::GetNodalValues(int inode, std::vector<double>& values)
{
    assert(values.size() == 3);
    
    const FENodeSet& nset = *GetNodeSet();
    const FENode& node = *nset.Node(inode);
    
    FEMaterialPoint mp;
    mp.m_r0 = node.m_r0;
    mp.m_index = inode;
    
    vec3d u0 = m_u0(mp);
    
    values[0] = u0.x;
    values[1] = u0.y;
    values[2] = u0.z;
}
