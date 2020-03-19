/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include "FEThermoFluidPressureLoad.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEThermoFluidPressureLoad, FENodeConstraintSet)
    ADD_PARAMETER(m_matID , "material");
    ADD_PARAMETER(m_p0    , "pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEThermoFluidPressureLoad::FEThermoFluidPressureLoad(FEModel* pfem) : FENodeConstraintSet(pfem), m_nset(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEThermoFluidPressureLoad::Activate()
{
    // don't forget to call base class
    FENLConstraint::Activate();
}

//-----------------------------------------------------------------------------
bool FEThermoFluidPressureLoad::Init()
{
    FEModel& fem = *GetFEModel();
    DOFS& dofs = fem.GetDOFS();
    FEMaterial* pm = fem.GetMaterial(m_matID-1);
    m_tfluid = pm->ExtractProperty<FEThermoFluid>();
    if (m_tfluid == nullptr) return false;
    FEThermoFluidMaterialPoint tp = m_tfluid->CreateMaterialPointData();
    FEFluidMaterialPoint& fp = *tp.ExtractData<FEFluidMaterialPoint>();
    m_dofT = dofs.GetDOF("T");
    m_dofEF = dofs.GetDOF("ef");
    m_alpha = fem.GetTime().alphaf;
    m_rhs = m_p0;

    // initialize surface
    m_nset.Init();
    
    int N = m_nset.Size();

    // create pressure constraints. For the nonlinear constraint c = p - p0,
    // dc = (∂p/∂T) dT + (∂p/∂J) dJ = 0
    // Constraint on T is (lam + eps*c)(∂p/∂T)
    // Constraint on J is (lam + eps*c)(∂p/∂J)
    for (int i=0; i<N; ++i) {

        FEAugLagLinearConstraint* pLC = new FEAugLagLinearConstraint;
        FEAugLagLinearConstraint::DOF dofT;
        FENode* node = m_nset.Node(i);
        double T = node->get(m_dofT)*m_alpha + node->get_prev(m_dofT)*(1-m_alpha);
        double J = 1 + node->get(m_dofEF)*m_alpha + node->get_prev(m_dofEF)*(1-m_alpha);
        tp.m_T = T;
        fp.m_Jf = J;
        double dpT = m_tfluid->Tangent_Pressure_Temperature(tp);
        dofT.node = node->GetID() - 1;    // zero-based
        dofT.bc = m_dofT;
        dofT.val = dpT;
        pLC->m_dof.push_back(dofT);
        FEAugLagLinearConstraint::DOF dofJ;
        double dpJ = m_tfluid->Tangent_Pressure_Strain(tp);
        dofJ.node = node->GetID() - 1;    // zero-based
        dofJ.bc = m_dofEF;
        dofJ.val = dpJ;
        pLC->m_dof.push_back(dofJ);
        // add the linear constraint to the system
        add(pLC);
    }
    
    return true;
}

//-----------------------------------------------------------------------------
double FEThermoFluidPressureLoad::constraint(FEAugLagLinearConstraint& LC)
{
    FEThermoFluidMaterialPoint tp = m_tfluid->CreateMaterialPointData();
    FEFluidMaterialPoint& fp = *tp.ExtractData<FEFluidMaterialPoint>();
    
    int n = (int)LC.m_dof.size();
    double c = 0;
    list<FEAugLagLinearConstraint::DOF>::iterator it = LC.m_dof.begin();
    FEMesh& mesh = GetFEModel()->GetMesh();
    for (int i=0; i<n; ++i, ++it)
    {
        FENode& node = mesh.Node(it->node);
        double T = node.get(m_dofT)*m_alpha + node.get_prev(m_dofT)*(1-m_alpha);
        double J = 1 + node.get(m_dofEF)*m_alpha + node.get_prev(m_dofEF)*(1-m_alpha);
        tp.m_T = T;
        fp.m_Jf = J;
        // calculate pressure
        c = m_tfluid->Pressure(tp);
        // update tangents
        if (it->bc == m_dofT) {
            double dpT = m_tfluid->Tangent_Pressure_Temperature(tp);
            it->val = dpT;
        }
        else {
            double dpJ = m_tfluid->Tangent_Pressure_Strain(tp);
            it->val = dpJ;
        }
    }

    return c;
}

