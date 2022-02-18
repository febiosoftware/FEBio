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
#include "FEFSIErosionVolumeRatio.h"
#include "FEFluidFSIDomain3D.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/log.h>
#include "FEBioMech/FEElasticMaterial.h"
#include <FECore/FELinearConstraintManager.h>
#include <algorithm>

BEGIN_FECORE_CLASS(FEFSIErosionVolumeRatio, FEMeshAdaptor)
ADD_PARAMETER(m_minJ, FE_RANGE_GREATER(0.0), "min_volume_ratio");
ADD_PARAMETER(m_maxElems, "max_elems");
ADD_PARAMETER(m_maxIters, "max_iters");
ADD_PARAMETER(m_metric  , "metric");
END_FECORE_CLASS();

FEFSIErosionVolumeRatio::FEFSIErosionVolumeRatio(FEModel* fem) : FEMeshAdaptor(fem)
{
    m_minJ = 0.0;
    m_maxElems = 0;
    m_maxIters = -1;
    m_metric = 0;
}

bool FEFSIErosionVolumeRatio::Apply(int iteration)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    if ((m_maxIters >= 0) && (iteration >= m_maxIters))
    {
        feLog("Max iterations reached.");
        return false;
    }
    
    int deactiveElems = 0;
    int activeElems = 0;
    
    map<int, bool> melem;
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        FEFluidFSIDomain* fsidom = dynamic_cast<FEFluidFSIDomain*>(&dom);
        if (fsidom) {
            int NE = dom.Elements();
            for (int j = 0; j < NE; ++j)
            {
                FEElement& el = dom.ElementRef(j);
                if (el.isActive())
                {
                    bool bdeactivate = false;
                    int nint = el.GaussPoints();
                    double elemVal = 0;
                    for (int n = 0; n < nint; ++n)
                    {
                        FEMaterialPoint* mp = el.GetMaterialPoint(n);
                        FEElasticMaterialPoint* ep = mp->ExtractData<FEElasticMaterialPoint>();
                        if (ep)
                        {
                            mat3ds C = ep->RightCauchyGreen();
                            
                            switch (m_metric)
                            {
                                case 0: elemVal = ep->m_J; break;
                                case 1:
                                {
                                    double lam[3];
                                    C.eigen(lam);
                                    elemVal = sqrt(lam[2]);
                                }
                                    break;
                                default:
                                    break;
                            }

                            if (elemVal <= m_minJ)
                            {
                                bdeactivate = true;
                                break;
                            }
                        }
                    }
                    
                    if (bdeactivate)
                    {
                        melem[el.GetID()] = true;
                        deactiveElems++;
                    }
                }
                // if an element is inactive, check whether we need to activate it again
                else
                {
                    // nodal coordinates
                    const int NELN = FEElement::MAX_NODES;
                    vec3d r[NELN];
                    FESolidElement* sel = dynamic_cast<FESolidElement*>(&el);
                    if (sel) {
                        FEFluidFSIDomain3D* fsi3D = dynamic_cast<FEFluidFSIDomain3D*>(fsidom);
                        fsi3D->GetCurrentNodalCoordinates(*sel, r);
                        int nint = el.GaussPoints();
                        double minJ = 1;
                        for (int n = 0; n < nint; ++n)
                        {
                            FEMaterialPoint* mp = el.GetMaterialPoint(n);
                            FEElasticMaterialPoint* ep = mp->ExtractData<FEElasticMaterialPoint>();
                            
                            // get the deformation gradient and determinant at intermediate time
                            double Jt;
                            mat3d Ft, Fp;
                            Jt = fsi3D->defgrad(*sel, Ft, n, r);
                            ep->m_J = Jt;
                            
                            ep->m_F = Ft;

                            mat3ds C = ep->RightCauchyGreen();
                            
                            switch (m_metric)
                            {
                                case 0: minJ = min(ep->m_J,minJ); break;
                                case 1:
                                {
                                    double lam[3];
                                    C.eigen(lam);
                                    minJ = min(sqrt(lam[2]),minJ);
                                }
                                    break;
                                default:
                                    break;
                            }
                        }
                        
                        if (minJ > m_minJ) {
                            melem[el.GetID()] = false;
                            activeElems++;
                        }
                    }
                }
            }
        }
    }
    if ((deactiveElems == 0) && (activeElems == 0))
    {
        feLog("Nothing to do.");
        return false;
    }
    
    int melems = melem.size();
    if (melems > m_maxElems) melems = m_maxElems;
    int nel = 0;
    for (std::map<int,bool>::iterator it=melem.begin(); it!=melem.end(); ++it)
    {
        FEElement* pe = mesh.FindElementFromID(it->first); assert(pe);
        if (it->second) {
            pe->setInactive();
            if (++nel > m_maxElems) break;
        }
        else {
            pe->setActive();
        }
    }
    
    // if any nodes were orphaned, we need to deactivate them as well
    int NN = mesh.Nodes();
    vector<int> tag(NN, 0);
    for (int i = 0; i < mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        int NE = dom.Elements();
        for (int j = 0; j < NE; ++j)
        {
            FEElement& el = dom.ElementRef(j);
            if (el.isActive())
            {
                int neln = el.Nodes();
                for (int n = 0; n < neln; ++n) tag[el.m_node[n]] = 1;
            }
        }
    }
    
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        if (tag[i] == 0)
        {
            node.SetFlags(FENode::EXCLUDE);
            int ndofs = node.dofs();
            for (int j = 0; j < ndofs; ++j)
                node.set_inactive(j);
        }
    }
    
    // remove any linear constraints of exclude nodes
    FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
    for (int j = 0; j < LCM.LinearConstraints();)
    {
        FELinearConstraint& lc = LCM.LinearConstraint(j);
        if (mesh.Node(lc.GetParentNode()).HasFlags(FENode::EXCLUDE))
        {
            LCM.RemoveLinearConstraint(j);
        }
        else ++j;
    }
    
    // also remove any linear constraints that have excluded child nodes
    for (int j = 0; j < LCM.LinearConstraints();)
    {
        FELinearConstraint& lc = LCM.LinearConstraint(j);
        
        bool del = false;
        int n = lc.Size();
        for (int k = 0; k < n; ++k)
        {
            if (mesh.Node(lc.GetChildDof(k).node).HasFlags(FENode::EXCLUDE))
            {
                del = true;
                break;
            }
        }
        if (del) LCM.RemoveLinearConstraint(j); else ++j;
    }
    
    // reactivate the linear constraints
    LCM.Activate();

//    feLog("Deactivate elements: %d\n", elems);
//    return (elems == 0);
    feLog("Deactivate elements: %d\n", melems);
    feLog("Reactivate elements: %d\n", activeElems);
    return (melems != 0);
}
