//
//  FEFSIErosionVolumeRatio.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/15/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEFSIErosionVolumeRatio.h"
#include "FEFluidFSIDomain3D.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/log.h>
#include "FEBioMech/FEElasticMaterial.h"
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
        return true;
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
        return true;
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
    // if any nodes were orphaned, we need to deactivate them
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
            vector<int>& id = node.m_ID;
            for (int j = 0; j < id.size(); ++j)
            {
                int n = id[j];
                if (n >= 0) id[j] = -n - 2;
            }
        }
    }
    
//    feLog("Deactivate elements: %d\n", elems);
//    return (elems == 0);
    feLog("Deactivate elements: %d\n", melems);
    feLog("Reactivate elements: %d\n", activeElems);
    return (melems == 0);
}
