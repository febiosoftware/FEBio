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



#pragma once
#include <FECore/FEAugLagLinearConstraint.h>
#include <FECore/FESurface.h>
#include "febiofluid_api.h"
#include <FECore/FEMesh.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FESurface.h>
#include <FECore/MeshTools.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FEPrescribedBC.h>
#include <map>

//-----------------------------------------------------------------------------
//! The FEConstraintImmersedBody class implements a constraint for the continuity of fluid velocity
//! with an immersedbody

class FEBIOFLUID_API FEConstraintImmersedBody : public FESurfaceConstraint
{
public:
    //! constructor
    FEConstraintImmersedBody(FEModel* pfem);
    
    //! destructor
    ~FEConstraintImmersedBody() {}
    
    //! Activation
    void Activate() override;

    //! initialization
    bool Init() override;
    
    //! Get the surface
    FESurface* GetSurface() override { return &m_surf; }
    
public:
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;

    //! add the linear constraint contributions to the residual
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

    //! add the linear constraint contributions to the stiffness matrix
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

    //! do the augmentation
    bool Augment(int naug, const FETimeInfo& tp) override;

    //! build connectivity for matrix profile
    void BuildMatrixProfile(FEGlobalMatrix& M) override;
    
    //! prep step
    void PrepStep() override;
    
    //! create intersection surface
    FESurface* IntersectionSurface();
    
    void ElementForce(FEGlobalVector& R, const FETimeInfo& tp, const int n, int* nn);
    
private:
    void GetIntersectedEdges();
    
protected:
    FESurface               m_surf;
    FELinearConstraintSet   m_lc;
    vector<int>             m_nodetag;
    vector<int>             m_edgetag;
    vector< vector<int> >   m_nodal_elems;  // list of elements associated with each node of the fluid domain
    vector<vec3d>           m_vs;           // solid velocity at intersection of immersed body with fluid domain edge
    vector<double>          m_vsn;          // normal solid velocity at same intersection
    FEEdgeList              m_EL;           // list of intersected edges
    vector< vector<int>>    m_elem_edges;   // list of intersected edges in each element of the fluid domain
    bool                    m_breset;
    FEPrescribedNodeSet*    m_pbcwx;
    FEPrescribedNodeSet*    m_pbcwy;
    FEPrescribedNodeSet*    m_pbcwz;
    FEPrescribedNodeSet*    m_pbcef;

    FEDofList   m_dofW;
    FEDofList   m_dofEF;
    FEDofList   m_dofV;
    
    DECLARE_FECORE_CLASS();
};
