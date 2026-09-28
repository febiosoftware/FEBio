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
#include <FEBioMech/FETiedElasticInterface.h>
#include <FECore/FESolidElement.h>
#include "FEFluidFSI.h"
#include <FEBioMech/FESolidMaterial.h>

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FETiedFluidFSISolidTraction :    public FEContactInterface
{
public:
    //! constructor
    FETiedFluidFSISolidTraction(FEModel* pfem);
    
    //! destructor
    ~FETiedFluidFSISolidTraction();
    
    //! initialization
    bool Init() override;
    
    //! interface activation
    void Activate() override;
    
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;
    
    //! return the primary and secondary surfaces
	FESurface* GetPrimarySurface() override { return &m_s1; }
	FESurface* GetSecondarySurface() override { return &m_s2; }
    
    //! return integration rule class
    bool UseNodalIntegration() override { return false; }
    
    //! build the matrix profile for use in the stiffness matrix
    void BuildMatrixProfile(FEGlobalMatrix& K) override;
    
public:
    //! calculate contact forces
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! calculate contact stiffness
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
    
    //! calculate Lagrangian augmentations
    bool Augment(int naug, const FETimeInfo& tp) override;
    
    //! update
    void Update() override;
    
protected:
    void InitialProjection(FETiedElasticSurface& s1, FETiedElasticSurface& s2);
    void ProjectSurface(FETiedElasticSurface& s1, FETiedElasticSurface& s2);
    
    //! return the fluid material shared by all elements attached to this surface
    //! (returns nullptr if the surface is not backed by a single fluid material)
    FEFluidFSI* GetFluidFSIMaterial(FETiedElasticSurface& s);
    FESolidMaterial* GetSolidMaterial(FETiedElasticSurface& s);
    
    void CalcAutoPenalty(FETiedElasticSurface& s);
    void SetAutoPenalty(FETiedElasticSurface& s1, FETiedElasticSurface& s2);

public:
    FETiedElasticSurface    m_s1;    //!< primary surface
    FETiedElasticSurface    m_s2;    //!< secondary surface
    
    double          m_atol;         //!< augmentation tolerance
    double          m_gtol;         //!< velocity gap tolerance
    double          m_stol;         //!< search tolerance
    double          m_srad;         //!< contact search radius
    int             m_naugmax;      //!< maximum nr of augmentations
    int             m_naugmin;      //!< minimum nr of augmentations
    
    double          m_epsn;         //!< traction penalty factor
    bool            m_bautopen;     //!< use autopenalty factor
    
    FEFluidFSI* m_pfluid = nullptr; //!< fluid-FSI pointer (set in Init)
    FESolidMaterial* m_psolid = nullptr; //!< solid pointer (set in Init)
    
    int             m_solid;        //!< integer pointer to solid surface (primary = 1, secondaru = 2);

	FEDofList		m_dofU;
   
    DECLARE_FECORE_CLASS();
};
