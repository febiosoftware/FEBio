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
#include <FEBioMech/FEContactInterface.h>
#include "FEFluidMaterial.h"

class FETiedSolidInterface;

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FETiedFluidFSISurface : public FESurfaceConstraint
{
public:
    //! constructor
    FETiedFluidFSISurface(FEModel* pfem);
    
    //! destructor
    ~FETiedFluidFSISurface() {}
    
    //! Activation
    void Activate() override;
    
    // allocate equations
    int InitEquations(int neq) override;

    //! initialization
    bool Init() override;
   
    //! Get the surface
    FESurface* GetSurface() override { return &m_surf; }
    
public:
	FEDofList	m_dofU;
    
public:
    std::vector<FESurfaceElement*>  m_pme;
    std::vector< std::vector<double>> m_rs;
    
    std::vector<int> m_tag;
    
    //! Set the sibling of this contact surface
    void SetSibling(FETiedFluidFSISurface* ps) { m_pSibling = ps; }
    FETiedFluidFSISurface* GetSibling() { return m_pSibling; }

    //! Set the parent of this contact surface
    void SetContactInterface(FETiedSolidInterface* ps) { m_pTiedSolidInterface = ps; }
    
    //! Get the parent of this contact surface
    FETiedSolidInterface* GetContactInterface() { return m_pTiedSolidInterface; }
    
public:
    void Update(const std::vector<double>& Ui, const std::vector<double>& ui) override;
    void UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) override;
    
    void PrepStep() override;
    
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

protected:
    FESurface               m_surf;
    FELinearConstraintSet   m_lcu;
    FETiedSolidInterface*   m_pTiedSolidInterface;
    FETiedFluidFSISurface*  m_pSibling;
    bool                    m_binit;

};

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FETiedSolidInterface : public FEContactInterface
{
public:
    //! constructor
    FETiedSolidInterface(FEModel* pfem);
    
    //! destructor
    ~FETiedSolidInterface();
    
    //! initialization
    bool Init() override;
    
    //! interface activation
    void Activate() override;
    
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;
    
    //! return the primary and secondary surfaces
    FESurface* GetPrimarySurface() override { return m_ss.GetSurface(); }
	FESurface* GetSecondarySurface() override { return m_ms.GetSurface(); }
    
    //! return integration rule class
    bool UseNodalIntegration() override { return false; }
    
    //! build the matrix profile for use in the stiffness matrix
    void BuildMatrixProfile(FEGlobalMatrix& K) override;

    int InitEquations(int neq) override;

public:
    //! calculate contact forces
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! calculate contact stiffness
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
    
    //! calculate Lagrangian augmentations
    bool Augment(int naug, const FETimeInfo& tp) override;
    
    // called at start of time step
    void PrepStep() override;
    
    //! update
    void Update(const std::vector<double>& Ui, const std::vector<double>& ui) override;
    void UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) override;

protected:
    void InitialNodalProjection(FETiedFluidFSISurface& ss, FETiedFluidFSISurface& ms);

public:
	FETiedFluidFSISurface    m_ss;    //!< primary surface
	FETiedFluidFSISurface    m_ms;    //!< secondary surface
    
    double          m_tol;          //!< augmentation tolerance
    bool            m_btwopass;     //!< flag for two pass analysis
    double          m_stol;         //!< search tolerance
    double          m_srad;         //!< contact search radius
    int             m_naugmax;      //!< maximum nr of augmentations
    int             m_naugmin;      //!< minimum nr of augmentations
    
    double          m_epsu;          //!< penalty factor for velocity

    FEFluidMaterial* m_pfluid;       //!< fluid pointer

	FEDofList		m_dofU;
   
    DECLARE_FECORE_CLASS();
};
