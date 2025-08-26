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

class FETiedFluidInterface;

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FETiedFluidSurface : public FESurface
{
public:
    //! constructor
    FETiedFluidSurface(FEModel* pfem);
    
    //! destructor
    ~FETiedFluidSurface() {}
    
    //! initialization
    bool Init() override;
   
    //! Get the surface
    FESurface* GetSurface() { return dynamic_cast<FESurface*>(this); }
    
public:
    std::vector<FESurfaceElement*>  m_pme;
    std::vector< std::vector<double>> m_rs;
    std::vector<double> m_gap;

    std::vector<int>m_tag;
};

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FETiedFluidInterface : public FEContactInterface
{
public:
    //! constructor
    FETiedFluidInterface(FEModel* pfem);
    
    //! destructor
    ~FETiedFluidInterface();
    
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
    void InitialNodalProjection(FETiedFluidSurface& ss, FETiedFluidSurface& ms);
    bool    m_binit;

public:
	FETiedFluidSurface    m_ss;    //!< primary surface
	FETiedFluidSurface    m_ms;    //!< secondary surface
    FELinearConstraintSet   m_lcv;
    FELinearConstraintSet   m_lcp;

    bool            m_btwopass;     //!< flag for two pass analysis
    double          m_stol;         //!< search tolerance
    double          m_srad;         //!< contact search radius
    bool            m_bfreedofs;    //!< flag to free DOFs inside tied interface

    double          m_epsv;          //!< penalty factor for velocity
    double          m_epsp;          //!< penalty factor for pressure

	FEDofList		m_dofWE;
   
    DECLARE_FECORE_CLASS();
};
