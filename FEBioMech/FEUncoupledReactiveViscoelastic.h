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
#include "FEUncoupledMaterial.h"
#include "FEBondRelaxation.h"
#include "FEReactiveVEMaterialPoint.h"
#include "FEDamageMaterialUC.h"
#include "FEUncoupledReactiveFatigue.h"
#include <FECore/FEFunction1D.h>

//-----------------------------------------------------------------------------
//! This class implements a large deformation reactive viscoelastic material
//! with uncoupled strain energy density formulation
//
class FEUncoupledReactiveViscoelasticMaterial :	public FEUncoupledMaterial
{
public:
    //! default constructor
    FEUncoupledReactiveViscoelasticMaterial(FEModel* pfem);
        
    //! get the elastic base material
	FEUncoupledMaterial* GetBaseMaterial() { return m_pBase; }
    
    //! Set the base material
    void SetBaseMaterial(FEUncoupledMaterial* pbase) { m_pBase = pbase; }
    
    //! get the elastic bond material
	FEUncoupledMaterial* GetBondMaterial() { return m_pBond; }
    
    //! Set the base material
    void SetBondMaterial(FEUncoupledMaterial* pbond) { m_pBond = pbond; }
    
public:
    //! data initialization
    bool Init() override;

    //! stress function
    mat3ds DevStress(FEMaterialPoint& pt) override;
    mat3ds DevStressStrongBonds(FEMaterialPoint& pt);
    mat3ds DevStressWeakBonds(FEMaterialPoint& pt);

    //! tangent function
    tens4ds DevTangent(FEMaterialPoint& pt) override;
    tens4ds DevTangentStrongBonds(FEMaterialPoint& pt);
    tens4ds DevTangentWeakBonds(FEMaterialPoint& pt);

    //! strain energy density function
    double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    double StrongBondDevSED(FEMaterialPoint& pt) override;
    double WeakBondDevSED(FEMaterialPoint& pt) override;

    //! cull generations
    void CullGenerations(FEMaterialPoint& pt);
    
    //! evaluate bond mass fraction for a given generation
    double BreakingBondMassFraction(FEMaterialPoint& pt, const int ig, const mat3ds D);
    
    //! evaluate bond mass fraction of reforming generation
    double ReformingBondMassFraction(FEMaterialPoint& pt);
    
    //! detect new generation
    bool NewGeneration(FEMaterialPoint& pt);
    
    //! return number of generations
    int RVEGenerations(FEMaterialPoint& pt);
    
    //! returns a pointer to a new material point object
    FEMaterialPointData* CreateMaterialPointData() override;

    //! specialized material points
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;

    //! get base material point
    FEMaterialPoint* GetBaseMaterialPoint(FEMaterialPoint& mp);
    
    //! get bond material point
    FEMaterialPoint* GetBondMaterialPoint(FEMaterialPoint& mp);
    
    //! evaluate scalar strain measure (same type as trigger strain for bond breaking)
    double ScalarStrain(FEMaterialPoint& mp);
    
private:
    FEUncoupledMaterial*	m_pBase;	//!< pointer to elastic solid material for strong bonds
	FEUncoupledMaterial*	m_pBond;	//!< pointer to elastic solid material for reactive bonds
	FEBondRelaxation*		m_pRelx;    //!< pointer to bond relaxation material for reactive bonds
    FEDamageCDF*            m_pWCDF;    //!< pointer to weak bond recruitment CDF
    
private:
    FEDamageMaterialUC*                 m_pDmg; //!< pointer to base material if it is a FEDamageMaterialUC
    FEUncoupledReactiveFatigue*         m_pFtg; //!< pointer to base material if it is a FEUncoupledReactiveFatigue
    double Damage(FEMaterialPoint& mp);         //!< return damage in the base material
    
public:
    double	m_wmin;		//!< minimum value of relaxation
    int     m_btype;    //!< bond kinetics type
    int     m_ttype;    //!< bond breaking trigger type
    double  m_emin;     //!< strain threshold for triggering new generation

    int     m_nmax;     //!< highest number of generations achieved in analysis
    
    DECLARE_FECORE_CLASS();
};
