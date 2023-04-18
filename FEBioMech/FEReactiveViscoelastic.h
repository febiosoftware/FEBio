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
#include "FEElasticMaterial.h"
#include "FEBondRelaxation.h"
#include "FEReactiveVEMaterialPoint.h"
#include "FEDamageMaterial.h"
#include "FEReactiveFatigue.h"
#include <FECore/FEFunction1D.h>

//-----------------------------------------------------------------------------
//! This class implements a large deformation reactive viscoelastic material
//
class FEReactiveViscoelasticMaterial :	public FEElasticMaterial
{
public:
	//! default constructor
	FEReactiveViscoelasticMaterial(FEModel* pfem);
       
	//! get the elastic base material
	FEElasticMaterial* GetBaseMaterial() { return m_pBase; }
    
	//! Set the base material
	void SetBaseMaterial(FEElasticMaterial* pbase) { m_pBase = pbase; }
    
	//! get the elastic bond material
	FEElasticMaterial* GetBondMaterial() { return m_pBond; }
    
	//! Set the base material
	void SetBondMaterial(FEElasticMaterial* pbond) { m_pBond = pbond; }
    
public:
    //! data initialization
    bool Init() override;
    
	//! stress function
	mat3ds Stress(FEMaterialPoint& pt) override;
    mat3ds StressStrongBonds(FEMaterialPoint& pt);
    mat3ds StressWeakBonds(FEMaterialPoint& pt);

	//! tangent function
	tens4ds Tangent(FEMaterialPoint& pt) override;
    tens4ds TangentStrongBonds(FEMaterialPoint& pt);
    tens4ds TangentWeakBonds(FEMaterialPoint& pt);

    //! strain energy density function
    double StrainEnergyDensity(FEMaterialPoint& pt) override;
    double StrongBondSED(FEMaterialPoint& pt) override;
    double WeakBondSED(FEMaterialPoint& pt) override;

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
	FEElasticMaterial*	m_pBase;	//!< pointer to elastic solid material for strong bonds
	FEElasticMaterial*	m_pBond;	//!< pointer to elastic solid material for reactive bonds
	FEBondRelaxation*   m_pRelx;    //!< pointer to bond relaxation material for reactive bonds
    FEDamageCDF*        m_pWCDF;    //!< pointer to weak bond recruitment CDF

private:
    FEDamageMaterial*   m_pDmg;     //!< pointer to base material if it is a FEDamageMaterial
    FEReactiveFatigue*  m_pFtg;     //!< pointer to base material if it is a FEReactiveFatigue
    double Damage(FEMaterialPoint& mp); //!< return damage in this material

public:
    double	m_wmin;		//!< minimum value of relaxation
    int     m_btype;    //!< bond kinetics type
    int     m_ttype;    //!< bond breaking trigger type
    double  m_emin;     //!< strain threshold for triggering new generation
    
    int     m_nmax;     //!< highest number of generations achieved in analysis
    
    DECLARE_FECORE_CLASS();
};
