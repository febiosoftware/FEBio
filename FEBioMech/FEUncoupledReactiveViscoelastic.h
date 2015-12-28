//
//  FEUncoupledReactiveViscoelastic.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 12/12/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEUncoupledReactiveViscoelastic__
#define __FEBioMech__FEUncoupledReactiveViscoelastic__

#include "FEUncoupledMaterial.h"
#include "FEBondRelaxation.h"
#include "FEReactiveVEMaterialPoint.h"

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
    FEElasticMaterial* GetBaseMaterial() { return m_pBase; }
    
    //! Set the base material
    void SetBaseMaterial(FEUncoupledMaterial* pbase) { m_pBase = pbase; }
    
    //! get the elastic bond material
    FEElasticMaterial* GetBondMaterial() { return m_pBond; }
    
    //! Set the base material
    void SetBondMaterial(FEUncoupledMaterial* pbond) { m_pBond = pbond; }
    
    //! Set the local coordinate system for a material point (overridden from FEMaterial)
    void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);
        
public:
    //! data initialization
    bool Init();
    
    //! stress function
    mat3ds DevStress(FEMaterialPoint& pt);
    
    //! tangent function
    tens4ds DevTangent(FEMaterialPoint& pt);
    
    //! strain energy density function
    double DevStrainEnergyDensity(FEMaterialPoint& pt);
    
    //! cull generations
    void CullGenerations(FEMaterialPoint& pt);
    
    //! evaluate bond mass fraction for a given generation
    double BreakingBondMassFraction(FEMaterialPoint& pt, const int ig);
    
    //! evaluate bond mass fraction of reforming generation
    double ReformingBondMassFraction(FEMaterialPoint& pt);
    
    //! detect new generation
    bool NewGeneration(FEMaterialPoint& pt);
    
    //! returns a pointer to a new material point object
    FEMaterialPoint* CreateMaterialPointData();
    
private:
    FEPropertyT<FEUncoupledMaterial>	m_pBase;	//!< pointer to elastic solid material for polymer base
	FEPropertyT<FEUncoupledMaterial>	m_pBond;	//!< pointer to elastic solid material for reactive bonds
	FEPropertyT<FEBondRelaxation>		m_pRelx;    //!< pointer to bond relaxation material for reactive bonds
    
public:
    double	m_wmin;		//!< minimum value of relaxation
    int     m_btype;    //!< bond kinetics type
    int     m_ttype;    //!< bond breaking trigger type
    
    DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEUncoupledReactiveViscoelastic__) */
