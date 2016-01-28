//
//  FEReactiveVEMaterialPoint.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 11/30/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEReactiveVEMaterialPoint__
#define __FEBioMech__FEReactiveVEMaterialPoint__

#include "FECore/FEMaterialPoint.h"
#include "FEReactiveViscoelastic.h"
#include "FEUncoupledReactiveViscoelastic.h"
#include <deque>

class FEReactiveViscoelasticMaterial;
class FEUncoupledReactiveViscoelasticMaterial;

//-----------------------------------------------------------------------------
//! Material point data for reactive viscoelastic materials
class FEReactiveVEMaterialPoint : public FEMaterialPoint
{
public:
    //! olverloaded constructors
    FEReactiveVEMaterialPoint(FEMaterialPoint *pt, FEReactiveViscoelasticMaterial *pe) : FEMaterialPoint(pt) { m_pRve = pe; m_pRuc = 0; }
    FEReactiveVEMaterialPoint(FEMaterialPoint *pt, FEUncoupledReactiveViscoelasticMaterial *pe) : FEMaterialPoint(pt) { m_pRve = 0; m_pRuc = pe; }
    
    //! copy material point data
    FEMaterialPoint* Copy();
    
    //! Initialize material point data
    void Init(bool bflag);
    
    //! Serialize data to archive
    void Serialize(DumpStream& ar);
    
public:
    // multigenerational material data
    deque <mat3d>  m_Fi;	//!< inverse of relative deformation gradient
    deque <double> m_Ji;	//!< determinant of Fi (store for efficiency)
    deque <double> m_v;     //!< time when generation starts breaking
    deque <double> m_w;     //!< mass fraction when generation starts breaking
    FEReactiveViscoelasticMaterial*  m_pRve; //!< pointer to parent material
    FEUncoupledReactiveViscoelasticMaterial*  m_pRuc; //!< pointer to parent material
};


#endif /* defined(__FEBioMech__FEReactiveVEMaterialPoint__) */
