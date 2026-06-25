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

//-----------------------------------------------------------------------------
//! Material point data for visco-elastic materials
class FESSVQLVMaterialPoint : public FEMaterialPointData
{
public:
    enum { MAX_TERMS = 6 };
    
public:
    //! constructor
    FESSVQLVMaterialPoint(FEMaterialPointData* mp = nullptr);
    
    //! copy material point data
    FEMaterialPointData* Copy();
    
    //! Update material point data
    void Update(const FETimeInfo& timeInfo);
    
    //! Initialize material point data
    void Init();
    
    //! Serialize data to archive
    void Serialize(DumpStream& ar);

    //! Evaluate right stretch tensor from right Cauchy-Green tensor
    mat3ds CtoU(mat3ds C);

public:
    mat3ds  m_E;                    //!< parallel spring Green-Lagrange strain tensor at current time
    mat3ds  m_Ep;                   //!< parallel spring Green-Lagrange strain tensor at previous time
    mat3ds  m_Es;                   //!< Maxwell spring Green-Lagrange strain tensor at current time
    mat3ds  m_Esp;                  //!< Maxwell spring Green-Lagrange strain tensor at previous time
    double  m_sed;                  //!< elastic strain energy density
    double  m_sedp;                 //!< sed at previous time step
    double  m_rd;                   //!< residual dissipation at current time
    mat3ds  m_Sm;                   //!< PK2 stess in the Maxwell element
    mat3ds  m_Edot;                 //!< time derivative of m_E at current time
    mat3ds  m_Edotp;                //!< time derivative of m_E at previous time
    mat3ds  m_Esdot;                //!< time derivative of m_Es at current time
    mat3ds  m_Esdotp;               //!< //!< time derivative of m_Es at previous time
};


//-----------------------------------------------------------------------------
//! This class implements a large deformation visco-elastic material
//
class FESSVQLV :    public FEElasticMaterial
{
public:
    //! default constructor
    FESSVQLV(FEModel* pfem);

    //! get the elastic base material
    FEElasticMaterial* GetBaseMaterial() { return m_Base; }

    //! get the elastic Maxwell material
    FEElasticMaterial* GetMxwlMaterial() { return m_Mxwl; }

    //! Set the base material
    void SetBaseMaterial(FEElasticMaterial* pbase) { m_Base = pbase; }

    //! Set the Maxwell material
    void SetMxwlMaterial(FEElasticMaterial* pmxwl) { m_Mxwl = pmxwl; }

public:
    //! stress function
    mat3ds Stress(FEMaterialPoint& pt) override;

    //! tangent function
    tens4ds Tangent(FEMaterialPoint& pt) override;

    //! strain energy density
    double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    // returns a pointer to a new material point object
    FEMaterialPointData* CreateMaterialPointData() override;

    //! update specialize material point data
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;
    
    bool UseSecantTangent() override { return m_secant_tangent; }
    
public:
    // material parameters
    FEParamDouble   m_eta;          //!< dashpot shear vicosity
    FEParamDouble   m_kappa;        //!< dashpot bulk vicosity
    int             m_itmin;        //!< minimum number of iterations
    int             m_itmax;        //!< maximum number of iterations

private:
    FEElasticMaterial*    m_Base;   //!< pointer to parallel elastic solid material
    FEElasticMaterial*    m_Mxwl;   //!< pointer to Maxwell elastic solid material
    bool    m_secant_tangent;

public:
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
