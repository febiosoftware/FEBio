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
class FESIVViscoelasticMaterialPoint : public FEMaterialPointData
{
public:
    enum { MAX_TERMS = 6 };

public:
    //! constructor
    FESIVViscoelasticMaterialPoint(FEMaterialPointData* mp = nullptr);

    //! copy material point data
    FEMaterialPointData* Copy();

    //! Initialize material point data
    void Init();

    //! Update material point data
    void Update(const FETimeInfo& timeInfo);

    //! Serialize data to archive
    void Serialize(DumpStream& ar);

public:
    mat3ds  m_H;                    //!< measure of Hencky strain
    mat3ds  m_Hp;                   //!< measure of Hencky strain at previous time step
    mat3ds  m_dHp;                  //!< measure of time derivative of Hencky strain at previous time step
    mat3d   m_R;                    //!< rotation matrix

    double  m_alpha[MAX_TERMS];     //!< exponent of right-stretch tensor in series spring
    double  m_alphap[MAX_TERMS];    //!< alpha at previoust time step
    double  m_dalphap[MAX_TERMS];   //!< alpha dot at  previoust time step
    double  m_mumr[MAX_TERMS];      //!< shear modulus calculated from reference state
    mat3ds  m_Hs[MAX_TERMS];        //!< Hs = alpha*H at current time step
    mat3ds  m_Hsp[MAX_TERMS];       //!< Hs at previous timestep
    mat3ds  m_dHsp[MAX_TERMS];      //!< Hs dot at previous timestep
    double  m_sed;                  //!< elastic strain energy density
    double  m_sedp;                 //!< sed at previous time step
};


//-----------------------------------------------------------------------------
//! This class implements a large deformation visco-elastic material
//
class FESIVViscoelastic :    public FEElasticMaterial
{
public:
    // NOTE: make sure that this parameter is the
    //       same as the MAX_TERMS in the FESIVViscoelasticMaterialPoint class
    enum { MAX_TERMS = FESIVViscoelasticMaterialPoint::MAX_TERMS };

public:
    //! default constructor
    FESIVViscoelastic(FEModel* pfem);

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
    
    //! convert right Hencky strain H to right stretch tensor U
    mat3ds H2U(mat3ds& H);
    
    //! solve a system of equations based on a tensorial nonlinear equation
    bool SolvedY(tens4ds& dTdY, mat3ds& Y, mat3ds& dY);
    
public:
    // material parameters
    double  m_g0;               //!< intitial visco-elastic coefficient
    double  m_g[MAX_TERMS];     //!< visco-elastic coefficients
    double  m_t[MAX_TERMS];     //!< relaxation times in reference configuration
    int     m_ttype;            //!< strain trigger type

private:
    FEElasticMaterial*    m_Base;    //!< pointer to parallel elastic solid material
    FEElasticMaterial*    m_Mxwl;    //!< pointer to Maxwell elastic solid material

public:
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
