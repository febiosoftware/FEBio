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
#include "FEActiveContractionMaterial.h"

//-----------------------------------------------------------------------------
class FEForceVelocityMaterialPoint : public FEMaterialPointData
{
public:
    enum { MAX_TERMS = 3 };
    
public:
    
	FEForceVelocityMaterialPoint();
    
	FEMaterialPointData* Copy();
    
    void Init();
    
    void Serialize(DumpStream& ar);
    
public:
    
    double m_lambdap;
    
    double m_H[MAX_TERMS];
    double m_Hp[MAX_TERMS];
    
};

//-----------------------------------------------------------------------------
//! A material class describing the active fiber contraction model of Estrada et al. in doi: 10.1115/1.4044030
class FEForceVelocityContraction : public FEActiveContractionMaterial
{
public:
    enum { MAX_TERMS = 3 };
    
public:
    FEForceVelocityContraction(FEModel* pfem);
    
    //! initialization
    bool Init() override;
    
    //! calculate the active stress
    mat3ds ActiveStress(FEMaterialPoint& mp, const vec3d& a0) override;
    
    //! active contraction stiffness contribution
    tens4ds ActiveStiffness(FEMaterialPoint& mp, const vec3d& a0) override;
    
    //! create material point data
    FEMaterialPointData* CreateMaterialPointData() override { return new FEForceVelocityMaterialPoint; }
    
public:
    //! update force-velocity material point
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp, const vec3d& a0) override;
    
protected:
    double  m_ascl;             //!< activation scale factor
    double  m_Tmax;             //!< activation scale factor
    double  m_ca0;              //!< intracellular calcium concentration
    double  m_camax;            //!< peak calcium concentration
    double  m_beta;             //!< shape of peak isometric tension-sarcomere length relation
    double  m_l0;               //!< unloaded length
    double  m_refl;             //!< sarcomere length
    double  m_alpha[MAX_TERMS];
    double  m_A[MAX_TERMS];
    double  m_at;
    bool    m_bfvel;            //!< flag for calculating force-velocity contraction response

    DECLARE_FECORE_CLASS();
};
