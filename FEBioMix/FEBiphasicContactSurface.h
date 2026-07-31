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
#include <FEBioMech/FEContactSurface.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
class FEBIOMIX_API FEBiphasicContactPoint : public FEContactMaterialPoint
{
public:
    vec3d   m_dg;       //!< vector gap
    double  m_Lmd;      //!< Lagrange multipliers for normal traction
    vec3d   m_Lmt;      //!< Lagrange multipliers for vector traction
    double  m_epsn;     //!< penalty factor
    double  m_epsp;     //!< pressure penalty factor
    double  m_p1;       //!< fluid pressure
    vec3d   m_nu;       //!< normal at integration points
    vec3d   m_s1;       //!< tangent along slip direction
    vec3d   m_tr;       //!< contact traction
    vec2d   m_rs;       //!< natural coordinates of projection
    vec2d   m_rsp;      //!< m_rs at the previous time step
    bool    m_bstick;   //!< stick flag
	double  m_Lmp;      //!< lagrange multipliers for fluid pressures
	double	m_pg;       //!< pressure "gap" for biphasic contact
    double  m_mueff;    //!< effective friction coefficient
    double  m_fls;      //!< local fluid load support

	void Init() override
	{
		FEContactMaterialPoint::Init();
        m_Lmd   = 0.0;
        m_Lmt   = m_tr = vec3d(0,0,0);
        m_Lmp   = 0.0;
        m_epsn  = 1.0;
        m_epsp  = 1.0;
        m_pg    = 0.0;
        m_p1    = 0.0;
        m_mueff = 0.0;
        m_fls   = 0.0;
        m_nu    = m_s1 = m_dg = vec3d(0,0,0);
        m_rs    = m_rsp = vec2d(0,0);
        m_bstick = false;
		m_Lmp = 0.0;
		m_pg = 0.0;
		m_mueff = 0.0;
		m_fls = 0.0;
	}

	void Serialize(DumpStream& ar) override;
};

//-----------------------------------------------------------------------------
//! This class describes a contact surface used in a biphasic/multiphasic analysis.
class FEBIOMIX_API FEBiphasicContactSurface : public FEContactSurface
{
public:
	//! constructor
	FEBiphasicContactSurface(FEModel* pfem);

	//! destructor
	~FEBiphasicContactSurface();

	//! initialization
	bool Init();

	//! serialization
	void Serialize(DumpStream& ar);

	void UnpackLM(FEElement& el, vector<int>& lm);

public:
	//! Get the total force exerted by the fluid
    virtual vec3d GetFluidForce();

    //! Get the total force exerted by the fluid
    virtual double GetFluidLoadSupport();
    
    //! Get the effective friction coefficient
    virtual void GetMuEffective    (int nface, double& pg);
    
    //! Get the local fluid load support
    virtual void GetLocalFLS       (int nface, double& pg);
    
    //! Get the local fluid load support projected from the element to the surface Gauss points
    void GetGPLocalFLS(int nface, double* pt, double ambp = 0);

protected:
	int	m_dofP;
};
