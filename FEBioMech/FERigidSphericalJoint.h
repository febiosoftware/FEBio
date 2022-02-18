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
#include "FECore/vec3d.h"
#include "FERigidConnector.h"
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
//! The FERigidSphericalJoint class implements a rigid spherical joint that
//! functions under static and dynamic conditions. This joint allows the
//! user to connect two rigid bodies at a point in space.

class FEBIOMECH_API FERigidSphericalJoint : public FERigidConnector
{
public:
    //! constructor
    FERigidSphericalJoint(FEModel* pfem);
    
    //! destructor
    virtual ~FERigidSphericalJoint() {}
    
    //! initialization
    bool Init() override;
    
    //! calculates the joint forces
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! calculates the joint stiffness
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
    
    //! calculate Lagrangian augmentation
    bool Augment(int naug, const FETimeInfo& tp) override;
    
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;
    
    //! update state
	void Update() override;
    
    //! Reset data
    void Reset() override;
    
    //! evaluate relative translation
    vec3d RelativeTranslation(const bool global = false) override;
    
    //! evaluate relative rotation
    vec3d RelativeRotation(const bool global = false) override;
    
    //! initial position
    vec3d InitialPosition() const;

    //! current position
    vec3d Position() const;

public: // parameters
    double	m_atol;     //! augmented Lagrangian tolerance
    double  m_gtol;     //! augmented Lagrangian gap tolerance
    double  m_qtol;     //! augmented Lagrangian angular gap tolerance
    double	m_eps;      //! penalty factor for constraining force
    double	m_ups;      //! penalty factor for constraining moment
    vec3d	m_q0;       //! initial position of joint
    int     m_naugmin;  //! minimum number of augmentations
    int     m_naugmax;  //! maximum number of augmentations
    bool    m_bq;           //! flag for prescribing rotation
    double  m_qpx;          //! prescribed rotation x-component
    double  m_qpy;          //! prescribed rotation y-component
    double  m_qpz;          //! prescribed rotation z-component
    double  m_Mpx;          //! prescribed moment along x
    double  m_Mpy;          //! prescribed moment along y
    double  m_Mpz;          //! prescribed moment along z
	bool	m_bautopen;	//!< auto-penalty for gap and ang tolerance

protected:
    vec3d	m_qa0;      //! initial relative position vector of joint w.r.t. A
    vec3d	m_qb0;      //! initial relative position vector of joint w.r.t. B
   
    vec3d	m_e0[3];        //! initial joint basis
    vec3d	m_ea0[3];       //! initial joint basis w.r.t. A
    vec3d	m_eb0[3];       //! initial joint basis w.r.t. B
    
    vec3d	m_L;        //! Lagrange multiplier for constraining force
    vec3d	m_U;        //! Lagrange multiplier for constraining moment
    
    DECLARE_FECORE_CLASS();
};
