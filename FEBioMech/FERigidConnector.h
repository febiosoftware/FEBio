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
#include <FECore/FENLConstraint.h>
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
class FERigidBody;

//-----------------------------------------------------------------------------
//! This is a virtual class for all rigid connectors, including
//! spherical, revolute, prismatic and cylindrical joints, as well
//! as springs and dampers that connect rigid bodies.

class FEBIOMECH_API FERigidConnector : public FENLConstraint
{
    FECORE_BASE_CLASS(FERigidConnector)

public:
    //! constructor
    FERigidConnector(FEModel* pfem);
    
    //! destructor
    virtual ~FERigidConnector();

	//! initialization
	bool Init() override;
    
    int GetConnectorID() { return m_nID; }

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	//! serialization
	void Serialize(DumpStream& ar) override;
    
    //! evaluate relative translation
    virtual vec3d RelativeTranslation(const bool global = false) { return vec3d(0,0,0); }
    
    //! evaluate relative rotation
    virtual vec3d RelativeRotation(const bool global = false) { return vec3d(0,0,0); }
    
public:
    int	m_nRBa;		//!< rigid body A that the connector connects
    int	m_nRBb;		//!< rigid body B that the connector connects
    
    vec3d	m_F;	//! constraining force
    vec3d	m_M;	//! constraining moment
    
protected:
    int		m_nID;		//!< ID of rigid connector
	bool	m_binit;	//!< initialization flag

	FERigidBody*	m_rbA;
	FERigidBody*	m_rbB;
    
    static int	m_ncount;	//!< used to create unique ID's for the nonlinear constraints
    
    DECLARE_FECORE_CLASS();
};
