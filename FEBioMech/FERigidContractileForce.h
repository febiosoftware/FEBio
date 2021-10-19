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

//-----------------------------------------------------------------------------
//! The FERigidContractileForce class implements a contractile force between
//! arbitrary points (not necessarily nodes) on two rigid bodies.

class FERigidContractileForce : public FERigidConnector
{
public:
    //! constructor
    FERigidContractileForce(FEModel* pfem);
    
    //! destructor
    ~FERigidContractileForce() {}
    
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
    vec3d RelativeTranslation(const bool global) override;
    
    //! evaluate relative rotation
    vec3d RelativeRotation(const bool global) override;

public: // parameters
    double	m_f0;       //! contractile force
    vec3d	m_a0;       //! initial absolute position vector of insertion on body A
    vec3d	m_b0;       //! initial absolute position vector of insertion on body B

	// output parameters
	vec3d	m_at;		//!< current absolute position of spring on body A
	vec3d	m_bt;		//!< current absolute position of spring on body B

protected:
    vec3d	m_qa0;      //! initial relative position vector of insertion on body A
    vec3d	m_qb0;      //! initial relative position vector of insertion on body B
    
    DECLARE_FECORE_CLASS();
};
