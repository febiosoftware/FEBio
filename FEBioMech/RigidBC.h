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
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FEInitialCondition.h>
#include "febiomech_api.h"

class FERigidBody;

//-----------------------------------------------------------------------------
class FEBIOMECH_API FERigidBC : public FEBoundaryCondition
{
	FECORE_BASE_CLASS(FERigidBC)

public:
	FERigidBC(FEModel* fem);

	bool Init() override;

	void Serialize(DumpStream& ar) override;

	FERigidBody& GetRigidBody();

	void CopyFrom(FEBoundaryCondition* pbc) override;

private:
	int		m_rb;			// rigid body ID

protected:
	int		m_rigidMat;		// rigid material ID
	bool	m_binit;
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API FERigidIC : public FEInitialCondition
{
	FECORE_BASE_CLASS(FERigidIC)

public:
	FERigidIC(FEModel* fem);

	bool Init() override;

	void Activate() override;

	void Serialize(DumpStream& ar) override;

	FERigidBody& GetRigidBody();

private:
	int		m_rigidMat;		// rigid material ID
	int		m_rb;			// rigid body ID

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! fixed rigid body constraint
class FEBIOMECH_API FERigidFixedBC : public FERigidBC
{
public:
	FERigidFixedBC(FEModel* pfem);
};

//-----------------------------------------------------------------------------
//! fixed rigid body constraint
class FEBIOMECH_API FERigidFixedBCNew : public FERigidFixedBC
{
public:
	FERigidFixedBCNew(FEModel* pfem);

	void Activate() override;

	void Deactivate() override;

public:
	bool	m_dof[6];		//!< constrained dof list

private:
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! fixed rigid body constraint
class FEBIOMECH_API FERigidFixedBCOld : public FERigidFixedBC
{
public:
	FERigidFixedBCOld(FEModel* pfem);

	bool Init() override;

	void Activate() override;

	void Deactivate() override;

public:
	vector<int>	m_dofs;		//!< constrained dof list

	DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
//! rigid body displacement

class FEBIOMECH_API FERigidPrescribedBC : public FERigidBC
{
public:
	FERigidPrescribedBC(FEModel* pfem);

	bool Init() override;

	double Value();

	void InitTimeStep();

	void Serialize(DumpStream& ar) override;

	void Activate() override;

	void Deactivate() override;

	void SetBC(int bc) { m_dof = bc; }
	int GetBC() const { return m_dof; }

	void SetRelativeFlag(bool b) { m_brel = b; }
	bool GetRelativeFlag() const { return m_brel; }

	void SetValue(double v) { m_val = v; }

protected:
	int		m_dof;		//!< displacement direction
	double	m_val;	//!< displacement value
	double	m_ref;	//!< reference value for relative displacement
	bool	m_brel;	//!< relative displacement flag

	bool	m_binit;	//!init flag
};

//-----------------------------------------------------------------------------
//! prescribed rigid body rotation

class FEBIOMECH_API FERigidDisplacement : public FERigidPrescribedBC
{
public:
	FERigidDisplacement(FEModel* pfem);

	bool Init() override;

private:
	int	m_bc;
	DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
//! prescribed rigid body rotation

class FEBIOMECH_API FERigidRotation : public FERigidPrescribedBC
{
public:
	FERigidRotation(FEModel* pfem);

	bool Init() override;

private:
	int	m_bc;
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Obsolete rigid prescribed bc class. Used only for backward compatibility.
class FEBIOMECH_API FERigidPrescribedOld : public FERigidPrescribedBC
{
public:
	FERigidPrescribedOld(FEModel* pfem);
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! rigid body initial velocity
class FEBIOMECH_API FERigidBodyVelocity : public FERigidIC
{
public:
	FERigidBodyVelocity(FEModel* pfem);

	void Activate() override;

public:
	vec3d	m_vel;	//!< initial velocity

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! rigid body initial angular velocity
class FEBIOMECH_API FERigidBodyAngularVelocity : public FERigidIC
{
public:
	FERigidBodyAngularVelocity(FEModel* pfem);

	void Activate() override;

public:
	vec3d	m_w;	//!< value

	DECLARE_FECORE_CLASS();
};
