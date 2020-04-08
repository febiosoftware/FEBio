/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "febiomech_api.h"

class FENodeSet;

//-----------------------------------------------------------------------------
class FEBIOMECH_API FERigidBC : public FEModelComponent
{
	FECORE_SUPER_CLASS

public:
	FERigidBC(FEModel* fem) : FEModelComponent(fem) {}
};

//-----------------------------------------------------------------------------
//! rigid node set
class FEBIOMECH_API FERigidNodeSet : public FERigidBC
{
public:
	enum SHELL_BC {
		HINGED_SHELL,
		CLAMPED_SHELL
	};

public:
	FERigidNodeSet(FEModel* pfem);
	FERigidNodeSet(const FERigidNodeSet& rs);
	void operator = (const FERigidNodeSet& rs);

	bool Init() override;

	void Serialize(DumpStream& ar) override;

	void Activate() override;
	void Deactivate() override;

	void SetRigidMaterialID(int rid);

	void SetNodeSet(FENodeSet& ns);

	void AddNode(int nid);

	void SetShellBC(SHELL_BC bc);

private: // parameters
	int			m_rigidMat;		//!< rigid body's material
	int			m_nshellBC;		//!< flag defining how shells are attached (0=hinged, 1=clamped)

private:
	vector<int>		m_node;	// node number

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! fixed rigid body constraint
class FEBIOMECH_API FERigidBodyFixedBC : public FERigidBC
{
public:
	FERigidBodyFixedBC(FEModel* pfem);

	bool Init();

	void Serialize(DumpStream& ar);

	void Activate();

	void Deactivate();

public:
	int				m_rigidMat;	//!< rigid body material ID
	vector<int>		m_dofs;		//!< constrained dof list

private:
	bool	m_binit;
	int		m_rb;		//!< rigid body's ID

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! rigid body displacement

class FEBIOMECH_API FERigidBodyDisplacement : public FERigidBC
{
public:
	FERigidBodyDisplacement(FEModel* pfem);

	bool Init() override;

	double Value();

	void Serialize(DumpStream& ar) override;

	void Activate() override;

	void Deactivate() override;

	void SetID(int id) { m_rigidMat = id; }
	int GetID() const { return m_rigidMat; }

	void SetBC(int bc) { m_dof = bc; }
	int GetBC() const { return m_dof; }

	void SetRelativeFlag(bool b) { m_brel = b; }
	void SetValue(double v) { m_val = v; }

private:
	int		m_rigidMat;		//!< rigid body material id
	int		m_dof;		//!< displacement direction
	double	m_val;	//!< displacement value
	double	m_ref;	//!< reference value for relative displacement
	bool	m_brel;	//!< relative displacement flag

private:
	bool	m_binit;	//!init flag

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! rigid body initial velocity
class FEBIOMECH_API FERigidBodyVelocity : public FERigidBC
{
public:
	FERigidBodyVelocity(FEModel* pfem);

	bool Init() override;

	void Activate() override;

public:
	int		m_rid;	//!< rigid body ID
	vec3d	m_vel;	//!< initial velocity

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! rigid body initial angular velocity
class FEBIOMECH_API FERigidBodyAngularVelocity : public FERigidBC
{
public:
	FERigidBodyAngularVelocity(FEModel* pfem);

	bool Init() override;

	void Activate() override;

public:
	int		m_rid;	//!< rigid body ID
	vec3d	m_w;	//!< value

	DECLARE_FECORE_CLASS();
};
