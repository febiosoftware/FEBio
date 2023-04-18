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
#include <vector>
#include <string>
#include "fecore_api.h"

class DumpStream;

//-----------------------------------------------------------------------------
// Degree of freedom types
#define DOF_OPEN		 0		// the dof is open and will be given an equation number
#define DOF_FIXED		 1		// the dof is fixed and will not be given an equation number
#define DOF_PRESCRIBED	 2		// the dof is prescribed. It will be given a negative equation number (equation = - index - 2)

//-----------------------------------------------------------------------------
// Types of variables
#define VAR_SCALAR	0		// scalar variable
#define VAR_VEC2	1		// 2D vector
#define VAR_VEC3	2		// 3D vector
#define VAR_ARRAY	3		// a variable array of scalars

//-----------------------------------------------------------------------------
//! Class that manages the variables and degrees of freedoms.
class FECORE_API DOFS
{
	// Class representing an individual degree of freedom
	class DOF_ITEM
	{
	public:
		enum { MAX_DOF_NAME = 8 };

	public:
		DOF_ITEM();
		DOF_ITEM(const char* sz);
		~DOF_ITEM();
		DOF_ITEM(const DOF_ITEM& d);
		void operator = (const DOF_ITEM& d);

		void SetName(const char* szdof);

		void Serialize(DumpStream& ar);

	public:
		char	sz[MAX_DOF_NAME];	//!< symbol of variable
		int		ndof;				//!< index of degree of freedom
		int		nvar;				//!< variable associated with this dof
	};

public:
	// A Variable is a logical grouping of degrees of freedoms.
	// (e.g. a displacement variable in 3D has 3 dofs.)
	// The order variable is the interpolation order that should be
	// assumed for this variable. By default, it is set to -1, which 
	// should be interpreted as the interpolation order implied by the 
	// nodes of the element type.
	class Var
	{
	public:
		Var();
		Var(const Var& v);
		void operator = (const Var& v);

		void Serialize(DumpStream& ar);

	public:
		int		m_ntype;				// type of variable
		int		m_order;				// interpolation order
		std::string				m_name;	// variable name
		std::vector<DOF_ITEM>	m_dof;	// list of dofs for this variable
	};

public:
	// constructors
	DOFS();
	DOFS(const DOFS& dofs);
	DOFS& operator = (const DOFS& dofs);
    
	//! destructor
	~DOFS();

	//! Reset dofs
	void Reset();

public:
	//! Add a (empty) variable
	//! returns the new variable's index (or -1 if the variable already exists)
	int AddVariable(const char* szname, int ntype = VAR_SCALAR);

	//! Get number of variables
	int Variables() const;

	//! Get the index of a variable (returns -1 if the variable does not exist)
	int GetVariableIndex(const char* szname);

	//! Get the variable name
	std::string GetVariableName(int n) const;

	// Add a degree of freedom
	// Can only be used on array variables.
	// returns >= 0 on success and -1 on failure (e.g. if the dof is already defined)
	int AddDOF(const char* szvar, const char* sz);

	// Add a degree of freedom
	// Can only be used on array variables.
	// returns >= 0 on success and -1 on failure (e.g. if the dof is already defined)
	int AddDOF(int nvar, const char* sz);

	//! return a dof index from the dof symbol
	//! this function returns -1 if the symbol is not recognized
	int GetDOF(const char* szdof, const char* varName = nullptr);

	//! return a list of dofs for a variable
	void GetDOFList(const char* varName, std::vector<int>& dofs);

	//! return a list of dofs for a variable
	void GetDOFList(int nvar, std::vector<int>& dofs);

	//! return a list of dofs from comma seperated list dof symbols
	bool ParseDOFString(const char* sz, std::vector<int>& dofs, const char* varName = nullptr);

	//! return a dof index from a variable
	//! this function returns -1 if the symbol is not recognized
	int GetDOF(const char* szvar, int n);

	//! return the index into the variable's dof list
	int GetIndex(const char* varName, const char* szdof);

	//! return a dof index from a variable index
	//! this function returns -1 if the symbol is not recognized
	int GetDOF(int nvar, int n);

	//! return the size of a variable
	//! (returns -1 if the variable is not defined)
	int GetVariableSize(const char* sz);

	//! return the size of a variable
	//! (returns -1 if the variable is not defined)
	int GetVariableSize(int nvar);

	//! return the type of a variable
	int GetVariableType(int nvar);

	//! return the total number of degrees of freedom
	int GetTotalDOFS() const;

	//! Set the name of a DOF
	void SetDOFName(const char* szvar, int n, const char* szname);
	void SetDOFName(int nvar, int n, const char* szname);

	//! Get a dof name
	const char* GetDOFName(int nvar, int n);
	const char* GetDOFName(int ndof);

	//! serialization
	void Serialize(DumpStream& ar);

	//! set the interpolation order for a variable
	void SetVariableInterpolationOrder(int nvar, int order);

	// return the interpolation order of a variable
	int GetVariableInterpolationOrder(int nvar);

	// Find the variable from a dof
	int FindVariableFromDOF(int ndof);

	// return the interpolation order for a degree of freedom
	int GetDOFInterpolationOrder(int ndof);

private:
	void Update();
	Var* GetVariable(const char* szvar);
	DOF_ITEM* GetDOFPtr(const char* szdof);
    
private:
	std::vector<Var>	m_var;		//!< array of variables
    int					m_maxdofs;  //!< total number of nodal DOFS (i.e. size of FENode::m_val array)
};

