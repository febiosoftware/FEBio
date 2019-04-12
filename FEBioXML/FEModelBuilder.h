/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include <FECore/FEModel.h>
#include <FEBioMech/RigidBC.h>
#include "febioxml_api.h"
#include <string>

class FESolver;

// This is a helper class for building the FEModel from file input. 
class FEBIOXML_API FEModelBuilder
{
public:
	struct FEBIOXML_API NodeSetPair
	{
		char		szname[256];
		FENodeSet*	pmaster;
		FENodeSet*	pslave;
	};

	struct FEBIOXML_API NodeSetSet
	{
		NodeSetSet() { count = 0; }
		enum { MAX_SETS = 32 };
		char		szname[256];
		FENodeSet*	set[MAX_SETS];
		int			count;

		void add(FENodeSet* ps) { set[count++] = ps; }
	};

public:
	//! constructor
	FEModelBuilder(FEModel& fem);

	//! set the module name
	void SetModuleName(const std::string& moduleName);

	//! Get the module name
	std::string GetModuleName() const;

	// create a new analysis step
	FEAnalysis* CreateNewStep();

	// create a material
	FEMaterial* CreateMaterial(const char* sztype);

	// get the current step (will create a new one if no step was defined yet)
	FEAnalysis*	GetStep();

	// add module component to current step
	void AddComponent(FEModelComponent* mc);

	// reset some data for reading next step
	void NextStep();

	//! Create a domain
	FEDomain* CreateDomain(FE_Element_Spec espec, FEMaterial* mat);

public:
	bool BuildSurface(FESurface& s, FEFacetSet& f, bool bnodal = false);

	bool BuildEdge(FEEdge& s, FESegmentSet& f);

	FE_Element_Spec ElementSpec(const char* sz);

	// Call this to initialize default variables when reading older files.
	void SetDefaultVariables();

public:
	void AddFixedBC(FEFixedBC* pbc);
	void AddPrescribedBC(FEPrescribedBC* pbc);
	void AddNodalLoad(FENodalLoad* pfc);
	void AddEdgeLoad(FEEdgeLoad* pel);
	void AddSurfaceLoad(FESurfaceLoad* psl);
	void AddInitialCondition(FEInitialCondition* pic);
	void AddContactInterface(FESurfacePairConstraint* pci);
	void AddModelLoad(FEModelLoad* pml);
	void AddNonlinearConstraint(FENLConstraint* pnc);

	void AddRigidFixedBC            (FERigidBodyFixedBC* prc);
	void AddRigidPrescribedBC       (FERigidBodyDisplacement* prc);
	void AddRigidBodyVelocity       (FERigidBodyVelocity* prv);
	void AddRigidBodyAngularVelocity(FERigidBodyAngularVelocity* prv);
	void AddRigidNodeSet            (FERigidNodeSet* rs);

public:
	void AddNodeSetPair(NodeSetPair& p) { m_nsetPair.push_back(p); }
	NodeSetPair* FindNodeSetPair(const char* szname);

	void AddNodeSetSet(NodeSetSet& p) { m_nsetSet.push_back(p); }
	NodeSetSet* FindNodeSetSet(const char* szname);

protected:
	FESolver* BuildSolver(FEModel& fem);

public:
	// Build the node ID table
	void BuildNodeList();

	// find a node index from its ID
	int FindNodeFromID(int nid);

	// convert an array of nodal ID to nodal indices
	void GlobalToLocalID(int* l, int n, vector<int>& m);

private:
	FEModel&		m_fem;				//!< model that is being constructed
	FEAnalysis*		m_pStep;			//!< pointer to current analysis step
	int				m_nsteps;			//!< nr of step sections read

public:
	int		m_maxid;		//!< max element ID

	bool	m_b3field_hex;	    //!< three-field element flag for hex (and wedge elements)
	bool	m_b3field_tet;	    //!< three-field element flag for quadratic tets
    bool    m_b3field_shell;    //!< three-field element flag for shells
    bool    m_b3field_quad;     //!< three-field element flag for quad shells
    bool    m_b3field_tri;      //!< three-field element flag for tri shells
	bool	m_but4;				//!< use UT4 formulation flag
	int		m_default_shell;	//!< shell formulation
	double	m_ut4_alpha;		//!< UT4 integration alpha value
	bool	m_ut4_bdev;			//!< UT4 integration deviatoric formulation flag
	double	m_udghex_hg;		//!< hourglass parameter for UDGhex integration
	FE_Element_Type		m_nhex8;	//!< hex integration rule
	FE_Element_Type		m_ntet4;	//!< tet4 integration rule
	FE_Element_Type		m_ntet10;	//!< tet10 integration rule
	FE_Element_Type		m_ntet15;	//!< tet15 integration rule
	FE_Element_Type		m_ntet20;	//!< tet20 integration rule
	FE_Element_Type		m_ntri3;	//!< tri3 integration rule
	FE_Element_Type		m_ntri6;	//!< tri6 integration rule
	FE_Element_Type		m_ntri7;	//!< tri7 integration rule
	FE_Element_Type		m_ntri10;	//!< tri10 integration rule
	FE_Element_Type		m_nquad4;	//!< quad4 integration rule
	FE_Element_Type		m_nquad8;	//!< quad8 integration rule
	FE_Element_Type		m_nquad9;	//!< quad9 integration rule

protected:
	vector<NodeSetPair>	m_nsetPair;
	vector<NodeSetSet>	m_nsetSet;

protected:
	int			m_node_off;		//!< node offset (i.e. lowest node ID)
	vector<int>	m_node_list;	//!< map node ID's to their nodes.
};
