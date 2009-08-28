// fem.h: interface for the FEM class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _FEM_H_07012006_
#define _FEM_H_07012006_

#include "Logfile.h"
#include "PlotFile.h"
#include "LoadCurve.h"
#include "FEElementLibrary.h"
#include "FECore/matrix.h"
#include "quatd.h"
#include "Archive.h"
#include "FEMesh.h"
#include "FESurface.h"
#include "FEContactInterface.h"
#include "FESlidingInterface.h"
#include "FERigidWallInterface.h"
#include "FETiedInterface.h"
#include "FEMaterial.h"
#include "FERigidBody.h"
#include "stack.h"
#include "version.h"
#include "FESolver.h"
#include "DataStore.h"
#include "FERigidJoint.h"
#include "FEAnalysis.h"
#include "FEAugLagLinearConstraint.h"

#include <list>
using namespace std;

#define MAX_STRING	256


// Macauley bracket
#define MBRACKET(x) ((x)>=0? (x): 0)

// Heavyside function
#define HEAVYSIDE(x) ((x)>=0?1:0)

// Heavyside function
inline double Heavyside(double x)
{
	const double eps = -1e-07;
	return (x >= eps ? 1 : 0);
}

// smooth step function
inline double smooth_step(double x, double a, double b)
{
	if (b==a) return Heavyside(x);

	double f = (x-a)/(b-a);
	f = (f<0?0:(f>1?1:f));
	return f*f*(3.0-2.0*f);
}

// ramp function
inline double ramp(double x, double a, double b)
{
	if (b==a) return Heavyside(x);

	double f = (x-a)/(b-a);
	return (f<0?0:(f>1?1:f));
}


//-----------------------------------------------------------------------------
//! A degree of freedom structure

class DOF
{
public:
	DOF() { node = bc = neq = -1; }
public:
	int	node;	// the node to which this dof belongs to
	int	bc;		// the degree of freedom
	int	neq;	// the equation number (or -1 if none)
};

//-----------------------------------------------------------------------------
//! linear constraint

class FELinearConstraint
{
public:
	class SlaveDOF : public DOF
	{
	public:
		SlaveDOF() : val(0){}
		double	val;	// coefficient value
	};

public:
	FELinearConstraint(){}
	FELinearConstraint(const FELinearConstraint& LC)
	{
		master = LC.master;
		int n = (int) LC.slave.size();
		list<SlaveDOF>::const_iterator it = LC.slave.begin();
		for (int i=0; i<n; ++i) slave.push_back(*it);
	}

	double FindDOF(int n)
	{
		int N = slave.size();
		list<SlaveDOF>::iterator it = slave.begin();
		for (int i=0; i<N; ++i, ++it) if (it->neq == n) return it->val;

		return 0;
	}

	void Serialize(Archive& ar)
	{
		if (ar.IsSaving())
		{
			ar.write(&master, sizeof(DOF), 1);
			int n = (int) slave.size();
			ar << n;
			list<SlaveDOF>::iterator it = slave.begin();
			for (int i=0; i<n; ++i, ++it) ar << it->val << it->node << it->bc << it->neq;
		}
		else
		{
			slave.clear();
			ar.read(&master, sizeof(DOF), 1);
			int n;
			ar >> n;
			for (int i=0; i<n; ++i)
			{
				SlaveDOF dof;
				ar >> dof.val >> dof.node >> dof.bc >> dof.neq;
				slave.push_back(dof);
			}
		}
	}

public:
	DOF			master;	// master degree of freedom
	list<SlaveDOF>	slave;	// list of slave nodes
};

//-----------------------------------------------------------------------------
//! concentrated nodal force boundary condition

class FENodalForce : public FEBoundaryCondition
{
public:
	double	s;		// scale factor
	int		node;	// node number
	int		bc;		// force direction
	int		lc;		// load curve
};

//-----------------------------------------------------------------------------
//! prescribed nodal displacement data

class FENodalDisplacement : public FEBoundaryCondition
{
public:
	double	s;		// scale factor
	int		node;	// node number
	int		bc;		// displacement direction
	int		lc;		// load curve
};

//-----------------------------------------------------------------------------
//! rigid node

class FERigidNode : public FEBoundaryCondition
{
public:
	int	nid;	// node number
	int	rid;	// rigid body number
};

//-----------------------------------------------------------------------------
//! This class describes a pressure load on a surface element

class FEPressureLoad : public FEBoundaryCondition
{
public:
	FEPressureLoad() { s[0] = s[1] = s[2] = s[3] = 1.0; blinear = false; }

public:
	double	s[4];		// nodal scale factors
	int		face;		// face number
	int		lc;			// load curve
	bool	blinear;	// linear or not (true is non-follower, false is follower)
};

///////////////////////////////////////////////////////////////////////////////
// STRUCT: FE_BODYFORCE
// body force data
//

struct FE_BODYFORCE
{
	double	s;		// scale factor
	int		lc;		// load curve number
};

// forward declaration of FESolver class
class FESolver;

// forward declaration of the FEM class
class FEM;

// define the FEBIO callback function
typedef void (*FEBIO_CALLBACK)(FEM*,void*);

//-----------------------------------------------------------------------------
//! Discrete element

struct FE_DISCRETE_ELEMENT
{
	int			n1, n2;	//!< nodes that this spring connects
	double		E;		//!< spring constant
};

//-----------------------------------------------------------------------------
//! The Finite Element Model class. 

//! This class stores solver parameters, geometry data, material data, and 
//! other data that is needed to solve the FE problem.
//! FEBio is designed to solve finite element problems. All the finite element
//! data is collected here in this class. This class also provides
//! routines to initalize, input, output and update the FE data. Although this
//! class provides the main solve routine it does not really solve anything.
//! The actual solving is done by the FESolver class.

class FEM
{
public:
	//! constructor - sets default variables
	FEM();

	//! destructor
	virtual ~FEM();

	//! read the configuration file
	bool Configure(const char* szfile);

	//! Restart from restart point
	bool Restart(const char* szfile);

	//! Initializes data structures
	bool Init();

	//! check data
	bool Check();

	//! Resets data structures
	bool Reset();

	//! Solves the problem
	bool Solve();

	//! Serialize the current state to/from restart file
	bool Serialize(Archive& ar);

	//! input data from file
	bool Input(const char* szfile);

	//! echoes the input data to the logfile
	void EchoInput();

	//! set callback function
	void SetCallback(FEBIO_CALLBACK pcb, void* pd) { m_pcb = pcb; m_pcbd = pd; }

	//! Add a material to the model
	void AddMaterial(FEMaterial* pm) { m_MAT.add(pm); }

	//! Add a parameter list
	void AddParameterList(FEParameterList* pl) { m_MPL.add(pl); }

	//! get the number of materials
	int Materials() { return m_MAT.size(); }

	//! return a pointer to a material
	FEMaterial* GetMaterial(int id) { return &m_MAT[id]; }

	//! return the elastic material
	FEElasticMaterial* GetElasticMaterial(int id)
	{
		FEMaterial* pm = &m_MAT[id];
		if (dynamic_cast<FENestedMaterial*>(pm)) pm = (dynamic_cast<FENestedMaterial*>(pm))->m_pBase;
		FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm);
		assert(pme);
		return pme;
	}

	//! Add a loadcurve to the model
	void AddLoadCurve(FELoadCurve* plc) { m_LC.add(plc); }

	//! get a loadcurve
	FELoadCurve* GetLoadCurve(int i) { return &m_LC[i]; }

	//! get the number of loadcurves
	int LoadCurves() { return m_LC.size(); }

	//! set the debug level
	void SetDebugFlag(bool b) { m_debug = b; }

	//! get the debug level
	bool GetDebugFlag() { return m_debug; }


	// set the i/o files
	void SetInputFilename(const char* szfile)
	{ 
		strcpy(m_szfile, szfile); 
		m_szfile_title = strrchr(m_szfile, '/');
		if (m_szfile_title == 0) 
		{
			m_szfile_title = strchr(m_szfile, '\\'); 
			if (m_szfile_title == 0) m_szfile_title = m_szfile; else ++m_szfile_title;
		}
		else ++m_szfile_title;
	}
	void SetLogFilename  (const char* szfile) { strcpy(m_szlog , szfile); }
	void SetPlotFilename (const char* szfile) { strcpy(m_szplot, szfile); }
	void SetDumpFilename (const char* szfile) { strcpy(m_szdump, szfile); }

	const char* GetLogfileName () { return m_szlog; }
	const char* GetPlotFileName() { return m_szplot; }

	//! set the problem title
	void SetTitle(const char* sz) { strcpy(m_sztitle, sz); }

	//! get the problem title
	const char* GetTitle() { return m_sztitle; }

	//! return a pointer to the named variable
	double* FindParameter(const char* szname);

	//! return number of contact interfaces
	int ContactInterfaces() { return m_CI.size(); } 

	//! find a boundary condition from the ID
	FEBoundaryCondition* FindBC(int nid);

	//! Set the sparse matrix symmetry flag
	void SetSymmetryFlag(bool bsymm) { m_bsymm = bsymm; }

protected:
	// copy constructor and assignment operator are protected since they
	// are reserved for a special purpose and cannot be used in a way
	// that the user may think.
	
	//! copy constructor
	FEM(FEM& fem);

	//! assignment operator
	void operator = (FEM& fem);

protected:
	void SerializeMaterials   (Archive& ar);
	void SerializeAnalysisData(Archive& ar);
	void SerializeGeometry    (Archive& ar);
	void SerializeContactData (Archive& ar);
	void SerializeBoundaryData(Archive& ar);
	void SerializeIOData      (Archive& ar);
	void SerializeLoadData    (Archive& ar);

public:
	//{ --- Initialization routines ---

		//! initialze equation numbering
		bool InitEquations();

		//! Initialize rigid bodies
		bool InitRigidBodies();

		//! Initialize poro-elastic data
		bool InitPoro();

		//! Initializes contact data
		bool InitContact();

		//! Iniatialize linear constraint data
		bool InitConstraints();
	//}


	//{ --- Update routines ---

		//! Update contact data
		void UpdateContact();
	//}

	//{ --- Miscellaneous routines ---

		//! see if body forces are present
		bool UseBodyForces() { return ((m_BF[0].lc>=0) || (m_BF[1].lc>=0) || (m_BF[2].lc>=0)); }

		//! call the callback function
		void DoCallback();

	//}

public:
	// --- Analysis Data ---
	//{
		ptr_vector<FEAnalysis>	m_Step;		//!< array of analysis steps
		int						m_nStep;	//!< current analysis step
		FEAnalysis*				m_pStep;	//!< pointer to current analysis step
		double					m_ftime;	//!< current time value
		int	m_nhex8;						//!< element type for hex8
		bool	m_bsym_poro;		//!< symmetric (old) poro-elastic flag
		bool	m_bplane_strain;	//!< run analysis in plain strain mode

		// body force loads
		FE_BODYFORCE	m_BF[3];		//!< body force data
		vec3d			m_acc;			//!< acceleration due to bodyforces
	//}

	// --- Geometry Data ---
	//{
		FEMesh	m_mesh;	//!< the FE mesh

		FESurface*	m_psurf;	//!< surface element array

		// rigid body data
		int						m_nreq;	//!< start of rigid body equations
		int						m_nrm;	//!< nr of rigid materials
		int						m_nrb;	//!< nr of rigid bodies in problem
		vector<FERigidBody>		m_RB;	//!< rigid body array

		// rigid joints
		int							m_nrj;	//!< nr of rigid joints
		ptr_vector<FERigidJoint>	m_RJ;	//!< rigid joint array

		// discrete elements
		vector<FE_DISCRETE_ELEMENT>	m_DE;	//!< discrete elements
	//}

	//{ --- Contact Data --
		bool								m_bcontact;	//!< contact flag
		ptr_vector<FEContactInterface>		m_CI;		//!< contact interface array
	//}

protected:
	// --- Material Data ---
	//{
		ptr_vector<FEMaterial>	m_MAT;		//!< array of materials
		ptr_vector<FEParameterList>	m_MPL;	//!< material parameter lists
	//}

	// --- Load Curve Data ---
	//{
		ptr_vector<FELoadCurve>	m_LC;	//!< load curve data
	//}

public:
	// --- Boundary Condition Data ---
	//{
		// displacement boundary data
		vector<FENodalDisplacement>		m_DC;	//!< prescribed displacement cards

		// concentrated nodal loads data
		vector<FENodalForce>	m_FC;		//!< concentrated nodal force cards

		// pressure boundary data
		vector<FEPressureLoad>	m_PC;		//!< pressure boundary cards

		// rigid displacements
		vector<FERigidBodyDisplacement>	m_RDC;	//!< rigid body displacements

		// rigid forces
		vector<FERigidBodyForce>	m_RFC;	//!< rigid body forces

		// rigid nodes
		vector<FERigidNode>		m_RN;		//!< rigid nodes

		// linear constraint data
		list<FELinearConstraint>	m_LinC;		//!< linear constraints data
		vector<int>					m_LCT;		//!< linear constraint table
		vector<FELinearConstraint*>	m_LCA;		//!< linear constraint array (temporary solution!)

		// Augmented Lagrangian linear constraint data
		list<FELinearConstraintSet*>	m_LCSet;	//!< aug lag linear constraint data
	//}

	// --- Direct Solver Data ---
	//{
		int		m_nsolver;	//!< type of solver selected
		int		m_neq;		//!< number of equations
		int		m_bwopt;	//!< bandwidth optimization flag
		bool	m_bsymm;	//!< symmetric flag
	//}
 
	// --- I/O-Data --- 
	//{
		Logfile		m_log;		//!< the log file
		PlotFile	m_plot;		//!< the plot file
		DataStore	m_Data;		//!< the data store used for data logging

protected:
		// file names
		char*	m_szfile_title;			//!< master input file title 
		char	m_szfile[MAX_STRING];	//!< master input file name (= path + title)
		char	m_szplot[MAX_STRING];	//!< plot output file name
		char	m_szlog [MAX_STRING];	//!< log output file name
		char	m_szdump[MAX_STRING];	//!< dump file name

		char	m_sztitle[MAX_STRING];	//!< problem title

		bool	m_debug;	//!< debug flag

		FEBIO_CALLBACK		m_pcb;	//!< pointer to callback function
		void*				m_pcbd;	//!< pointer to callback data
	//}

	// some friends of this class
	friend class FEAnalysis;
	friend class FESolver;
	friend class stack<FEM>;
};

#endif // _FEM_H_07012006_
