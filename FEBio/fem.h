// fem.h: interface for the FEM class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _FEM_H_07012006_
#define _FEM_H_07012006_

#include "PlotFile.h"
#include "FECore/LoadCurve.h"
#include "FECore/DumpFile.h"
#include "FEMesh.h"
#include "FEContactInterface.h"
#include "FEMaterial.h"
#include "FERigidBody.h"
#include "DataStore.h"
#include "FERigidJoint.h"
#include "FEAnalysis.h"
#include "FELinearConstraint.h"
#include "FEAugLagLinearConstraint.h"
#include "Timer.h"
#include "FESurfaceLoad.h"
#include "FEBodyForce.h"
#include "FEPointConstraint.h"
#include "FECore/FEModel.h"

#include <stack>
#include <list>
using namespace std;

#define MAX_STRING	256

//-----------------------------------------------------------------------------
// forward declaration of the FEM class
class FEM;

//-----------------------------------------------------------------------------
// FEBIO callback structure
typedef void (*FEBIO_CB_FNC)(FEM*,void*);
struct FEBIO_CALLBACK {
	FEBIO_CB_FNC	m_pcb;
	void*	m_pd;
};

//-----------------------------------------------------------------------------
//! The Finite Element Model class. 

//! This class stores solver parameters, geometry data, material data, and 
//! other data that is needed to solve the FE problem.
//! FEBio is designed to solve finite element problems. All the finite element
//! data is collected here in this class. This class also provides
//! routines to initalize, input, output and update the FE data. Although this
//! class provides the main solve routine it does not really solve anything.
//! The actual solving is done by one of the classes derived from the FESolver class.

class FEM : public FEModel
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

	//! Resets data structures
	bool Reset();

	//! Solves the problem
	bool Solve();

	//! Serialize the current state to/from restart file
	bool Serialize(DumpFile& ar);

	//! input data from file
	bool Input(const char* szfile);

	//! Add a material to the model
	void AddMaterial(FEMaterial* pm) { m_MAT.push_back(pm); }

	//! get the number of materials
	int Materials() { return m_MAT.size(); }

	//! return a pointer to a material
	FEMaterial* GetMaterial(int id) { return m_MAT[id]; }

	//! return the elastic material
	FEElasticMaterial* GetElasticMaterial(int id);

	//! return the elastic material
	static FEElasticMaterial* GetElasticMaterial(FEMaterial* pm);

	//! set the debug level
	void SetDebugFlag(bool b) { m_debug = b; }

	//! get the debug level
	bool GetDebugFlag() { return m_debug; }

	// set the i/o files
	void SetInputFilename(const char* szfile);
	void SetLogFilename  (const char* szfile) { strcpy(m_szlog , szfile); }
	void SetPlotFilename (const char* szfile) { strcpy(m_szplot, szfile); }
	void SetDumpFilename (const char* szfile) { strcpy(m_szdump, szfile); }

	void SetPlotFileNameExtension(const char* szext);

	// Get the I/O files
	const char* GetInputFileName () { return m_szfile; }
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

	static map<string, double> m_Const;
	static void SetGlobalConstant(const string& s, double v);
	static double GetGlobalConstant(const string& s);
	
public:
	//! copy constructor
	FEM(const FEM& fem);

	//! assignment operator
	void operator = (const FEM& fem);

protected:
	void ShallowCopy(FEM& fem);


protected:
	void SerializeMaterials   (DumpFile& ar);
	void SerializeAnalysisData(DumpFile& ar);
	void SerializeGeometry    (DumpFile& ar);
	void SerializeContactData (DumpFile& ar);
	void SerializeBoundaryData(DumpFile& ar);
	void SerializeIOData      (DumpFile& ar);
	void SerializeLoadData    (DumpFile& ar);
	void SerializeConstants   (DumpFile& ar);

public:
	//{ --- Initialization routines ---

		//! initialze equation numbering
		bool InitEquations();

		//! Initialize rigid bodies
		bool InitRigidBodies();

		//! Initialize poroelastic/biphasic and solute data
		bool InitPoroSolute();

		//! Initializes contact data
		bool InitContact();

		//! Iniatialize linear constraint data
		bool InitConstraints();

		//! Initialize material data
		bool InitMaterials();
	//}

	//{ --- Update routines ---

		//! Update contact data
		void UpdateContact();
	//}

	//{ --- Miscellaneous routines ---

		//! set callback function
		void AddCallback(FEBIO_CB_FNC pcb, void* pd);

		//! call the callback function
		void DoCallback();
	//}

public:
	// --- Analysis Data ---
	//{
		vector<FEAnalysis*>		m_Step;		//!< array of analysis steps
		int						m_nStep;	//!< current analysis step
		FEAnalysis*				m_pStep;	//!< pointer to current analysis step
		double					m_ftime;	//!< current time value
		double					m_ftime0;	//!< start time of current step
		bool	m_bsym_poro;		//!< symmetric (old) poro-elastic flag
		int		m_nplane_strain;	//!< run analysis in plain strain mode

		// body force loads
		vector<FEBodyForce*>	m_BF;		//!< body force data

		// Create timer to track total running time
		Timer	m_TotalTime;
	//}

	// --- Geometry Data ---
	//{
		FEMesh	m_mesh;	//!< the FE mesh

		// rigid body data
		int						m_nreq;	//!< start of rigid body equations
		int						m_nrm;	//!< nr of rigid materials
		int						m_nrb;	//!< nr of rigid bodies in problem
		vector<FERigidBody>		m_RB;	//!< rigid body array

		// rigid joints
		vector<FERigidJoint*>		m_RJ;	//!< rigid joint array
	//}

	//{ --- Contact Data --
		vector<FEContactInterface*>		m_CI;		//!< contact interface array
	//}

protected:
	// --- Material Data ---
	//{
		vector<FEMaterial*>			m_MAT;	//!< array of materials
	//}

public:
	// --- Boundary Condition Data ---
	//{
		// surface loads
		vector<FESurfaceLoad*>	m_SL;		//!< surface loads

		// rigid displacements
		vector<FERigidBodyDisplacement*>	m_RDC;	//!< rigid body displacements

		// rigid forces
		vector<FERigidBodyForce*>	m_RFC;	//!< rigid body forces

		// linear constraint data
		list<FELinearConstraint>	m_LinC;		//!< linear constraints data
		vector<int>					m_LCT;		//!< linear constraint table
		vector<FELinearConstraint*>	m_LCA;		//!< linear constraint array (temporary solution!)

		// Augmented Lagrangian linear constraint data
		list<FELinearConstraintSet*>	m_LCSet;	//!< aug lag linear constraint data

		// Point constriant data
		vector<FEPointConstraint>	m_PC;		//!< point constraint data
	//}

	// --- Direct Solver Data ---
	//{
		int		m_nsolver;	//!< type of solver selected
		int		m_neq;		//!< number of equations
		int		m_npeq;		//!< number of equations related to pressure dofs
		int		m_nceq;		//!< number of equations related to concentration dofs
		int		m_bwopt;	//!< bandwidth optimization flag
		bool	m_bsymm;	//!< symmetric flag
	//}
 
	// --- I/O-Data --- 
	//{
		PlotFile*	m_plot;		//!< the plot file
		DataStore	m_Data;		//!< the data store used for data logging
		
		bool		m_bInterruptable;	//!< true if this model can be interrupted with ctrl+c

protected:
		// file names
		char*	m_szfile_title;			//!< master input file title 
		char	m_szfile[MAX_STRING];	//!< master input file name (= path + title)
		char	m_szplot[MAX_STRING];	//!< plot output file name
		char	m_szlog [MAX_STRING];	//!< log output file name
		char	m_szdump[MAX_STRING];	//!< dump file name

		char	m_sztitle[MAX_STRING];	//!< problem title

		bool	m_debug;			//!< debug flag

		list<FEBIO_CALLBACK>	m_pcb;	//!< pointer to callback function
	//}

	// some friends of this class
	friend class FEAnalysis;
	friend class stack<FEM>;
};

#endif // _FEM_H_07012006_
