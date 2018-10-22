#pragma once

//-----------------------------------------------------------------------------
// Element Class:
// Defines the general category of element.
enum FE_Element_Class {
	FE_ELEM_INVALID_CLASS,

	FE_ELEM_SOLID,
	FE_ELEM_SHELL,
	FE_ELEM_BEAM,
	FE_ELEM_SURFACE,
	FE_ELEM_TRUSS,
	FE_ELEM_DISCRETE,
	FE_ELEM_2D,
	FE_ELEM_EDGE,

	FE_ELEM_WIRE = 100	// temporary. Can change.
};

//-----------------------------------------------------------------------------
// Element shapes:
// This defines the general element shape classes. This classification differs from the
// element types below, in that the latter is defined by a shape and integration rule.
enum FE_Element_Shape {
	FE_ELEM_INVALID_SHAPE,

	ET_HEX8,
	ET_HEX20,
	ET_HEX27,
	ET_PENTA6,
	ET_TET4,
	ET_TET10,
	ET_TET15,
	ET_TET20,
	ET_QUAD4,
    ET_QUAD8,
    ET_QUAD9,
	ET_TRI3,
    ET_TRI6,
	ET_TRI7,
	ET_TRI10,
	ET_TRUSS2,
	ET_LINE2,
	ET_DISCRETE,
    ET_PENTA15,
	ET_PYRA5
};

//-----------------------------------------------------------------------------
// Element types:
//  Note that these numbers are actually indices into the m_Traits array
//  of the ElementLibrary class so make sure the numbers correspond
//  with the entries into this array
//

enum FE_Element_Type {
	// 3D solid elements
	FE_HEX8G8,	
	FE_HEX8RI,
	FE_HEX8G1,	
	FE_TET4G1,	
	FE_TET4G4,		
	FE_PENTA6G6,	
	FE_TET10G1,
	FE_TET10G4,
	FE_TET10G8,
	FE_TET10GL11,
	FE_TET10G4RI1,
	FE_TET10G8RI4,
	FE_TET15G4,
	FE_TET15G8,
	FE_TET15G11,
	FE_TET15G15,
	FE_TET15G15RI4,
	FE_TET20G15,
    FE_HEX20G8,
	FE_HEX20G27,
	FE_HEX27G27,
    FE_PENTA15G8,
    FE_PENTA15G21,
	FE_PYRA5G8,

	// 2.5D surface elements
	FE_QUAD4G4,
	FE_QUAD4NI,
	FE_TRI3G1,
	FE_TRI3G3,
	FE_TRI3G7,
	FE_TRI3NI,
	FE_TRI6G3,
	FE_TRI6G4,
	FE_TRI6G7,
	FE_TRI6MG7,
	FE_TRI6GL7,
	FE_TRI6NI,
	FE_TRI7G3,
	FE_TRI7G4,
	FE_TRI7G7,
	FE_TRI7GL7,
	FE_TRI10G7,
	FE_TRI10G12,
	FE_QUAD8G9,
    FE_QUAD8NI,
	FE_QUAD9G9,
    FE_QUAD9NI,

	// shell elements
    FE_SHELL_QUAD4G8,
    FE_SHELL_QUAD4G12,
    FE_SHELL_QUAD8G18,
    FE_SHELL_QUAD8G27,
    FE_SHELL_TRI3G6,
    FE_SHELL_TRI3G9,
    FE_SHELL_TRI6G14,
    FE_SHELL_TRI6G21,

	// truss elements
	FE_TRUSS,

	// discrete elements
	FE_DISCRETE,

	// 2D elements
	FE2D_TRI3G1,
	FE2D_TRI6G3,
	FE2D_QUAD4G4,
	FE2D_QUAD8G9,
	FE2D_QUAD9G9,

	// line elements
	FE_LINE2G1,

	// unspecified
	FE_ELEM_INVALID_TYPE = 0xFFFF
};

//-----------------------------------------------------------------------------
// Shell formulations
enum SHELL_FORMULATION {
	NEW_SHELL,
	OLD_SHELL,
	EAS_SHELL,
	ANS_SHELL
};

//-----------------------------------------------------------------------------
//! Helper class for creating domain classes.
struct FE_Element_Spec
{
	FE_Element_Class    eclass;
	FE_Element_Shape	eshape;
	FE_Element_Type		etype;
	bool				m_bthree_field_hex;
	bool				m_bthree_field_tet;
    bool                m_bthree_field_shell;
	bool				m_but4;
	int					m_shell_formulation;

	FE_Element_Spec()
	{
		eclass = FE_ELEM_INVALID_CLASS;
		eshape = FE_ELEM_INVALID_SHAPE;
		etype  = FE_ELEM_INVALID_TYPE;
		m_bthree_field_hex = false;
		m_bthree_field_tet = false;
        m_bthree_field_shell = false;
		m_shell_formulation = NEW_SHELL;
		m_but4 = false;
	}

	bool operator == (const FE_Element_Spec& s)
	{
		if ((eclass == s.eclass) &&
			(eshape == s.eshape) &&
			(etype  == s.etype )) return true;
		return false;
	}
};

//-----------------------------------------------------------------------------
//! This lists the super-class id's that can be used to register new classes
//! with the kernel. It effectively defines the base class that a class
//! is derived from.
enum SUPER_CLASS_ID {
	FEINVALID_ID,					// an invalid ID
	FEOBJECT_ID,					// derived from FECoreBase (TODO: work in progress)
	FETASK_ID,                   	// derived from FECoreTask
	FESOLVER_ID,                 	// derived from FESolver
	FEMATERIAL_ID,               	// derived from FEMaterial
	FEBODYLOAD_ID,               	// derived from FEBodyLoad
	FESURFACELOAD_ID,            	// derived from FESurfaceLoad
	FENLCONSTRAINT_ID,           	// derived from FENLConstraint
	FECOORDSYSMAP_ID,            	// derived from FECoordSysMap
	FEPLOTDATA_ID,               	// derived from FEPlotData
	FEANALYSIS_ID,               	// derived from FEAnalysis
	FESURFACEPAIRINTERACTION_ID, 	// derived from FESurfacePairInteraction
	FENODELOGDATA_ID,            	// derived from FENodeLogData
	FEELEMLOGDATA_ID,            	// derived from FELogElemata
	FEOBJLOGDATA_ID,            	// derived from FELogObjectData
	FEBC_ID,						// derived from FEBoundaryCondition (TODO: This does not work yet)
	FEGLOBALDATA_ID,				// derived from FEGlobalData
	FERIGIDOBJECT_ID,				// derived from FECoreBase (TODO: work in progress)
	FENLCLOGDATA_ID,             	// derived from FELogNLConstraintData
	FECALLBACK_ID,					// derived from FECallBack
	FEDOMAIN_ID,					// derived from FEDomain (TODO: work in progress)
	FEIC_ID,						// derived from initial condition
	FEEDGELOAD_ID,					// derived from FEEdgeLoad
	FEDATAGENERATOR_ID,				// derived from FEDataGenerator
	FELOADCONTROLLER_ID,			// derived from FELoadContoller (TODO: work in progress)
	FEMODEL_ID,						// derived from FEModel (TODO: work in progress)
	FEMODELDATA_ID,					// derived from FEModelData (TODO: work in progress)
	FEVECTORGENERATOR_ID,			// derived from FEVectorGenerator (NOTE: work in progress!)
	FEFUNCTION1D_ID					// derived from FEFunction1D (TODO: work in progress)
};

/////////////////////////////////////////////////////////////////////////////
// ENUM: Linear solvers
//  Defines the supported linear solvers. Note that some of these solvers
//  are only available on certain platforms
//

enum FE_Linear_Solver_Type {
	SKYLINE_SOLVER,
	PSLDLT_SOLVER,		// use only where available
	SUPERLU_SOLVER,		// use only where available
	SUPERLU_MT_SOLVER,	// use only where available
	PARDISO_SOLVER, 	// use only where available
	LU_SOLVER,
	WSMP_SOLVER,		// use only where available
	CG_ITERATIVE_SOLVER,
	RCICG_SOLVER,			// use only where available
	FGMRES_SOLVER,			// use only where available
	FGMRES_ILUT_SOLVER,		// use only where available
	FGMRES_ILU0_SOLVER,		// use only where available
	BIPN_SOLVER,
	HYPRE_GMRES,
	STOKES_SOLVER,
	CG_STOKES_SOLVER,
	SCHUR_SOLVER
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: Analysis types
//  Types of analysis that can be performed
// TODO: Make this a FESolver attribute
enum FE_Analysis_Type {
	FE_STATIC		= 0,
	FE_DYNAMIC		= 1,
	FE_STEADY_STATE	= 2
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: rigid surfaces

enum FE_Rigid_Surface_Type {
	FE_RIGID_PLANE,
	FE_RIGID_SPHERE
};

//-----------------------------------------------------------------------------
// Plot level sets the frequency of writes to the plot file.
enum FE_Plot_Level {
	FE_PLOT_NEVER,			// don't output anything
	FE_PLOT_MAJOR_ITRS,		// only output major iterations (i.e. converged time steps)
	FE_PLOT_MINOR_ITRS,		// output minor iterations (i.e. every Newton iteration)
	FE_PLOT_MUST_POINTS,	// output only on must-points
	FE_PLOT_FINAL,			// only output final converged state
	FE_PLOT_AUGMENTATIONS,	// plot state before augmentations
	FE_PLOT_STEP_FINAL		// output the final step of a step
};

//-----------------------------------------------------------------------------
// Print level sets the frequency of output to the screen and log file
enum FE_Print_Level {
	FE_PRINT_NEVER,
	FE_PRINT_PROGRESS,
	FE_PRINT_MAJOR_ITRS,
	FE_PRINT_MINOR_ITRS,
	FE_PRINT_MINOR_ITRS_EXP,
};

//-----------------------------------------------------------------------------
// Output level sets the frequency of data output is written to the log or data files.
enum FE_Output_Level {
	FE_OUTPUT_NEVER,
	FE_OUTPUT_MAJOR_ITRS,
	FE_OUTPUT_MINOR_ITRS,
	FE_OUTPUT_MUST_POINTS,
	FE_OUTPUT_FINAL
};

//-----------------------------------------------------------------------------
// Dump level determines the times the restart file is written
enum FE_Dump_Level {
	FE_DUMP_NEVER,			// never write a dump file
	FE_DUMP_MAJOR_ITRS,		// create a dump file at the end of each converged time step
	FE_DUMP_STEP			// create a dump file at the end of an analysis step
};

//-----------------------------------------------------------------------------
//! Domain classes
//! The domain class defines the general catergory of element types
//! NOTE: beams are not supported yet.
#define	FE_DOMAIN_SOLID		1
#define	FE_DOMAIN_SHELL		2
#define	FE_DOMAIN_BEAM		3
#define	FE_DOMAIN_SURFACE	4
#define	FE_DOMAIN_TRUSS		5
#define	FE_DOMAIN_DISCRETE	6
#define	FE_DOMAIN_2D		7
#define FE_DOMAIN_EDGE		8

// --- data types ---
enum Var_Type { 
	PLT_FLOAT,		// scalar             : single fp
	PLT_VEC3F,		// 3D vector          : 3 fps
	PLT_MAT3FS,		// symm 2o tensor     : 6 fps
	PLT_MAT3FD,		// diagonal 2o tensor : 3 fps
	PLT_TENS4FS,	// symm 4o tensor     : 21 fps
	PLT_MAT3F,		// 2o tensor          : 9 fps
	PLT_ARRAY,		// variable array (see dictionary for size)
	PLT_ARRAY_VEC3F	// array of vec3f (see dictionary for size)
};

// --- storage format ---
// FMT_NODE : one value stored for each node of a region
// FMT_ITEM : one value stored for each item (e.g. element) of a region
// FMT_MULT : one value for each node of each item of a region
// FMT_REGION: one value per region (surface, domain)
enum Storage_Fmt { FMT_NODE, FMT_ITEM, FMT_MULT, FMT_REGION };

//-----------------------------------------------------------------------------
enum FEDataType {
	FE_DOUBLE,
	FE_VEC2D,
	FE_VEC3D,
	FE_MAT3D
};

//-----------------------------------------------------------------------------
//! Different matrix types. This is used when requesting a sparse matrix format
//! from a linear solver. 
//! \sa LinearSolver::CreateSparseMatrix.
enum Matrix_Type {
	REAL_SYMMETRIC,
	REAL_UNSYMMETRIC,
	COMPLEX_SYMMETRIC,
	COMPLEX_UNSYMMETRIC
};
