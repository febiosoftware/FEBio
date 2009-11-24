// FE_enum.h: defines enumerations used in FEBio
//
//////////////////////////////////////////////////////////////////////

#ifndef _FE_ENUM_H_05132007_
#define _FE_ENUM_H_05132007_

///////////////////////////////////////////////////////////////////////////////
// ENUM: Element types
//  Note that these numbers are actually indices into the m_Traits array
//  of the ElementLibrary class so make sure the numbers correspond
//  with the entries into this array
//

enum FE_Element_Type {
	// 3D solid elements
	FE_HEX,			// = 0
	FE_RIHEX,		// = 1
	FE_UDGHEX,		// = 2
	FE_TET,			// = 3
	FE_PENTA,		// = 4

	// 2.5D surface elements
	FE_QUAD,		// = 5
	FE_NIQUAD,		// = 6
	FE_TRI,			// = 7
	FE_NITRI,		// = 8

	// shell elements
	FE_SHELL_QUAD,	// = 9
	FE_SHELL_TRI	// = 10
};

/////////////////////////////////////////////////////////////////////////////
// ENUM: Fiber distribution generation functions
//  Each element can be associated with a fiber direction which is generated
//  from a generator
//

enum FE_Fiber_Type {
	FE_FIBER_LOCAL		= 0,
	FE_FIBER_SPHERICAL	= 1,
	FE_FIBER_VECTOR		= 2,
	FE_FIBER_RANDOM2D	= 3,
	FE_FIBER_USER		= 4
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
	CG_ITERATIVE_SOLVER
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: Module types
//
enum FE_Module_Type {
	FE_SOLID,
	FE_POROELASTIC,
	FE_HEAT
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: Analysis types
//  Types of analysis that can be performed with FEBio.
//
enum FE_Analysis_Type {
	FE_STATIC		= 0,
	FE_DYNAMIC		= 1,
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: Rigid node types
//  Designates the type of displacement dofs
//  free      = node belongs to a deformable mesh
//  rigid     = node belongs to a rigid body
//  interface = node is on the face between deformable and rigid mesh
//

enum FE_Node_Type {
	FE_NODE_FREE,
	FE_NODE_RIGID,
	FE_NODE_INTERFACE
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: Unpack flags

enum FE_Unpack_Flags {
	FE_UNPACK_JAC0		= 1,
	FE_UNPACK_JACT		= 2,
	FE_UNPACK_DEFGRAD	= 4,
	FE_UNPACK_LM		= 8,
	FE_UNPACK_ALL		= 15
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: Contact interface types

enum FE_Contact_Types {
	FE_CONTACT_SLIDING		= 1,
	FE_CONTACT_RIGIDWALL	= 2,
	FE_CONTACT_TIED			= 3,
	FE_FACET2FACET_SLIDING	= 4,
	FE_CONTACT_SLIDING2		= 5,
	FE_PERIODIC_BOUNDARY	= 6,
	FE_SURFACE_CONSTRAINT	= 7
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: rigid surfaces

enum FE_Rigid_Surface_Type {
	FE_RIGID_PLANE,
	FE_RIGID_SPHERE
};

///////////////////////////////////////////////////////////////////////////////

enum FE_Plot_Level {
	FE_PLOT_NEVER,
	FE_PLOT_MAJOR_ITRS,
	FE_PLOT_MINOR_ITRS,
	FE_PLOT_MUST_POINTS
};

///////////////////////////////////////////////////////////////////////////////

enum FE_Print_Level {
	FE_PRINT_NEVER,
	FE_PRINT_PROGRESS,
	FE_PRINT_MAJOR_ITRS,
	FE_PRINT_MINOR_ITRS,
	FE_PRINT_MINOR_ITRS_EXP,
};

#endif // _FE_ENUM_H_05132007_
