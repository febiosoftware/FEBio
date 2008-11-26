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
	FE_TET,			// = 2
	FE_PENTA,		// = 3

	// 2.5D surface elements
	FE_QUAD,		// = 4
	FE_NIQUAD,		// = 5
	FE_TRI,			// = 6
	FE_NITRI,		// = 7

	// shell elements
	FE_SHELL_QUAD,	// = 8
	FE_SHELL_TRI	// = 9
};

/////////////////////////////////////////////////////////////////////////////
// ENUM: Material models
//  Numbers are chosen to be consistent with NIKE3D's material numbers
//  The materials that are not supported by NIKE3D have material
//  numbers >= 100.
//

enum FE_Material_Type {
	FE_NEOHOOKEAN				= 1,
	FE_MOONEY_RIVLIN			= 15,
	FE_TISO_MOONEY_RIVLIN		= 18,
	FE_RIGID					= 20,
	FE_OGDEN_MATERIAL			= 63,
	FE_LINEAR_ELASTIC			= 100,
	FE_STVENANT_KIRCHHOFF		= 101,
	FE_INCOMP_NEOHOOKEAN		= 102,
	FE_PORO_ELASTIC				= 103,
	FE_LINEAR_ORTHOTROPIC		= 104,
	FE_VERONDA_WESTMANN			= 115,
	FE_TISO_VERONDA_WESTMANN    = 118,
	FE_TC_NONLINEAR_ORTHOTROPIC = 215,
	FE_MUSCLE_MATERIAL			= 300,
	FE_TENDON_MATERIAL			= 301,
	FE_TISO2D_MOONEY_RIVLIN		= 315,
	FE_RAND_FIBER_MATERIAL		= 316,
	FE_VISCO_ELASTIC			= 415,
	FE_NEOHOOKEAN_TRANSISO		= 501		// material added by Shawn Reese
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
	LU_SOLVER,
	WSMP_SOLVER,			// use only where available
	CG_ITERATIVE_SOLVER
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: Analysis types
//  Types of analysis that can be performed with FEBio. 
//

enum FE_Analysis_Type {
	FE_STATIC		= 0,
	FE_DYNAMIC		= 1,
	FE_STATIC_PORO	= 4
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
	FE_CONTACT_TIED			= 3
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
