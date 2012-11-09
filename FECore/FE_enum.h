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
	FE_HEX8G8,	
	FE_HEX8RI,
	FE_HEX8G1,	
	FE_TET4G1,	
	FE_TET4G4,		
	FE_PENTA6G6,	
	FE_TET10G4,
	FE_TET10G8,
	FE_HEX20G27,

	// 2.5D surface elements
	FE_QUAD4G4,
	FE_QUAD4NI,
	FE_TRI3G1,
	FE_TRI3G3,
	FE_TRI3NI,
	FE_TRI6G3,
	FE_TRI6G4,
	FE_TRI6G7,
	FE_TRI6GL7,
	FE_TRI6NI,
	FE_QUAD8G9,

	// shell elements
	FE_SHELL_QUAD,
	FE_SHELL_TRI,

	// truss elements
	FE_TRUSS,

	// discrete elements
	FE_DISCRETE,
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
	FE_FIBER_USER		= 3
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
	RCICG_SOLVER		// use only where available
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: Step types
//
enum FE_Step_Type {
	FE_SOLID,
	FE_BIPHASIC,
	FE_HEAT,
	FE_POROSOLUTE,
	FE_LINEAR_SOLID,
	FE_HEAT_SOLID,
	FE_EXPLICIT_SOLID,
};

///////////////////////////////////////////////////////////////////////////////
// ENUM: Analysis types
//  Types of analysis that can be performed with FEBio.
//
enum FE_Analysis_Type {
	FE_STATIC		= 0,
	FE_DYNAMIC		= 1,
	FE_STEADY_STATE	= 2
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
// ENUM: Contact interface types

enum FE_Contact_Types {
	FE_CONTACT_SLIDING		 = 1,
	FE_CONTACT_RIGIDWALL	 = 2,
	FE_CONTACT_TIED			 = 3,
	FE_FACET2FACET_SLIDING	 = 4,
	FE_CONTACT_SLIDING2		 = 5,
	FE_PERIODIC_BOUNDARY	 = 6,
	FE_SURFACE_CONSTRAINT	 = 7,
	FE_CONTACT_SLIDING3		 = 8,
	FE_CONTACT_TIED_BIPHASIC = 9,
	FE_CONTACT_SLIDINGBW	 = 10
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
	FE_PLOT_MUST_POINTS,
	FE_PLOT_FINAL
};

///////////////////////////////////////////////////////////////////////////////

enum FE_Print_Level {
	FE_PRINT_NEVER,
	FE_PRINT_PROGRESS,
	FE_PRINT_MAJOR_ITRS,
	FE_PRINT_MINOR_ITRS,
	FE_PRINT_MINOR_ITRS_EXP,
};

//-----------------------------------------------------------------------------
//! Domain Types
#define FE_SOLID_DOMAIN				1
#define FE_SHELL_DOMAIN				2
#define FE_SURFACE_DOMAIN			3
#define FE_TRUSS_DOMAIN				4
#define FE_RIGID_SOLID_DOMAIN		5
#define FE_RIGID_SHELL_DOMAIN		6
#define FE_UDGHEX_DOMAIN			7
#define FE_UT4_DOMAIN				8
#define FE_HEAT_SOLID_DOMAIN		9
#define FE_DISCRETE_DOMAIN			10
#define FE_3F_SOLID_DOMAIN			11
#define FE_BIPHASIC_DOMAIN			12
#define FE_BIPHASIC_SOLUTE_DOMAIN	13
#define FE_LINEAR_SOLID_DOMAIN		14
#define FE_TRIPHASIC_DOMAIN			15
#define FE_MULTIPHASIC_DOMAIN		16

//-----------------------------------------------------------------------------
//! surface load types
#define FE_PRESSURE_LOAD		1
#define FE_TRACTION_LOAD		2
#define FE_FLUID_FLUX			3
#define FE_PORO_TRACTION		4
#define FE_SOLUTE_FLUX			5
#define FE_HEAT_FLUX			6

//-----------------------------------------------------------------------------
//! body loads
#define FE_CONST_BODY_FORCE			1
#define FE_NONCONST_BODY_FORCE		2
#define FE_CENTRIFUGAL_BODY_FORCE	3

#endif // _FE_ENUM_H_05132007_
