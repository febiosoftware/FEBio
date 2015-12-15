#ifndef __FECore__DOFS__
#define __FECore__DOFS__

//-----------------------------------------------------------------------------
// Degree of freedom types
#define DOF_OPEN		 0		// the dof is open and will be given an equation number
#define DOF_FIXED		-1		// the dof is fixed and will not be given an equation number
#define DOF_PRESCRIBED	 2		// the dof is prescribed. It will be given a negative equation number (equation = - index - 2)

//-----------------------------------------------------------------------------
// status of a dof
#define DOF_INACTIVE	-1		// the dof is inactive and should not be assigned an equation (regardless of its type)
#define DOF_ACTIVE		 0		// the dof is active and an equation can be assigned (depending on its type)

//-----------------------------------------------------------------------------
// Max nr of nodal degrees of freedom

// At this point the nodal dofs are used as follows:
//
#define DOF_X			0		// x-displacement
#define DOF_Y			1		// y-displacement
#define DOF_Z			2		// z-displacement
#define DOF_U			3		// x-rotation
#define DOF_V			4		// y-rotation
#define DOF_W			5		// z-rotation
#define DOF_P			6		// fluid pressure
#define DOF_RU			7		// rigid x-rotation
#define DOF_RV			8		// rigid y-rotation
#define DOF_RW			9		// rigid z-rotation
#define DOF_T			10		// temperature
#define DOF_VX			11		// x-fluid velocity
#define DOF_VY			12		// y-fluid velocity
#define DOF_VZ			13		// z-fluid velocity
#define DOF_E           14      // fluid dilatation
#define DOF_C			15		// solute concentration
//
// The rotational degrees of freedom are only used for rigid nodes and shells.
// The fluid pressure is only used for poroelastic problems.
// The rigid rotational degrees of freedom are only used for rigid nodes and only during the creation of the stiffenss matrix
// The temperature is only used during heat-conduction problems
// The solute concentration is only used in solute transport problems.

//-----------------------------------------------------------------------------
//! Class that is used for setting number of nodal degrees of freedom

//! Note that this class is implemented as a singleton, in other words, only one
//! instance can be created.

class DOFS
{
public:
	//! obtain a pointer to the DOFS
	static DOFS* GetInstance();
    
	//! destructor
	~DOFS();

	//! Reset dofs
	void Reset();

	//! return a dof index from the dof symbol
	//! this function returns -1 if the symbol is not recognized
	int GetDOF(const char* sz);
    
private:
	//! constructor is private so that you cannot create it directly
	DOFS();
	DOFS(const DOFS& dofs) {}
    
protected:
	static DOFS* m_pdofs;	//!< the one and only DOFS
    
private:
    int     MAX_NDOFS;      //!< total number of nodal DOFS
    int     MAX_CDOFS;      //!< number of solute DOFS

public:
    void    SetNDOFS(int ndofs);
    int     GetNDOFS();
    void    SetCDOFS(int cdofs);
    int     GetCDOFS();
};
#endif /* defined(__FECore__DOFS__) */
