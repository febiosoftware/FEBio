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
