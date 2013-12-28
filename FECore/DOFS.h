#ifndef __FECore__DOFS__
#define __FECore__DOFS__

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
