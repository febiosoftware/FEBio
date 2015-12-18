#ifndef __FECore__DOFS__
#define __FECore__DOFS__
#include <vector>

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
	struct DOF_ITEM
	{
		const char*		sz;		//!< symbol of variable
		int				nsize;	//!< number of degrees of freedom
		int				ndof;	//!< start index of degree of freedom arrays
	};

public:
	// constructors
	DOFS();
	DOFS(const DOFS& dofs);
	DOFS& operator = (const DOFS& dofs);
    
	//! destructor
	~DOFS();

	//! Reset dofs
	void Reset();

public:
	//! return the total number of degrees of freedom
	int GetNDOFS() const { return m_maxdofs; }

	//! return a dof index from the dof symbol
	//! this function returns -1 if the symbol is not recognized
	int GetDOF(const char* sz, int n = 0);

	// Add a degree of freedom
	// returns >= 0 on success of -1 on failure (e.g. if the dof is already defined)
	int AddDOF(const char* sz, int n = 1);
 
	// get the symbol of a dof
	const char* GetDOFSymbol(int ndof);

	//! change the size of a DOF variable
	bool ChangeDOFSize(const char* sz, int n);

	//! return the size of a variable
	//! (returns -1 if the symbol is not defined)
	int GetDOFSize(const char* sz);

private:
	void Update();
    
private:
	std::vector<DOF_ITEM>		m_dof;		//!< array of dof symbols
    int							m_maxdofs;  //!< total number of nodal DOFS (i.e. size of FENode::m_val array)
};
#endif /* defined(__FECore__DOFS__) */
