// FENodeReorder.h: interface for the FENodeReorder class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FENODEREORDER_H__53098BEC_6B39_412C_B20B_14C73B5A5C95__INCLUDED_)
#define AFX_FENODEREORDER_H__53098BEC_6B39_412C_B20B_14C73B5A5C95__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FELevelStructure.h"

//-----------------------------------------------------------------------------
//! This class implements an algoritm that calculates a permutation of 
//! the node numbering in order to obtain a bandwidth reduced stiffness matrix

//! The algorithm comes from "An algorithm for reducing the bandwidth and 
//! profile of a sparse matrix", by N.E.Gibbs e.a. It applies the algorithm
//! on the node numberings in stead of the actual sparse matrix since that
//! was easier to implement :). In the future I would like to extend it to
//! work with the actual sparse matrix. 

class FECORE_API FENodeReorder
{

public:
	//! default constructor
	FENodeReorder();

	//! destructor
	virtual ~FENodeReorder();

	//! calculates the permutation vector
	void Apply(FEMesh& m, vector<int>& P);
};

#endif // !defined(AFX_FENODEREORDER_H__53098BEC_6B39_412C_B20B_14C73B5A5C95__INCLUDED_)
