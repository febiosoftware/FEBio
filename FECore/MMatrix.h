#pragma once
#include "MItem.h"

//-----------------------------------------------------------------------------
class MMatrix;

typedef MItem* (*FUNCMATPTR )(MMatrix& A);
typedef MItem* (*FUNCMAT2PTR)(MMatrix& A, MMatrix& B);

//-----------------------------------------------------------------------------
// Matrix class
class MMatrix : public MItem
{
public:
	MMatrix();
	~MMatrix();
	void Create(int nrow, int ncol);

public:
	MItem* copy();

	MItem* Item(int i, int j) { return m_d[i][j]; }
	std::vector<MItem*>& operator [] (int n) { return m_d[n]; }

	MItem* operator () (int i, int j) { return m_d[i][j]->copy(); }

	int rows   () { return m_nrow; }
	int columns() { return m_ncol; }

protected:
	std::vector< std::vector<MItem*>	>	m_d;
	int	m_nrow;
	int	m_ncol;
};

//-----------------------------------------------------------------------------
// function of matrices
class MFuncMat : public MUnary
{
public: 
	MFuncMat(FUNCMATPTR pf, const std::string& s, MItem* pi) : MUnary(pi, MFMAT), m_name(s), m_pf(pf) {}
	MItem* copy() { return new MFuncMat(m_pf, m_name, m_pi->copy()); }
	const std::string& Name() { return m_name; }
	FUNCMATPTR	funcptr() { return m_pf; }

protected:
	std::string	m_name;
	FUNCMATPTR	m_pf;
};

//-----------------------------------------------------------------------------
// function of 2 matrices
class MFuncMat2 : public MBinary
{
public: 
	MFuncMat2(FUNCMAT2PTR pf, const std::string& s, MItem* pl, MItem* pr) : MBinary(pl, pr, MFMAT2), m_name(s), m_pf(pf) {}
	MItem* copy() { return new MFuncMat2(m_pf, m_name, m_pleft->copy(), m_pright->copy()); }
	const std::string& Name() { return m_name; }
	FUNCMAT2PTR	funcptr() { return m_pf; }

protected:
	std::string	m_name;
	FUNCMAT2PTR	m_pf;
};

inline MMatrix*	  mmatrix (MItem* pi) { return static_cast<MMatrix  *>(pi); }
inline MFuncMat*  mfncmat (MItem* pi) { return static_cast<MFuncMat *>(pi); }
inline MFuncMat2* mfncmat2(MItem* pi) { return static_cast<MFuncMat2*>(pi); }

inline MMatrix*   mmatrix (MITEM& i) { return static_cast<MMatrix  *>(i.ItemPtr()); }
inline MFuncMat*  mfncmat (MITEM& i) { return static_cast<MFuncMat *>(i.ItemPtr()); }
inline MFuncMat2* mfncmat2(MITEM& i) { return static_cast<MFuncMat2*>(i.ItemPtr()); }

//-----------------------------------------------------------------------------
// Matrix operations
MMatrix* operator + (MMatrix& A, MMatrix& B);
MMatrix* operator - (MMatrix& A, MMatrix& B);
MMatrix* operator * (MMatrix& A, MITEM& n);
MMatrix* operator / (MMatrix& A, MITEM& n);
MMatrix* operator * (MMatrix& A, MMatrix& B);

MItem* matrix_transpose  (MMatrix& A);
MItem* matrix_trace      (MMatrix& A);
MItem* matrix_determinant(MMatrix& A);
MItem* matrix_inverse    (MMatrix& A);

MItem* matrix_contract(MMatrix& A, MMatrix& B);

MMatrix* matrix_identity(int n);
