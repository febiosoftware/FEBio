#pragma once
#include "MItem.h"

//-----------------------------------------------------------------------------
class MMatrix;

typedef MItem* (*FUNCMATPTR )(const MMatrix& A);
typedef MItem* (*FUNCMAT2PTR)(const MMatrix& A, const MMatrix& B);

//-----------------------------------------------------------------------------
// Matrix class
class MMatrix : public MItem
{
public:
	MMatrix();
	~MMatrix();
	void Create(int nrow, int ncol);

public:
	MItem* copy() const override;

	MItem* Item(int i, int j) { return m_d[i][j]; }
	const MItem* Item(int i, int j) const { return m_d[i][j]; }

	std::vector<MItem*>& operator [] (int n) { return m_d[n]; }
	const std::vector<MItem*>& operator [] (int n) const { return m_d[n]; }

	MItem* operator () (int i, int j) { return m_d[i][j]->copy(); }
	const MItem* operator () (int i, int j) const { return m_d[i][j]->copy(); }

	int rows   () const { return m_nrow; }
	int columns() const { return m_ncol; }

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
	MItem* copy() const override { return new MFuncMat(m_pf, m_name, m_pi->copy()); }
	const std::string& Name() { return m_name; }
	FUNCMATPTR	funcptr() const { return m_pf; }

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
	MItem* copy() const override { return new MFuncMat2(m_pf, m_name, m_pleft->copy(), m_pright->copy()); }
	const std::string& Name() { return m_name; }
	FUNCMAT2PTR	funcptr() const { return m_pf; }

protected:
	std::string	m_name;
	FUNCMAT2PTR	m_pf;
};

inline MMatrix*	mmatrix (MItem* pi) { return static_cast<MMatrix  *>(pi); }
inline MFuncMat*  mfncmat (MItem* pi) { return static_cast<MFuncMat *>(pi); }
inline MFuncMat2* mfncmat2(MItem* pi) { return static_cast<MFuncMat2*>(pi); }

inline const MMatrix*	mmatrix (const MItem* pi) { return static_cast<const MMatrix  *>(pi); }
inline const MFuncMat*  mfncmat (const MItem* pi) { return static_cast<const MFuncMat *>(pi); }
inline const MFuncMat2* mfncmat2(const MItem* pi) { return static_cast<const MFuncMat2*>(pi); }

inline const MMatrix*   mmatrix (const MITEM& i) { return static_cast<const MMatrix  *>(i.ItemPtr()); }
inline const MFuncMat*  mfncmat (const MITEM& i) { return static_cast<const MFuncMat *>(i.ItemPtr()); }
inline const MFuncMat2* mfncmat2(const MITEM& i) { return static_cast<const MFuncMat2*>(i.ItemPtr()); }

//-----------------------------------------------------------------------------
// Matrix operations
MMatrix* operator + (const MMatrix& A, const MMatrix& B);
MMatrix* operator - (const MMatrix& A, const MMatrix& B);
MMatrix* operator * (const MMatrix& A, const MITEM& n);
MMatrix* operator / (const MMatrix& A, const MITEM& n);
MMatrix* operator * (const MMatrix& A, const MMatrix& B);

MItem* matrix_transpose  (const MMatrix& A);
MItem* matrix_trace      (const MMatrix& A);
MItem* matrix_determinant(const MMatrix& A);
MItem* matrix_inverse    (const MMatrix& A);

MItem* matrix_contract(const MMatrix& A, const MMatrix& B);

MMatrix* matrix_identity(int n);
