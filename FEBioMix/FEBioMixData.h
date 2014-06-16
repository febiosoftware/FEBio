#pragma once
#include "FECore/DataStore.h"
#include "FECore/NodeDataRecord.h"
#include "FECore/ElementDataRecord.h"

//=============================================================================
// N O D E  D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FENodePressure : public FENodeLogData
{ 
public: 
	FENodePressure(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeConcentration : public FENodeLogData
{ 
public: 
	FENodeConcentration(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeConcentration_ : public FENodeLogData
{ 
protected: 
	FENodeConcentration_(FEModel* pfem, int nsol) : FENodeLogData(pfem), m_nsol(nsol) {} 
	double value(int node); 
private:
	int	m_nsol;
};

//-----------------------------------------------------------------------------
template <int N> class FENodeConcentration_T : public FENodeConcentration_
{ 
public: 
	FENodeConcentration_T(FEModel* pfem) : FENodeConcentration_(pfem, N) {} 
	double value(int node) { return FENodeConcentration_::value(node); }
};

//=============================================================================
// E L E M E N T   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FELogElemFluidPressure : public FELogElemData
{
public:
	FELogElemFluidPressure(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemFluidFluxX : public FELogElemData
{
public:
	FELogElemFluidFluxX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemFluidFluxY : public FELogElemData
{
public:
	FELogElemFluidFluxY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemFluidFluxZ : public FELogElemData
{
public:
	FELogElemFluidFluxZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemSoluteConcentration : public FELogElemData
{
public:
	FELogElemSoluteConcentration(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemSoluteFluxX : public FELogElemData
{
public:
	FELogElemSoluteFluxX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemSoluteFluxY : public FELogElemData
{
public:
	FELogElemSoluteFluxY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemSoluteFluxZ : public FELogElemData
{
public:
	FELogElemSoluteFluxZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemSoluteRefConcentration : public FELogElemData
{
public:
	FELogElemSoluteRefConcentration(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemSoluteConcentration_ : public FELogElemData
{
protected:
	FELogElemSoluteConcentration_(FEModel* pfem, int nsol) : FELogElemData(pfem), m_nsol(nsol) {}
	double value(FEElement& el);
private:
	int m_nsol;
};

template <int N> class FELogElemSoluteConcentration_T : public FELogElemSoluteConcentration_
{
public:
	FELogElemSoluteConcentration_T(FEModel* pfem) : FELogElemSoluteConcentration_(pfem, N) {}
	double value(FEElement& el) { return FELogElemSoluteConcentration_::value(el); }
};

//-----------------------------------------------------------------------------
class FELogElemSoluteFluxX_ : public FELogElemData
{
protected:
	FELogElemSoluteFluxX_(FEModel* pfem, int nsol) : FELogElemData(pfem), m_nsol(nsol) {}
	double value(FEElement& el);
private:
	int	m_nsol;
};

template <int N> class FELogElemSoluteFluxX_T : public FELogElemSoluteFluxX_
{
public:
	FELogElemSoluteFluxX_T(FEModel* pfem) : FELogElemSoluteFluxX_(pfem, N) {}
	double value(FEElement& el) { return FELogElemSoluteFluxX_::value(el); }
};

//-----------------------------------------------------------------------------
class FELogElemSoluteFluxY_ : public FELogElemData
{
protected:
	FELogElemSoluteFluxY_(FEModel* pfem, int nsol) : FELogElemData(pfem), m_nsol(nsol) {}
	double value(FEElement& el);
private:
	int	m_nsol;
};

template <int N> class FELogElemSoluteFluxY_T : public FELogElemSoluteFluxY_
{
public:
	FELogElemSoluteFluxY_T(FEModel* pfem) : FELogElemSoluteFluxY_(pfem, N) {}
	double value(FEElement& el) { return FELogElemSoluteFluxY_::value(el); }
};

//-----------------------------------------------------------------------------
class FELogElemSoluteFluxZ_ : public FELogElemData
{
protected:
	FELogElemSoluteFluxZ_(FEModel* pfem, int nsol) : FELogElemData(pfem), m_nsol(nsol) {}
	double value(FEElement& el);
private:
	int	m_nsol;
};

template <int N> class FELogElemSoluteFluxZ_T : public FELogElemSoluteFluxZ_
{
public:
	FELogElemSoluteFluxZ_T(FEModel* pfem) : FELogElemSoluteFluxZ_(pfem, N) {}
	double value(FEElement& el) { return FELogElemSoluteFluxZ_::value(el); }
};

//-----------------------------------------------------------------------------
class FELogElemElectricPotential : public FELogElemData
{
public:
	FELogElemElectricPotential(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemCurrentDensityX : public FELogElemData
{
public:
	FELogElemCurrentDensityX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemCurrentDensityY : public FELogElemData
{
public:
	FELogElemCurrentDensityY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemCurrentDensityZ : public FELogElemData
{
public:
	FELogElemCurrentDensityZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemSBMConcentration_: public FELogElemData
{
protected:
	FELogElemSBMConcentration_(FEModel* pfem, int nsol) : FELogElemData(pfem), m_nsol(nsol) {}
	double value(FEElement& el);
private:
	int	m_nsol;
};

template <int N> class FELogElemSBMConcentration_T: public FELogElemSBMConcentration_
{
public:
	FELogElemSBMConcentration_T(FEModel* pfem) : FELogElemSBMConcentration_(pfem, N) {}
	double value(FEElement& el) { return FELogElemSBMConcentration_::value(el); }
};
