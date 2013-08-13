#pragma once
#include "FECore/DataStore.h"

//-----------------------------------------------------------------------------
class FELogElemData
{
public:
	FELogElemData(FEModel* pfem) : m_pfem(pfem){}
	virtual ~FELogElemData(){}
	virtual double value(FEElement& el) = 0;
protected:
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
class FELogElemPosX : public FELogElemData
{
public:
	FELogElemPosX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemPosY : public FELogElemData
{
public:
	FELogElemPosY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemPosZ : public FELogElemData
{
public:
	FELogElemPosZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemJacobian : public FELogElemData
{
public:
	FELogElemJacobian(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStrainX : public FELogElemData
{
public:
	FELogElemStrainX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStrainY : public FELogElemData
{
public:
	FELogElemStrainY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStrainZ : public FELogElemData
{
public:
	FELogElemStrainZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStrainXY : public FELogElemData
{
public:
	FELogElemStrainXY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStrainYZ : public FELogElemData
{
public:
	FELogElemStrainYZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStrainXZ : public FELogElemData
{
public:
	FELogElemStrainXZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStrain1 : public FELogElemData
{
public:
	FELogElemStrain1(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStrain2 : public FELogElemData
{
public:
	FELogElemStrain2(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStrain3 : public FELogElemData
{
public:
	FELogElemStrain3(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStressX : public FELogElemData
{
public:
	FELogElemStressX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStressY : public FELogElemData
{
public:
	FELogElemStressY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStressZ : public FELogElemData
{
public:
	FELogElemStressZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStressXY : public FELogElemData
{
public:
	FELogElemStressXY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStressYZ : public FELogElemData
{
public:
	FELogElemStressYZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStressXZ : public FELogElemData
{
public:
	FELogElemStressXZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStress1 : public FELogElemData
{
public:
	FELogElemStress1(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStress2 : public FELogElemData
{
public:
	FELogElemStress2(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStress3 : public FELogElemData
{
public:
	FELogElemStress3(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemDeformationGradientXX : public FELogElemData
{
public:
	FELogElemDeformationGradientXX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemDeformationGradientXY : public FELogElemData
{
public:
	FELogElemDeformationGradientXY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemDeformationGradientXZ : public FELogElemData
{
public:
	FELogElemDeformationGradientXZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemDeformationGradientYX : public FELogElemData
{
public:
	FELogElemDeformationGradientYX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemDeformationGradientYY : public FELogElemData
{
public:
	FELogElemDeformationGradientYY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemDeformationGradientYZ : public FELogElemData
{
public:
	FELogElemDeformationGradientYZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemDeformationGradientZX : public FELogElemData
{
public:
	FELogElemDeformationGradientZX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemDeformationGradientZY : public FELogElemData
{
public:
	FELogElemDeformationGradientZY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemDeformationGradientZZ : public FELogElemData
{
public:
	FELogElemDeformationGradientZZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

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

//-----------------------------------------------------------------------------
class ElementDataRecord : public DataRecord
{
	struct ELEMREF
	{
		int	ndom;
		int	nid;
	};

public:
	ElementDataRecord(FEModel* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
	int Size() { return (int) m_Data.size(); }

protected:
	void BuildELT();

protected:
	vector<ELEMREF>	m_ELT;
	vector<FELogElemData*>	m_Data;
};
