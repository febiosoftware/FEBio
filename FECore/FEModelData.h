#pragma once
#include "FECoreBase.h"
#include "ElementDataRecord.h"


//-------------------------------------------------------------------
// Class that can used to query model data
class FECORE_API FEModelData : public FECoreBase
{
	FECORE_SUPER_CLASS

	struct ELEMREF
	{
		int	ndom;
		int	nid;
	};

public:
	// constructor
	FEModelData(FEModel* fem, FELogElemData* eval, vector<int>& items);

	// update model data (i.e. evaluate it based on current state of model)
	void Update();

private:
	FELogElemData*		m_eval;		//!< class that evaluates data
	std::vector<int>	m_item;		//!< item list
	std::vector<double>	m_data;		//!< model data values

protected:
	void BuildELT();
	vector<ELEMREF>	m_ELT;
	int	m_offset;

	DECLARE_FECORE_CLASS();
};
