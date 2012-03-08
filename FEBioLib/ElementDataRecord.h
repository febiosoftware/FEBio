#pragma once
#include "DataStore.h"

//-----------------------------------------------------------------------------
// TODO: should I create a different class for each data record? Like for the plot file?
class ElementDataRecord : public DataRecord
{
	enum {X, Y, Z, J, 
		EX, EY, EZ, EXY, EYZ, EXZ, 
		E1, E2, E3,
		SX, SY, SZ, SXY, SYZ, SXZ,
		S1, S2, S3,
		FX, FY, FZ, FYZ, FZX, FXY, FYX, FXZ, FZY, 
		P, WX, WY, WZ, C, JX, JY, JZ, CRC,
		C1, J1X, J1Y, J1Z, C2, J2X, J2Y, J2Z,
		C3, J3X, J3Y, J3Z, C4, J4X, J4Y, J4Z,
		C5, J5X, J5Y, J5Z, C6, J6X, J6Y, J6Z,
		PSI, IEX, IEY, IEZ
	};

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

protected:
	void BuildELT();

protected:
	vector<ELEMREF>	m_ELT;
};
