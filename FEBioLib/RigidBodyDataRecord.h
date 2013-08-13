#pragma once
#include "FECore/DataStore.h"

class FERigidBody;

//-----------------------------------------------------------------------------
class FELogRigidBodyData
{
public:
	FELogRigidBodyData(FEModel* pfem) : m_pfem(pfem) {}
	virtual ~FELogRigidBodyData(){}
	virtual double value(FERigidBody& rb) = 0;
private:
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyPosX : public FELogRigidBodyData
{
public:
	FELogRigidBodyPosX(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyPosY : public FELogRigidBodyData
{
public:
	FELogRigidBodyPosY(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyPosZ : public FELogRigidBodyData
{
public:
	FELogRigidBodyPosZ(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyQuatX : public FELogRigidBodyData
{
public:
	FELogRigidBodyQuatX(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyQuatY : public FELogRigidBodyData
{
public:
	FELogRigidBodyQuatY(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyQuatZ : public FELogRigidBodyData
{
public:
	FELogRigidBodyQuatZ(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyQuatW : public FELogRigidBodyData
{
public:
	FELogRigidBodyQuatW(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyForceX : public FELogRigidBodyData
{
public:
	FELogRigidBodyForceX(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyForceY : public FELogRigidBodyData
{
public:
	FELogRigidBodyForceY(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyForceZ : public FELogRigidBodyData
{
public:
	FELogRigidBodyForceZ(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyTorqueX : public FELogRigidBodyData
{
public:
	FELogRigidBodyTorqueX(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyTorqueY : public FELogRigidBodyData
{
public:
	FELogRigidBodyTorqueY(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class FELogRigidBodyTorqueZ : public FELogRigidBodyData
{
public:
	FELogRigidBodyTorqueZ(FEModel* pfem) : FELogRigidBodyData(pfem){}
	double value(FERigidBody& rb);
};

//-----------------------------------------------------------------------------
class RigidBodyDataRecord : public DataRecord
{
public:
	RigidBodyDataRecord(FEModel* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
	int Size() { return (int) m_Data.size(); }

private:
	vector<FELogRigidBodyData*>	m_Data;
};
