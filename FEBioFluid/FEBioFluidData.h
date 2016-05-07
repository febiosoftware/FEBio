#pragma once
#include "FECore/NodeDataRecord.h"
#include "FECore/ElementDataRecord.h"

//=============================================================================
// N O D E  D A T A
//=============================================================================

//=============================================================================
// E L E M E N T   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FELogElemFluidPosX : public FELogElemData
{
public:
    FELogElemFluidPosX(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemFluidPosY : public FELogElemData
{
public:
    FELogElemFluidPosY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemFluidPosZ : public FELogElemData
{
public:
    FELogElemFluidPosZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElasticFluidPressure : public FELogElemData
{
public:
	FELogElasticFluidPressure(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidVolumeRatio : public FELogElemData
{
public:
	FELogFluidVolumeRatio(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidDensity : public FELogElemData
{
public:
    FELogFluidDensity(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidStressPower : public FELogElemData
{
public:
    FELogFluidStressPower(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidVelocityX : public FELogElemData
{
public:
    FELogFluidVelocityX(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidVelocityY : public FELogElemData
{
public:
    FELogFluidVelocityY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidVelocityZ : public FELogElemData
{
public:
    FELogFluidVelocityZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidAccelerationX : public FELogElemData
{
public:
	FELogFluidAccelerationX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidAccelerationY : public FELogElemData
{
public:
    FELogFluidAccelerationY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidAccelerationZ : public FELogElemData
{
public:
    FELogFluidAccelerationZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidVorticityX : public FELogElemData
{
public:
    FELogFluidVorticityX(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidVorticityY : public FELogElemData
{
public:
    FELogFluidVorticityY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidVorticityZ : public FELogElemData
{
public:
    FELogFluidVorticityZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidStressXX : public FELogElemData
{
public:
    FELogFluidStressXX(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidStressYY : public FELogElemData
{
public:
    FELogFluidStressYY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidStressZZ : public FELogElemData
{
public:
    FELogFluidStressZZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidStressXY : public FELogElemData
{
public:
    FELogFluidStressXY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidStressYZ : public FELogElemData
{
public:
    FELogFluidStressYZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidStressXZ : public FELogElemData
{
public:
    FELogFluidStressXZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidRateOfDefXX : public FELogElemData
{
public:
    FELogFluidRateOfDefXX(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidRateOfDefYY : public FELogElemData
{
public:
    FELogFluidRateOfDefYY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidRateOfDefZZ : public FELogElemData
{
public:
    FELogFluidRateOfDefZZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidRateOfDefXY : public FELogElemData
{
public:
    FELogFluidRateOfDefXY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidRateOfDefYZ : public FELogElemData
{
public:
    FELogFluidRateOfDefYZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogFluidRateOfDefXZ : public FELogElemData
{
public:
    FELogFluidRateOfDefXZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

