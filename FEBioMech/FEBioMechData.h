/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#pragma once
#include "FECore/NodeDataRecord.h"
#include "FECore/ElementDataRecord.h"
#include "ObjectDataRecord.h"
#include "FECore/NLConstraintDataRecord.h"
#include "FECore/FENLConstraint.h"

//=============================================================================
// N O D E  D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FENodeXPos : public FENodeLogData
{ 
public: 
	FENodeXPos(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeYPos : public FENodeLogData 
{ 
public: 
	FENodeYPos(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeZPos : public FENodeLogData
{ 
public: 
	FENodeZPos(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeXDisp : public FENodeLogData
{ 
public: 
	FENodeXDisp(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeYDisp : public FENodeLogData
{ 
public: 
	FENodeYDisp(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeZDisp : public FENodeLogData
{ 
public: 
	FENodeZDisp(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeXVel : public FENodeLogData
{ 
public: 
	FENodeXVel(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeYVel : public FENodeLogData
{ 
public: 
	FENodeYVel(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeZVel : public FENodeLogData
{ 
public: 
	FENodeZVel(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeXAcc : public FENodeLogData
{ 
public: 
	FENodeXAcc(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeYAcc : public FENodeLogData
{ 
public: 
	FENodeYAcc(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeZAcc : public FENodeLogData
{ 
public: 
	FENodeZAcc(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeForceX: public FENodeLogData
{ 
public: 
	FENodeForceX(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeForceY: public FENodeLogData
{ 
public: 
	FENodeForceY(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//-----------------------------------------------------------------------------
class FENodeForceZ: public FENodeLogData
{ 
public: 
	FENodeForceZ(FEModel* pfem) : FENodeLogData(pfem){} 
	double value(int node); 
};

//=============================================================================
// E L E M E N T   D A T A
//=============================================================================

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
class FELogElemInfStrainX : public FELogElemData
{
public:
	FELogElemInfStrainX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemInfStrainY : public FELogElemData
{
public:
	FELogElemInfStrainY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemInfStrainZ : public FELogElemData
{
public:
	FELogElemInfStrainZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemInfStrainXY : public FELogElemData
{
public:
	FELogElemInfStrainXY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemInfStrainYZ : public FELogElemData
{
public:
	FELogElemInfStrainYZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemInfStrainXZ : public FELogElemData
{
public:
	FELogElemInfStrainXZ(FEModel* pfem) : FELogElemData(pfem){}
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
// Base class for elasticity tensor output
class FELogElemElasticity_ : public FELogElemData
{
protected:
	FELogElemElasticity_(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el, int n);
};

//-----------------------------------------------------------------------------
// template class for generating the different tensor components
template <int n> class FELogElemElasticity : public FELogElemElasticity_
{
public:
	FELogElemElasticity(FEModel* pfem) : FELogElemElasticity_(pfem){}
	double value(FEElement& el) { return FELogElemElasticity_::value(el, n); }
};

//-----------------------------------------------------------------------------
class FELogElemStrainEnergyDensity : public FELogElemData
{
public:
	FELogElemStrainEnergyDensity(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemDevStrainEnergyDensity : public FELogElemData
{
public:
	FELogElemDevStrainEnergyDensity(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemFiberStretch : public FELogElemData
{
public:
	FELogElemFiberStretch(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemFiberVectorX : public FELogElemData
{
public:
	FELogElemFiberVectorX(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemFiberVectorY : public FELogElemData
{
public:
	FELogElemFiberVectorY(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemFiberVectorZ : public FELogElemData
{
public:
	FELogElemFiberVectorZ(FEModel* pfem) : FELogElemData(pfem){}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
//! Damage reduction factor
class FELogDamage : public FELogElemData
{
public:
    FELogDamage(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//=============================================================================
// R I G I D   B O D Y    D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FELogRigidBodyPosX : public FELogObjectData
{
public:
	FELogRigidBodyPosX(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyPosY : public FELogObjectData
{
public:
	FELogRigidBodyPosY(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyPosZ : public FELogObjectData
{
public:
	FELogRigidBodyPosZ(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyVelX : public FELogObjectData
{
public:
	FELogRigidBodyVelX(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyVelY : public FELogObjectData
{
public:
	FELogRigidBodyVelY(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyVelZ : public FELogObjectData
{
public:
	FELogRigidBodyVelZ(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAccX : public FELogObjectData
{
public:
	FELogRigidBodyAccX(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAccY : public FELogObjectData
{
public:
	FELogRigidBodyAccY(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAccZ : public FELogObjectData
{
public:
	FELogRigidBodyAccZ(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAngPosX : public FELogObjectData
{
public:
	FELogRigidBodyAngPosX(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAngPosY : public FELogObjectData
{
public:
	FELogRigidBodyAngPosY(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAngPosZ : public FELogObjectData
{
public:
	FELogRigidBodyAngPosZ(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAngVelX : public FELogObjectData
{
public:
	FELogRigidBodyAngVelX(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAngVelY : public FELogObjectData
{
public:
	FELogRigidBodyAngVelY(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAngVelZ : public FELogObjectData
{
public:
	FELogRigidBodyAngVelZ(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAngAccX : public FELogObjectData
{
public:
	FELogRigidBodyAngAccX(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAngAccY : public FELogObjectData
{
public:
	FELogRigidBodyAngAccY(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyAngAccZ : public FELogObjectData
{
public:
	FELogRigidBodyAngAccZ(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyQuatX : public FELogObjectData
{
public:
	FELogRigidBodyQuatX(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyQuatY : public FELogObjectData
{
public:
	FELogRigidBodyQuatY(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyQuatZ : public FELogObjectData
{
public:
	FELogRigidBodyQuatZ(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyQuatW : public FELogObjectData
{
public:
	FELogRigidBodyQuatW(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyR11 : public FELogObjectData
{
public:
	FELogRigidBodyR11(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyR12 : public FELogObjectData
{
public:
	FELogRigidBodyR12(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyR13 : public FELogObjectData
{
public:
	FELogRigidBodyR13(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyR21 : public FELogObjectData
{
public:
	FELogRigidBodyR21(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyR22 : public FELogObjectData
{
public:
	FELogRigidBodyR22(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyR23 : public FELogObjectData
{
public:
	FELogRigidBodyR23(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyR31 : public FELogObjectData
{
public:
	FELogRigidBodyR31(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyR32 : public FELogObjectData
{
public:
	FELogRigidBodyR32(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyR33 : public FELogObjectData
{
public:
	FELogRigidBodyR33(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyForceX : public FELogObjectData
{
public:
	FELogRigidBodyForceX(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyForceY : public FELogObjectData
{
public:
	FELogRigidBodyForceY(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyForceZ : public FELogObjectData
{
public:
	FELogRigidBodyForceZ(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyTorqueX : public FELogObjectData
{
public:
	FELogRigidBodyTorqueX(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyTorqueY : public FELogObjectData
{
public:
	FELogRigidBodyTorqueY(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyTorqueZ : public FELogObjectData
{
public:
	FELogRigidBodyTorqueZ(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidBodyKineticEnergy : public FELogObjectData
{
public:
    FELogRigidBodyKineticEnergy(FEModel* pfem) : FELogObjectData(pfem){}
	double value(FERigidBody& rb) override;
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorForceX : public FELogNLConstraintData
{
public:
    FELogRigidConnectorForceX(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorForceY : public FELogNLConstraintData
{
public:
    FELogRigidConnectorForceY(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorForceZ : public FELogNLConstraintData
{
public:
    FELogRigidConnectorForceZ(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorMomentX : public FELogNLConstraintData
{
public:
    FELogRigidConnectorMomentX(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorMomentY : public FELogNLConstraintData
{
public:
    FELogRigidConnectorMomentY(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorMomentZ : public FELogNLConstraintData
{
public:
    FELogRigidConnectorMomentZ(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogVolumeConstraint : public FELogNLConstraintData
{
public:
    FELogVolumeConstraint(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogVolumePressure : public FELogNLConstraintData
{
public:
    FELogVolumePressure(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

