/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

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
#include <FECore/NodeDataRecord.h>
#include <FECore/ElementDataRecord.h>
#include <FECore/FaceDataRecord.h>
#include "ObjectDataRecord.h"
#include <FECore/NLConstraintDataRecord.h>
#include <FECore/FENLConstraint.h>
#include <FECore/SurfaceDataRecord.h>

//=============================================================================
// N O D E  D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FENodeXPos : public FELogNodeData
{ 
public: 
	FENodeXPos(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeYPos : public FELogNodeData 
{ 
public: 
	FENodeYPos(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeZPos : public FELogNodeData
{ 
public: 
	FENodeZPos(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeXDisp : public FELogNodeData
{ 
public: 
	FENodeXDisp(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeYDisp : public FELogNodeData
{ 
public: 
	FENodeYDisp(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeZDisp : public FELogNodeData
{ 
public: 
	FENodeZDisp(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeXVel : public FELogNodeData
{ 
public: 
	FENodeXVel(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeYVel : public FELogNodeData
{ 
public: 
	FENodeYVel(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeZVel : public FELogNodeData
{ 
public: 
	FENodeZVel(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeXAcc : public FELogNodeData
{ 
public: 
	FENodeXAcc(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeYAcc : public FELogNodeData
{ 
public: 
	FENodeYAcc(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeZAcc : public FELogNodeData
{ 
public: 
	FENodeZAcc(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeForceX: public FELogNodeData
{ 
public: 
	FENodeForceX(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeForceY: public FELogNodeData
{ 
public: 
	FENodeForceY(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//-----------------------------------------------------------------------------
class FENodeForceZ: public FELogNodeData
{ 
public: 
	FENodeForceZ(FEModel* pfem) : FELogNodeData(pfem){} 
	double value(const FENode& node) override;
};

//=============================================================================
// F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
class FELogContactGap : public FELogFaceData
{
public:
	FELogContactGap(FEModel* fem) : FELogFaceData(fem) {}
	double value(FESurfaceElement& el) override;
};

//-----------------------------------------------------------------------------
class FELogContactPressure : public FELogFaceData
{
public:
	FELogContactPressure(FEModel* fem) : FELogFaceData(fem) {}
	double value(FESurfaceElement& el) override;
};

//-----------------------------------------------------------------------------
class FELogContactTractionX : public FELogFaceData
{
public:
	FELogContactTractionX(FEModel* fem) : FELogFaceData(fem) {}
	double value(FESurfaceElement& el) override;
};

//-----------------------------------------------------------------------------
class FELogContactTractionY : public FELogFaceData
{
public:
	FELogContactTractionY(FEModel* fem) : FELogFaceData(fem) {}
	double value(FESurfaceElement& el) override;
};

//-----------------------------------------------------------------------------
class FELogContactTractionZ : public FELogFaceData
{
public:
	FELogContactTractionZ(FEModel* fem) : FELogFaceData(fem) {}
	double value(FESurfaceElement& el) override;
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
class FELogElemStrainEffective : public FELogElemData
{
public:
	FELogElemStrainEffective(FEModel* pfem) : FELogElemData(pfem) {}
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
class FELogElemRightStretchX : public FELogElemData
{
public:
    FELogElemRightStretchX(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightStretchY : public FELogElemData
{
public:
    FELogElemRightStretchY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightStretchZ : public FELogElemData
{
public:
    FELogElemRightStretchZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightStretchXY : public FELogElemData
{
public:
    FELogElemRightStretchXY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightStretchYZ : public FELogElemData
{
public:
    FELogElemRightStretchYZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightStretchXZ : public FELogElemData
{
public:
    FELogElemRightStretchXZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightStretch1 : public FELogElemData
{
public:
    FELogElemRightStretch1(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightStretch2 : public FELogElemData
{
public:
    FELogElemRightStretch2(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightStretch3 : public FELogElemData
{
public:
    FELogElemRightStretch3(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightStretchEffective : public FELogElemData
{
public:
    FELogElemRightStretchEffective(FEModel* pfem) : FELogElemData(pfem) {}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftStretchX : public FELogElemData
{
public:
    FELogElemLeftStretchX(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftStretchY : public FELogElemData
{
public:
    FELogElemLeftStretchY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftStretchZ : public FELogElemData
{
public:
    FELogElemLeftStretchZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftStretchXY : public FELogElemData
{
public:
    FELogElemLeftStretchXY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftStretchYZ : public FELogElemData
{
public:
    FELogElemLeftStretchYZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftStretchXZ : public FELogElemData
{
public:
    FELogElemLeftStretchXZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftStretch1 : public FELogElemData
{
public:
    FELogElemLeftStretch1(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftStretch2 : public FELogElemData
{
public:
    FELogElemLeftStretch2(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftStretch3 : public FELogElemData
{
public:
    FELogElemLeftStretch3(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftStretchEffective : public FELogElemData
{
public:
    FELogElemLeftStretchEffective(FEModel* pfem) : FELogElemData(pfem) {}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightHenckyX : public FELogElemData
{
public:
    FELogElemRightHenckyX(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightHenckyY : public FELogElemData
{
public:
    FELogElemRightHenckyY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightHenckyZ : public FELogElemData
{
public:
    FELogElemRightHenckyZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightHenckyXY : public FELogElemData
{
public:
    FELogElemRightHenckyXY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightHenckyYZ : public FELogElemData
{
public:
    FELogElemRightHenckyYZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightHenckyXZ : public FELogElemData
{
public:
    FELogElemRightHenckyXZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightHencky1 : public FELogElemData
{
public:
    FELogElemRightHencky1(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightHencky2 : public FELogElemData
{
public:
    FELogElemRightHencky2(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightHencky3 : public FELogElemData
{
public:
    FELogElemRightHencky3(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemRightHenckyEffective : public FELogElemData
{
public:
    FELogElemRightHenckyEffective(FEModel* pfem) : FELogElemData(pfem) {}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftHenckyX : public FELogElemData
{
public:
    FELogElemLeftHenckyX(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftHenckyY : public FELogElemData
{
public:
    FELogElemLeftHenckyY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftHenckyZ : public FELogElemData
{
public:
    FELogElemLeftHenckyZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftHenckyXY : public FELogElemData
{
public:
    FELogElemLeftHenckyXY(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftHenckyYZ : public FELogElemData
{
public:
    FELogElemLeftHenckyYZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftHenckyXZ : public FELogElemData
{
public:
    FELogElemLeftHenckyXZ(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftHencky1 : public FELogElemData
{
public:
    FELogElemLeftHencky1(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftHencky2 : public FELogElemData
{
public:
    FELogElemLeftHencky2(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftHencky3 : public FELogElemData
{
public:
    FELogElemLeftHencky3(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemLeftHenckyEffective : public FELogElemData
{
public:
    FELogElemLeftHenckyEffective(FEModel* pfem) : FELogElemData(pfem) {}
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
class FELogElemStressEffective : public FELogElemData
{
public:
	FELogElemStressEffective(FEModel* pfem) : FELogElemData(pfem) {}
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
class FELogElemPK2StressX : public FELogElemData
{
public:
	FELogElemPK2StressX(FEModel* pfem) : FELogElemData(pfem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemPK2StressY : public FELogElemData
{
public:
	FELogElemPK2StressY(FEModel* pfem) : FELogElemData(pfem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemPK2StressZ : public FELogElemData
{
public:
	FELogElemPK2StressZ(FEModel* pfem) : FELogElemData(pfem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemPK2StressXY : public FELogElemData
{
public:
	FELogElemPK2StressXY(FEModel* pfem) : FELogElemData(pfem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemPK2StressYZ : public FELogElemData
{
public:
	FELogElemPK2StressYZ(FEModel* pfem) : FELogElemData(pfem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemPK2StressXZ : public FELogElemData
{
public:
	FELogElemPK2StressXZ(FEModel* pfem) : FELogElemData(pfem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElemStressEigenVector : public FELogElemData
{
public:
	FELogElemStressEigenVector(FEModel* pfem) : FELogElemData(pfem), m_eigenVector(-1), m_component(-1) {}
	double value(FEElement& el);

protected:
	int	m_eigenVector;
	int	m_component;
};

template <int eigenVector, int component> class FELogElemStressEigenVector_T : public FELogElemStressEigenVector
{
public:
	FELogElemStressEigenVector_T(FEModel* pfem) : FELogElemStressEigenVector(pfem)
	{
		m_eigenVector = eigenVector;
		m_component = component;
	}
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

//-----------------------------------------------------------------------------
//! Damage reduction factor
class FELogOctahedralPlasticStrain : public FELogElemData
{
public:
    FELogOctahedralPlasticStrain(FEModel* pfem) : FELogElemData(pfem){}
    double value(FEElement& el);
};

//-----------------------------------------------------------------------------
//! Discrete element stretch
class FELogDiscreteElementStretch : public FELogElemData
{
public:
	FELogDiscreteElementStretch(FEModel* fem) : FELogElemData(fem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
//! Discrete element elongation
class FELogDiscreteElementElongation : public FELogElemData
{
public:
	FELogDiscreteElementElongation(FEModel* fem) : FELogElemData(fem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
//! Discrete element force
class FELogDiscreteElementForce : public FELogElemData
{
public:
	FELogDiscreteElementForce(FEModel* fem) : FELogElemData(fem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
//! Discrete element force
class FELogDiscreteElementForceX : public FELogElemData
{
public:
	FELogDiscreteElementForceX(FEModel* fem) : FELogElemData(fem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
//! Discrete element force
class FELogDiscreteElementForceY : public FELogElemData
{
public:
	FELogDiscreteElementForceY(FEModel* fem) : FELogElemData(fem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
//! Discrete element force
class FELogDiscreteElementForceZ : public FELogElemData
{
public:
	FELogDiscreteElementForceZ(FEModel* fem) : FELogElemData(fem) {}
	double value(FEElement& el);
};

//-----------------------------------------------------------------------------
class FELogElementMixtureStress : public FELogElemData
{
public:
	FELogElementMixtureStress(FEModel* fem, int n, int m) : FELogElemData(fem), m_comp(n), m_metric(m) {}
	double value(FEElement& el);

private:
	int	m_comp;
	int	m_metric;
};

template <int N, int M> class FELogElementMixtureStress_T : public FELogElementMixtureStress
{
public:
	FELogElementMixtureStress_T(FEModel* fem) : FELogElementMixtureStress(fem, N, M) {}
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
class FELogRigidConnectorTranslationX : public FELogNLConstraintData
{
public:
    FELogRigidConnectorTranslationX(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorTranslationY : public FELogNLConstraintData
{
public:
    FELogRigidConnectorTranslationY(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorTranslationZ : public FELogNLConstraintData
{
public:
    FELogRigidConnectorTranslationZ(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorRotationX : public FELogNLConstraintData
{
public:
    FELogRigidConnectorRotationX(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorRotationY : public FELogNLConstraintData
{
public:
    FELogRigidConnectorRotationY(FEModel* pfem) : FELogNLConstraintData(pfem){}
    double value(FENLConstraint& rc);
};

//-----------------------------------------------------------------------------
class FELogRigidConnectorRotationZ : public FELogNLConstraintData
{
public:
    FELogRigidConnectorRotationZ(FEModel* pfem) : FELogNLConstraintData(pfem){}
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

//=============================================================================
// S U R F A C E   D A T A
//=============================================================================

class FELogContactArea : public FELogSurfaceData
{
public:
	FELogContactArea(FEModel* fem) : FELogSurfaceData(fem) {}
	double value(FESurface& surface) override;
};
