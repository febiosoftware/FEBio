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
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
// Virtual base class for damage criterion

class FEDamageCriterion : public FEMaterial
{
public:
	FEDamageCriterion(FEModel* pfem) : FEMaterial(pfem) {}
    
	//! damage
	virtual double DamageCriterion(FEMaterialPoint& pt) = 0;
    
};

//-----------------------------------------------------------------------------
// Simo's damage criterion

class FEDamageCriterionSimo : public FEDamageCriterion
{
public:
	FEDamageCriterionSimo(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// Strain energy density as damage criterion

class FEDamageCriterionSED : public FEDamageCriterion
{
public:
	FEDamageCriterionSED(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// Specific strain energy as damage criterion

class FEDamageCriterionSSE : public FEDamageCriterion
{
public:
    FEDamageCriterionSSE(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage
    double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// von Mises stress as damage criterion

class FEDamageCriterionVMS : public FEDamageCriterion
{
public:
	FEDamageCriterionVMS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// max shear stress as damage criterion

class FEDamageCriterionMSS : public FEDamageCriterion
{
public:
	FEDamageCriterionMSS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// max normal stress as damage criterion

class FEDamageCriterionMNS : public FEDamageCriterion
{
public:
	FEDamageCriterionMNS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// max normal Lagrange strain as damage criterion

class FEDamageCriterionMNLS : public FEDamageCriterion
{
public:
	FEDamageCriterionMNLS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};
