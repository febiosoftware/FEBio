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

#include "stdafx.h"
#include "FEElasticDomain.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEElasticDomain::FEElasticDomain(FEModel* pfem)
{
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");

    m_dofRU = pfem->GetDOFIndex("Ru");
    m_dofRV = pfem->GetDOFIndex("Rv");
    m_dofRW = pfem->GetDOFIndex("Rw");
    
    m_dofSX = pfem->GetDOFIndex("sx");
    m_dofSY = pfem->GetDOFIndex("sy");
    m_dofSZ = pfem->GetDOFIndex("sz");
    
    m_dofSXP = pfem->GetDOFIndex("sxp");
    m_dofSYP = pfem->GetDOFIndex("syp");
    m_dofSZP = pfem->GetDOFIndex("szp");
    
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
    
    m_dofSVX  = pfem->GetDOFIndex("svx");
    m_dofSVY  = pfem->GetDOFIndex("svy");
    m_dofSVZ  = pfem->GetDOFIndex("svz");
    m_dofSVXP = pfem->GetDOFIndex("svxp");
    m_dofSVYP = pfem->GetDOFIndex("svyp");
    m_dofSVZP = pfem->GetDOFIndex("svzp");
    
    m_dofSAX  = pfem->GetDOFIndex("sax");
    m_dofSAY  = pfem->GetDOFIndex("say");
    m_dofSAZ  = pfem->GetDOFIndex("saz");
    m_dofSAXP = pfem->GetDOFIndex("saxp");
    m_dofSAYP = pfem->GetDOFIndex("sayp");
    m_dofSAZP = pfem->GetDOFIndex("sazp");
}
