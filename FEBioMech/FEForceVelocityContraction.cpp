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



#include "stdafx.h"
#include "FEForceVelocityContraction.h"
#include "FEElasticMaterial.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
FEForceVelocityMaterialPoint::FEForceVelocityMaterialPoint()
{

}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEForceVelocityMaterialPoint::Copy()
{
    FEForceVelocityMaterialPoint* pt = new FEForceVelocityMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEForceVelocityMaterialPoint::Init()
{
	FEMaterialPointData::Init();
    
    m_lambdap = 1;
    
    for (int i=0; i<MAX_TERMS; ++i) {
        m_H[i] = 0;
        m_Hp[i] = 0;
    };
}

//-----------------------------------------------------------------------------
void FEForceVelocityMaterialPoint::Serialize(DumpStream& ar)
{
    
    if (ar.IsSaving())
    {
        ar << m_lambdap;
        for (int i=0; i<MAX_TERMS; ++i) ar << m_H[i] << m_Hp[i];
    }
    else
    {
        ar >> m_lambdap;
        for (int i=0; i<MAX_TERMS; ++i) ar >> m_H[i] >> m_Hp[i];
    }
    
	FEMaterialPointData::Serialize(ar);
}

//////////////////////////////////////////////////////////////////////
// FEForceVelocityContraction
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEForceVelocityContraction, FEActiveContractionMaterial);
    ADD_PARAMETER(m_ascl , "ascl");
    ADD_PARAMETER(m_Tmax , "Tmax");
    ADD_PARAMETER(m_ca0  , "ca0");
    ADD_PARAMETER(m_camax, "camax");
    ADD_PARAMETER(m_beta , "beta");
    ADD_PARAMETER(m_l0   , "l0");
    ADD_PARAMETER(m_refl , "refl");
    ADD_PARAMETER(m_alpha[0] , "alpha1");
    ADD_PARAMETER(m_alpha[1] , "alpha2");
    ADD_PARAMETER(m_alpha[2] , "alpha3");
    ADD_PARAMETER(m_A[0] , "A1");
    ADD_PARAMETER(m_A[1] , "A2");
    ADD_PARAMETER(m_A[2] , "A3");
    ADD_PARAMETER(m_at   , "a_t");
    ADD_PARAMETER(m_bfvel, "force_velocity" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEForceVelocityContraction::FEForceVelocityContraction(FEModel* pfem) : FEActiveContractionMaterial(pfem)
{
    m_ascl = 0;
    m_Tmax = 1.0;
    m_ca0 = 1.0;
    m_camax = 0.0;
    m_l0 = m_refl = 0;
    m_beta = 0;
    m_A[0] = m_A[1] = m_A[2] = 0;
    m_alpha[0] = m_alpha[1] = m_alpha[2] = 0;
    m_at = 0;
    m_bfvel = true;
}

//-----------------------------------------------------------------------------
bool FEForceVelocityContraction::Init()
{
    if (FEActiveContractionMaterial::Init() == false) return false;
    
    // for backward compatibility we set m_camax to m_ca0 if it is not defined
    if (m_camax == 0.0) m_camax = m_ca0;
    if (m_camax <= 0.0) { feLogError("camax must be larger than zero"); return false; }
    
    return true;
}

//-----------------------------------------------------------------------------
mat3ds FEForceVelocityContraction::ActiveStress(FEMaterialPoint& mp, const vec3d& a0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the deformation gradient
    mat3d F = pt.m_F;
    double J = pt.m_J;
    double Jm13 = pow(J, -1.0 / 3.0);
    
    // calculate the current material axis lam*a = F*a0;
    vec3d a = F*a0;
    
    double lam, lamd;
    lam = a.unit();
    lamd = lam*Jm13; // i.e. lambda tilde
    
    mat3ds AxA = dyad(a);
    
    // get the activation
    double saf = 0.0;
    double FVstress = 0.0;
    double dt = GetTimeInfo().timeIncrement;
    if (m_ascl > 0)
    {
        double ctenslm = m_ascl;
        
        // current sarcomere length
        double strl = m_refl*lamd;
        
        // sarcomere length change
        double dl = strl - m_l0;
        
        if (dl >= 0)
        {
            // calcium sensitivity
            double eca50i = (exp(m_beta*dl) - 1);
            
            // ratio of Camax/Ca0
            double rca = m_camax/m_ca0;
            
            // active fiber stress
            saf = m_Tmax*(eca50i / ( eca50i + rca*rca ))*ctenslm;
            
            if (m_bfvel)
            {
                FEForceVelocityMaterialPoint& vpt = *mp.ExtractData<FEForceVelocityMaterialPoint>();
                
                double lambdap = vpt.m_lambdap;
                
                double Qac = 0;
                double Qac2 = 0;
                
                // double dt = GetFEModel()->GetTime().timeIncrement;
                // double dt = 1;
                double g;
                double h;
                // double Hnew;
                // double Hnew = 10;
                
                for (int i=0; i<MAX_TERMS; ++i)
                {
                    g = exp(-dt/m_alpha[i]);
                    h = (1 - g)/(dt/m_alpha[i]);
                    
                    if ((lamd - lambdap) <= 0)
                    {
                        vpt.m_H[i] = vpt.m_Hp[i]*g + (lamd - lambdap)*h;
                    }
                    
                    if ((lamd - lambdap) == 0)
                    {
                        vpt.m_H[i] = vpt.m_Hp[i]*g + (lamd - lambdap)*h;
                    }
                    
                    if ((lamd - lambdap) > 0)
                    {
                        vpt.m_H[i] = vpt.m_Hp[i]*g; // POSITIVE VELOCITY SET TO 0
                    }

                    Qac += vpt.m_H[i]*m_A[i];
                    Qac2 = abs(Qac);
                }
                
                // return the total Cauchy stress,
                // which is the push-forward of sarcomere
                
                // FVstress = saf*(1 - m_at*Qac2)/(1+Qac2); //NEGATIVE ABSOLUTE VALUE OF Q
                FVstress = saf*(1 + m_at*Qac)/(1-Qac);
            }
            else
            {
                FVstress = saf;
            }
        }
    }
    return AxA*FVstress;

}

//-----------------------------------------------------------------------------
tens4ds FEForceVelocityContraction::ActiveStiffness(FEMaterialPoint& mp, const vec3d& a0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // get the deformation gradient
    mat3d F = pt.m_F;
    double J = pt.m_J;
    double Jm13 = pow(J, -1.0 / 3.0);
    
    // calculate the current material axis lam*a = F*a0;
    vec3d a = F*a0;
    
    // normalize material axis and store fiber stretch
    double lam, lamd;
    lam = a.unit();
    lamd = lam*Jm13; // i.e. lambda tilde
    
    // calculate dyad of a: AxA = (a x a)
    mat3ds AxA = dyad(a);
    tens4ds AxAxAxA = dyad1s(AxA);
    
    double c = 0;
    if (m_ascl > 0)
    {
        // current sarcomere length
        double strl = m_refl*lamd;
        
        // sarcomere length change
        double dl = strl - m_l0;
        
        if (dl >= 0)
        {
            // calcium sensitivity
            double eca50i = (exp(m_beta*dl) - 1);
            
            double decl = m_beta*m_refl*exp(m_beta*dl);
            
            // ratio of Camax/Ca0
            double rca = m_camax / m_ca0;
            
            double d = eca50i + rca*rca;
            
            // active fiber stress
            double saf = m_Tmax*(eca50i / d)*m_ascl;
            
            double dsf = m_Tmax*m_ascl*decl*(1.0/ d - eca50i / (d*d));
            
            c = (lamd*dsf - 2.0*saf);
        }
    }
    
    return AxAxAxA*c;
}

//-----------------------------------------------------------------------------
// update force-velocity material point
void FEForceVelocityContraction::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp, const vec3d& a0)
{
    FEForceVelocityMaterialPoint& pt = *mp.ExtractData<FEForceVelocityMaterialPoint>();
    FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
    double Jm13 = pow(pe.m_J, -1.0/3.0);
    
    vec3d a = pe.m_F*a0;
    double lambda1 = a.unit();
    pt.m_lambdap = lambda1*Jm13;
    
    for (int i=0; i<FEForceVelocityMaterialPoint::MAX_TERMS; ++i) {
        pt.m_Hp[i] = pt.m_H[i];
    }
}
