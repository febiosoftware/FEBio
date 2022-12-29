//
//  FEPlasticFlowCurve.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/6/21.
//  Copyright © 2021 febio.org. All rights reserved.
//

#include "FEPlasticFlowCurve.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
//        F E P L A S T I C F L O W C U R V E M A T E R I A L P O I N T
//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPointData* FEPlasticFlowCurveMaterialPoint::Copy()
{
    FEPlasticFlowCurveMaterialPoint* pt = new FEPlasticFlowCurveMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FEPlasticFlowCurveMaterialPoint::Init()
{
    // don't forget to initialize the base class
    FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
void FEPlasticFlowCurveMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
    
    ar & m_Ky & m_w;
    ar & m_binit;
}

//-----------------------------------------------------------------------------
//              F E P L A S T I C F L O W C U R V E
//-----------------------------------------------------------------------------
vector<double> FEPlasticFlowCurve::BondYieldMeasures(FEMaterialPoint& mp)
{
    FEPlasticFlowCurveMaterialPoint& fp = *mp.ExtractData<FEPlasticFlowCurveMaterialPoint>();
    return fp.m_Ky;
}

//-----------------------------------------------------------------------------
vector<double> FEPlasticFlowCurve::BondMassFractions(FEMaterialPoint& mp)
{
    FEPlasticFlowCurveMaterialPoint& fp = *mp.ExtractData<FEPlasticFlowCurveMaterialPoint>();
    return fp.m_w;
}

//-----------------------------------------------------------------------------
size_t FEPlasticFlowCurve::BondFamilies(FEMaterialPoint& mp)
{
    FEPlasticFlowCurveMaterialPoint& fp = *mp.ExtractData<FEPlasticFlowCurveMaterialPoint>();
    return fp.m_Ky.size();
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPointData* FEPlasticFlowCurve::CreateMaterialPointData()
{
    return new FEPlasticFlowCurveMaterialPoint(new FEMaterialPointData());
}

//-----------------------------------------------------------------------------
//              F E P L A S T I C F L O W C U R V E P A P E R
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPlasticFlowCurvePaper, FEPlasticFlowCurve)
    ADD_PARAMETER(m_Ymin   , FE_RANGE_GREATER_OR_EQUAL(0.0), "Y0"  );
    ADD_PARAMETER(m_Ymax   , FE_RANGE_GREATER_OR_EQUAL(0.0), "Ymax");
    ADD_PARAMETER(m_wmin   , FE_RANGE_CLOSED(0.0, 1.0)     , "w0"  );
    ADD_PARAMETER(m_we     , FE_RANGE_CLOSED(0.0, 1.0)     , "we"  );
    ADD_PARAMETER(m_n      , FE_RANGE_GREATER(0)           , "nf"  );
    ADD_PARAMETER(m_bias   , FE_RANGE_LEFT_OPEN(0.0, 1.0)  , "r"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEPlasticFlowCurvePaper::FEPlasticFlowCurvePaper(FEModel* pfem) : FEPlasticFlowCurve(pfem)
{
    m_n = 1;
    m_wmin = 1;
    m_wmax = 1;
    m_we = 0;
    m_Ymin = 0;
    m_Ymax = 0;
    m_bias = 0.9;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEPlasticFlowCurvePaper::InitFlowCurve(FEMaterialPoint& mp)
{
    FEPlasticFlowCurveMaterialPoint& fp = *mp.ExtractData<FEPlasticFlowCurveMaterialPoint>();
    
    if (fp.m_binit == false) {
        double wmin = m_wmin(mp);
        double wmax = m_wmax(mp);
        double we = m_we(mp);
        double Ymin = m_Ymin(mp);
        double Ymax = m_Ymax(mp);
        double bias = m_bias(mp);
        wmax = 1 - we;
        if (wmax < wmin) {
            if (m_n ==1)
                feLogError("w0 + we = 1 must be satisfied");
            else
                feLogError("w0 + we < 1 must be satisfied");
            return false;
        }
        
        fp.m_Ky.assign(m_n,0);
        fp.m_w.assign(m_n+1,0);
        vector<double> Kp(m_n,0);
        
        if (m_n == 1) {
            fp.m_Ky[0] = Ymin;
            fp.m_w[0] = wmin;
        }
        else {
            // use bias r to reduce intervals in Ky and w as they increase proportionally
            double r = bias;
            // r= 1 uses uniform intervals
            if (r == 1) {
                fp.m_w[0] = wmin;
                Kp[0] = Ymin;
                fp.m_Ky[0] = Kp[0];
                double sw = fp.m_w[0];
                for (int i=1; i<m_n; ++i) {
                    fp.m_w[i] = (wmax - wmin)/(m_n-1);
                    Kp[i] = Ymin + (Ymax - Ymin)*i/(m_n-1);
                    fp.m_Ky[i] = fp.m_Ky[i-1] + (Kp[i]-Kp[i-1])/(1-sw);
                    sw += fp.m_w[i];
                }
            }
            else {
                double c = (1-r)/(1-pow(r, m_n-1));
                fp.m_w[0] = wmin;
                fp.m_w[1] = c*(wmax-wmin);
                Kp[0] = Ymin;
                Kp[1] = Kp[0] + c*(Ymax - Ymin);
                double sw = fp.m_w[0];
                fp.m_Ky[0] = Kp[0];
                fp.m_Ky[1] = fp.m_Ky[0] + (Kp[1]-Kp[0])/(1-sw);
                sw += fp.m_w[1];
                for (int i=2; i<m_n; ++i) {
                    fp.m_w[i] = fp.m_w[i-1]*r;
                    Kp[i] = Kp[i-1] + (Kp[i-1]-Kp[i-2])*r;
                    fp.m_Ky[i] = fp.m_Ky[i-1] + (Kp[i]-Kp[i-1])/(1-sw);
                    sw += fp.m_w[i];
                }
            }
        }
        fp.m_w[m_n] = we;
        
        fp.m_binit = true;
        
        return true;
    }

    return false;
}

//-----------------------------------------------------------------------------
//              F E P L A S T I C F L O W C U R V E U S E R
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPlasticFlowCurveUser, FEPlasticFlowCurve)
    ADD_PROPERTY(m_Y  , "plastic_response");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEPlasticFlowCurveUser::FEPlasticFlowCurveUser(FEModel* pfem) : FEPlasticFlowCurve(pfem)
{
    m_Y = nullptr;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEPlasticFlowCurveUser::InitFlowCurve(FEMaterialPoint& mp)
{
    FEPlasticFlowCurveMaterialPoint& fp = *mp.ExtractData<FEPlasticFlowCurveMaterialPoint>();
    
    if (fp.m_binit == false) {
        m_Y->Init();
        
        // get number of points on flow curve
        int m_n = m_Y->Points();
        
        fp.m_Ky.assign(m_n,0);
        fp.m_w.assign(m_n+1,0);
        vector<double> Kp(m_n,0);
        // get first point on assumption that it represents the yield point
        // and extract Young's modulus
        LOADPOINT p = m_Y->LoadPoint(0);
        double Ey = p.value/p.time;
        for (int i=0; i<m_n; ++i) {
            LOADPOINT p = m_Y->LoadPoint(i);
            Kp[i] = p.value;
            fp.m_Ky[i] = Ey*p.time;
        }
        
        double sw = 0;
        if (m_n == 1) {
            fp.m_w[0] = 1;
            fp.m_w[m_n] = 0;
        }
        else {
            for (int i=0; i<m_n-1; ++i) {
                fp.m_w[i] = 1 - (Kp[i+1]-Kp[i])/(fp.m_Ky[i+1]-fp.m_Ky[i]) - sw;
                sw += fp.m_w[i];
            }
            fp.m_w[m_n-1] = 1 - sw;
            fp.m_w[m_n] = 0;
        }
        fp.m_binit = true;
        return true;
    }

    return false;
}

//-----------------------------------------------------------------------------
//              F E P L A S T I C F L O W C U R V E M A T H
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPlasticFlowCurveMath, FEPlasticFlowCurve)
    ADD_PARAMETER(m_n      , FE_RANGE_GREATER(0)           , "nf"  );
    ADD_PARAMETER(m_emin   , FE_RANGE_GREATER(0)           , "e0"  );
    ADD_PARAMETER(m_emax   , FE_RANGE_GREATER(0)           , "emax");
    ADD_PARAMETER(m_Ymath, "plastic_response")->setLongName("plastic flow curve");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEPlasticFlowCurveMath::FEPlasticFlowCurveMath(FEModel* pfem) : FEPlasticFlowCurve(pfem)
{
    m_n = 1;
    m_emin = 0;
    m_emax = 1;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEPlasticFlowCurveMath::InitFlowCurve(FEMaterialPoint& mp)
{
    FEPlasticFlowCurveMaterialPoint& fp = *mp.ExtractData<FEPlasticFlowCurveMaterialPoint>();
    
    if (fp.m_binit == false) {

        FEMathFunction Y(GetFEModel());
        Y.SetMathString(m_Ymath);
        if (Y.Init() == false) return false;
        
        fp.m_Ky.assign(m_n,0);
        fp.m_w.assign(m_n+1,0);
        vector<double> Kp(m_n,0);
        if (m_n == 1) {
            fp.m_w[0] = 1;
            fp.m_w[m_n] = 0;
            Kp[0] = fp.m_Ky[0] = Y.value(m_emin);
        }
        else {
            // set uniform increments in Kp and find corresponding strains
            // then evaluate Ky at those strains
            Kp[0] = fp.m_Ky[0] = Y.value(m_emin);
            Kp[m_n-1] = Y.value(m_emax);
            // Extract Young's modulus
            double Ey = Kp[0]/m_emin;
            double dKp = (Kp[m_n-1] - Kp[0])/(m_n-1);
            double e = m_emin;
            vector<double> enat(m_n,0);
            enat[0] = m_emin;
            for (int i=1; i<m_n; ++i) {
                Kp[i] = Kp[0] + i*dKp;
                if (Y.invert(Kp[i], e) == false) return false;
                fp.m_Ky[i] = Ey*e;
                enat[i] = e;
            }
            // evaluate bond mass fractions
            fp.m_w[0] = (fp.m_Ky[1] - Kp[1])/(fp.m_Ky[1] - Kp[0]);
            double sw = fp.m_w[0];
            for (int i=2; i<m_n; ++i) {
                fp.m_w[i-1] = 1 - sw - (Kp[i] - Kp[i-1])/(fp.m_Ky[i] - fp.m_Ky[i-1]);
                sw += fp.m_w[i-1];
            }
            fp.m_w[m_n-1] = 1 - sw;
            fp.m_w[m_n] = 0;
        }
        
        fp.m_binit = true;
        return FEPlasticFlowCurve::Init();
    }
    return false;
}

