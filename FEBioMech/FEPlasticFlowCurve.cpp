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
    m_wmin = m_wmax = 1;
    m_we = 0;
    m_Ymin = m_Ymax = 0;
    m_bias = 0.9;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEPlasticFlowCurvePaper::Init()
{
    if (m_binit == false) {
        m_wmax = 1 - m_we;
        if (m_wmax < m_wmin) {
            if (m_n ==1)
                feLogError("w0 + we = 1 must be satisfied");
            else
                feLogError("w0 + we < 1 must be satisfied");
            return false;
        }
        
        Ky.assign(m_n,0);
        w.assign(m_n+1,0);
        vector<double> Kp(m_n,0);
        
        if (m_n == 1) {
            Ky[0] = m_Ymin;
            w[0] = m_wmin;
        }
        else {
            // use bias r to reduce intervals in Ky and w as they increase proportionally
            double r = m_bias;
            // r= 1 uses uniform intervals
            if (r == 1) {
                w[0] = m_wmin;
                Kp[0] = m_Ymin;
                Ky[0] = Kp[0];
                double sw = w[0];
                for (int i=1; i<m_n; ++i) {
                    w[i] = (m_wmax - m_wmin)/(m_n-1);
                    Kp[i] = m_Ymin + (m_Ymax - m_Ymin)*i/(m_n-1);
                    Ky[i] = Ky[i-1] + (Kp[i]-Kp[i-1])/(1-sw);
                    sw += w[i];
                }
            }
            else {
                double c = (1-r)/(1-pow(r, m_n-1));
                w[0] = m_wmin;
                w[1] = c*(m_wmax-m_wmin);
                Kp[0] = m_Ymin;
                Kp[1] = Kp[0] + c*(m_Ymax - m_Ymin);
                double sw = w[0];
                Ky[0] = Kp[0];
                Ky[1] = Ky[0] + (Kp[1]-Kp[0])/(1-sw);
                sw += w[1];
                for (int i=2; i<m_n; ++i) {
                    w[i] = w[i-1]*r;
                    Kp[i] = Kp[i-1] + (Kp[i-1]-Kp[i-2])*r;
                    Ky[i] = Ky[i-1] + (Kp[i]-Kp[i-1])/(1-sw);
                    sw += w[i];
                }
            }
        }
        w[m_n] = m_we;
        
        m_binit = true;
    }

    return true;
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
bool FEPlasticFlowCurveUser::Init()
{
    if (m_binit == false) {
        m_Y->Init();
        
        // get number of points on flow curve
        int m_n = m_Y->Points();
        
        Ky.assign(m_n,0);
        w.assign(m_n+1,0);
        vector<double> Kp(m_n,0);
        // get first point on assumption that it represents the yield point
        // and extract Young's modulus
        LOADPOINT p = m_Y->LoadPoint(0);
        double Ey = p.value/p.time;
        for (int i=0; i<m_n; ++i) {
            LOADPOINT p = m_Y->LoadPoint(i);
            Kp[i] = p.value;
            Ky[i] = Ey*p.time;
        }
        
        double sw = 0;
        if (m_n == 1) {
            w[0] = 1;
            w[m_n] = 0;
        }
        else {
            for (int i=0; i<m_n-1; ++i) {
                w[i] = 1 - (Kp[i+1]-Kp[i])/(Ky[i+1]-Ky[i]) - sw;
                sw += w[i];
            }
            w[m_n] = 1 - sw;
        }
        m_binit = true;
    }

    return true;
}

//-----------------------------------------------------------------------------
//              F E P L A S T I C F L O W C U R V E M A T H
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPlasticFlowCurveMath, FEPlasticFlowCurve)
    ADD_PARAMETER(m_n      , FE_RANGE_GREATER(0)           , "nf"  );
    ADD_PARAMETER(m_emin   , FE_RANGE_GREATER(0)           , "emin");
    ADD_PARAMETER(m_emax   , FE_RANGE_GREATER(0)           , "emax");
    ADD_PROPERTY(m_Y  , "plastic_response");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEPlasticFlowCurveMath::FEPlasticFlowCurveMath(FEModel* pfem) : FEPlasticFlowCurve(pfem)
{
    m_n = 1;
    m_emin = 0;
    m_emax = 1;
    m_Y = nullptr;
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEPlasticFlowCurveMath::Init()
{
    if (m_binit == false) {
        m_Y->Init();
        
        Ky.assign(m_n,0);
        w.assign(m_n+1,0);
        vector<double> Kp(m_n,0);
        if (m_n == 1) {
            w[0] = 1;
            w[m_n] = 0;
            Ky[0] = m_Y->value(m_emin);
        }
        else {
            // set uniform increments in Kp and find corresponding strains
            // then evaluate Ky at those strains
            Kp[0] = Ky[0] = m_Y->value(m_emin);
            Kp[m_n-1] = m_Y->value(m_emax);
            // Extract Young's modulus
            double Ey = Kp[0]/m_emin;
            double dKp = (Kp[m_n-1] - Kp[0])/(m_n-1);
            double e = m_emin;
            vector<double> enat(m_n,0);
            enat[0] = m_emin;
            for (int i=1; i<m_n; ++i) {
                Kp[i] = Kp[0] + i*dKp;
                if (m_Y->invert(Kp[i], e) == false) return false;
                Ky[i] = Ey*e;
                enat[i] = e;
            }
            // evaluate bond mass fractions
            w[0] = (Ky[1] - Kp[1])/(Ky[1] - Kp[0]);
            double sw = w[0];
            for (int i=2; i<m_n; ++i) {
                w[i-1] = 1 - sw - (Kp[i] - Kp[i-1])/(Ky[i] - Ky[i-1]);
                sw += w[i-1];
            }
            w[m_n-1] = 1 - sw;
            w[m_n] = 0;
        }
        
        m_binit = true;
    }
    
    return true;
}

