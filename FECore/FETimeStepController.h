#pragma once
#include "FEParameterList.h"
class FEAnalysis;
class DumpStream;

//-------------------------------------------------------------------
// Class to control the time step
class FECORE_API FETimeStepController : public FEParamContainer
{
public:
	FETimeStepController(FEAnalysis* step);

	// initialization
	bool Init();

	//! reset
	void Reset();

	//! serialize
	void Serialize(DumpStream& ar) override;

	//! copy from
	void CopyFrom(FETimeStepController* tc);

public:
	//! Do a running restart
	void Retry();

	//! Update Time step
	void AutoTimeStep(int niter);

	//! Adjust for must points
	double CheckMustPoints(double t, double dt);

private:
	FEAnalysis*	m_step;

public:
	int		m_nretries;		//!< nr of retries tried so far
	int		m_maxretries;	//!< max nr of retries allowed per time step
	int		m_naggr;		//!< aggressivness parameter
	int		m_nmplc;		//!< must point load curve number
	int		m_nmust;		//!< current must-point
	int		m_next_must;	//!< next must-point to visit
	int		m_iteopt;		//!< optimum nr of iterations
	double	m_dtmin;		//!< min time step size
	double	m_dtmax;		//!< max time step size

private:
	double	m_ddt;			//!< used by auto-time stepper
	double	m_dtp;			//!< previous time step size

	bool	m_dtforce;		//!< force max time step

	DECLARE_FECORE_CLASS();
};
