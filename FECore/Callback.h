#pragma once
#include <list>
#include "fecore_export.h"

//-----------------------------------------------------------------------------
// Forward declarations.
class FEModel;

//-----------------------------------------------------------------------------
// callback events
#define CB_ALWAYS		0xFFFFFFFF		//!< Call for all reasons
#define CB_INIT			0x00000001		//!< Call after model initialization (i.e. FEModel::Init())
#define CB_STEP_ACTIVE	0x00000002		//!< call after step was activated (i.e. 
#define CB_MAJOR_ITERS	0x00000004		//!< Call at the end of each major converged iteration
#define CB_MINOR_ITERS	0x00000008		//!< Call for each minor iteration
#define CB_SOLVED		0x00000010		//!< Call at the end of FEModel::Solve
#define CB_UPDATE_TIME	0x00000020		//!< Call when time is updated and right before time step is solved (in FEAnalysis::Solve)
#define CB_AUGMENT		0x00000040		//!< The model is entering augmentations (called before Augment)
#define CB_STEP_SOLVED	0x00000080		//!< The step was solved

typedef unsigned int FECORE_CB_WHEN;
typedef bool(*FECORE_CB_FNC)(FEModel*, unsigned int, void*);

//-----------------------------------------------------------------------------
// callback structure
struct FECORE_CALLBACK {
	FECORE_CB_FNC	m_pcb;		// pointer to callback function
	void*			m_pd;		// pointer to user data
	FECORE_CB_WHEN	m_nwhen;	// when to call function
};

//-----------------------------------------------------------------------------
// class that handles callbacks
class FECOREDLL_EXPORT CallbackHandler
{
public:
	CallbackHandler();
	virtual ~CallbackHandler();

	//! set callback function
	void AddCallback(FECORE_CB_FNC pcb, unsigned int nwhen, void* pd);

	//! call the callback function
	//! This function returns false if the run is to be aborted
	bool DoCallback(FEModel* fem, unsigned int nevent);

private:
	std::list<FECORE_CALLBACK>	m_pcb;	//!< pointer to callback function
};
