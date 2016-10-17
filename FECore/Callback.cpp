#include "stdafx.h"
#include "Callback.h"

//-----------------------------------------------------------------------------
CallbackHandler::CallbackHandler()
{
}

//-----------------------------------------------------------------------------
CallbackHandler::~CallbackHandler()
{
}

//-----------------------------------------------------------------------------
//! This function adds a callback routine
//
void CallbackHandler::AddCallback(FECORE_CB_FNC pcb, unsigned int nwhen, void *pd)
{
	FECORE_CALLBACK cb;
	cb.m_pcb = pcb;
	cb.m_pd = pd;
	cb.m_nwhen = nwhen;

	m_pcb.push_back(cb);
}

//-----------------------------------------------------------------------------
//! Call the callback functions.
//
bool CallbackHandler::DoCallback(FEModel* fem, unsigned int nevent)
{
	std::list<FECORE_CALLBACK>::iterator it = m_pcb.begin();
	for (int i = 0; i<(int)m_pcb.size(); ++i, ++it)
	{
		// call the callback function
		if (it->m_nwhen & nevent)
		{
			bool bret = (it->m_pcb)(fem, nevent, it->m_pd);
			if (bret == false) return false;
		}
	}

	return true;
}
