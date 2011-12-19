#include "stdafx.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
FEModel::FEModel(void)
{
}

//-----------------------------------------------------------------------------
FEModel::~FEModel(void)
{
}

//-----------------------------------------------------------------------------
// This function adds a callback routine
//
void FEModel::AddCallback(FEBIO_CB_FNC pcb, void *pd)
{
	FEBIO_CALLBACK cb;
	cb.m_pcb = pcb;
	cb.m_pd = pd;

	m_pcb.push_back(cb);
}

//-----------------------------------------------------------------------------
// Call the callback function if there is one defined
//
void FEModel::DoCallback()
{
	list<FEBIO_CALLBACK>::iterator it = m_pcb.begin();
	for (int i=0; i<(int) m_pcb.size(); ++i, ++it)
	{
		// call the callback function
		(it->m_pcb)(this, it->m_pd);
	}
}
