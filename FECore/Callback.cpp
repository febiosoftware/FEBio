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
#include "Callback.h"

//-----------------------------------------------------------------------------
CallbackHandler::CallbackHandler()
{
	m_event = 0;
}

//-----------------------------------------------------------------------------
CallbackHandler::~CallbackHandler()
{
}

//-----------------------------------------------------------------------------
//! This function adds a callback routine
//
void CallbackHandler::AddCallback(FECORE_CB_FNC pcb, unsigned int nwhen, void *pd, CBInsertPolicy insert)
{
	FECORE_CALLBACK cb;
	cb.m_pcb = pcb;
	cb.m_pd = pd;
	cb.m_nwhen = nwhen;

	if (insert == CBInsertPolicy::CB_ADD_END)
		m_pcb.push_back(cb);
	else
		m_pcb.push_front(cb);
}

//-----------------------------------------------------------------------------
//! Call the callback functions.
//
bool CallbackHandler::DoCallback(FEModel* fem, unsigned int nevent)
{
	m_event = nevent;
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
	m_event = 0;

	return true;
}

//! Get the current callback reason (or zero if not inside DoCallback)
unsigned int CallbackHandler::CurrentEvent() const
{
	return m_event;
}
