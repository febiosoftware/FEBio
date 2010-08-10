// FEElement.cpp: implementation of the FEElement class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElement.h"
#include <math.h>

//-----------------------------------------------------------------------------
FEElementState::FEElementState(const FEElementState& s)
{
	m_data.resize( s.m_data.size() );
	for (size_t i=0; i<m_data.size(); ++i) m_data[i] = s.m_data[i]->Copy();
}

FEElementState& FEElementState::operator = (const FEElementState& s)
{
	Clear();
	m_data.resize( s.m_data.size() );
	for (size_t i=0; i<m_data.size(); ++i) m_data[i] = s.m_data[i]->Copy();
	return (*this);
}

//-----------------------------------------------------------------------------
FESolidElement::FESolidElement(const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
	m_gid = el.m_gid;

	// copy solid element data
	m_eJ = el.m_eJ;
	m_ep = el.m_ep;
	m_Lk = el.m_Lk;
}

FESolidElement& FESolidElement::operator = (const FESolidElement& el)
{
	// set the traits of the element
	if (el.m_pT) { SetTraits(el.m_pT); m_State = el.m_State; }

	// copy base class data
	m_mat = el.m_mat;
	m_nrigid = el.m_nrigid;
	m_node = el.m_node;
	m_nID = el.m_nID;
	m_gid = el.m_gid;

	// copy solid element data
	m_eJ = el.m_eJ;
	m_ep = el.m_ep;
	m_Lk = el.m_Lk;

	return (*this);
}
