// PlotFile.cpp: implementation of the PlotFile class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "PlotFile.h"
#include "fem.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PlotFile::PlotFile()
{
}

PlotFile::~PlotFile()
{
	Close();
}

//-----------------------------------------------------------------------------
void PlotFile::Close()
{
	m_ar.Close();
}
