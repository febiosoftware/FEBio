#include "stdafx.h"
#include "FENewtonStrategy.h"

FENewtonStrategy::FENewtonStrategy()
{
	m_maxups = 10;
	m_max_buf_size = 0; // when zero, it should default to m_maxups
	m_cycle_buffer = true;
}

FENewtonStrategy::~FENewtonStrategy()
{
}
