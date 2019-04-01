#pragma once

class Interruption
{
public:
	Interruption();
	virtual ~Interruption();

	static void handler(int sig);
	static bool	m_bsig;
};
