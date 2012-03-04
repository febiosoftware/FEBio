#pragma once

#include <Flx_Dialog.h>

class CDlgEditFind : public Flx_Dialog
{
public:
	CDlgEditFind();

public:
	bool	m_bcase;
	char	m_sztxt[256];
};
