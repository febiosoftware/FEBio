#pragma once

#include <Flx_Dialog.h>

class CDlgEditFind : public Flx_Dialog
{
public:
	CDlgEditFind();

	int InitDialog();

public:
	bool	m_bcase;
	char	m_sztxt[256];
};

class CDlgEditGoToLine : public Flx_Dialog
{
public:
	CDlgEditGoToLine();

public:
	int	m_nline;
};
