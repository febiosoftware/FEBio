#pragma once

#include <Flx_Dialog.h>
#include <FL/Fl_Input_Choice.H>

class CDlgEditFind : public Flx_Dialog
{
public:
	CDlgEditFind();

	int InitDialog();

protected:
	void OnOk(Fl_Widget* pw, void* pd);

public:
	bool	m_bcase;
	char	m_sztxt[256];
	Fl_Input_Choice* m_pinput;
};

class CDlgEditGoToLine : public Flx_Dialog
{
public:
	CDlgEditGoToLine();

public:
	int	m_nline;
};

class CDlgEditFilter : public Flx_Dialog
{
public:
	CDlgEditFilter();

	int InitDialog();

protected:
	void OnOk(Fl_Widget* pw, void* pd);

public:
	bool	m_bcase;
	char	m_sztxt[256];
	Fl_Input_Choice* m_pinput;
};
