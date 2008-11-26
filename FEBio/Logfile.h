// Logfile.h: interface for the Logfile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LOGFILE_H__18090874_EE74_42E8_AB1B_1874D975D646__INCLUDED_)
#define AFX_LOGFILE_H__18090874_EE74_42E8_AB1B_1874D975D646__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//-----------------------------------------------------------------------------
//! Class that is used for logging purposes

//! This class can output to different 
//! files at the same time. In FEBio, it is used for logging purposes.
//! At this time it outputs data to the screen (stdout) and to an external text file.

class Logfile  
{
public:
	enum MODE { FILE_ONLY = 1, SCREEN_ONLY, FILE_AND_SCREEN };

public:
	Logfile();
	virtual ~Logfile();

	bool open(const char* szfile);
	bool append(const  char* szfile);

	void printf(const char* sz, ...);
	void printbox(const char* sztitle, const char* sz, ...);

	MODE SetMode(MODE mode);

	MODE GetMode() {return m_mode; }

	void flush() { if (m_fp) fflush(m_fp); }

	const char* FileName() { return m_szfile; }

	operator FILE* () { return m_fp; }

protected:
	FILE*	m_fp;	// the actual log file

	MODE	m_mode;	// mode of log file

	char	m_szfile[256];	// file name of logfile
};

#endif // !defined(AFX_LOGFILE_H__18090874_EE74_42E8_AB1B_1874D975D646__INCLUDED_)
