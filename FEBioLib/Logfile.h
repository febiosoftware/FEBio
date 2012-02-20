// Logfile.h: interface for the Logfile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LOGFILE_H__18090874_EE74_42E8_AB1B_1874D975D646__INCLUDED_)
#define AFX_LOGFILE_H__18090874_EE74_42E8_AB1B_1874D975D646__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <stdio.h>

//-----------------------------------------------------------------------------
// class used to create an abstract interface to a screen
class LogStream
{
public:
	LogStream() {}
	virtual ~LogStream() {}

	// override function to print
	virtual void print(const char* sz) = 0;
};

//-----------------------------------------------------------------------------
//! Class that is used for logging purposes

//! This class can output to different 
//! files at the same time. In FEBio, it is used for logging purposes.
//! At this time it outputs data to the screen (stdout) and to an external text file.
//! Note that this class is implemented as a singleton, in other words, only one
//! instance can be created.

class Logfile  
{
public:
	enum MODE { NEVER = 0, FILE_ONLY = 1, SCREEN_ONLY, FILE_AND_SCREEN };

public:

	//! obtain a pointer to the logfile
	static Logfile* GetInstance();

	//! destructor
	virtual ~Logfile();

	//! open a new logfile
	bool open(const char* szfile);

	//! append to existing file
	bool append(const  char* szfile);

	//! formatted printing
	void printf(const char* sz, ...);

	//! print a nice box
	void printbox(const char* sztitle, const char* sz, ...);

	//! set the loggin mode
	MODE SetMode(MODE mode);

	//! get the loggin mode
	MODE GetMode() {return m_mode; }

	//! flush the logfile
	void flush() { if (m_fp) fflush(m_fp); }

	//! return the file name
	const char* FileName() { return m_szfile; }

	//! return a file pointer
	operator FILE* () { return m_fp; }

	//! returns if the logfile is ready to be written to
	bool is_valid() { return (m_fp != 0); }

	// set the log stream
	void SetLogStream(LogStream* ps) { m_ps = ps; }

private:
	//! constructor is private so that you cannot create it directly
	Logfile();
	Logfile(const Logfile& log){}

protected:
	FILE*	m_fp;	//!< the actual log file

	LogStream*	m_ps;	//!< This stream is used to output to the screen

	MODE	m_mode;	//!< mode of log file

	char	m_szfile[256];	//!< file name of logfile

	static Logfile* m_plog;	//!< the one and only logfile
};


#endif // !defined(AFX_LOGFILE_H__18090874_EE74_42E8_AB1B_1874D975D646__INCLUDED_)
