// Logfile.h: interface for the Logfile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LOGFILE_H__18090874_EE74_42E8_AB1B_1874D975D646__INCLUDED_)
#define AFX_LOGFILE_H__18090874_EE74_42E8_AB1B_1874D975D646__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <stdio.h>
#include "fecore_export.h"

//-----------------------------------------------------------------------------
// class used to create an abstract interface to a screen
class FECOREDLL_EXPORT LogStream
{
public:
	LogStream() {}
	virtual ~LogStream() {}

	// override function to print
	virtual void print(const char* sz) = 0;

	// flush the stream
	virtual void flush() {}
};

//-----------------------------------------------------------------------------
// A stream that outputs to a file
class FECOREDLL_EXPORT LogFileStream : public LogStream
{
public:
	// constructor
	LogFileStream();

	// destructor
	~LogFileStream();

	// open the file
	bool open(const char* szfile);

	// open for appending
	bool append(const char* szfile);

	// close the file stream
	void close();

	// get the file handle
	FILE* GetFileHandle() { return m_fp; }

public:
	// print text to the file
	void print(const char* sz);

	// flush the stream
	void flush();

private:
	FILE*	m_fp;
};

//-----------------------------------------------------------------------------
//! Class that is used for logging purposes

//! This class can output to different 
//! files at the same time.
//! At this time it outputs data to the screen (stdout) and to an external text file.
//! Note that this class is implemented as a singleton, in other words, only one
//! instance can be created.

class FECOREDLL_EXPORT Logfile
{
public:
	enum MODE { LOG_NEVER = 0, LOG_FILE = 1, LOG_SCREEN, LOG_FILE_AND_SCREEN };

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
	void flush();

	//! close the logfile
	void close();

	//! return the file name
	const char* FileName() { return m_szfile; }

	//! returns if the logfile is ready to be written to
	bool is_valid() { return (m_fp != 0); }

	// set the log stream
	void SetLogStream(LogStream* ps) { m_ps = ps; }

	// return the file handle
	operator FILE* () { return (m_fp ? m_fp->GetFileHandle() : 0); }

private:
	//! constructor is private so that you cannot create it directly
	Logfile();
	Logfile(const Logfile& log){}

protected:
	LogFileStream*	m_fp;	//!< the actual log file

	LogStream*	m_ps;	//!< This stream is used to output to the screen

	MODE	m_mode;	//!< mode of log file

	char	m_szfile[256];	//!< file name of logfile

	static Logfile* m_plog;	//!< the one and only logfile
};


#endif // !defined(AFX_LOGFILE_H__18090874_EE74_42E8_AB1B_1874D975D646__INCLUDED_)
