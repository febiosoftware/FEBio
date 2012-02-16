#include "stdafx.h"
#include "FileBrowser.h"
#include <FL/filename.H>
#include <assert.h>
#include "Wnd.h"

//-----------------------------------------------------------------------------
// we had some problem figuring out how to access this function, so I've just
// copied it here.


size_t						/* O - Length of string */
strlcat(char	      *dst,	/* O - Destination string */
        const char *src,	/* I - Source string */
        size_t     size) /* I - Size of destination string buffer */
{
	size_t	srclen;		/* Length of source string */
	size_t	dstlen;		/* Length of destination string */
	dstlen = strlen(dst);
	size   -= dstlen + 1;
	if (!size) return (dstlen);	/* No room, return immediately... */
	srclen = strlen(src);
	if (srclen > size) srclen = size;
	memcpy(dst + dstlen, src, srclen);
	dst[dstlen + srclen] = '\0';
	return (dstlen + srclen);
}

//-----------------------------------------------------------------------------
size_t						/* O - Length of string */
strlcpy(char       *dst,	/* O - Destination string */
        const char *src,	/* I - Source string */
        size_t      size)	/* I - Size of destination string buffer */
{	
	size_t	srclen;		/* Length of source string */
	size --;
	srclen = strlen(src);
	if (srclen > size) srclen = size;
	memcpy(dst, src, srclen);
	dst[srclen] = '\0';
	return (srclen);
}

//-----------------------------------------------------------------------------
CFileBrowser::CFileBrowser(int x, int y, int w, int h, CWnd* pwnd) : Flx_Group(x, y, w, h), m_pWnd(pwnd)
{
	begin();
	{
		Fl_File_Browser* pf = m_pfile = new Fl_File_Browser(x, y, w, h);
		pf->type(FL_HOLD_BROWSER);
		pf->iconsize(16);
		AddCallback(pf, (FLX_CALLBACK) &CFileBrowser::OnSelectFile);
		pf->filter("*.feb");
	}
	end();
	box(FL_NO_BOX);

	// set the working directory
	set_dir(".");
}

//-----------------------------------------------------------------------------
CFileBrowser::~CFileBrowser()
{
}


//-----------------------------------------------------------------------------
void CFileBrowser::OnSelectFile(Fl_Widget* pw, void* pd)
{
	Fl_File_Browser* pf = dynamic_cast<Fl_File_Browser*>(pw);
	assert(pf);

	// get the selected filename
	const char* szfile = pf->text(pf->value());
	if (!szfile) return;

	// set the path name
	char szpath[1024];
	if (m_szdir[0] == 0) strcpy(szpath, szfile);
	else if (strcmp(m_szdir, "/") == 0) sprintf(szpath, "/%s", szfile);
	else sprintf(szpath, "%s/%s", m_szdir, szfile);

	// if the user double-clicks we either read the folder
	// or open the file
	if (Fl::event_clicks())
	{
		if (is_dir(szpath))
		{
			// user selected a directory so load that folder
			set_dir(szpath);

			// Reset the click count so that a click in the same spot won't
			// be treated as a triple-click.  We use a value of -1 because
			// the next click will increment click count to 0, which is what
			// we really want...
			Fl::event_clicks(-1);
		}
		else
		{
			m_pWnd->OpenFile(szpath);
			m_pfile->deselect();
		}
	}
}

//-----------------------------------------------------------------------------
bool CFileBrowser::is_dir(const char* szfile)
{
#if (defined(WIN32) && ! defined(__CYGWIN__)) || defined(__EMX__)
    return ((strlen(szfile) == 2 && szfile[1] == ':') || _fl_filename_isdir_quick(szfile));
#else
		// TODO: upgrade to 1.1.10 and call _fl_filename_isdir_quick instead
    return (fl_filename_isdir(szfile));
#endif /* WIN32 || __EMX__ */
}

//-----------------------------------------------------------------------------
void CFileBrowser::set_dir(const char* szpath)
{
	char	*dirptr;			// Pointer into directory

	// NULL == current directory
	if (szpath == NULL) szpath = ".";

#ifdef WIN32
	// See if the filename contains backslashes...
	char	*slash;				// Pointer to slashes
	char	fixpath[1024];		// Path with slashes converted
	if (strchr(szpath, '\\'))
	{
		// Convert backslashes to slashes...
		strlcpy(fixpath, szpath, sizeof(fixpath));

		for (slash = strchr(fixpath, '\\'); slash; slash = strchr(slash + 1, '\\'))
			*slash = '/';

		szpath = fixpath;
	}
#endif // WIN32

	if (szpath[0] != '\0')
	{
		// Make the directory absolute...
#if (defined(WIN32) && ! defined(__CYGWIN__))|| defined(__EMX__)
		if (szpath[0] != '/' && szpath[0] != '\\' && szpath[1] != ':')
#else
		if (szpath[0] != '/' && szpath[0] != '\\')
#endif /* WIN32 || __EMX__ */
			fl_filename_absolute(m_szdir, szpath);
		else
			strlcpy(m_szdir, szpath, sizeof(m_szdir));

		// Strip any trailing slash...
		dirptr = m_szdir + strlen(m_szdir) - 1;
		if ((*dirptr == '/' || *dirptr == '\\') && dirptr > m_szdir)	*dirptr = '\0';

		// See if we have a trailing .. or . in the filename...
		dirptr = m_szdir + strlen(m_szdir) - 3;
		if (dirptr >= m_szdir && strcmp(dirptr, "/..") == 0) 
		{
			// Yes, we have "..", so strip the trailing path...
			*dirptr = '\0';
			while (dirptr > m_szdir) 
			{
				if (*dirptr == '/') break;
				dirptr --;
			}

			if (dirptr >= m_szdir && *dirptr == '/') *dirptr = '\0';
		} 
		else if ((dirptr + 1) >= m_szdir && strcmp(dirptr + 1, "/.") == 0) 
		{
			// Strip trailing "."...
			dirptr[1] = '\0';
		}
	}
	else m_szdir[0] = '\0';

	// Rescan the directory...
	char pathname[1024];		// New pathname for filename field

	// Clear the current filename
	strlcpy(pathname, m_szdir, sizeof(pathname));
	if (pathname[0] && pathname[strlen(pathname) - 1] != '/') 
	{
		strlcat(pathname, "/", sizeof(pathname));
	}

	 // Build the file list...
	m_pfile->load(m_szdir);
}
