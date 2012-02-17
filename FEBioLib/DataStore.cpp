// DataStore.cpp: implementation of the DataStore class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "DataStore.h"
#include "log.h"

//-----------------------------------------------------------------------------
UnknownDataField::UnknownDataField(const char* sz)
{
	m_szdata[0] = 0;
	int l = strlen(sz);
	if (l > 63) l = 63;
	if (l>0) strncpy(m_szdata, sz, l);
}

//////////////////////////////////////////////////////////////////////
// DataStore
//////////////////////////////////////////////////////////////////////

DataStore::DataStore()
{
}

DataStore::~DataStore()
{
}

void DataStore::Clear()
{
	for (size_t i=0; i<m_data.size(); ++i) delete m_data[i];
	m_data.clear();
}

//-----------------------------------------------------------------------------

void DataStore::Write()
{
	for (size_t i=0; i<m_data.size(); ++i)
	{
		DataRecord& DR = *m_data[i];
		DR.Write();
	}
}

//-----------------------------------------------------------------------------

void DataStore::AddRecord(DataRecord* prec)
{
	static int nid = 1;
	prec->m_nid = nid++;
	m_data.push_back(prec);
}


//////////////////////////////////////////////////////////////////////
// DataRecord
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
DataRecord::DataRecord(FEModel* pfem, const char* szfile)
{
	m_pfem = pfem;
	m_nid = 0;
	m_szname[0] = 0;

	strcpy(m_szdelim, " ");
	
	m_bcomm = true;

	m_fp = 0;
	m_szfile[0] = 0;

	if (szfile)
	{
		strcpy(m_szfile, szfile);
		m_fp = fopen(szfile, "wt");
		fprintf(m_fp, "*Title:%s\n", m_pfem->GetTitle());
	}
}

//-----------------------------------------------------------------------------
DataRecord::~DataRecord()
{
	if (m_fp)
	{
		fclose(m_fp);
		m_fp = 0;
	}
}

//-----------------------------------------------------------------------------
void DataRecord::SetName(const char* sz)
{
	strcpy(m_szname, sz);
}

//-----------------------------------------------------------------------------
void DataRecord::SetDelim(const char* sz)
{
	strcpy(m_szdelim, sz);
}

//-----------------------------------------------------------------------------
bool DataRecord::Write()
{
	int nstep = m_pfem->GetCurrentStep()->m_ntimesteps;
	double ftime = m_pfem->m_ftime;
	double val;

	FILE* fplog = (FILE*) clog;

	// make a note in the log file
	fprintf(fplog, "\nData Record #%d\n", m_nid);
	fprintf(fplog, "===========================================================================\n");
	fprintf(fplog, "Step = %d\n", nstep);
	fprintf(fplog, "Time = %.9lg\n", ftime);
	fprintf(fplog, "Data = %s\n", m_szname);

	// see if we are saving the data to the logfile or to a 
	// seperate data file
	FILE* fp = m_fp;
	if (fp == 0)
	{
		// we store the data in the logfile
		fp = fplog;
	}
	else if (m_bcomm)
	{
		// we save the data in a seperate file
		fprintf(fplog, "File = %s\n", m_szfile);

		// make a note in the data file
		fprintf(fp,"*Step  = %d\n", nstep);
		fprintf(fp,"*Time  = %.9lg\n", ftime);
		fprintf(fp,"*Data  = %s\n", m_szname);
	}

	// save the data
	for (size_t i=0; i<m_item.size(); ++i)
	{
		fprintf(fp, "%d%s", m_item[i], m_szdelim);
		int nd = (int) m_data.size();
		for (int j=0; j<nd; ++j)
		{
			val = Evaluate(m_item[i], m_data[j]);
			fprintf(fp, "%lg", val);
			if (j!=nd-1) fprintf(fp, "%s", m_szdelim);
			else fprintf(fp, "\n");
		}
	}

	return true;
}

//-----------------------------------------------------------------------------

void DataRecord::SetItemList(const char* szlist)
{
	int i, n = 0, n0, n1, nn;
	char* ch;
	char* sz = (char*) szlist;
	int nread;
	do
	{
		ch = strchr(sz, ',');
		if (ch) *ch = 0;
		nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
		switch (nread)
		{
		case 1:
			n1 = n0;
			nn = 1;
			break;
		case 2:
			nn = 1;
			break;
		case 3:
			break;
		default:
			n0 = 0;
			n1 = -1;
			nn = 1;
		}

		for (i=n0; i<=n1; i += nn) ++n;

		if (ch) *ch = ',';
		sz = ch+1;
	}
	while (ch != 0);

	if (n != 0)
	{
		m_item.resize(n);

		sz = (char*) szlist;
		n = 0;
		do
		{
			ch = strchr(sz, ',');
			if (ch) *ch = 0;
			nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
			switch (nread)
			{
			case 1:
				n1 = n0;
				nn = 1;
				break;
			case 2:
				nn = 1;
			}

			for (i=n0; i<=n1; i += nn) m_item[n++] = i;
			assert(n <= (int) m_item.size());

			if (ch) *ch = ',';
			sz = ch+1;
		}
		while (ch != 0);
	}
	else SelectAllItems();
}

//-----------------------------------------------------------------------------

void DataRecord::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		ar << m_nid;
		ar << m_szname;
		ar << m_szdelim;
		ar << m_szfile;
		ar << m_bcomm;
		ar << m_item;
	}
	else
	{
		ar >> m_nid;
		ar >> m_szname;
		ar >> m_szdelim;
		ar >> m_szfile;
		ar >> m_bcomm;
		ar >> m_item;

		if (m_fp) fclose(m_fp);
		m_fp = 0;
		if (m_szfile[0] != 0)
		{
			// reopen data file for appending
			m_fp = fopen(m_szfile, "a+");
		}
	}
}
