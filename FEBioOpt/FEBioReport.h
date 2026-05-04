/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#pragma once
#include "febioopt_api.h"
#include <string>
#include <vector>
#include <memory>
#include <variant>

class FEBIOOPT_API FEReportItem
{
public:
	FEReportItem() = default;
	virtual ~FEReportItem() = default;

	virtual std::string Type() const = 0;
};

using FEReportItemPtr = std::unique_ptr<FEReportItem>;

class FEBIOOPT_API FEReportText : public FEReportItem
{
public:
	std::string text;

	FEReportText(const std::string& text) : text(text) {}

	std::string Type() const override { return "text"; }
};

class FEBIOOPT_API FEReportFile : public FEReportItem
{
public:
	std::string filename;
	std::string description;

	FEReportFile(const std::string& filename, const std::string& description = "") : filename(filename), description(description) {}

	std::string Type() const override { return "file"; }
};

class FEBIOOPT_API FEReportValue: public FEReportItem
{
public:
	std::string name;
	std::string value;
	std::string units;

	FEReportValue(const std::string& name, const std::string& value, const std::string& units = "") : name(name), value(value), units(units) {}

	std::string Type() const override { return "value"; }
};

class FEBIOOPT_API FEReportTable : public FEReportItem
{
public:
	using TableEntry = std::variant<double, std::string>;

	enum ColumnType { Numeric, Text };

	struct TableColumn
	{
		std::string name;
		std::string units;
		ColumnType type;
		std::vector<TableEntry> data;
		TableColumn(const std::string& name, ColumnType type, std::vector<TableEntry>&& data, const std::string& units = "") : name(name), type(type), data(std::move(data)), units(units) {}
	};

public:
	std::string id;
	std::vector<TableColumn> columns;

	FEReportTable() = default;

	std::string Type() const override { return "table"; }

	void AddColumn(const std::string& name, const std::vector<std::string>& data);
	void AddColumn(const std::string& name, const std::vector<double>& data, const std::string& units = "");

	size_t Columns() const { return columns.size(); }
	const TableColumn& GetColumn(size_t i) const;

	std::string GetColumnName(size_t i) const { return GetColumn(i).name; }

	size_t Rows() const { return columns.empty() ? 0 : columns[0].data.size(); }

	TableEntry GetEntry(size_t row, size_t column) const;
};

class FEBIOOPT_API FEReportTableView : public FEReportItem
{
public:
	std::string tableId;
	std::string tableTitle;
	std::string tableCaption;

	std::string Type() const override { return "table_view"; }

	FEReportTableView& SetTitle(const std::string& title) { tableTitle = title; return *this; }
	FEReportTableView& SetCaption(const std::string& caption) { tableCaption = caption; return *this; }
};

class FEBIOOPT_API FEReportChart : public FEReportItem
{
public:
	enum ChartType { Line, Bar, Pie };

	enum DataRole { X, Y, Label, Value };

	struct ChartData {
		DataRole role;
		std::string tableId; // reference to a table
		std::string columnName; // name of the column in the referenced table
	};

	struct ChartDataSeries
	{
		std::string name;
		std::vector<ChartData> data;

		ChartDataSeries& AddData(DataRole role, const std::string& tableId, const std::string& columnName) {
			data.push_back({ role, tableId, columnName });
			return *this;
		}

		ChartData FindDataByRole(DataRole role) const {
			for (const auto& d : data) {
				if (d.role == role) return d;
			}
			return ChartData{ role, "", "" }; // return an empty data if not found
		}

		std::vector<ChartData> FindAllDataByRole(DataRole role) const {
			std::vector<ChartData> results;
			for (const auto& d : data) {
				if (d.role == role) results.push_back(d);
			}
			return results;
		}
	};

public:
	ChartType chartType;
	std::string chartTitle;
	std::string chartCaption;
	std::vector<ChartDataSeries> dataSeries;

	std::string Type() const override { return "chart"; }
	
	FEReportChart& SetTitle(const std::string& title) { chartTitle = title; return *this; }
	FEReportChart& SetCaption(const std::string& caption) { chartCaption = caption; return *this; }

	ChartDataSeries& AddDataSeries(const std::string& name) {
		dataSeries.emplace_back(ChartDataSeries{name});
		return dataSeries.back();
	}
};

class FEBIOOPT_API FEReportSection
{
public:
	std::string name;
	std::vector<FEReportItemPtr> m_items;

	FEReportSection() = default;
	FEReportSection(const FEReportSection&) = delete;
	FEReportSection& operator=(const FEReportSection&) = delete;

	size_t Items() const { return m_items.size(); }
	const FEReportItem& GetItem(size_t i) const;

	FEReportText& AddText(const std::string& text);

	FEReportFile& AddFile(const std::string& filename, const std::string& description);

	FEReportValue& AddValue(const std::string& name, const std::string& value, const std::string& units = "");
	FEReportValue& AddValue(const std::string& name, int value, const std::string& units = "");
	FEReportValue& AddValue(const std::string& name, double value, const std::string& units = "");

	FEReportTable& AddTable();

	FEReportTableView& AddTableView(const FEReportTable& table);

	FEReportChart& AddChart(FEReportChart::ChartType chartType);
};

using FEReportSectionPtr = std::unique_ptr<FEReportSection>;

class FEBIOOPT_API FEBioReport
{
	struct Imp;

public:
	FEBioReport();
	~FEBioReport();

	FEBioReport(const FEBioReport& other) = delete;

	FEBioReport& operator=(const FEBioReport& other) = delete;

	// Write the report to a file
	bool Write(const std::string& filename) const;

	// load a report from a file
	bool Load(const std::string& filename);

public: // API for building reports

	// set the study title
	void SetTitle(const std::string& title);

	// set the options file
	void SetOptionsFile(const std::string& filename);

	// set the status of the study (0 = failed, 1 = success)
	void SetStatus(int status);

	// add a section to the report
	FEReportSection& AddSection(const std::string& name);


public: // API for reading reports

	std::string GetTitle() const;
	std::string GetOptionsFile() const;
	int GetStatus() const;

	size_t Sections() const;
	const FEReportSection& GetSection(size_t i) const;

	FEReportTable GetTable(const std::string& tableId) const;

	FEReportTable::TableColumn GetTableColumn(const std::string& tableId, const std::string& columnName) const;

private:
	Imp& m;
};
