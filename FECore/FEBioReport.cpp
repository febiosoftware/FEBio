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
#include "FEBioReport.h"
#include "XMLWriter.h"
#include "XMLReader.h"
#include <sstream>

using namespace std;

const FEReportItem& FEReportSection::GetItem(size_t i) const 
{ 
	if (i >= m_items.size()) throw std::out_of_range("Invalid item index"); 
	return *m_items[i]; 
}

void FEReportTable::AddColumn(const std::string& name, const std::vector<std::string>& data)
{
	columns.emplace_back(name, FEReportTable::Text, std::vector<TableEntry>(data.begin(), data.end()));
}

void FEReportTable::AddColumn(const std::string& name, const std::vector<double>& data, const std::string& units)
{
	std::vector<TableEntry> columnData;
	columnData.reserve(data.size());
	for (const auto& value : data)
	{
		columnData.push_back(value);
	}
	columns.emplace_back(name + (units.empty() ? "" : " (" + units + ")"), FEReportTable::Numeric, std::move(columnData), units);
}

const FEReportTable::TableColumn& FEReportTable::GetColumn(size_t i) const
{
	if (i >= columns.size()) throw std::out_of_range("Invalid column index");
	return columns[i];
}

FEReportTable::TableEntry FEReportTable::GetEntry(size_t row, size_t column) const
{
	FEReportTable::TableEntry entry;
	if (column < columns.size())
	{
		const auto& col = columns[column];
		if (row < col.data.size())
			entry = col.data[row];
	}
	return entry;
}

struct FEBioReport::Imp
{
	string title;
	string optionsFile;
	int m_status = 0;

	vector<FEReportSectionPtr> m_sections;
};

FEBioReport::FEBioReport() : m(*new Imp)
{
}

FEBioReport::~FEBioReport() { delete &m; }

void FEBioReport::SetTitle(const std::string& title)
{
	m.title = title;
}

void FEBioReport::SetOptionsFile(const std::string& filename)
{
	m.optionsFile = filename;
}

void FEBioReport::SetStatus(int status)
{
	m.m_status = status;
}

FEReportSection& FEBioReport::AddSection(const std::string& name)
{
	auto section = make_unique<FEReportSection>(this);
	section->name = name;
	m.m_sections.push_back(std::move(section));
	return *m.m_sections.back();
}

FEReportText& FEReportSection::AddText(const std::string& text)
{
	auto item = make_unique<FEReportText>(text);
	m_items.push_back(std::move(item));
	return *static_cast<FEReportText*>(m_items.back().get());
}

FEReportFile& FEReportSection::AddFile(const std::string& filename, const std::string& description)
{
	auto item = make_unique<FEReportFile>(filename, description);
	m_items.push_back(std::move(item));
	return *static_cast<FEReportFile*>(m_items.back().get());
}

FEReportValue& FEReportSection::AddValue(const std::string& name, const std::string& value, const std::string& units)
{
	auto item = make_unique<FEReportValue>(name, value, units);
	m_items.push_back(std::move(item));
	return *static_cast<FEReportValue*>(m_items.back().get());
}

FEReportValue& FEReportSection::AddValue(const std::string& name, int value, const std::string& units)
{
	auto item = make_unique<FEReportValue>(name, std::to_string(value), units);
	m_items.push_back(std::move(item));
	return *static_cast<FEReportValue*>(m_items.back().get());
}

FEReportValue& FEReportSection::AddValue(const std::string& name, double value, const std::string& units)
{
	auto item = make_unique<FEReportValue>(name, std::to_string(value), units);
	m_items.push_back(std::move(item));
	return *static_cast<FEReportValue*>(m_items.back().get());
}

void FEReportSection::AddValues(const FEParameterList& pl)
{
	if (pl.Parameters() == 0) return;
	FEParamIteratorConst it = pl.first();
	for (int i = 0; i < pl.Parameters(); ++i, ++it)
	{
		const FEParam& param = *it;
		switch (param.type())
		{
		case FE_PARAM_BOOL: AddValue(param.name(), param.value<bool>()); break;
		case FE_PARAM_INT:
		{
			if (param.enums() != nullptr)
				AddValue(param.name(), param.enumKey());
			else
				AddValue(param.name(), param.value<int>());
			break;
		}
		case FE_PARAM_DOUBLE: AddValue(param.name(), param.value<double>()); break;
		case FE_PARAM_STD_STRING: AddValue(param.name(), param.value<std::string>()); break;
		case FE_PARAM_STRING: AddValue(param.name(), param.cvalue()); break;
		default:
			assert(false);
		}
	}
}

FEReportTable& FEReportSection::AddTable()
{
	// count all tables in the current section to generate a unique ID for the new table
	size_t tableCount = report->TableCount();

	size_t tableIndex = tableCount + 1; // start from 1 for better readability

	auto item = make_unique<FEReportTable>();
	item->id = "table_" + std::to_string(tableIndex);
	m_items.push_back(std::move(item));
	return *static_cast<FEReportTable*>(m_items.back().get());
}

FEReportTableView& FEReportSection::AddTableView(const FEReportTable& table)
{
	auto item = make_unique<FEReportTableView>();
	item->tableId = table.id;
	m_items.push_back(std::move(item));
	return *static_cast<FEReportTableView*>(m_items.back().get());
}

FEReportChart& FEReportSection::AddChart(FEReportChart::ChartType chartType)
{
	auto item = make_unique<FEReportChart>();
	item->chartType = chartType;
	m_items.push_back(std::move(item));
	return *static_cast<FEReportChart*>(m_items.back().get());
}

std::string FEBioReport::GetTitle() const 
{ 
	return m.title; 
}

std::string FEBioReport::GetOptionsFile() const
{ 
	return m.optionsFile; 
}

int FEBioReport::GetStatus() const
{
	return m.m_status;
}

size_t FEBioReport::Sections() const
{
	return m.m_sections.size();
}

const FEReportSection& FEBioReport::GetSection(size_t i) const
{
	if (i >= m.m_sections.size()) throw std::out_of_range("Invalid section index");
	return *m.m_sections[i];
}

FEReportTable FEBioReport::GetTable(const std::string& tableId) const
{
	FEReportTable table;
	for (const auto& section : m.m_sections) {
		for (const auto& item : section->m_items) {
			if (auto tableItem = dynamic_cast<FEReportTable*>(item.get())) {
				if (tableItem->id == tableId) {
					table = *tableItem;
					return table;
				}
			}
		}
	}
	return table;
}

FEReportTable::TableColumn FEBioReport::GetTableColumn(const std::string& tableId, const std::string& columnName) const
{
	FEReportTable table = GetTable(tableId);
	for (const auto& column : table.columns) {
		if (column.name == columnName) {
			return column;
		}
	}
	return FEReportTable::TableColumn("", FEReportTable::Text, {});
}

size_t FEBioReport::TableCount() const
{
	size_t count = 0;
	for (const auto& section : m.m_sections) {
		for (const auto& item : section->m_items) {
			if (dynamic_cast<FEReportTable*>(item.get())) {
				++count;
			}
		}
	}
	return count;
}

bool FEBioReport::Write(const std::string& filename) const
{
	XMLWriter xml;
	if (!xml.open(filename.c_str())) return false;

	xml.add_branch("febio_report");
	{
		xml.add_branch("metadata");
		{
			xml.add_leaf("title", m.title);
			xml.add_leaf("options_file", m.optionsFile);
			xml.add_leaf("status", m.m_status);
		}
		xml.close_branch(); // metadata

		for (const auto& section : m.m_sections) 
		{
			XMLElement el("section");
			el.add_attribute("name", section->name);
			xml.add_branch(el);
			{
				for (const auto& item : section->m_items) 
				{
					XMLElement xmlItem("item");
					xmlItem.add_attribute("type", item->Type());
					if (auto textItem = dynamic_cast<FEReportText*>(item.get())) 
					{
						xml.add_branch(xmlItem);
						{
							xml.add_leaf("text", textItem->text);
						}
						xml.close_branch();
					}
					if (auto fileItem = dynamic_cast<FEReportFile*>(item.get())) 
					{
						xml.add_branch(xmlItem);
						{

							xml.add_leaf("filename", fileItem->filename);
							if (!fileItem->description.empty())
								xml.add_leaf("description", fileItem->description);
						}
						xml.close_branch();
					}
					if (auto valueItem = dynamic_cast<FEReportValue*>(item.get()))
					{
						xml.add_branch(xmlItem);
						{
							xml.add_leaf("name", valueItem->name);
							xml.add_leaf("value", valueItem->value);
							if (!valueItem->units.empty())
								xml.add_leaf("units", valueItem->units);
						}
						xml.close_branch();
					}
					if (auto tableItem = dynamic_cast<FEReportTable*>(item.get()))
					{
						xmlItem.add_attribute("id", tableItem->id);
						xml.add_branch(xmlItem);
						{
							for (const auto& column : tableItem->columns) {
								XMLElement xmlColumn("column");
								xmlColumn.add_attribute("name", column.name);
								xmlColumn.add_attribute("type", column.type == FEReportTable::Numeric ? "numeric" : "text");
								if (!column.units.empty())
									xmlColumn.add_attribute("units", column.units);
								xml.add_branch(xmlColumn);
								{
									std::stringstream colValues;
									size_t count = column.data.size();
									for (int i = 0; i < count; ++i) {
										const FEReportTable::TableEntry& entry = column.data[i];

										if (std::holds_alternative<double>(entry)) {
											colValues << std::get<double>(entry);
										}
										else if (std::holds_alternative<std::string>(entry)) {
											colValues << std::get<std::string>(entry);
										}
										if (i < count - 1) colValues << ",";
									}
									xml.add_leaf("values", colValues.str());
								}
								xml.close_branch(); // column
							}
						}
						xml.close_branch();
					}
					if (auto tableViewItem = dynamic_cast<FEReportTableView*>(item.get()))
					{
						xml.add_branch(xmlItem);
						{
							XMLElement srcTag("source");
							srcTag.add_attribute("ref", tableViewItem->tableId);
							xml.add_empty(srcTag);
							if (!tableViewItem->tableTitle.empty())
								xml.add_leaf("title", tableViewItem->tableTitle);
							if (!tableViewItem->tableCaption.empty())
								xml.add_leaf("caption", tableViewItem->tableCaption);
						}
						xml.close_branch();
					}
					if (auto chartItem = dynamic_cast<FEReportChart*>(item.get()))
					{
						switch (chartItem->chartType)
						{
						case FEReportChart::Line   : xmlItem.add_attribute("chart_type", "line"   ); break;
						case FEReportChart::Bar    : xmlItem.add_attribute("chart_type", "bar"    ); break;
						case FEReportChart::Pie    : xmlItem.add_attribute("chart_type", "pie"    ); break;
							default: break;
						}

						xml.add_branch(xmlItem);
						{
							if (!chartItem->chartTitle.empty())
								xml.add_leaf("title", chartItem->chartTitle);
							if (!chartItem->chartCaption.empty())
								xml.add_leaf("caption", chartItem->chartCaption);

							for (const auto& series : chartItem->dataSeries)
							{
								XMLElement seriesTag("series");
								if (!series.name.empty())
									seriesTag.add_attribute("name", series.name);
								xml.add_branch(seriesTag);
								{
									for (const auto& data : series.data) {
										XMLElement dataTag("data");
										switch (data.role)
										{
										case FEReportChart::X: dataTag.add_attribute("role", "x"); break;
										case FEReportChart::Y: dataTag.add_attribute("role", "y"); break;
										case FEReportChart::Label: dataTag.add_attribute("role", "label"); break;
										case FEReportChart::Value: dataTag.add_attribute("role", "value"); break;
										default: break;
										}
										dataTag.add_attribute("table", data.tableId);
										dataTag.add_attribute("column", data.columnName);
										xml.add_empty(dataTag);
									}
								}
								xml.close_branch(); // series
							}
						}
						xml.close_branch(); // chart item
					}
				}
			}
			xml.close_branch(); // section
		}
	}
	xml.close_branch(); // febio_report

	xml.close();
	
	return true;
}

bool FEBioReport::Load(const std::string& filename)
{
	XMLReader xml;
	if (!xml.Open(filename.c_str()))
		return false;

	XMLTag tag;
	if (!xml.FindTag("febio_report", tag))
	{
		xml.Close();
		return false;
	}

	++tag;
	while (!tag.isend())
	{
		if ((tag == "metadata") && !tag.isleaf())
		{
			++tag;
			while (!tag.isend())
			{
				if      (tag == "title"       ) tag.value(m.title);
				else if (tag == "options_file") tag.value(m.optionsFile);
				else if (tag == "status"      ) tag.value(m.m_status);
				++tag;
			}
		}
		else if (tag == "section" && !tag.isleaf())
		{
			string sectionName;
			XMLAtt* att = tag.AttributePtr("name");
			if (att) sectionName = att->m_val;
			FEReportSection& section = AddSection(sectionName);
			++tag;
			while (!tag.isend())
			{
				if ((tag == "item") && !tag.isleaf())
				{
					string itemType;
					XMLAtt* typeAtt = tag.AttributePtr("type");
					if (typeAtt) itemType = typeAtt->m_val;
					if (itemType == "text")
					{
						FEReportText& textItem = section.AddText("");
						++tag;
						while (!tag.isend())
						{
							if (tag == "text") tag.value(textItem.text);
							++tag;
						}
					}
					else if (itemType == "file")
					{
						FEReportFile& fileItem = section.AddFile("", "");
						++tag;
						while (!tag.isend())
						{
							if      (tag == "filename"   ) tag.value(fileItem.filename);
							else if (tag == "description") tag.value(fileItem.description);
							++tag;
						}
					}
					else if (itemType == "value")
					{
						FEReportValue& valueItem = section.AddValue("", "");
						++tag;
						while (!tag.isend())
						{
							if      (tag == "name" ) tag.value(valueItem.name);
							else if (tag == "value") tag.value(valueItem.value);
							else if (tag == "units") tag.value(valueItem.units);
							++tag;
						}
					}
					else if (itemType == "table")
					{
						FEReportTable& tableItem = section.AddTable();
						XMLAtt* idAtt = tag.AttributePtr("id");
						if (idAtt) tableItem.id = idAtt->m_val;
						++tag;
						while (!tag.isend())
						{
							if (tag == "column" && !tag.isleaf())
							{
								string columnName;
								string columnType;
								string columnUnits;
								XMLAtt* colNameAtt = tag.AttributePtr("name");
								XMLAtt* colTypeAtt = tag.AttributePtr("type");
								XMLAtt* colUnitsAtt = tag.AttributePtr("units");
								if (colNameAtt) columnName = colNameAtt->m_val;
								if (colTypeAtt) columnType = colTypeAtt->m_val;
								if (colUnitsAtt) columnUnits = colUnitsAtt->m_val;
								vector<string> values;
								++tag;
								while (!tag.isend())
								{
									if (tag == "values")
									{
										string valStr;
										tag.value(valStr);
										std::stringstream ss(valStr);
										string item;
										while (std::getline(ss, item, ',')) {
											values.push_back(item);
										}
									}
									++tag;
								}
								if (columnType == "numeric") {
									vector<double> numericValues;
									for (const auto& val : values) {
										numericValues.push_back(std::stod(val));
									}
									tableItem.AddColumn(columnName, numericValues, columnUnits);
								}
								else {
									tableItem.AddColumn(columnName, values);
								}
							}
							++tag;
						}
					}
					else if (itemType == "table_view")
					{
						FEReportTableView& tableViewItem = section.AddTableView(FEReportTable());
						++tag;
						while (!tag.isend())
						{
							if (tag == "source")
							{
								XMLAtt* refAtt = tag.AttributePtr("ref");
								if (refAtt) tableViewItem.tableId = refAtt->m_val;
							}
							else if (tag == "title")
							{
								tag.value(tableViewItem.tableTitle);
							}
							else if (tag == "caption")
							{
								tag.value(tableViewItem.tableCaption);
							}
							else
							{
								tag.skip();
							}
							++tag;
						}
					}
					else if (itemType == "chart")
					{
						FEReportChart& chartItem = section.AddChart(FEReportChart::Line);
						XMLAtt* chartTypeAtt = tag.AttributePtr("chart_type");
						if (chartTypeAtt) {
							string chartTypeStr = chartTypeAtt->m_val;
							if      (chartTypeStr == "line"   ) chartItem.chartType = FEReportChart::Line;
							else if (chartTypeStr == "bar"    ) chartItem.chartType = FEReportChart::Bar;
							else if (chartTypeStr == "pie"    ) chartItem.chartType = FEReportChart::Pie;
						}
						++tag;
						while (!tag.isend())
						{
							if (tag == "title")
							{
								tag.value(chartItem.chartTitle);
							}
							else if (tag == "caption")
							{
								tag.value(chartItem.chartCaption);
							}
							else if (tag == "series" && !tag.isleaf())
							{
								FEReportChart::ChartDataSeries& series = chartItem.AddDataSeries("");
								XMLAtt* seriesNameAtt = tag.AttributePtr("name");
								if (seriesNameAtt) series.name = seriesNameAtt->m_val;
								++tag;
								while (!tag.isend())
								{
									if (tag == "data")
									{
										string roleString, tableId, columnName;
										XMLAtt* roleAtt = tag.AttributePtr("role");
										XMLAtt* tableAtt = tag.AttributePtr("table");
										XMLAtt* columnAtt = tag.AttributePtr("column");

										FEReportChart::DataRole role = FEReportChart::X; // default role
										if (roleAtt) {
											roleString = roleAtt->m_val;
											if      (roleString == "x"    ) role = FEReportChart::X;
											else if (roleString == "y"    ) role = FEReportChart::Y;
											else if (roleString == "label") role = FEReportChart::Label;
											else if (roleString == "value") role = FEReportChart::Value;
										}

										if (tableAtt) tableId = tableAtt->m_val;
										if (columnAtt) columnName = columnAtt->m_val;
										FEReportChart::ChartData data{ role, tableId, columnName };
										series.data.emplace_back(data);
									}
									else
									{
										tag.skip();
									}
									++tag;
								}
							}
							else
							{
								tag.skip();
							}
							++tag;
						}
					}
					else
					{
						// Handle other item types here
						tag.skip();
					}
				}
				else
				{
					tag.skip();
				}
				++tag;
			}
		}
		else
		{
			tag.skip();
		}
		++tag;
	}

	xml.Close();

	return true;
}
