#include "CSVFile.h"
#include <fstream>
#include <sstream>

bool CSVFile::open(std::string FileName, std::string delim)
{
    this->delimiter = delim;
    std::ifstream file;
    file.open(FileName);
    if (file.is_open()) {
        std::vector<std::string> columnName = readNextRow(file);
        for (size_t i = 0; i < columnName.size(); i++) {
            this->columnName.insert(std::pair<std::string, size_t>(columnName[i], i));
        }
        this->sizeColumn_ = columnName.size();
        columnName.clear();

        while (!file.eof()) {
            auto buff = readNextRow(file);
            if(buff.size() == sizeColumn_)
                dataFile.push_back(buff);
        }

        this->sizeRow_ = dataFile.size();
        
    }
    else {
        return false;
    }
    return true;
}

size_t CSVFile::getSizeRow()
{
    return sizeRow_;
}

size_t CSVFile::getSizeColumn()
{
    return sizeColumn_;
}

size_t CSVFile::getColumnNumByName(std::string columnName)
{
    size_t buff;
    if (this->columnName.find(columnName) != this->columnName.end()) {
        return (*this->columnName.find(columnName)).second;
    }
    else {
        throw "This column not exists";
    }
}

std::string CSVFile::getColumnValueByNum(size_t rowNum, size_t columnNum)
{
    if (rowNum < sizeRow_ && columnNum < sizeColumn_) {
        return dataFile[rowNum][columnNum];
    }
    else {
        throw "Data not in file";
    }
}

std::string CSVFile::getColumnValueByName(size_t rowNum, std::string columnName)
{
    size_t columnNum = getColumnNumByName(columnName);

    return getColumnValueByNum(rowNum, columnNum);
}

std::vector<std::string> CSVFile::getColumnByNum(std::size_t columnNum)
{
    if (columnNum < sizeColumn_) {
        std::vector<std::string> buff;
        for (size_t i = 0; i < sizeRow_; i++) {
            buff.push_back(dataFile[i][columnNum]);
        }
        return buff;
    }
    else {
        throw "Exit from massiv data";
    }
}

std::vector<std::string> CSVFile::getColumnByName(std::string columnName)
{
    size_t columnNum = getColumnNumByName(columnName);
    
    return getColumnByNum(columnNum);
}

std::vector<std::string> CSVFile::readNextRow(std::ifstream& str)
{
    std::string line;
    std::vector<std::string> data;
    std::getline(str, line);

    data.clear();
    std::string::size_type start = 0;
    std::string::size_type end = 0;

    while ((end = line.find(this->delimiter, start)) != std::string::npos)
    {
        data.push_back(line.substr(start, end - start));
        start = end + this->delimiter.length();
    }

    data.push_back(line.substr(start, line.size()));

    return data;
}
