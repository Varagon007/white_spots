#pragma once
#include <string>
#include <vector>
#include <map>

class CSVFile
{
public:
	//открытие файла
	bool open(std::string FileName, std::string delim = ";");
	//сколько всего строчек
	size_t getSizeRow();
	//число столбцов
	size_t getSizeColumn();
	//номер номер столбца по названию
	size_t getColumnNumByName(std::string columnName);
	//вернуть €чейку по позиции
	std::string getColumnValueByNum(size_t rowNum, size_t columnNum);
	//вернуть €чейку по названию столца
	std::string getColumnValueByName(size_t rowNum, std::string columnName);
	//вернуть столбец по номеру
	std::vector<std::string> getColumnByNum(std::size_t columnNum);
	//вернуть столбец по названию
	std::vector<std::string> getColumnByName(std::string columnName);
private:
	//число строк
	size_t sizeRow_;
	//число столбцов
	size_t sizeColumn_;
	//разделитель
	std::string delimiter = ";";
	//название колонки, чтобы быстро вытащить идентификатор
	std::multimap<std::string, size_t> columnName;
	//содержание файла
	//первое строка , второе столбец
	std::vector<std::vector<std::string>> dataFile;

	std::vector<std::string> readNextRow(std::ifstream& str);
};

