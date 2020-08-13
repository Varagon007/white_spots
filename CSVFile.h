#pragma once
#include <string>
#include <vector>
#include <map>

class CSVFile
{
public:
	//�������� �����
	bool open(std::string FileName, std::string delim = ";");
	//������� ����� �������
	size_t getSizeRow();
	//����� ��������
	size_t getSizeColumn();
	//����� ����� ������� �� ��������
	size_t getColumnNumByName(std::string columnName);
	//������� ������ �� �������
	std::string getColumnValueByNum(size_t rowNum, size_t columnNum);
	//������� ������ �� �������� ������
	std::string getColumnValueByName(size_t rowNum, std::string columnName);
	//������� ������� �� ������
	std::vector<std::string> getColumnByNum(std::size_t columnNum);
	//������� ������� �� ��������
	std::vector<std::string> getColumnByName(std::string columnName);
private:
	//����� �����
	size_t sizeRow_;
	//����� ��������
	size_t sizeColumn_;
	//�����������
	std::string delimiter = ";";
	//�������� �������, ����� ������ �������� �������������
	std::multimap<std::string, size_t> columnName;
	//���������� �����
	//������ ������ , ������ �������
	std::vector<std::vector<std::string>> dataFile;

	std::vector<std::string> readNextRow(std::ifstream& str);
};

