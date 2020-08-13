#pragma once
#include "PoligonBuilding.h"
#include "BS.h"
class BsLinkToBuild
{
public:
	//� ������ ������� �������� ������ ������������ ��� � ����������� �� �������� ��
	double getWeight();

	double recalcWeight(BS* Bs);

	void changeBuildingWeight(double weightIn);

	BsLinkToBuild(PoligonBuilding* toBuild);

private:

	double	weight;

	PoligonBuilding* Build;
};