#pragma once
#include "PoligonBuilding.h"
#include "BS.h"
class BsLinkToBuild
{
public:
	//в данной функции заложена логика возвращающая вес в зависимости от покрытия БС
	double getWeight();

	double recalcWeight(BS* Bs);

	void changeBuildingWeight(double weightIn);

	BsLinkToBuild(PoligonBuilding* toBuild);

private:

	double	weight;

	PoligonBuilding* Build;
};