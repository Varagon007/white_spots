#include "BsLinkToBuild.h"

double BsLinkToBuild::getWeight()
{   
    return weight;
}

double BsLinkToBuild::recalcWeight(BS* Bs)
{
    double Rbuff = Bs->Coordinate.Distance(&Build->PointCenter);
    if (2 * Bs->getR() >= Rbuff) {
        if (Bs->getR() >= Rbuff) {
            return weight = Build->getWeight();
        }
        else {
            return weight = Build->getWeight() * pow(Bs->getR() / Rbuff, 10);
        }
    }
    else {
        return 0.0;
    }

    return 0.0;
}

void BsLinkToBuild::changeBuildingWeight(double weightIn)
{
    Build -> changeWeight(weightIn);
}

BsLinkToBuild::BsLinkToBuild(PoligonBuilding* toBuild)
{
    this->Build = toBuild;
}
