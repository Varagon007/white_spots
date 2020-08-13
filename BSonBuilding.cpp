#include "BSonBuilding.h"

std::multimap <int, int> BSonBuilding::BuildingtoBS;

void BSonBuilding::addLink(PoligonBuilding* Building)
{
    //LinkToBuilding.insert(BsLinkToBuild(Building));
    LinkToBuilding.push_back(BsLinkToBuild(Building));
}

void BSonBuilding::recalcWeight()
{
    this->resetWeight();
    for (auto iBuild = LinkToBuilding.begin(); iBuild != LinkToBuilding.end(); iBuild++) {
        Weight += (*iBuild).recalcWeight(this);
    }
}

void BSonBuilding::resetWeight()
{
    Weight = 0;
}

double BSonBuilding::getWeight()
{
    return Weight;
}

BSonBuilding::BSonBuilding()
{  
   // BS();
}

BSonBuilding::BSonBuilding(PoligonBuilding setupBuild)
{
    this->setupBuild = setupBuild;
    //OGRPoint Point_buff;
    //this->setupBuild.Polygon.Centroid(&Point_buff);
    this->Coordinate = setupBuild.PointCenter;
   // BS();
}

bool BSonBuilding::operator<(const BSonBuilding& vc) const
{
    return this->Weight < vc.Weight;
}

bool BSonBuilding::operator>(const BSonBuilding& vc) const
{
    return this->Weight > vc.Weight;
}
