#include "BSonBuilding.h"

std::multimap <int, int> BSonBuilding::BuildingtoBS;

void BSonBuilding::addWeight(double WeightIn, PoligonBuilding* Build_in, int IdBSCustom, int IdBuildCustom)
{
    BSWeight += WeightIn;
    if (Build_in->getWeight() == 0) {
        LinkToBuilding.insert(std::pair<double, PoligonBuilding*>(0, Build_in));
    }
    else {
        LinkToBuilding.insert(std::pair<double, PoligonBuilding*>(WeightIn / (Build_in->getWeight()), Build_in));
    }

    BSonBuilding::addLink(IdBSCustom, IdBuildCustom);
    
}

void BSonBuilding::changeWeight(double WeightIn)
{
    if (WeightIn < 0 && fabs(WeightIn) >= BSWeight) {
        BSWeight = 0;
    }
    else {
        BSWeight += WeightIn;
    }
}

void BSonBuilding::recalcWeight()
{
    this->resetWeight();
    for (auto iBuild = LinkToBuilding.begin(); iBuild != LinkToBuilding.end(); iBuild++) {
        BSWeight += (*iBuild).first * (*iBuild).second->getWeight();
    }
}

void BSonBuilding::resetWeight()
{
    BSWeight = 0;
}

double BSonBuilding::getWeight()
{
    return BSWeight;
}

BSonBuilding BSonBuilding::operator=(const PoligonBuilding& rv)
{
	this->ID_build = rv.ID_build;
    strcpy_s(this->Type, rv.Type);
    strcpy_s(this->Wall_Mat, rv.Wall_Mat);
    this->Building = rv.Building;
    this->Etage = rv.Etage;
    this->Area = rv.Area;
    this->People_D = rv.People_D;
    this->People_N = rv.People_N;
    strcpy_s(this->Volcano_type, rv.Volcano_type);
    this->Clutter_num = rv.Clutter_num;
    strcpy_s(this->Clutter_type, rv.Clutter_type);
    this->Lon = rv.Lon;
    this->Lat = rv.Lat;
    strcpy_s(this->Address_NP, rv.Address_NP);
    strcpy_s(this->Address_street, rv.Address_street);
    strcpy_s(this->Address_number, rv.Address_number);
    strcpy_s(this->FIAS,rv.FIAS);
    this->Polygon = rv.Polygon;
    this->PointCenter = rv.PointCenter;
	return (*this);
}
