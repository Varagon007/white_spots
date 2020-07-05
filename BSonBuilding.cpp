#include "BSonBuilding.h"

std::multimap <int, int> BSonBuilding::BuildingtoBS;

void BSonBuilding::addWeight(double WeightIn, int IdBSCustom, int IdBuildCustom)
{
    Weight += WeightIn;
    LinkToBuilding.insert(IdBuildCustom);
    BSonBuilding::addLink(IdBuildCustom, IdBSCustom);
}

void BSonBuilding::changeWeight(double WeightIn)
{
    if (WeightIn < 0 && fabs(WeightIn) >= Weight) {
        Weight = 0;
    }
    else {
        Weight += WeightIn;
    }
}

double BSonBuilding::getWeight()
{
    return Weight;
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
