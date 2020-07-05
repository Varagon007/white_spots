#pragma once
#include "PoligonBuilding.h"
#include "set"
class BSonBuilding :
    public PoligonBuilding
{
private:
    double Weight;
public:

    static std::multimap <int, int> BuildingtoBS;
    //на вход поступает id кастомный
    void addWeight(double WeightIn, int IdBSCustom, int IdBuildCustom);

    void changeWeight(double WeightIn);

    double getWeight();

    static void addLink(int IdBs, int IdBuilding) {
        BuildingtoBS.insert(std::pair<int, int>(IdBuilding, IdBs));
    }

    //вес базовой станции
    std::set <int> LinkToBuilding;

    BSonBuilding() {};

    BSonBuilding operator=(const PoligonBuilding& rv);

    bool operator < (const BSonBuilding& vc) const {
        return this->Weight < vc.Weight;
    }

    bool operator > (const BSonBuilding& vc) const {
        return this->Weight > vc.Weight;
    }
};

