#pragma once
#include "PoligonBuilding.h"
#include "set"
//класс в котором хранится информация о БС, запланированных на стройке на здании, дополнительно хранит в себе информацию о том с каких зданий берется покрытия, через указатель
class BSonBuilding :
    public PoligonBuilding
{
private:
    double BSWeight = 0;
public:

    static std::multimap <int, int> BuildingtoBS;
    //на вход поступает id кастомный
    void addWeight(double WeightIn, PoligonBuilding* Build_in, int IdBSCustom, int IdBuildCustom);

    void changeWeight(double WeightIn);

    void recalcWeight();

    void resetWeight();

    double getWeight();

    static void addLink(int IdBs, int IdBuilding) {
        BuildingtoBS.insert(std::pair<int, int>(IdBuilding, IdBs));
    }

    //вес базовой станции
    std::set <std::pair<double, PoligonBuilding*>> LinkToBuilding;

    BSonBuilding() {};

    BSonBuilding operator=(const PoligonBuilding& rv);

    bool operator < (const BSonBuilding& vc) const {
        return this->BSWeight < vc.BSWeight;
    }

    bool operator > (const BSonBuilding& vc) const {
        return this->BSWeight > vc.BSWeight;
    }
};

