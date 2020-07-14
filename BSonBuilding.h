#pragma once
#include "PoligonBuilding.h"
#include "set"
class BSonBuilding :
    public PoligonBuilding
{
private:
    double BSWeight = 0;
public:

    static std::multimap <int, int> BuildingtoBS;
    //�� ���� ��������� id ���������
    void addWeight(double WeightIn, PoligonBuilding* Build_in, int IdBSCustom, int IdBuildCustom);

    void changeWeight(double WeightIn);

    double getWeight();

    static void addLink(int IdBs, int IdBuilding) {
        BuildingtoBS.insert(std::pair<int, int>(IdBuilding, IdBs));
    }

    //��� ������� �������
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

