#pragma once
#include "PoligonBuilding.h"
#include "BsLinkToBuild.h"
#include "BS.h"
#include "set"
//����� � ������� �������� ���������� � ��, ��������������� �� ������� �� ������, ������������� ������ � ���� ���������� � ��� � ����� ������ ������� ��������, ����� ���������
class BSonBuilding :
    public BS
{
private:
    double Weight = 0;

    static std::multimap <int, int> BuildingtoBS;

    static void addLink(int IdBs, int IdBuilding) {
        BuildingtoBS.insert(std::pair<int, int>(IdBuilding, IdBs));
    }
    


public:

    PoligonBuilding setupBuild;

    void addLink(PoligonBuilding* Building);

    void recalcWeight();

    void resetWeight();

    double getWeight();

    //����� �� ������� �������
    std::vector <BsLinkToBuild> LinkToBuilding;

    BSonBuilding();

    BSonBuilding(PoligonBuilding setupBuild);
     
    bool operator < (const BSonBuilding& vc) const;

    bool operator > (const BSonBuilding& vc) const;
};

