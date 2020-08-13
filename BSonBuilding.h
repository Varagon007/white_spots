#pragma once
#include "PoligonBuilding.h"
#include "BsLinkToBuild.h"
#include "BS.h"
#include "set"
//класс в котором хранится информация о БС, запланированных на стройке на здании, дополнительно хранит в себе информацию о том с каких зданий берется покрытия, через указатель
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

    //линки до базовой станции
    std::vector <BsLinkToBuild> LinkToBuilding;

    BSonBuilding();

    BSonBuilding(PoligonBuilding setupBuild);
     
    bool operator < (const BSonBuilding& vc) const;

    bool operator > (const BSonBuilding& vc) const;
};

