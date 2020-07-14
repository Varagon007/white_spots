#include "PoligonBuilding.h"

std::set<double> PoligonBuilding::LevelList{};

int PoligonBuilding::ID_build_custom_gen = 0;

double PoligonBuilding::getAddArea(double RxLev, double Area, int k = 1) {
    double buff;
    if (LevelArea.find(RxLev) == LevelArea.end()) {
        buff = k * Area;
    }
    else {
        buff = LevelArea[RxLev] + k * Area;
    }
    return buff;
}

void PoligonBuilding::addArea(double RxLev, double Area, int Floor) {
    

    int FloorBuilding = trunc((Building - 1.5) / 3);

    switch (Floor) {
    case 1:
        switch (FloorBuilding) {
        case 1:
            LevelArea[RxLev] = getAddArea(RxLev, Area);
            break;
        case 2:
            LevelArea[RxLev] = getAddArea(RxLev, Area, 2);
            break;
        case 3:
            LevelArea[RxLev] = getAddArea(RxLev, Area, 3);
            break;
        default:
            LevelArea[RxLev] = getAddArea(RxLev, Area, 2);
            break;
        }
        break;
    case 4:
        if (FloorBuilding > 3) {
            switch (FloorBuilding) {
            case 4:
                LevelArea[RxLev] = getAddArea(RxLev, Area, 2);
                break;
            case 5:
                LevelArea[RxLev] = getAddArea(RxLev, Area, 3);
                break;
            case 6:
                LevelArea[RxLev] = getAddArea(RxLev, Area, 4);
                break;
            default:
                LevelArea[RxLev] = getAddArea(RxLev, Area, 3);
                break;

            }
        }
        break;
    case 7:
        if (FloorBuilding > 6) {
            switch (FloorBuilding) {
            case 7:
                LevelArea[RxLev] = getAddArea(RxLev, Area, 2);
                break;
            case 8:
                LevelArea[RxLev] = getAddArea(RxLev, Area, 3);
                break;
            case 9:
                LevelArea[RxLev] = getAddArea(RxLev, Area, 4);
                break;
            case 10:
                LevelArea[RxLev] = getAddArea(RxLev, Area, 5);
                break;
            default:
                LevelArea[RxLev] = getAddArea(RxLev, Area, 5);
                break;
            }
        }
        break;
    }

    PoligonBuilding::addLevel(RxLev);
}
void PoligonBuilding::createWeight()
{
    for (auto iLevelArea = LevelArea.begin(); iLevelArea != LevelArea.end(); iLevelArea++) {
        if ((*iLevelArea).first == -113){
            Weight += (*iLevelArea).second * 0.05;
        }
        if ((*iLevelArea).first == -120) {
            Weight += (*iLevelArea).second * 0.5;
        }
        if ((*iLevelArea).first == -200) {
            Weight += (*iLevelArea).second * 1;
        }
    }
}
void PoligonBuilding::changeWeight(double WeightIn)
{
    if (WeightIn < 0 && fabs(WeightIn) > Weight) {
        Weight = 0;
    }
    else {
        Weight += WeightIn;
    }
}
double PoligonBuilding::getWeight()
{
    return this->Weight;
}

