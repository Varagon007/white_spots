#include "PoligonBuilding.h"

std::set<double> PoligonBuilding::LevelList{};

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
};
