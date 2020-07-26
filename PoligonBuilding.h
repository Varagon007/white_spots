#pragma once

#include "ogrsf_frmts.h"
#include <set>
#include <map>
//данный класс существует для хранения полной информации о здании, в том числе о том сколько в здании плохого покрытия
class PoligonBuilding
{
private:
    double Weight = 0;
    int ID_build_custom;
    static int ID_build_custom_gen;

protected:
    static std::set<double> LevelList;
    double getAddArea(double RxLev, double Area, int k);

    static void addLevel(double RxLev) {
        LevelList.insert(RxLev);
    };

public:
    
    int         ID_build;
    char        Type[50];
    char        Wall_Mat[50];
    double      Building;
    int         Etage;
    double      Area;
    int         People_D;
    int         People_N;
    char        Volcano_type[50];
    int         Clutter_num;
    char        Clutter_type[50];
    double      Lon;
    double      Lat;
    char        Address_NP[150];
    char        Address_street[100];
    char        Address_number[10];
    char        FIAS[200];
    OGRPolygon  Polygon;
    OGRPoint    PointCenter;

    std::map<double, double> LevelArea;
        
    static std::set<double> getLevel() {
        return LevelList;
    };
    
    void addArea(double RxLev, double Area, int Floor);

    void createWeight();

    void changeWeight(double WeightIn);

    static int getId();

    double getWeight();

    PoligonBuilding() {};

    PoligonBuilding(
        int         ID_build,
        char* Type,
        char* Wall_Mat,
        double      Building,
        int         Etage,
        double      Area,
        int         People_D,
        int         People_N,
        char* Volcano_type,
        int         Clutter_num,
        char* Clutter_type,
        double      Lon,
        double      Lat,
        char* Address_NP,
        char* Address_street,
        char* Address_number,
        char* FIAS,
        OGRPolygon* Polygon
    ) {
        this->ID_build = ID_build;
        strcpy_s(this->Type, Type);
        strcpy_s(this->Wall_Mat, Wall_Mat);
        this->Building = Building;
        this->Etage = Etage;
        this->Area = Area;
        this->People_D = People_D;
        this->People_N = People_N;
        strcpy_s(this->Volcano_type, Volcano_type);
        this->Clutter_num = Clutter_num;
        strcpy_s(this->Clutter_type, Clutter_type);
        this->Lon = Lon;
        this->Lat = Lat;
        strcpy_s(this->Address_NP, Address_NP);
        strcpy_s(this->Address_street, Address_street);
        strcpy_s(this->Address_number, Address_number);
        strcpy_s(this->FIAS, FIAS);
        this->Polygon = *Polygon;
        this->Polygon.Centroid(&this->PointCenter);
    };

    PoligonBuilding(int         ID_build,
        char* Type,
        char* Wall_Mat,
        double      Building,
        int         Etage,
        double      Area,
        int         People_D,
        int         People_N,
        char* Volcano_type,
        int         Clutter_num,
        char* Clutter_type,
        double      Lon,
        double      Lat,
        char* Address_NP,
        char* Address_street,
        char* Address_number,
        char* FIAS,
        OGRGeometry* Geometry) {
        this->ID_build = ID_build;
        strcpy_s(this->Type, Type);
        strcpy_s(this->Wall_Mat, Wall_Mat);
        this->Building = Building;
        this->Etage = Etage;
        this->Area = Area;
        this->People_D = People_D;
        this->People_N = People_N;
        strcpy_s(this->Volcano_type, Volcano_type);
        this->Clutter_num = Clutter_num;
        strcpy_s(this->Clutter_type, Clutter_type);
        this->Lon = Lon;
        this->Lat = Lat;
        strcpy_s(this->Address_NP, Address_NP);
        strcpy_s(this->Address_street, Address_street);
        strcpy_s(this->Address_number, Address_number);
        strcpy_s(this->FIAS, FIAS);
        this->Polygon = *(OGRPolygon*)Geometry;
        this->Polygon.Centroid(&this->PointCenter);
    };
};

