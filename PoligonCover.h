#pragma once

#include "ogrsf_frmts.h"

class PoligonCover
{
    static int ID_poly_gen;
public:
    int             ID_poly;
    double          RxLev;
    OGRPolygon      Polygon;

    static int get_id() { return ID_poly_gen++; }

    static int get_now_id() { return ID_poly_gen; }

    PoligonCover(double RxLev, OGRPolygon* Polygon) {
        this->ID_poly = get_id();
        this->RxLev = RxLev;
        this->Polygon = *Polygon;
    };

    PoligonCover(double RxLev, OGRPolygon* Polygon, int id) {
        this->ID_poly = id;
        this->RxLev = RxLev;
        this->Polygon = *Polygon;
    };

    PoligonCover(double RxLev, OGRPolygon Polygon) {
        this->ID_poly = get_id();
        this->RxLev = RxLev;
        this->Polygon = Polygon;
    };

    PoligonCover(double RxLev, OGRPolygon Polygon, int id) {
        this->ID_poly = id;
        this->RxLev = RxLev;
        this->Polygon = Polygon;
    };

    PoligonCover(double RxLev, OGRGeometry* Geometry) {
        this->ID_poly = get_id();
        this->RxLev = RxLev;
        this->Polygon = *(OGRPolygon*)Geometry;
    };

    PoligonCover(double RxLev, OGRGeometry* Geometry, int id) {
        this->ID_poly = id;
        this->RxLev = RxLev;
        this->Polygon = *(OGRPolygon*)Geometry;
    };

};

