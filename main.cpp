#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include "ogrsf_frmts.h"
#include <omp.h>
#include <ctime>

class PoligonAsset {
    static int ID_poly_gen;
public:
    int             ID_poly;
    double          RxLev;
    OGRPolygon      Polygon;

    static int get_id() { return ID_poly_gen++; }

    static int get_now_id() { return ID_poly_gen; }

    PoligonAsset(double RxLev, OGRPolygon* Polygon) {
        this->ID_poly = get_id();
        this->RxLev = RxLev;
        this->Polygon = *Polygon;
    };

    PoligonAsset(double RxLev, OGRPolygon* Polygon, int id) {
        this->ID_poly = id;
        this->RxLev = RxLev;
        this->Polygon = *Polygon;
    };

    PoligonAsset(double RxLev, OGRPolygon Polygon) {
        this->ID_poly = get_id();
        this->RxLev = RxLev;
        this->Polygon = Polygon;
    };

    PoligonAsset(double RxLev, OGRPolygon Polygon, int id) {
        this->ID_poly = id;
        this->RxLev = RxLev;
        this->Polygon = Polygon;
    };

    PoligonAsset(double RxLev, OGRGeometry* Geometry) {
        this->ID_poly = get_id();
        this->RxLev = RxLev;
        this->Polygon = *(OGRPolygon*)Geometry;
    };

    PoligonAsset(double RxLev, OGRGeometry* Geometry, int id) {
        this->ID_poly = id;
        this->RxLev = RxLev;
        this->Polygon = *(OGRPolygon*)Geometry;
    };
    
};

int PoligonAsset::ID_poly_gen = 0;

class PoligonBuilding {
    static int ID_gen;
    
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
    OGRPolygon    Polygon;
    std::map<double, double> LevelArea;

    static int get_id() { return ID_gen++; }

    void addArea(double RxLev, double Area) {
        if (LevelArea.find(RxLev) == LevelArea.end()) {
            LevelArea[RxLev] = Area;
        }
        else {
            double buff = LevelArea[RxLev] + Area;
            LevelArea[RxLev] = buff;
        }
    }

    PoligonBuilding(
                    int         ID_build,
                    char*       Type,
                    char*       Wall_Mat,
                    double      Building,
                    int         Etage,
                    double      Area,
                    int         People_D,
                    int         People_N,
                    char*       Volcano_type,
                    int         Clutter_num,
                    char*       Clutter_type,
                    double      Lon,
                    double      Lat,
                    char*       Address_NP,
                    char*       Address_street,
                    char*       Address_number,
                    char*       FIAS,
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
    };
    
};

int PoligonBuilding::ID_gen = 0;



int main() {

    setlocale(LC_ALL, "Russian");
    double start = clock();
    std::map<int, PoligonAsset> AssetDate;
    std::map<int, PoligonBuilding> BuildingDate;
    int max_thread = omp_get_max_threads() - 1;
    omp_set_num_threads(max_thread);
    int Asset_step_m = 10;

    std::ofstream   write_kml_poly, write_kml_build, linkfile, polyfile, buildfile;

    write_kml_poly.open("respoly.kml");
    write_kml_build.open("resbuild.kml");
    linkfile.open("link.csv");
    polyfile.open("poly.csv");
    buildfile.open("build.csv");

    double rxlev;
    
    GDALAllRegister();

    GDALDataset* poDS;

    poDS = (GDALDataset*)GDALOpenEx("1.MIF", GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (poDS == NULL)
    {
        std::cout << "Open failed" << std::endl;
        return 0;
    }

    poDS->ResetReading();

    std::cout << "Всего слоев: " << poDS->GetLayerCount() << std::endl;
    
    for (auto i = 0; i < poDS->GetLayerCount(); i++) {

        OGRLayer* poLayer = poDS->GetLayer(i);
        OGRFeature* poFeature;
        OGRFeatureDefn* poFDefn = poLayer->GetLayerDefn();


        while ((poFeature = poLayer->GetNextFeature()) != NULL) {

            for (int iField = 0; iField < poFDefn->GetFieldCount(); iField++)
            {
                OGRFieldDefn* poFieldDefn = poFDefn->GetFieldDefn(iField);
               
                switch (poFieldDefn->GetType())
                {
                case OFTReal:
                    rxlev = poFeature->GetFieldAsDouble(iField);
                    break;
                }
            }

            OGRGeometry* poGeometry = poFeature->GetGeometryRef();
            
            if (poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon) {

                OGRMultiPolygon* Mpol = (OGRMultiPolygon*)poGeometry;

                std::cout << "Rxlev: " << rxlev;
                std::cout << " Total Area: " << Mpol->get_Area() << std::endl;

                for (auto j = 0; j < Mpol->getNumGeometries(); j++) {

                    if (Mpol->getGeometryRef(j) != NULL && 
                        wkbFlatten(Mpol->getGeometryRef(j)->getGeometryType()) == wkbPolygon) {

                        int id = PoligonAsset::get_id();

                        AssetDate.insert(std::pair<int, PoligonAsset>(id, PoligonAsset(rxlev, Mpol->getGeometryRef(j), id)));
                    }
                    else {
                        std::cout << "Неопознанный тип геометрии" << std::endl;
                        std::cout << Mpol->getGeometryRef(j)->getGeometryName() << std::endl;

                    }
                }

            }
            
            OGRFeature::DestroyFeature(poFeature);

        }
    }
   
    GDALClose(poDS);
   
    
    poDS = (GDALDataset*)GDALOpenEx("B.MIF", GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (poDS == NULL)
    {
        std::cout << "Open failed" << std::endl;
        return 0;
    }
    
    poDS->ResetReading();

    std::cout << "Всего слоев: " << poDS->GetLayerCount() << std::endl;

    for (auto i = 0; i < poDS->GetLayerCount(); i++) {
        
        OGRLayer* poLayer = poDS->GetLayer(i);
        OGRFeature* poFeature;
        OGRFeatureDefn* poFDefn = poLayer->GetLayerDefn();

        while ((poFeature = poLayer->GetNextFeature()) != NULL) {
            
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

            for (int iField = 0; iField < poFDefn->GetFieldCount(); iField++)
            {
                OGRFieldDefn* poFieldDefn = poFDefn->GetFieldDefn(iField);
               
                if (strcmp(poFieldDefn->GetNameRef(),"ID_build")==0) {
                    ID_build = poFeature->GetFieldAsInteger(iField);
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Type") == 0) {
                    strcpy_s(Type, poFeature->GetFieldAsString(iField));
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Wall_Mat") == 0) {
                    strcpy_s(Wall_Mat, poFeature->GetFieldAsString(iField));
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Building") == 0) {
                    Building = poFeature->GetFieldAsDouble(iField);
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Etage") == 0) {
                    Etage = poFeature->GetFieldAsInteger(iField);
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Area") == 0) {
                    Area = poFeature->GetFieldAsDouble(iField);
                }

                if (strcmp(poFieldDefn->GetNameRef(), "People_D") == 0) {
                    People_D = poFeature->GetFieldAsInteger(iField);
                }

                if (strcmp(poFieldDefn->GetNameRef(), "People_N") == 0) {
                    People_N = poFeature->GetFieldAsInteger(iField);
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Volcano_type") == 0) {
                    strcpy_s(Volcano_type, poFeature->GetFieldAsString(iField));
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Clutter_num") == 0) {
                    Clutter_num = poFeature->GetFieldAsInteger(iField);
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Clutter_type") == 0) {
                    strcpy_s(Clutter_type, poFeature->GetFieldAsString(iField));
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Долгота") == 0) {
                    Lon = poFeature->GetFieldAsDouble(iField);
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Широта") == 0) {
                    Lat = poFeature->GetFieldAsDouble(iField);
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Address_NP") == 0) {
                    strcpy_s(Address_NP, poFeature->GetFieldAsString(iField));
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Address_street") == 0) {
                    strcpy_s(Address_street, poFeature->GetFieldAsString(iField));
                }

                if (strcmp(poFieldDefn->GetNameRef(), "Address_number") == 0) {
                    strcpy_s(Address_number, poFeature->GetFieldAsString(iField));
                }

                if (strcmp(poFieldDefn->GetNameRef(), "FIAS") == 0) {
                    strcpy_s(FIAS, poFeature->GetFieldAsString(iField));
                }

            }
            
            
            if (poFeature->GetGeometryRef() != NULL && wkbFlatten(poFeature->GetGeometryRef()->getGeometryType()) == wkbPolygon) {
                
                int id = PoligonBuilding::get_id();
                
                BuildingDate.insert(std::pair<int, PoligonBuilding> (id, PoligonBuilding(ID_build, Type, Wall_Mat, Building, Etage, Area, People_D, People_N, Volcano_type, Clutter_num, Clutter_type, Lon, Lat, Address_NP, Address_street, Address_number, FIAS, poFeature->GetGeometryRef())));
                
            }
            
            OGRFeature::DestroyFeature(poFeature);
            
        }

    }


    GDALClose(poDS);
    
    //вот тут начинается работа
    
    std::cout << "Всего полигонов из Ассета: ";
    std::cout << AssetDate.size() << std::endl;
    
    std::cout << "Всего зданий: ";
    std::cout << BuildingDate.size() << std::endl;

    std::multimap <int, int> LinkPolyToBuild;
    
#pragma omp parallel for 
    for (int i = 0; i < PoligonAsset::get_now_id(); i++){

        for (auto iBuildingDate = BuildingDate.begin(); iBuildingDate != BuildingDate.end(); iBuildingDate++) {
            
            auto poly = AssetDate.find(i);
            
            if ((*iBuildingDate).second.Polygon.Contains(&((*poly).second.Polygon)) || (*iBuildingDate).second.Polygon.Overlaps(&((*poly).second.Polygon)) || (*poly).second.Polygon.Contains(&((*iBuildingDate).second.Polygon))) {
            
#pragma omp critical
                {
                    LinkPolyToBuild.insert(std::pair<int, int>(i, (*iBuildingDate).first));
                }
            }
        }
    }
    
    std::cout << "Всего созданных линков: ";
    std::cout << LinkPolyToBuild.size() << std::endl;

    for (auto i = LinkPolyToBuild.begin(); i != LinkPolyToBuild.end(); i++) {
        linkfile << (*i).first << ";" << (*i).second << std::endl;
    }


    linkfile << "Poly_id;Building_id" << std::endl;
//#pragma omp parallel for
    for (int i = 0; i < PoligonAsset::get_now_id(); i++) {

        if (LinkPolyToBuild.count(i) == 1) {
            auto PolyId = LinkPolyToBuild.find(i);
            auto Build = BuildingDate.find((*PolyId).second);
            auto PolyAss = AssetDate.find(i);
//#pragma omp critical
            {
                (*Build).second.addArea((*PolyAss).second.RxLev, (*PolyAss).second.Polygon.get_Area());
            }
        }
        else {
            //если линкуемся к нескольким зданиям, то разбиваем полигон на части
            OGRPolygon Polygon = (*(AssetDate.find(i))).second.Polygon;
            double RxLev = (*(AssetDate.find(i))).second.RxLev;

            double MaxX = NULL;
            double MaxY = NULL;

            double MinX = NULL;
            double MinY = NULL;

            OGRLineString* Border;
            OGRMultiLineString* Borders;


            switch (wkbFlatten(Polygon.getBoundary()->getGeometryType())) {
            case wkbLineString:

                Border = (OGRLineString*) Polygon.getBoundary();

                for (int k = 0; k < Border->getNumPoints(); k++) {
                    if (MaxX < Border->getX(k) || MaxX == NULL) {
                        MaxX = Border->getX(k);
                    }
                    if (MinX > Border->getX(k) || MinX == NULL) {
                        MinX = Border->getX(k);
                    }
                    if (MaxY < Border->getY(k) || MaxY == NULL) {
                        MaxY = Border->getY(k);
                    }
                    if (MinY > Border->getY(k) || MinY == NULL) {
                        MinY = Border->getY(k);
                    }
                }

                break;

            case wkbMultiLineString:

                Borders = (OGRMultiLineString*)Polygon.getBoundary();

                for (int l = 0; l < Borders->getNumGeometries(); l++) {
                    switch (wkbFlatten(Borders->getGeometryRef(l)->getGeometryType())) {
                    case wkbLineString:

                        Border = (OGRLineString*)Borders->getGeometryRef(l);

                        for (int k = 0; k < Border->getNumPoints(); k++) {
                            if (MaxX < Border->getX(k) || MaxX == NULL) {
                                MaxX = Border->getX(k);
                            }
                            if (MinX > Border->getX(k) || MinX == NULL) {
                                MinX = Border->getX(k);
                            }
                            if (MaxY < Border->getY(k) || MaxY == NULL) {
                                MaxY = Border->getY(k);
                            }
                            if (MinY > Border->getY(k) || MinY == NULL) {
                                MinY = Border->getY(k);
                            }
                        }
                        break;
                    default:
//#pragma omp critical
                        {
                            std::cout << "Неопознанный тип геометрии" << std::endl;
                            std::cout << Borders->getGeometryRef(l)->getGeometryName() << std::endl;
                        }
                        break;
                    }
                }

                break;
            default:
//#pragma omp critical
                {
                    std::cout << "Неопознанный тип геометрии" << std::endl;
                    std::cout << Polygon.getBoundary()->getGeometryName() << std::endl;
                }
                break;
            }

            for (double x = MinX; x < MaxX; x += Asset_step_m) {
                for (double y = MinY; y < MaxY; y += Asset_step_m) {

                    OGRLinearRing buffBorder;

                    buffBorder.assignSpatialReference(Polygon.getSpatialReference());

                    buffBorder.addPoint(x, y);
                    buffBorder.addPoint(x, y + Asset_step_m);
                    buffBorder.addPoint(x + Asset_step_m, y + Asset_step_m);
                    buffBorder.addPoint(x + Asset_step_m, y);
                    buffBorder.addPoint(x, y);

                    OGRPolygon buff;
                    buff.addRing(buffBorder.toCurve());
                    buff.assignSpatialReference(Polygon.getSpatialReference());
                  

                    //если вырезанный полигон попадает в настоящий
                    if (Polygon.Contains(&buff))
                    {

                        std::map<double, int> AreaBuildCover;

                        for (auto j = LinkPolyToBuild.lower_bound(i); j != LinkPolyToBuild.upper_bound(i); j++) {

                            double S = Asset_step_m * Asset_step_m;

                            auto BuildPoly = (*(BuildingDate.find((*j).second))).second.Polygon;
                            int BuildID = (*(BuildingDate.find((*j).second))).first;

                            if (BuildPoly.Contains(&buff) || buff.Contains(&BuildPoly)) {
                                AreaBuildCover[S] = BuildID;
                                continue;
                            }

                            if (BuildPoly.Overlaps(&buff)) {
                                try {
                                    OGRGeometry* bb = BuildPoly.Intersection(&buff);

                                    if (bb != NULL) {
                                        switch (wkbFlatten(bb->getGeometryType())) {
                                        case wkbPolygon:
                                            S = (double)((OGRPolygon*)bb)->get_Area() / S;
                                            break;
                                        case wkbMultiPolygon:
                                            S = (double)((OGRMultiPolygon*)bb)->get_Area() / S;
                                            break;
                                        }
                                    }
                                }
                                catch (...) {

                                }

                                AreaBuildCover[S] = BuildID;
                                continue;
                            }

                            if (BuildPoly.Touches(&buff)) {
                                AreaBuildCover[-1] = BuildID;
                            }
                        }

                        if (AreaBuildCover.size() == 1) {
//#pragma omp critical
                            {
                                (*(BuildingDate.find((*(AreaBuildCover.begin())).second))).second.addArea(RxLev, Asset_step_m * Asset_step_m);
                            }
                        }

                        if (AreaBuildCover.size() > 1) {
                            
                            auto j = AreaBuildCover.end();
                            j--;
//#pragma omp critical
                            {
                                (*(BuildingDate.find((*j).second))).second.addArea(RxLev, Asset_step_m * Asset_step_m);
                            }
                            
                        }

                        if (AreaBuildCover.size() == 0) {
                            std::cout << "ALARM" << std::endl;
                            return 0;
                        }
                    }

                }
            }


        }

    }
   
    polyfile << "Poly_id;RxLev;Area" << std::endl;
    
    for (auto iAssetDate = AssetDate.begin(); iAssetDate != AssetDate.end(); iAssetDate++) {
        polyfile << (*iAssetDate).first << ";" << (*iAssetDate).second.RxLev << ";" << (*iAssetDate).second.Polygon.get_Area() << std::endl;
    }
    
    buildfile << "Building_id_cust;Building_id;RxLev_100;RxLev_113;RxLev_120" << std::endl;
    
    for (auto iBuildingDate = BuildingDate.begin(); iBuildingDate != BuildingDate.end(); iBuildingDate++) {
        buildfile << (*iBuildingDate).first << ";" << (*iBuildingDate).second.ID_build;
        if ((*iBuildingDate).second.LevelArea.find(-100) != (*iBuildingDate).second.LevelArea.end())
        {
            buildfile << ";" << (*iBuildingDate).second.LevelArea[-100];
        }
        else {
            buildfile << ";";
        }

        if ((*iBuildingDate).second.LevelArea.find(-113) != (*iBuildingDate).second.LevelArea.end())
        {
            buildfile << ";" << (*iBuildingDate).second.LevelArea[-113];
        }
        else {
            buildfile << ";";
        }

        if ((*iBuildingDate).second.LevelArea.find(-120) != (*iBuildingDate).second.LevelArea.end())
        {
            buildfile << ";" << (*iBuildingDate).second.LevelArea[-120];
        }
        else {
            buildfile << ";";
        }

        buildfile << std::endl;
    }
    
    
    OGRSpatialReference srTo;

    srTo.SetGeogCS("My geographic CRS",
       "World Geodetic System 1984",
       "My WGS84 Spheroid",
       SRS_WGS84_SEMIMAJOR, SRS_WGS84_INVFLATTENING,
       "Greenwich", 0.0,
       "degree", 0.0174532925199433);
   
   
    write_kml_poly << "<?xml version = \"1.0\" encoding = \"UTF-8\"?>";
    write_kml_build << "<?xml version = \"1.0\" encoding = \"UTF-8\"?>";

    write_kml_poly << "<kml xmlns = \"http://www.opengis.net/kml/2.2\" xmlns:gx = \"http://www.google.com/kml/ext/2.2\" xmlns:kml = \"http://www.opengis.net/kml/2.2\" xmlns:atom = \"http://www.w3.org/2005/Atom\"><Document>";
    write_kml_build << "<kml xmlns = \"http://www.opengis.net/kml/2.2\" xmlns:gx = \"http://www.google.com/kml/ext/2.2\" xmlns:kml = \"http://www.opengis.net/kml/2.2\" xmlns:atom = \"http://www.w3.org/2005/Atom\"><Document>";

    for (auto iAssetDate = AssetDate.begin(); iAssetDate != AssetDate.end(); iAssetDate++) {
        (*iAssetDate).second.Polygon.transformTo(&srTo);
        write_kml_poly << "<Folder><name>" << (*iAssetDate).second.ID_poly << "</name>" << std::endl;
        write_kml_poly << "<Placemark>";
        (*iAssetDate).second.Polygon.swapXY();
        write_kml_poly << (*iAssetDate).second.Polygon.exportToKML() << std::endl;
        write_kml_poly << "</Placemark>";
        write_kml_poly << "</Folder>";
    }

    for (auto iBuildingDate = BuildingDate.begin(); iBuildingDate != BuildingDate.end(); iBuildingDate++) {
        (*iBuildingDate).second.Polygon.transformTo(&srTo);
        write_kml_build << "<Folder><name>" << (*iBuildingDate).first << "</name>" << std::endl;
        write_kml_build << "<Placemark>";
        (*iBuildingDate).second.Polygon.swapXY();
        write_kml_build << (*iBuildingDate).second.Polygon.exportToKML() << std::endl;
        write_kml_build << "</Placemark>";
        write_kml_build << "</Folder>";
    }

    write_kml_poly << "</Document></kml>";
    write_kml_build << "</Document></kml>";
    
    write_kml_poly.close();
    write_kml_build.close();
    linkfile.close();
    polyfile.close();
    buildfile.close();

    std::cout << "Time in seconds: " << ((clock() - start) / CLOCKS_PER_SEC) << std::endl;

    system("pause");

    return 0;
}
