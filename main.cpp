#include <iostream>
#include <fstream>
#include <string.h>
#include <map>
#include <vector>
#include "ogrsf_frmts.h"
#include <omp.h>
#include <set>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <PoligonCover.h>
#include <PoligonBuilding.h>


void Read_Building(const char* InputFile,std::vector<PoligonBuilding> &BuildingDate) {
    
    GDALDataset* poDS;

    poDS = (GDALDataset*)GDALOpenEx(InputFile, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (poDS == NULL)
    {
        throw "Файл не найден";
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

                if (strcmp(poFieldDefn->GetNameRef(), "ID_build") == 0) {
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

                BuildingDate.push_back(PoligonBuilding(ID_build, Type, Wall_Mat, Building, Etage, Area, People_D, People_N, Volcano_type, Clutter_num, Clutter_type, Lon, Lat, Address_NP, Address_street, Address_number, FIAS, poFeature->GetGeometryRef()));

            }

            OGRFeature::DestroyFeature(poFeature);

        }

    }


    GDALClose(poDS);
}

void Read_Cover(const char* InputFile, std::vector<PoligonCover>& CoverDate) {

    GDALDataset* poDS;
    double rxlev;

    poDS = (GDALDataset*)GDALOpenEx(InputFile, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (poDS == NULL)
    {
        throw "Файл не найден";
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

                //std::cout << "Rxlev: " << rxlev;
                //std::cout << " Total Area: " << Mpol->get_Area() << std::endl;

                for (auto j = 0; j < Mpol->getNumGeometries(); j++) {

                    if (Mpol->getGeometryRef(j) != NULL &&
                        wkbFlatten(Mpol->getGeometryRef(j)->getGeometryType()) == wkbPolygon) {

                        int id = PoligonCover::get_id();

                        CoverDate.push_back(PoligonCover(rxlev, Mpol->getGeometryRef(j), id));
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
}

void AddCoverToBuilding(int Floor , int Asset_step_m, std::vector<PoligonBuilding> &BuildingDate, std::vector<PoligonCover> &CoverDate, std::multimap <int, int> &LinkPolyToBuild) {
    for (int i = 0; i < CoverDate.size(); i++) {

        if (LinkPolyToBuild.count(i) == 1) {
            auto PolyId = LinkPolyToBuild.find(i);
            BuildingDate[(*PolyId).second].addArea(CoverDate[i].RxLev, CoverDate[i].Polygon.get_Area(), Floor);
        }
        else {
            //если линкуемся к нескольким зданиям, то разбиваем полигон на части

            OGRPolygon Polygon = CoverDate[i].Polygon;
            double RxLev = CoverDate[i].RxLev;

            double MaxX = NULL;
            double MaxY = NULL;

            double MinX = NULL;
            double MinY = NULL;

            OGRLineString* Border;
            OGRMultiLineString* Borders;


            switch (wkbFlatten(Polygon.getBoundary()->getGeometryType())) {
            case wkbLineString:

                Border = (OGRLineString*)Polygon.getBoundary();

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
                        std::cout << "Неопознанный тип геометрии" << std::endl;
                        std::cout << Borders->getGeometryRef(l)->getGeometryName() << std::endl;
                        break;
                    }
                }

                break;
            default:
                std::cout << "Неопознанный тип геометрии" << std::endl;
                std::cout << Polygon.getBoundary()->getGeometryName() << std::endl;
                break;
            }

#pragma omp parallel for 
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

                            auto BuildPoly = BuildingDate[(*j).second].Polygon;
                            int BuildID = (*j).second;

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
#pragma omp critical
                            {
                                BuildingDate[(*(AreaBuildCover.begin())).second].addArea(RxLev, Asset_step_m * Asset_step_m, Floor);
                            }
                        }
                        if (AreaBuildCover.size() > 1) {
#pragma omp critical
                            {
                                auto j = AreaBuildCover.end();
                                j--;
                                BuildingDate[(*j).second].addArea(RxLev, Asset_step_m * Asset_step_m, Floor); 
                            }

                        }

                        if (AreaBuildCover.size() == 0) {
                            std::cout << "ALARM" << std::endl;
                            throw "ALARM";
                        }
                    }

                }
            }


        }

    }

}

void OutKml(std::vector<PoligonCover> &CoverDate,std::string FileName) {
    
    OGRSpatialReference srTo;

    srTo.SetGeogCS("My geographic CRS",
        "World Geodetic System 1984",
        "My WGS84 Spheroid",
        SRS_WGS84_SEMIMAJOR, SRS_WGS84_INVFLATTENING,
        "Greenwich", 0.0,
        "degree", 0.0174532925199433);

    std::ofstream   write;

    std::string file = FileName;
    file += ".kml";

    write.open(file);

    write << "<?xml version = \"1.0\" encoding = \"UTF-8\"?>";
    write << "<kml xmlns = \"http://www.opengis.net/kml/2.2\" xmlns:gx = \"http://www.google.com/kml/ext/2.2\" xmlns:kml = \"http://www.opengis.net/kml/2.2\" xmlns:atom = \"http://www.w3.org/2005/Atom\"><Document>";

    for (size_t i = 0; i < CoverDate.size(); i++) {
        CoverDate[i].Polygon.transformTo(&srTo);
        write << "<Folder><name>" << CoverDate[i].ID_poly << "</name>" << std::endl;
        write << "<Placemark>";
        CoverDate[i].Polygon.swapXY();
        write << CoverDate[i].Polygon.exportToKML() << std::endl;
        write << "</Placemark>";
        write << "</Folder>";
    }

    write << "</Document></kml>";

    write.close();
}

void OutKml(std::vector<PoligonBuilding> &BuildingDate, const char* FileName) {
    
    if (BuildingDate.size() == 0) {
        return;
    }

    OGRSpatialReference srTo, basic;

    srTo.SetGeogCS("My geographic CRS",
        "World Geodetic System 1984",
        "My WGS84 Spheroid",
        SRS_WGS84_SEMIMAJOR, SRS_WGS84_INVFLATTENING,
        "Greenwich", 0.0,
        "degree", 0.0174532925199433);

    std::ofstream   write;

    std::string file = FileName;
    file += ".kml";

    write.open(file);

    write << "<?xml version = \"1.0\" encoding = \"UTF-8\"?>";
    write << "<kml xmlns = \"http://www.opengis.net/kml/2.2\" xmlns:gx = \"http://www.google.com/kml/ext/2.2\" xmlns:kml = \"http://www.opengis.net/kml/2.2\" xmlns:atom = \"http://www.w3.org/2005/Atom\"><Document>";

    for (size_t i = 0; i < BuildingDate.size(); i++) {
        basic.SetGeocCS(BuildingDate[i].Polygon.getSpatialReference()->GetName());
        BuildingDate[i].Polygon.transformTo(&srTo);
        write << "<Folder><name>" << i << "</name>" << std::endl;
        write << "<Placemark>";
        BuildingDate[i].Polygon.swapXY();
        write << BuildingDate[i].Polygon.exportToKML() << std::endl;
        write << "</Placemark>";
        write << "</Folder>";
        //BuildingDate[i].Polygon.transformTo(&basic);
    }

    write << "</Document></kml>";
    write.close();
}

void calcBadAreaByCover(std::vector<PoligonBuilding> &BuildingDate, const char* filename, int Floor, int Asset_step_m = 10) {

    std::vector<PoligonCover> CoverDate;
    std::ofstream   linkfile, polyfile;

    std::string buff = filename;

    linkfile.open(buff + "_link.csv");
    polyfile.open(buff + ".csv");

    try {
        Read_Cover(filename, CoverDate);
    }
    catch (const char* a) {
        std::cout << a << std::endl;
    }

    //вот тут начинается работа

    std::cout << "Всего полигонов на " << Floor << " этаже: ";
    std::cout << CoverDate.size() << std::endl;

    polyfile << "Poly_id;RxLev;Area" << std::endl;

    for (size_t i = 0; i < CoverDate.size(); i++) {
        polyfile << i << ";" << CoverDate[i].RxLev << ";" << CoverDate[i].Polygon.get_Area() << std::endl;
    }

    polyfile.close();

    std::multimap <int, int> LinkPolyToBuild;

#pragma omp parallel for 
    for (int i = 0; i < CoverDate.size(); i++) {

        for (int j = 0; j < BuildingDate.size(); j++) {

            if (BuildingDate[j].Polygon.Intersect(&(CoverDate[i].Polygon)) || BuildingDate[j].Polygon.Contains(&(CoverDate[i].Polygon)) || BuildingDate[j].Polygon.Overlaps(&(CoverDate[i].Polygon)) || CoverDate[i].Polygon.Contains(&(BuildingDate[j].Polygon)) || CoverDate[i].Polygon.Overlaps(&(BuildingDate[j].Polygon))) {

#pragma omp critical
                {
                    LinkPolyToBuild.insert(std::pair<int, int>(i, j));
                }
            }
        }
    }

    std::cout << "Всего созданных линков на " << Floor << " этаже: ";
    std::cout << LinkPolyToBuild.size() << std::endl;

    linkfile << "Poly_id;Building_id" << std::endl;

    for (auto i = LinkPolyToBuild.begin(); i != LinkPolyToBuild.end(); i++) {
        linkfile << (*i).first << ";" << (*i).second << std::endl;
    }

    linkfile.close();

    try {
        AddCoverToBuilding(Floor, Asset_step_m, BuildingDate, CoverDate, LinkPolyToBuild);
        std::string Out_file;
        Out_file = "poly_floor_" + std::to_string(Floor);
        OutKml(CoverDate, Out_file);
        CoverDate.erase(CoverDate.begin(), CoverDate.end());
        LinkPolyToBuild.erase(LinkPolyToBuild.begin(), LinkPolyToBuild.end());
    }
    catch (const char* a) {
        std::cout << a << std::endl;
    }

}

int main() {

    setlocale(LC_ALL, "Russian");
    time_t rawtime;
    struct tm* timeinfo;
    
    time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Старт: " << std::put_time(localtime(&now), "%F %T") << std::endl;

    double start = clock();
    

    std::vector<PoligonBuilding> BuildingDate;
    int max_thread = omp_get_max_threads() - 1;
    omp_set_num_threads(max_thread);
    int Asset_step_m = 10;
    
    std::ofstream   buildfile;

    buildfile.open("build.csv");

    GDALAllRegister();

    try {
        Read_Building("B.MIF", BuildingDate);
    }
    catch (const char* a) {
        std::cout << a << std::endl;
    }
    
    std::cout << "Всего зданий: ";
    std::cout << BuildingDate.size() << std::endl;


    calcBadAreaByCover(BuildingDate, "1.MIF", 1, Asset_step_m);
    calcBadAreaByCover(BuildingDate, "4.MIF", 4, Asset_step_m);
    calcBadAreaByCover(BuildingDate, "7.MIF", 7, Asset_step_m);


    buildfile << "Building_id_cust;Building_id;Lat;Lon;height";
    
    std::set<double> LevelList = PoligonBuilding::getLevel();

    for (auto iLevelList = LevelList.begin(); iLevelList != LevelList.end(); iLevelList++) {
        buildfile << ";RXLev_" << (*iLevelList);
    }
    
    buildfile << std::endl;

    for (size_t i = 0; i < BuildingDate.size(); i++) {

        buildfile << i << ";" << BuildingDate[i].ID_build << ";" << BuildingDate[i].Lat << ";" << BuildingDate[i].Lon << ";" << BuildingDate[i].Building;

        for (auto iLevelList = LevelList.begin(); iLevelList != LevelList.end(); iLevelList++) {

            if (BuildingDate[i].LevelArea.find((*iLevelList)) != BuildingDate[i].LevelArea.end())
            {
                buildfile << ";" << BuildingDate[i].LevelArea[(*iLevelList)];
            }
            else {
                buildfile << ";";
            }

        }
        buildfile << std::endl;
    }
    
    buildfile.close();

    std::cout << BuildingDate[0].Polygon.get_Area() << std::endl;

    OutKml(BuildingDate, "building");

    std::cout << BuildingDate[0].Polygon.get_Area() << std::endl;

    now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Финиш: " << std::put_time(localtime(&now), "%F %T") << std::endl;

    double tt = ((clock() - start) / CLOCKS_PER_SEC);
    int hh = trunc(tt / 3600);
    int mi = trunc((tt - hh * 3600)/60);
    int ss = tt - 3600 * hh - mi * 60;
    std::cout << "Затраченное время: " << hh << " часов " << mi << " минут " << ss << " секунд" << std::endl;
    std::cout << "Затраченное время (в секундах): " << tt << std::endl;

    system("pause");

    return 0;
}