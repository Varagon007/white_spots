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
#include <algorithm>
#include <execution>
#include <PoligonCover.h>
#include <PoligonBuilding.h>
#include <BSonBuilding.h>
#include <CSVFile.h>
 
//функция считавыющая информацию о границах с MID/MIF файла
OGRPolygon Read_Border(const char* InputFile, OGRSpatialReference* RefGeo) {

    GDALDataset* poDS;

    poDS = (GDALDataset*)GDALOpenEx(InputFile, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (poDS == NULL)
    {
        throw "Файл не найден";
    }

    poDS->ResetReading();

    if (poDS->GetLayerCount() > 1) {
        throw "Несколько слоев в границе";
    }

    OGRLayer* poLayer = poDS->GetLayer(0);
    OGRFeature* poFeature;
    OGRPolygon res;

    while ((poFeature = poLayer->GetNextFeature()) != NULL) {


        if (poFeature->GetGeometryRef() != NULL && wkbFlatten(poFeature->GetGeometryRef()->getGeometryType()) == wkbPolygon) {

            res = *(OGRPolygon *) poFeature->GetGeometryRef();
            
            res.transformTo(RefGeo);
        }

        OGRFeature::DestroyFeature(poFeature);

    }

    GDALClose(poDS);

    return res;
}

//функция считавыющая информацию о здании с MID/MIF файла
void Read_Building(const char* InputFile,std::vector<PoligonBuilding> &BuildingDate) {
    
    GDALDataset* poDS;

    poDS = (GDALDataset*)GDALOpenEx(InputFile, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (poDS == NULL)
    {
        throw "Файл не найден";
    }

    poDS->ResetReading();

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
//функция считывающая информацию о покрытии
void Read_Cover(const char* InputFile, std::vector<PoligonCover>& CoverDate) {

    GDALDataset* poDS;
    double rxlev;

    poDS = (GDALDataset*)GDALOpenEx(InputFile, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (poDS == NULL)
    {
        throw "Файл не найден";
    }

    poDS->ResetReading();

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
//функция возвращающая границы в виде прямоугольника
void getBorder(OGRLineString* Border, double* minX, double* minY, double* maxX, double* maxY) {
    for (int k = 0; k < Border->getNumPoints(); k++) {
        if (*maxX < Border->getX(k) || *maxX == NULL) {
            *maxX = Border->getX(k);
        }
        if (*minX > Border->getX(k) || *minX == NULL) {
            *minX = Border->getX(k);
        }
        if (*maxY < Border->getY(k) || *maxY == NULL) {
            *maxY = Border->getY(k);
        }
        if (*minY > Border->getY(k) || *minY == NULL) {
            *minY = Border->getY(k);
        }
        
    }
}

//функция которая позволяет посчитать сколько в здании какого покрытия
void AddCoverToBuilding(int Floor , int Asset_step_m, std::vector<PoligonBuilding> &BuildingDate, std::vector<PoligonCover> &CoverDate, std::multimap <int, int> &LinkPolyToBuild) {
    
    std::ofstream errorLink;

    OGRSpatialReference srTo;

    srTo.SetGeogCS("My geographic CRS",
        "World Geodetic System 1984",
        "My WGS84 Spheroid",
        SRS_WGS84_SEMIMAJOR, SRS_WGS84_INVFLATTENING,
        "Greenwich", 0.0,
        "degree", 0.0174532925199433);

    errorLink.open("errorLink.txt", std::fstream::out | std::fstream::app);

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

            //OGRLineString* Border;
            OGRMultiLineString* Borders;


            switch (wkbFlatten(Polygon.getBoundary()->getGeometryType())) {
            case wkbLineString:

                 getBorder((OGRLineString*)Polygon.getBoundary(), &MinX, &MinY, &MaxX, &MaxY);

                break;

            case wkbMultiLineString:

                Borders = (OGRMultiLineString*)Polygon.getBoundary();

                for (int l = 0; l < Borders->getNumGeometries(); l++) {
                    switch (wkbFlatten(Borders->getGeometryRef(l)->getGeometryType())) {
                    case wkbLineString:

                        getBorder((OGRLineString*)Borders->getGeometryRef(l), &MinX, &MinY, &MaxX, &MaxY);

                        break;
                    default:
                        std::cout << "Неопознанный тип геометрии" << std::endl;
                        std::cout << Borders->getGeometryRef(l)->getGeometryName() << std::endl;

                        std::clog << "Неопознанный тип геометрии" << std::endl;
                        std::clog << Borders->getGeometryRef(l)->getGeometryName() << std::endl;
                        throw "Mistake";
                        break;
                    }
                }

                break;
            default:
                std::cout << "Неопознанный тип геометрии" << std::endl;
                std::cout << Polygon.getBoundary()->getGeometryName() << std::endl;

                std::clog << "Неопознанный тип геометрии" << std::endl;
                std::clog << Polygon.getBoundary()->getGeometryName() << std::endl;
                throw "Mistake";
                break;
            }

#pragma omp parallel for 
            for (int x = trunc(MinX); x < (int)ceil(MaxX); x += Asset_step_m) {
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

                            if (BuildPoly.Overlaps(&buff)|| buff.Overlaps(&BuildPoly)) {
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
#pragma omp critical
                            {

                                int reLink = -1;
                                double reDist = NULL;

                                for (auto j = LinkPolyToBuild.lower_bound(i); j != LinkPolyToBuild.upper_bound(i); j++) {

                                    double S = Asset_step_m * Asset_step_m;

                                    auto BuildPoly = BuildingDate[(*j).second].Polygon;
                                    
                                    int BuildID = (*j).second;

                                    if (reDist == NULL || buff.Distance(&BuildPoly) < reDist) {
                                        reDist = buff.Distance(&BuildPoly);
                                        reLink = BuildID;
                                    }
                                }

                                if (reLink != -1)
                                    BuildingDate[reLink].addArea(RxLev, Asset_step_m* Asset_step_m, Floor);

                                std::cout << "ALARM" << std::endl;
                                std::clog << "Присутствует непривязанный полигон" << std::endl;
                                errorLink << "<Placemark>";
                                buff.transformTo(&srTo);
                                buff.swapXY();
                                errorLink << buff.exportToKML();
                                errorLink << "</Placemark>" << std::endl;
                            }

                        }
                    }

                }
            }


        }

    }

    errorLink.close();

}
//построение KML c заранее известными цветами зданий
void OutKml(std::vector<PoligonCover> &CoverDateIn,std::string FileName) {
    
    std::vector<PoligonCover> CoverDate(CoverDateIn);

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
    write << "<Style id=\"NSff50b000\"><IconStyle><color>ff00ff00</color><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/polygon.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><colorMode>normal</colorMode><scale>0.8</scale></LabelStyle><LineStyle><color>ff50b000</color><colorMode>normal</colorMode><width>1</width></LineStyle><PolyStyle><color>ff50b000</color><outline>1</outline></PolyStyle></Style>"
        << "<Style id=\"HSff50b000\"><IconStyle><color>ff00ff00</color><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/polygon.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><colorMode>normal</colorMode><scale>0.8</scale></LabelStyle><LineStyle><color>ff00ffff</color><colorMode>normal</colorMode><width>2</width></LineStyle><PolyStyle><color>ff50b000</color><outline>1</outline></PolyStyle></Style>"
        << "<Style id=\"NSff00ffff\"><IconStyle><color>ff00ff00</color><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/polygon.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><colorMode>normal</colorMode><scale>0.8</scale></LabelStyle><LineStyle><color>ff00ffff</color><colorMode>normal</colorMode><width>1</width></LineStyle><PolyStyle><color>ff00ffff</color><outline>1</outline></PolyStyle></Style>"
        << "<Style id=\"HSff00ffff\"><IconStyle><color>ff00ff00</color><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/polygon.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><colorMode>normal</colorMode><scale>0.8</scale></LabelStyle><LineStyle><color>ff00ffff</color><colorMode>normal</colorMode><width>2</width></LineStyle><PolyStyle><color>ff00ffff</color><outline>1</outline></PolyStyle></Style>"
        << "<Style id=\"NSff262626\"><IconStyle><color>ff00ff00</color><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/polygon.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><colorMode>normal</colorMode><scale>0.8</scale></LabelStyle><LineStyle><color>ff262626</color><colorMode>normal</colorMode><width>1</width></LineStyle><PolyStyle><color>ff262626</color><outline>1</outline></PolyStyle></Style>"
        << "<Style id=\"HSff262626\"><IconStyle><color>ff00ff00</color><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/polygon.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><colorMode>normal</colorMode><scale>0.8</scale></LabelStyle><LineStyle><color>ff00ffff</color><colorMode>normal</colorMode><width>2</width></LineStyle><PolyStyle><color>ff262626</color><outline>1</outline></PolyStyle></Style>"
        << "<Style id=\"NSff0000ff\"><IconStyle><color>ff00ff00</color><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/polygon.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><colorMode>normal</colorMode><scale>0.8</scale></LabelStyle><LineStyle><color>ff0000ff</color><colorMode>normal</colorMode><width>1</width></LineStyle><PolyStyle><color>ff0000ff</color><outline>1</outline></PolyStyle></Style>"
        << "<Style id=\"HSff0000ff\"><IconStyle><color>ff00ff00</color><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/polygon.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><colorMode>normal</colorMode><scale>0.8</scale></LabelStyle><LineStyle><color>ff00ffff</color><colorMode>normal</colorMode><width>2</width></LineStyle><PolyStyle><color>ff0000ff</color><outline>1</outline></PolyStyle></Style>"
        << "<Style id=\"NSffffffff\"><IconStyle><color>ff00ff00</color><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/polygon.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><colorMode>normal</colorMode><scale>0.8</scale></LabelStyle><LineStyle><color>ffffffff</color><colorMode>normal</colorMode><width>1</width></LineStyle><PolyStyle><color>ffffffff</color><outline>1</outline></PolyStyle></Style>"
        << "<Style id=\"HSffffffff\"><IconStyle><color>ff00ff00</color><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/polygon.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><colorMode>normal</colorMode><scale>0.8</scale></LabelStyle><LineStyle><color>ff00ffff</color><colorMode>normal</colorMode><width>2</width></LineStyle><PolyStyle><color>ffffffff</color><outline>1</outline></PolyStyle></Style>"
        << "<StyleMap id=\"SMff50b000\"><Pair><key>normal</key><styleUrl>#NSff50b000</styleUrl></Pair><Pair><key>highlight</key><styleUrl>#HSff50b000</styleUrl></Pair></StyleMap>"
        << "<StyleMap id=\"SMff00ffff\"><Pair><key>normal</key><styleUrl>#NSff00ffff</styleUrl></Pair><Pair><key>highlight</key><styleUrl>#HSff00ffff</styleUrl></Pair></StyleMap>"
        << "<StyleMap id=\"SMff262626\"><Pair><key>normal</key><styleUrl>#NSff262626</styleUrl></Pair><Pair><key>highlight</key><styleUrl>#HSff262626</styleUrl></Pair></StyleMap>"
        << "<StyleMap id=\"SMff0000ff\"><Pair><key>normal</key><styleUrl>#NSff0000ff</styleUrl></Pair><Pair><key>highlight</key><styleUrl>#HSff0000ff</styleUrl></Pair></StyleMap>"
        << "<StyleMap id=\"SMffffffff\"><Pair><key>normal</key><styleUrl>#NSffffffff</styleUrl></Pair><Pair><key>highlight</key><styleUrl>#HSffffffff</styleUrl></Pair></StyleMap>" 
        << std::endl;

    for (size_t i = 0; i < CoverDate.size(); i++) {
        CoverDate[i].Polygon.transformTo(&srTo);
        write << "<Folder><name>" << CoverDate[i].ID_poly << "</name>" << std::endl;
        write << "<Placemark>";
        if (CoverDate[i].RxLev == -100) {
            write << "<styleUrl>#SMff50b000</styleUrl>";
        }
        if (CoverDate[i].RxLev == -113) {
            write << "<styleUrl>#SMff00ffff</styleUrl>";
        }
        if (CoverDate[i].RxLev == -120) {
            write << "<styleUrl>#SMff0000ff</styleUrl>";
        }
        if (CoverDate[i].RxLev == -200) {
            write << "<styleUrl>#SMff262626</styleUrl>";
        }
        CoverDate[i].Polygon.swapXY();
        write << CoverDate[i].Polygon.exportToKML() << std::endl;
        write << "</Placemark>";
        write << "</Folder>";
    }

    write << "</Document></kml>";

    write.close();
}
//построение kml по зданиям
void OutKml(std::vector<PoligonBuilding> &BuildingDateIn, const char* FileName) {
    
    if (BuildingDateIn.size() == 0) {
        return;
    }

    std::vector<PoligonBuilding> BuildingDate(BuildingDateIn);

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
//функция в которой маппится плохое покрытие к зданиям
void calcBadAreaByCover(std::vector<PoligonBuilding> &BuildingDate, const char* filename, int Floor, int Asset_step_m = 10) {


    time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    double start = clock();

    std::cout << "Загрузка порытия " << Floor << " этажа: " << std::put_time(localtime(&now), "%F %T") << std::endl;
    std::clog << "Загрузка порытия " << Floor << " этажа: " << std::put_time(localtime(&now), "%F %T") << std::endl;

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

    std::clog << "Всего полигонов на " << Floor << " этаже: ";
    std::clog << CoverDate.size() << std::endl;

    polyfile << "Poly_id;RxLev;Area" << std::endl;

    for (size_t i = 0; i < CoverDate.size(); i++) {
        polyfile << i << ";" << CoverDate[i].RxLev << ";" << CoverDate[i].Polygon.get_Area() << std::endl;
    }

    polyfile.close();

    std::multimap <int, int> LinkPolyToBuild;

#pragma omp parallel for 
    for (int i = 0; i < CoverDate.size(); i++) {

        for (int j = 0; j < BuildingDate.size(); j++) {

            if (BuildingDate[j].Polygon.Distance(&(CoverDate[i].Polygon)) < 1000) {

                if (BuildingDate[j].Polygon.Intersect(&(CoverDate[i].Polygon)) || BuildingDate[j].Polygon.Contains(&(CoverDate[i].Polygon)) || BuildingDate[j].Polygon.Overlaps(&(CoverDate[i].Polygon)) || CoverDate[i].Polygon.Contains(&(BuildingDate[j].Polygon)) || CoverDate[i].Polygon.Overlaps(&(BuildingDate[j].Polygon))) {

#pragma omp critical
                    {
                        LinkPolyToBuild.insert(std::pair<int, int>(i, j));
                    }
                }
            }
        }
    }

    std::cout << "Всего созданных линков на " << Floor << " этаже: ";
    std::cout << LinkPolyToBuild.size() << std::endl;

    std::clog << "Всего созданных линков на " << Floor << " этаже: ";
    std::clog << LinkPolyToBuild.size() << std::endl;

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

    now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Загрузка порытия " << Floor << " этажа (Завершено): " << std::put_time(localtime(&now), "%F %T") << std::endl;
    std::clog << "Загрузка порытия " << Floor << " этажа (Завершено): " << std::put_time(localtime(&now), "%F %T") << std::endl;

    double tt = ((clock() - start) / CLOCKS_PER_SEC);
    int hh = trunc(tt / 3600);
    int mi = trunc((tt - hh * 3600) / 60);
    int ss = tt - 3600 * hh - mi * 60;
    std::cout << "Затраченное время: " << hh << " часов " << mi << " минут " << ss << " секунд" << std::endl;
    std::cout << "Затраченное время (в секундах): " << tt << std::endl;

    std::clog << "Затраченное время: " << hh << " часов " << mi << " минут " << ss << " секунд" << std::endl;
    std::clog << "Затраченное время (в секундах): " << tt << std::endl;
}

void OutFile(std::vector<PoligonBuilding> &BuildingDate, std::string OutName = "build.csv") {

    std::ofstream buildfile;

    buildfile.open(OutName);

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
}

void OutFile(std::vector<BSonBuilding> &BSDate, std::string OutName = "bslist.csv") {

    std::ofstream bsfile;

    bsfile.open(OutName);

    bsfile << "ID_bs;Lat;Lon;Weight;Build_list" << std::endl;

    for (size_t i = 0; i < BSDate.size(); i++) {
        bsfile << BSDate[i].getBsId() << ";" << BSDate[i].setupBuild.Lat << ";" << BSDate[i].setupBuild.Lon << ";" << BSDate[i].getWeight();
        /*
        for (auto j = BSDate[i].LinkToBuilding.begin(); j != BSDate[i].LinkToBuilding.end(); j++) {
            bsfile << ";" << (*j).first << "/" << (*j).second->ID_build;
        }
        */
        bsfile << std::endl;
    }

    bsfile.close();

}
//возвращает true если файл есть в наличии
bool checkFile(std::string FileIn) {
    std::ifstream ifs;
    bool res = true;

    ifs.open(FileIn);

    if (!ifs.is_open()) {
        res = false;
    }

    ifs.close();

    return res;
}

void topByCircle(std::vector<BSonBuilding> &BSDate, std::vector<PoligonBuilding> &BuildingDate) {

    std::ofstream bsfile;
    size_t i = 0;

#pragma omp parallel for
    for (int i = 0; i < BSDate.size(); i++) {
        double R = BSDate[i].getR();
        for (size_t j = 0; j < BuildingDate.size(); j++) {
            double Rbuff = BSDate[i].Coordinate.Distance(&BuildingDate[j].PointCenter);
#pragma omp critical   
            {
                if (Rbuff <= 2 * R) {
                    BSDate[i].addLink(&BuildingDate[j]);
                }
            }

        }
        BSDate[i].recalcWeight();
    }

    OutFile(BSDate);

    std::cout << "Выборка топ: Старт" << std::endl;
    std::clog << "Выборка топ: Старт" << std::endl;

    bsfile.open("bslist_filter.csv");

    bsfile << "ID_bs;Lat;Lon;Weight;Rmax" << std::endl;

    for (auto iBSDate = BSDate.begin(); iBSDate != BSDate.end(); iBSDate++, i++) {

        std::cout << "Опредедение Топ " << i << std::endl;
        std::clog << "Опредедение Топ " << i << std::endl;

        auto TopI = std::max_element(std::execution::par, iBSDate, BSDate.end());

        std::iter_swap(TopI, iBSDate);

        bsfile << (*iBSDate).getBsId() << ";" << (*iBSDate).setupBuild.Lat << ";" << (*iBSDate).setupBuild.Lon << ";" << (*iBSDate).getWeight() << ";" << (*iBSDate).getR();

        bsfile << std::endl;

        if ((*iBSDate).getWeight() == 0) {
            break;
        }
        
        for (auto iLinkToBuilding = (*iBSDate).LinkToBuilding.begin(); iLinkToBuilding != (*iBSDate).LinkToBuilding.end(); iLinkToBuilding++) {
            (*iLinkToBuilding).changeBuildingWeight(-((*iLinkToBuilding).getWeight()));
        }
        
        std::cout << "Пересчет весов" << std::endl;
        std::clog << "Пересчет весов" << std::endl;

        auto jBSDate = iBSDate;
        jBSDate++;

        for (; jBSDate != BSDate.end(); jBSDate++) {
            (*jBSDate).setR((*iBSDate));
            (*jBSDate).recalcWeight();
        }

    }

    bsfile.close();
}

int main(int argc, char* argv[]) {

    setlocale(LC_ALL, "Russian");

    double x = NULL;

    int max_thread = omp_get_max_threads() - 1;
    
    //шаг ассета
    int Asset_step_m = 10;
    //уровень с которого строится БС
    int RxLevBad = -113;
    //радиус шайбы
    double R = 1000;
    //файлы с покрытием, сейчас принимает 1/4/7 этажи
    std::vector<std::pair<const char*, int>> CoverFile = {
        std::pair<const char*, int>("1.MIF", 1)
        ,std::pair<const char*, int>("4.MIF", 4)
        ,std::pair<const char*, int>("7.MIF", 7)
    };

    auto Building = "B.MIF";
    auto BSCurrent = "BS.csv";

    std::vector<PoligonBuilding> BuildingDate;
    std::vector<BSonBuilding> BSDate;
    std::vector<BS> BSList;

    if (argc > 1) {
        if (argc == 2 && strcmp(argv[1], "-h") == 0) {
            std::cout << "Примеры вызова:" << std::endl;
            std::cout << "white_vc.exe -f1 1.MIF -f4 4.MIF -f7 7.MIF" << std::endl;
            std::cout << "white_vc.exe -f1 1.MIF" << std::endl;
            std::cout << "white_vc.exe -f1 1.MIF -f4 4.MIF -f7 7.MIF -th 2 -as 10 -r 250 -rx -113" << std::endl;
            std::cout << std::endl;
            std::cout << "Список параметров для настройки:" << std::endl;
            std::cout << "Пар.  def     описание" << std::endl;
            std::cout << "-as   10      шаг в системе планирования, влияет на нарезку полигонов" << std::endl;
            std::cout << "-th   max-1   число ядер ЦП, которое будет использовано для работы" << std::endl;
            std::cout << "-rx   -113    уровень с которого начинает строиться здание" << std::endl;
            std::cout << "-r    1000    радиус 100% закрытия планируемой БС" << std::endl;
            std::cout << "-f1   -       название файла с покрытием для 1 этажа (регистр важен)" << std::endl;
            std::cout << "-f4   -       название файла с покрытием для 4 этажа (регистр важен)" << std::endl;
            std::cout << "-f7   -       название файла с покрытием для 7 этажа (регистр важен)" << std::endl;
            std::cout << "-b    B.MIF   название файла со зданиями (регистр важен)" << std::endl;
            std::cout << "-bs   BS.csv  название файла с базовыми станциями" << std::endl;
            std::cout << std::endl;
            std::cout << "Без списка файлов с этажами программа не запустится, предварительной проверки на наличии файла нет" << std::endl;
            std::cout << "Файл со зданиями должен называться B.MIF" << std::endl;
            return 0;
        }
        if (argc % 2 == 0) {
            std::cout << "Неверное число аргументов "<< argv[1] << std::endl;
            std::cout << "Для вызова справки вызовите из командной строки с параметром -h (например, white_vc.exe -h)" << std::endl;
            return 0;
        }
        for (size_t i = 1; i < argc; i+=2) {
            if (strcmp(argv[i], "-as")==0) {
                try {
                    int bb = std::atoi(argv[i + 1]);
                    Asset_step_m = bb;
                    std::cout << "Установлен шаг системы планирования: " << bb << std::endl;
                }
                catch (...) {
                    std::cout << "Невозможно выполнить конвертацию -as" << std::endl;
                    return 0;
                }
                continue;
            }
            if (strcmp(argv[i], "-th")==0) {
                try {
                    int bb = std::atoi(argv[i+1]);
                    if (bb > omp_get_max_threads()) {
                        throw;
                    }
                    max_thread = bb;
                    std::cout << "Установлено максимальное количество потоков: " << bb << std::endl;
                }
                catch (...) {
                    std::cout << "Невозможно выполнить конвертацию -th или столько ресурсов на ПК нет" << std::endl;
                    return 0;
                }
                continue;
            }
            if (strcmp(argv[i], "-rx")==0) {
                try {
                    int bb = std::atoi(argv[i + 1]);
                    RxLevBad = bb;
                    std::cout << "Установлен уровень плохого покрытия: " << bb << std::endl;
                }
                catch (...) {
                    std::cout << "Невозможно выполнить конвертацию -rx" << std::endl;
                    return 0;
                }
                continue;
            }
            if (strcmp(argv[i], "-r")==0) {
                try {
                    int bb = std::atoi(argv[i + 1]);
                    R = bb;
                    BS::setRMax(bb);
                    std::cout << "Установлен размер 100% закрытия БС: " << bb << std::endl;
                }
                catch (...) {
                    std::cout << "Невозможно выполнить конвертацию -r" << std::endl;
                    return 0;
                }
                continue;
            }
            if (strcmp(argv[i], "-f1") == 0) {
                CoverFile.push_back(std::pair<const char*, int>(argv[i+1], 1));
                std::cout << "Выбран файл для 1 этажа: " << argv[i + 1] << std::endl;
                continue;
            }
            if (strcmp(argv[i], "-f4") == 0) {
                CoverFile.push_back(std::pair<const char*, int>(argv[i + 1], 4));
                std::cout << "Выбран файл для 4 этажа: " << argv[i + 1] << std::endl;
                continue;
            }
            if (strcmp(argv[i], "-f7") == 0) {
                CoverFile.push_back(std::pair<const char*, int>(argv[i + 1], 7));
                std::cout << "Выбран файл для 7 этажа: " << argv[i + 1] << std::endl;
                continue;
            }
            if (strcmp(argv[i], "-b") == 0) {
                Building = argv[i + 1];
                std::cout << "Выбран файл для зданий: " << Building << std::endl;
                continue;
            }
            if (strcmp(argv[i], "-bs") == 0) {
                BSCurrent = argv[i + 1];
                std::cout << "Выбран файл для зданий: " << BSCurrent << std::endl;
                continue;
            }
           
            std::cout << "Параметр не известен" << std::endl;
            return 0;
        }
    }

    if (CoverFile.size() == 0) {
        std::cout << "Не добавлен ни 1 файл с покрытием для анализа" << std::endl;
        return 0;
    }

    for (size_t i = 0; i < CoverFile.size(); i++) {
        if (!checkFile(CoverFile[i].first)) {
            std::cout << "Не найден файл c покрытием " << CoverFile[i].first << " ВЫХОД" << std::endl;
            return 0;
        }
    }

    if (!checkFile(Building)) {
        std::cout << "Не найден файл со зданиями.... ВЫХОД" << std::endl;
        return 0;
    }

    if (!checkFile(BSCurrent)) {
        std::cout << "Не найден файл со зданиями.... ВЫХОД" << std::endl;
        return 0;
    }
    
    omp_set_num_threads(max_thread);

    std::ofstream ofs("logfile.txt");
    std::clog.rdbuf(ofs.rdbuf());

    time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Старт: " << std::put_time(localtime(&now), "%F %T") << std::endl;
    std::clog << "Старт: " << std::put_time(localtime(&now), "%F %T") << std::endl;

    double start = clock();
    
    GDALAllRegister();

    try {
        Read_Building(Building, BuildingDate);
    }
    catch (const char* a) {
        std::cout << a << std::endl;
        return 0;
    }
      
    std::cout << "Всего зданий: ";
    std::cout << BuildingDate.size() << std::endl;

    std::clog << "Всего зданий: ";
    std::clog << BuildingDate.size() << std::endl;

    
    for (auto iCoverFile = 0; iCoverFile < CoverFile.size(); iCoverFile++) {
        calcBadAreaByCover(BuildingDate, CoverFile[iCoverFile].first, CoverFile[iCoverFile].second, Asset_step_m);
    }
    
    OutFile(BuildingDate);
    
    OutKml(BuildingDate, "building");

    //создаем файл с текущими БС
    
    std::cout << "Читаем файл с базовыми станциями" << std::endl;
    std::clog << "Читаем файл с базовыми станциями" << std::endl;

    try {
        CSVFile BSfile;
        BSfile.open(BSCurrent, "	");

        size_t LatNum = BSfile.getColumnNumByName("LATITUDE");
        size_t LonNum = BSfile.getColumnNumByName("LONGITUDE");

        for (size_t i = 0; i < BSfile.getSizeRow(); i++) {
            BSList.push_back(BS(std::stod(BSfile.getColumnValueByNum(i, LatNum)), std::stod(BSfile.getColumnValueByNum(i, LonNum))));
        }

    }
    catch (const char* a) {
        std::cout << a << std::endl;
        return 0;
    }
    
    std::cout << "На вход подано Базовых станций: " << BSList.size() << std::endl;
    std::clog << "На вход подано Базовых станций: " << BSList.size() << std::endl;



    //------------------------------------
    //создаем массив 
    std::cout << "Заполняем вес исходных зданий и формируем список БС кандидатов" << BSDate.size() << std::endl;
    std::clog << "Заполняем вес исходных зданий и формируем список БС кандидатов" << BSDate.size() << std::endl;
    
    for (size_t i = 0; i < BuildingDate.size(); i++) {
        auto LevelList = PoligonBuilding::getLevel();
        for (auto iLevelList = LevelList.begin(); iLevelList != LevelList.end(); iLevelList++) {
            if ((*iLevelList) <= RxLevBad) {
                if (BuildingDate[i].LevelArea.find((*iLevelList)) != BuildingDate[i].LevelArea.end()) {
                    BSonBuilding buff(BuildingDate[i]);
                    BSDate.push_back(buff);
                    break;
                }
            }
        }
        BuildingDate[i].createWeight();
    }

#pragma omp parallel for
    for (int i = 0; i < BSList.size(); i++) {
        BSList[i].Coordinate.transformTo(BSDate[0].Coordinate.getSpatialReference());
    }

#pragma omp parallel for
    for (int i = 0; i < BSDate.size(); i++) {
        BSDate[i].setR(BSList);
    }
    
    now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::cout << "Зданий кандидатов: " << BSDate.size() << std::endl;
    std::clog << "Зданий кандидатов: " << BSDate.size() << std::endl;
        
    topByCircle(BSDate, BuildingDate);
   
    now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "Финиш: " << std::put_time(localtime(&now), "%F %T") << std::endl;
    std::clog << "Финиш: " << std::put_time(localtime(&now), "%F %T") << std::endl;

    double tt = ((clock() - start) / CLOCKS_PER_SEC);
    int hh = trunc(tt / 3600);
    int mi = trunc((tt - hh * 3600)/60);
    int ss = tt - 3600 * hh - mi * 60;

    std::cout << "Затраченное время: " << hh << " часов " << mi << " минут " << ss << " секунд" << std::endl;
    std::cout << "Затраченное время (в секундах): " << tt << std::endl;

    std::clog << "Затраченное время: " << hh << " часов " << mi << " минут " << ss << " секунд" << std::endl;
    std::clog << "Затраченное время (в секундах): " << tt << std::endl;

    return 0;
}