#include "BS.h"

int	BS::IdNow = 0;

double	BS::RMax = 1000;

OGRSpatialReference sr;

double BS::getR()
{
	return R;
}

double BS::getRMax()
{
	return RMax;
}

void BS::setRMax(double R)
{
	if (R < 0) {
		RMax = 0;
	}
	else {
		RMax = R;
	}
}

void BS::setR(double R)
{
	if (R < 0) {
		this->R = 0;
	}
	else {
		if (R < RMax) {
			this->R = R;
		}
		else {
			this->R = RMax;
		}
	}
}

void BS::setR(BS BorderBSList)
{
	double RRes = this->R, Rbuff;

	BorderBSList.Coordinate.transformTo(Coordinate.getSpatialReference());

	Rbuff = Coordinate.Distance(&BorderBSList.Coordinate);

	if (Rbuff < RRes) {
		RRes = Rbuff;
	}
	
	setR(RRes);
}

void BS::setR(std::vector<BS> BorderBSList)
{
	double RRes = RMax, Rbuff;

	for (int i = 0; i < BorderBSList.size(); i++) {

		Rbuff = Coordinate.Distance(&BorderBSList[i].Coordinate);

		if (Rbuff < RRes) {
			RRes = Rbuff;
		}
	}

	setR(RRes);
			
}

void BS::checkR(BS BS)
{
}

int BS::getBsId()
{
	return BsId;
}

BS::BS(double R)
{
	setR(R);
	BsId = getNext();
}

BS::BS(double X, double Y, double R)
{

	sr.SetGeogCS("My geographic CRS",
		"World Geodetic System 1984",
		"My WGS84 Spheroid",
		SRS_WGS84_SEMIMAJOR, SRS_WGS84_INVFLATTENING,
		"Greenwich", 0.0,
		"degree", 0.0174532925199433);

	setR(R);

	Coordinate.setX(X);
	Coordinate.setY(Y);

	Coordinate.assignSpatialReference(&sr);

	BsId = getNext();

}

BS::BS(double X, double Y, OGRSpatialReference srNow, double R)
{
	setR(R);

	Coordinate.setX(X);
	Coordinate.setY(Y);

	Coordinate.assignSpatialReference(&srNow);

}

int BS::getNext()
{
	return IdNow++;
}
