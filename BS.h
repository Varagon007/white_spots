#pragma once
#include "ogrsf_frmts.h"
#include <string>

// Basic class of Base station of mobile operator 
class BS
{
public:
	//coordinate base station
	OGRPoint	Coordinate;
	//name NS
	std::string BsName;

	double		getR();
	
	void		setR(double R);

	void		setR(BS BorderBS);

	void		setR(std::vector<BS> BorderBSList);

	int			getBsId();
	
	static double	getRMax();

	static void		setRMax(double R);

	BS(double R = 0);

	BS(double X, double Y, double R = 0);

	BS(double X, double Y, OGRSpatialReference srNow, double R = 0);


private:
	//unique id BS
	int				BsId;
	//radius of coverage BS
	double			R;
	//maximum BS radius
	static double	RMax;
	//increment
	static int		IdNow;

	static int		getNext();
};

