#include "phantom1.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

Phantom1::Phantom1() {
}
Phantom1::~Phantom1() {
}

double get_density(double x,double y,double z) {
	double rad=sqrt(x*x+y*y);
	if ((rad>10)||(rad<6)) return 0;
	if ((z<-3)||(z>3)) return 0;
	return 1;
}

SpinProperties Phantom1::spinPropertiesAt(double x,double y,double z) const {
	SpinProperties ret;
	ret.T2star=30;
	ret.chemical_shift=0;
	ret.density=get_density(x,y,z);
	return ret;
}

void Phantom1::getBoundingBox(double &xmin,double &xmax,double &ymin,double &ymax,double &zmin,double &zmax) const {
	xmin=-10; xmax=10;
	ymin=-10; ymax=10;
	zmin=-3; zmax=3;
}

double Phantom1::T1() const {
	return 100;
}
double Phantom1::T2() const {
	return 40;
}
