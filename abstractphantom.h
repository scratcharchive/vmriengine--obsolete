#ifndef abstractphantom_H
#define abstractphantom_H

#include <QList>

struct SpinProperties {
	double T2star; //ms
	double chemical_shift; //ppm
	double density;
};

class AbstractPhantom {
public:
	AbstractPhantom() {}
	virtual ~AbstractPhantom() {}
	
	double computeVolume() const;
	void generateRandomLocations(QList<double> &x,QList<double> &y,QList<double> &z,int N) const;
	
	virtual double T1() const=0;
	virtual double T2() const=0;
	virtual SpinProperties spinPropertiesAt(double x,double y,double z) const=0;
	virtual void getBoundingBox(double &xmin,double &xmax,double &ymin,double &ymax,double &zmin,double &zmax) const=0;
};

#endif
