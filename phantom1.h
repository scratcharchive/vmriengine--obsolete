#ifndef phantom1_H
#define phantom1_H

#include "abstractphantom.h"

class Phantom1 : public AbstractPhantom {
public:
	Phantom1();
	virtual ~Phantom1();
	
	virtual double T1() const;
	virtual double T2() const;
	virtual SpinProperties spinPropertiesAt(double x,double y,double z) const;
	virtual void getBoundingBox(double &xmin,double &xmax,double &ymin,double &ymax,double &zmin,double &zmax) const;
};

#endif
