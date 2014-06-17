#include "abstractphantom.h"
#include <stdlib.h>
#include <stdio.h>

void make_rand(double &x,double xmin,double xmax) {
	double ret=(rand()%100000)*1.0/100000;
	x=ret*(xmax-xmin)+xmin;
}

void make_rand(double &x,double &y,double &z,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax) {
	make_rand(x,xmin,xmax);
	make_rand(y,ymin,ymax);
	make_rand(z,zmin,zmax);
}



double AbstractPhantom::computeVolume() const {
	double xmin,xmax,ymin,ymax,zmin,zmax;
	getBoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);
	
	long num_hits=0;
	long num_trials=0;
	long num_trials_per_region=100;
	for (int iz=0; iz<10; iz++)
	for (int iy=0; iy<10; iy++)
	for (int ix=0; ix<10; ix++) {
		double xmin0=xmin+ix*(xmax-xmin)/10;
		double xmax0=xmin+(ix+1)*(xmax-xmin)/10;
		double ymin0=ymin+iy*(ymax-ymin)/10;
		double ymax0=ymin+(iy+1)*(ymax-ymin)/10;
		double zmin0=zmin+iz*(zmax-zmin)/10;
		double zmax0=zmin+(iz+1)*(zmax-zmin)/10;
		for (long i=0; i<num_trials_per_region; i++) {
			double x0,y0,z0;
			make_rand(x0,y0,z0,xmin0,xmax0,ymin0,ymax0,zmin0,zmax0);
			SpinProperties PP=spinPropertiesAt(x0,y0,z0);
			if (PP.density) {
				num_hits++;
			}
			num_trials++;
		}
	}
	return num_hits*1.0/num_trials*(xmax-xmin)*(ymax-ymin)*(zmax-zmin);
}

void AbstractPhantom::generateRandomLocations(QList<double> &x,QList<double> &y,QList<double> &z,int N) const {
	//TO DO: Figure out how to make this more uniform
	
	double xmin,xmax,ymin,ymax,zmin,zmax;
	getBoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);
	
	long ct=0;
	while (x.count()<N) {
		double x0,y0,z0;
		make_rand(x0,y0,z0,xmin,xmax,ymin,ymax,zmin,zmax);
		SpinProperties PP=spinPropertiesAt(x0,y0,z0);
		if (PP.density) {
			x << x0;
			y << y0;
			z << z0;
			ct=0;
		}
		else ct++;
		if (ct>10000) {
			printf("Unexpected problem in generateRandomLocations\n");
			return;
		}
	}
	
}
