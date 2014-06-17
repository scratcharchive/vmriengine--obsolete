#ifndef vmriinterface_H
#define vmriinterface_H

#include <QList>

class VmriInterfacePrivate;
class VmriInterface {
public:
	friend class VmriInterfacePrivate;
	
	VmriInterface();
	virtual ~VmriInterface();
	
	void setT1(double T1);
	void setT2(double T2);
	void setIsochromats(const QList<double> &x,const QList<double> &y,const QList<double> &z,const QList<double> &f,const QList<double> &d);
	void addRFWaveform(const QList<double> &data_real,const QList<double> &data_imag,double dt);
	void initialize();
	
	void evolve(double A[3],double t,double E1,double E2);
	void excite(int rf_waveform_index,double A[3],double dt,double phase,double frequency);
	void readout(double *ret_real,double *ret_imag,int num_readout_points,double A[3],double dt,double frequency);
	
private:
	VmriInterfacePrivate *d;
};

#endif
