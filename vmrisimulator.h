#ifndef vmrisimulator_h
#define vmrisimulator_cpp

#include <QString>
#include <QMap>
#include <QVariant>

struct VmriReadout {
	QList<double> real;
	QList<double> imag;
	int readout_index;
	float dwell_time;
};

class VmriSimulatorPrivate;
class VmriSimulator {
public:
	friend class VmriSimulatorPrivate;
	
	VmriSimulator();
	virtual ~VmriSimulator();
	
	void setT1(double T1);
	void setT2(double T2);
	void setIsochromats(const QMap<QString,QVariant> &isochromats);
	void setSimBlocks(const QMap<QString,QVariant> &simblocks);
	void simulate();
	QList<VmriReadout> getReadouts();
	
private:
	VmriSimulatorPrivate *d;
};

#endif