#include "vmrisimulator.h"
#include "vmriinterface.h"
#include <stdio.h>
#include <math.h>
#include <QDebug>
#include <ctime>

#define PI 3.141592653589793


class VmriSimulatorPrivate {
public:
	VmriInterface m_interface;
	QMap<QString,QVariant> m_isochromats;
	QMap<QString,QVariant> m_simblocks;
	double m_T1,m_T2;
	QList<VmriReadout> m_readouts;
	
	void do_evolve(const QList<double> &moment,double duration);
};

VmriSimulator::VmriSimulator() {
	d=new VmriSimulatorPrivate;
	d->m_T1=100;
	d->m_T2=40;
}

VmriSimulator::~VmriSimulator() {
	delete d;
}

void VmriSimulator::setT1(double T1) {
	d->m_T1=T1;
}
void VmriSimulator::setT2(double T2) {
	d->m_T2=T2;
}

void VmriSimulator::setIsochromats(const QMap<QString,QVariant> &isochromats) {
	d->m_isochromats=isochromats;
}

void VmriSimulator::setSimBlocks(const QMap<QString,QVariant> &simblocks) {
	d->m_simblocks=simblocks;
}

QList<double> to_double_list(const QList<QVariant> &L) {
	QList<double> ret;
	for (long i=0; i<L.count(); i++) {
		ret << L[i].toDouble();
	}
	return ret;
}

void VmriSimulator::simulate() {
	qDebug()  << "Simulating";
	
	d->m_interface.setT1(d->m_T1);
	d->m_interface.setT2(d->m_T2);
	
	QList<double> xx=to_double_list(d->m_isochromats.value("x").toList());
	QList<double> yy=to_double_list(d->m_isochromats.value("y").toList());
	QList<double> zz=to_double_list(d->m_isochromats.value("z").toList());
	QList<double> ff=to_double_list(d->m_isochromats.value("f").toList());
	QList<double> dd=to_double_list(d->m_isochromats.value("d").toList());
	d->m_interface.setIsochromats(xx,yy,zz,ff,dd);
	
	QList<QVariant> rf_waveforms=d->m_simblocks["rf_waveforms"].toList();
	QList<QVariant> blocks=d->m_simblocks["blocks"].toList();
	
	qDebug()  << "Adding RF waveforms...";
	for (int i=0; i<rf_waveforms.count(); i++) {
		QMap<QString,QVariant> WW=rf_waveforms[i].toMap();
		QList<double> data_real=to_double_list(WW["data_real"].toList());
		QList<double> data_imag=to_double_list(WW["data_imag"].toList());
		double dt=WW["dt"].toDouble();
		d->m_interface.addRFWaveform(data_real,data_imag,dt);
	}
	
	d->m_interface.initialize();
	
	qDebug()  << "Executing blocks...";
	
	std::clock_t start;
	start=std::clock();
	
	for (int i=0; i<blocks.count(); i++) {
		QMap<QString,QVariant> BB=blocks[i].toMap();
		QString block_type=BB["block_type"].toString();
		double duration=BB["duration"].toDouble();
		if (block_type=="evolve") {
		    printf("Evolve (%g ms)...",duration);
		    start=std::clock();
			QList<double> gradient_moment=to_double_list(BB["gradient_moment"].toList());
			d->do_evolve(gradient_moment,duration);
			printf("(elapsed=%g). ", (std::clock() - start ) / (double) CLOCKS_PER_SEC*1000);
		}
		else if (block_type=="rf_pulse") {
		    printf("RF (%g ms)...",duration);
		    start=std::clock();
			QList<double> gradient_amplitude=to_double_list(BB["gradient_amplitude"].toList());
			double A[3];
			double factor=42.57; //uT/mm * Hz/uT = Hz/mm
			A[0]=gradient_amplitude.value(0)*factor;
			A[1]=gradient_amplitude.value(1)*factor;
			A[2]=gradient_amplitude.value(2)*factor;
			int rf_waveform_index=BB["rf_waveform_index"].toInt();
			double dt=BB["dt"].toDouble();
			double phase=BB["phase"].toDouble();
			double frequency=BB["frequency"].toDouble();
			d->m_interface.excite(rf_waveform_index,A,dt,phase,frequency);
			printf("(elapsed=%g). ", (std::clock() - start ) / (double) CLOCKS_PER_SEC*1000);
		}
		else if (block_type=="readout") {
		    printf("Readout (%g ms).\n",duration);
		    start=std::clock();
			double T1=d->m_T1;
			double T2=d->m_T2;
			
			QList<double> gradient_amplitude=to_double_list(BB["gradient_amplitude"].toList());
			double A[3];
			double factor=42.57; //uT/mm * Hz/uT = Hz/mm
			A[0]=gradient_amplitude.value(0)*factor;
			A[1]=gradient_amplitude.value(1)*factor;
			A[2]=gradient_amplitude.value(2)*factor;
			double dt=BB["dt"].toDouble();
			double phase=BB["phase"].toDouble();
			double frequency=BB["frequency"].toDouble();
			int N=BB["N"].toInt();
			int readout_index=BB["readout_index"].toInt();
	
			QList<double> output_real;
			QList<double> output_imag;
	
			double *data_real=(double *)malloc(N*sizeof(double));
			double *data_imag=(double *)malloc(N*sizeof(double));
			
			d->m_interface.readout(data_real,data_imag,N,A,dt,frequency);
			for (int i=0; i<N; i++) {
				output_real << data_real[i];
				output_imag << data_imag[i];
			}
			free(data_real);
			free(data_imag);
	
			//finally, apply global phase and decay
			//evolve by half of a step
			double cosphase=cos(phase*PI/180);
			double sinphase=sin(phase*PI/180);
			for (int t=0; t<N; t++) {
				double re1=output_real[t]*cosphase-output_imag[t]*sinphase;
				double im1=output_real[t]*sinphase+output_imag[t]*cosphase;
				re1*=exp(-(t+0.5)*dt/(1000)/T2);
				im1*=exp(-(t+0.5)*dt/(1000)/T2);
				output_real[t]=re1;
				output_imag[t]=im1;
			}
	
			double dur0=dt*N;
			QList<double> moment0; moment0 << 0 << 0 << 0;
			moment0[0]=gradient_amplitude[0]*dur0;
			moment0[1]=gradient_amplitude[1]*dur0;
			moment0[2]=gradient_amplitude[2]*dur0;
			d->do_evolve(moment0,dur0);

			VmriReadout readout;
			readout.real=output_real;
			readout.imag=output_imag;
			readout.readout_index=readout_index; //this is set outside this function
			readout.dwell_time=dt;
			
			d->m_readouts << readout;
			printf("(elapsed=%g). ", (std::clock() - start ) / (double) CLOCKS_PER_SEC*1000);
		}
	}
}

void VmriSimulatorPrivate::do_evolve(const QList<double> &moment,double duration) {
	double A[3];
	double factor=42.57/(1000*1000)*2*PI; //uT/mm.us * Hz/uT / 10^6 * 2pi = rad/mm
	A[0]=moment.value(0)*factor;
	A[1]=moment.value(1)*factor;
	A[2]=moment.value(2)*factor;
	double E1=exp(-duration/m_T1);
	double E2=exp(-duration/m_T2);
	m_interface.evolve(A,duration,E1,E2);
}

QList<VmriReadout> VmriSimulator::getReadouts() {
	return d->m_readouts;
}

