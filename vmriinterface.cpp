#include "vmriinterface.h"
#include <QList>
#include <QDebug>

#define PI 3.141592653589793
#include <math.h>

struct Isochromat {
	double x,y,z,f;
	double Mx,My,Mz;
	double density;
};
struct Matrix44 {
	double data[4][4];
};
struct RF_operator {
	QList<Matrix44> matrices;
};

class VmriInterfacePrivate {
public:
	double m_T1,m_T2;
	double m_gamma;
	QList<Isochromat> m_isochromats;
	QList<RF_operator> m_RF_operators;
	
	void rotateZ(double *M,double theta);
	void rotateY(double *M,double theta);
	void decayMag(double *M,double E1,double E2);
	double dot_product(double *V1,double *V2);
	void apply_matrix(double *Mag,const Matrix44 &MM);
	Matrix44 compute_RF_matrix(const QList<double> &data_real,const QList<double> &data_imag,double dt,double dtheta);
	Matrix44 get_RF_operator_matrix(int waveform_index,double dtheta);
};

VmriInterface::VmriInterface() {
	d=new VmriInterfacePrivate;
	d->m_T1=100;
	d->m_T2=40;
	d->m_gamma=42.57;
}

VmriInterface::~VmriInterface() {
	delete d;
}

void VmriInterface::setT1(double T1) {
	d->m_T1=T1;
}
void VmriInterface::setT2(double T2) {
	d->m_T2=T2;
}

void VmriInterface::setIsochromats(const QList<double> &x,const QList<double> &y,const QList<double> &z,const QList<double> &f,const QList<double> &density) {
	d->m_isochromats.clear();
	int N=x.count();
	for (int i=0; i<N; i++) {
		Isochromat II;
		II.x=x[i];
		II.y=y[i];
		II.z=z[i];
		II.f=f[i];
		II.density=density[i];
		II.Mx=0;
		II.My=0;
		II.Mz=1;
		d->m_isochromats << II;
	}
}
void VmriInterface::addRFWaveform(const QList<double> &data_real,const QList<double> &data_imag,double dt) {
	RF_operator XX;
	
	int N=data_real.count();
	
	for (int i=0; i<N; i++) {
		double dtheta=-PI+i*1.0/N*(2*PI);
		Matrix44 MM=d->compute_RF_matrix(data_real,data_imag,dt,dtheta);
		XX.matrices << MM;
	}
	
	d->m_RF_operators << XX;
}
void VmriInterface::initialize() {
}

void VmriInterface::evolve(double A[3],double t,double E1,double E2) {;
	for (int i=0; i<d->m_isochromats.count(); i++) {
		Isochromat *II=&d->m_isochromats[i];
		double pos0[3];
		pos0[0]=II->x;
		pos0[1]=II->y;
		pos0[2]=II->z;
		double theta=d->dot_product(A,pos0); //radians
		theta+=(II->f)*t*2*PI/(1000*1000);
		double Mag[3];
		Mag[0]=II->Mx;
		Mag[1]=II->My;
		Mag[2]=II->Mz;
		d->rotateZ(Mag,theta);
		d->decayMag(Mag,E1,E2);
		II->Mx=Mag[0];
		II->My=Mag[1];
		II->Mz=Mag[2];
	}
}

Matrix44 VmriInterfacePrivate::get_RF_operator_matrix(int waveform_index,double dtheta) {
	RF_operator *OO=&m_RF_operators[waveform_index];
	double frac=(dtheta-(-PI))/(2*PI);
	int ind=round(frac*OO->matrices.count());
	if ((ind<0)||(ind>=OO->matrices.count())) {
		Matrix44 MM;
		for (int i=0; i<4; i++)
		for (int j=0; j<4; j++) {
			if (i==j) MM.data[i][j]=1;
			else MM.data[i][j]=0;
		}
		return MM;
	}
	return OO->matrices[ind];
}

void VmriInterface::excite(int rf_waveform_index,double A[3],double dt,double phase,double frequency) {
	for (int i=0; i<d->m_isochromats.count(); i++) {
		Isochromat *II=&d->m_isochromats[i];
		double pos0[3];
		pos0[0]=II->x;
		pos0[1]=II->y;
		pos0[2]=II->z;
		double freq=II->f + d->dot_product(A,pos0) -frequency; //Hz
		double dtheta=freq*dt/(1000*1000);
		Matrix44 MM=d->get_RF_operator_matrix(rf_waveform_index,dtheta);
		
		double Mag[3];
		Mag[0]=II->Mx;
		Mag[1]=II->My;
		Mag[2]=II->Mz;
		d->rotateZ(Mag,phase*PI/180);
		d->apply_matrix(Mag,MM);
		d->rotateZ(Mag,-phase*PI/180);
		II->Mx=Mag[0];
		II->My=Mag[1];
		II->Mz=Mag[2];
	}
}
void VmriInterface::readout(double *output_real,double *output_imag,int N,double A[3],double dt,double frequency) {
	//we are going to accumulate the spins in the Fourier domain
	//this has a couple advantages... first it may be faster, second we can filter out the high freqs
	
	QList<double> frequencies; //Hz
	QList<double> signals_real;
	QList<double> signals_imag;
	
	for (int ii=0; ii<d->m_isochromats.count(); ii++) {
		Isochromat *II=&d->m_isochromats[ii];
		double pos0[3];
		pos0[0]=II->x;
		pos0[1]=II->y;
		pos0[2]=II->z;
		frequencies << II->f+d->dot_product(A,pos0)-frequency; //this will be in Hz
		signals_real << II->Mx*II->density;
		signals_imag << II->My*II->density;
	}
	
	int oversamp=4; //oversampling factor
	QList<double> X_real;
	QList<double> X_imag;
	for (int i=0; i<N*oversamp; i++) {
		X_real << 0;
		X_imag << 0;
	}
	//The sampling time is dt
	//so the sampling bandwidth is 1/dt
	double sampling_bandwidth=1.0/dt*(1000*1000); //Hz
	double freq_step=sampling_bandwidth/(N*oversamp);
	double min_freq=-freq_step*(N*oversamp/2);
	double max_freq=freq_step*(N*oversamp/2-1);
	for (int i=0; i<frequencies.count(); i++) {
		double f0=frequencies[i];
		if ((min_freq<=f0)&&(f0<=max_freq)) {
			int ind0=round((f0-min_freq)/freq_step);
			if ((0<=ind0)&&(ind0<N*oversamp)) {
				//in the future, consider spreading this using sinc interpolation
				X_real[ind0]+=signals_real[i];
				X_imag[ind0]+=signals_imag[i];
			}
		}
	}
	
	for (int t=0; t<N; t++) {
		output_real[t]=0;
		output_imag[t]=0;
	}
	
	for (int i=0; i<N*oversamp; i++) {
		double re0=X_real[i];
		double im0=X_imag[i];
		if ((re0)||(im0)) {
			double freq=min_freq+freq_step*i;
			double dtheta=freq*dt/(1000*1000)*2*PI;
			double cosdtheta=cos(dtheta);
			double sindtheta=sin(dtheta);
			double cosdtheta2=cos(dtheta/2);
			double sindtheta2=sin(dtheta/2);
		
			//evolve by half of a step
			double re1=cosdtheta2*re0-sindtheta2*im0;
			double im1=sindtheta2*re0+cosdtheta2*im0;
			re0=re1;
			im0=im1;

			for (int t=0; t<N; t++) {
				output_real[t]+=re0;
				output_imag[t]+=im0;
				
				//evolve by full step
				re1=cosdtheta*re0-sindtheta*im0;
				im1=sindtheta*re0+cosdtheta*im0;
				re0=re1;
				im0=im1;
				
			}
		}
	}
	
	//finally, apply global phase and decay
	//evolve by half of a step
	//this will need to be done elsewhere!!!
	/*
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
	*/
	
	//this should be done elsewhere
	/*
	double dur0=dt*N;
	double moment0[3];
	moment0[0]=gradient_amplitude[0]*dur0;
	moment0[1]=gradient_amplitude[1]*dur0;
	moment0[2]=gradient_amplitude[2]*dur0;
	evolve(dur0,moment0);
	*/
}

void VmriInterfacePrivate::rotateZ(double *M,double theta) {
	double costheta=cos(theta);
	double sintheta=sin(theta);
	double ret0=costheta*M[0]-sintheta*M[1];
	double ret1=sintheta*M[0]+costheta*M[1];
	M[0]=ret0;
	M[1]=ret1;
}
void VmriInterfacePrivate::rotateY(double *M,double theta) {
	double costheta=cos(theta);
	double sintheta=sin(theta);
	double ret0=costheta*M[0]+sintheta*M[2];
	double ret2=-sintheta*M[0]+costheta*M[2];
	M[0]=ret0;
	M[2]=ret2;
}
void VmriInterfacePrivate::decayMag(double *M,double E1,double E2) {
	M[0]*=E2;
	M[1]*=E2;
	M[2]=1-(1-M[2])*E1;
}

double VmriInterfacePrivate::dot_product(double *V1,double *V2) {
	return V1[0]*V2[0]+V1[1]*V2[1]+V1[2]*V2[2];
}
void VmriInterfacePrivate::apply_matrix(double *Mag,const Matrix44 &MM) {
	double ret0=Mag[0]*MM.data[0][0]+Mag[1]*MM.data[0][1]+Mag[2]*MM.data[0][2]+MM.data[0][3];
	double ret1=Mag[1]*MM.data[1][0]+Mag[1]*MM.data[1][1]+Mag[2]*MM.data[1][2]+MM.data[1][3];
	double ret2=Mag[2]*MM.data[2][0]+Mag[1]*MM.data[2][1]+Mag[2]*MM.data[2][2]+MM.data[2][3];
	Mag[0]=ret0;
	Mag[1]=ret1;
	Mag[2]=ret2;
}
Matrix44 VmriInterfacePrivate::compute_RF_matrix(const QList<double> &data_real,const QList<double> &data_imag,double dt,double dtheta) {
	double T1=m_T1;
	double T2=m_T2;
	
	double E1=exp(-(dt/1000)/T1);
	double E2=exp(-(dt/1000)/T2);
	double E1b=exp(-(dt/2/1000)/T1);
	double E2b=exp(-(dt/2/1000)/T2);
	
	double M0[3]; M0[0]=0; M0[1]=0; M0[2]=0;
	double M1[3]; M1[0]=1; M1[1]=0; M1[2]=0;
	double M2[3]; M2[0]=0; M2[1]=1; M2[2]=0;
	double M3[3]; M3[0]=0; M3[1]=0; M3[2]=1;
	
	rotateZ(M0,dtheta/2);
	rotateZ(M1,dtheta/2);
	rotateZ(M2,dtheta/2);
	rotateZ(M3,dtheta/2);
	
	decayMag(M0,E1b,E2b);
	decayMag(M1,E1b,E2b);
	decayMag(M2,E1b,E2b);
	decayMag(M3,E1b,E2b);
	
	int N=data_real.count();
	
	for (int i=0; i<N; i++) {
		if (i>0) {
			rotateZ(M0,dtheta);
			rotateZ(M1,dtheta);
			rotateZ(M2,dtheta);
			rotateZ(M3,dtheta);
			
			decayMag(M0,E1,E2);
			decayMag(M1,E1,E2);
			decayMag(M2,E1,E2);
			decayMag(M3,E1,E2);
		}
		double re0=data_real[i]; //uT
		double im0=data_imag[i]; //uT
		double mag0=sqrt(re0*re0+im0*im0); //uT
		double phi=atan2(im0,re0); //radians
		double theta=mag0*m_gamma*dt/(1000*1000)*2*PI; //radians
		
		rotateZ(M0,phi);
		rotateZ(M1,phi);
		rotateZ(M2,phi);
		rotateZ(M3,phi);
		
		rotateY(M0,theta);
		rotateY(M1,theta);
		rotateY(M2,theta);
		rotateY(M3,theta);
		
		rotateZ(M0,-phi);
		rotateZ(M1,-phi);
		rotateZ(M2,-phi);
		rotateZ(M3,-phi);
	}
	
	rotateZ(M0,dtheta/2);
	rotateZ(M1,dtheta/2);
	rotateZ(M2,dtheta/2);
	rotateZ(M3,dtheta/2);
	
	decayMag(M0,E1b,E2b);
	decayMag(M1,E1b,E2b);
	decayMag(M2,E1b,E2b);
	decayMag(M3,E1b,E2b);
	
	Matrix44 matrix;
	matrix.data[0][0]=M1[0]; matrix.data[0][1]=M2[0]; matrix.data[0][2]=M3[0]; matrix.data[0][3]=M0[0];
	matrix.data[1][0]=M1[1]; matrix.data[1][1]=M2[1]; matrix.data[1][2]=M3[1]; matrix.data[1][3]=M0[1];
	matrix.data[2][0]=M1[2]; matrix.data[2][1]=M2[2]; matrix.data[2][2]=M3[2]; matrix.data[2][3]=M0[2];
	matrix.data[3][0]=0; matrix.data[3][1]=0; matrix.data[3][2]=0; matrix.data[3][3]=1;
	
	return matrix;
}

