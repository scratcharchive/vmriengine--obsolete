
#include "vmrisimulator.h"
#include <stdio.h>
#include "qjson.h"
#include "textfile.h"
#include <QMap>
#include <QStringList>
#include "phantom1.h"
#include <QDebug>

void usage() {
	printf("usage: vmriengine --phantom=[phantom1|other] --simblocks=[simblocks.json] --output=[output.raw] --num_iso=[10000]\n");
}

int main(int argc,char **argv) {
	
	QMap<QString,QString> parameters;
	
	for (int i=1; i<argc; i++) {
		QString str=argv[i];
		QStringList vals=str.split("=");
		if (vals.count()==2) {
			if (vals[0].mid(0,2)=="--") {
				parameters[vals[0].mid(2)]=vals[1];
			}
		}
	}
	
	QString phantom_name=parameters["phantom"];
	QString simblocks_file=parameters["simblocks"];
	QString output_file=parameters["output"];
	int num_isochromats=parameters.value("num_iso","10000").toInt();
	
	
	if ((phantom_name.isEmpty())||(simblocks_file.isEmpty())||(output_file.isEmpty())) {
		usage();
		return 1;
	}
	
	if (phantom_name!="phantom1") {
		qWarning() << "Unrecognized phantom: "+phantom_name;
		usage();
		return 1;
	}
	
	Phantom1 PP;
	QList<double> x,y,z,f,d;
	PP.generateRandomLocations(x,y,z,num_isochromats);
	for (int i=0; i<num_isochromats; i++) {
		SpinProperties SP=PP.spinPropertiesAt(x[i],y[i],z[i]);
		f << SP.chemical_shift*42.57*1.5;
		d << SP.density;
	}
	
	QMap<QString,QVariant> simblocks0=parseJSON(read_text_file(simblocks_file)).toMap();
	
	VmriSimulator VS;
	VS.setT1(PP.T1());
	VS.setT2(PP.T2());
	VS.setIsochromats(x,y,z,f,d);
	VS.setSimBlocks(simblocks0);
	VS.simulate();
	
	QList<VmriReadout> readouts=VS.getReadouts();
	
	FILE *outf=fopen(output_file.toAscii().data(),"wb");
	for (int i=0; i<readouts.count(); i++) {
		VmriReadout RO=readouts[i];
		
		quint32 N0=RO.real.count();
		quint32 readout_index0=RO.readout_index;
		float dwell_time0=RO.dwell_time;
		fwrite(&N0,sizeof(quint32),1,outf);
		fwrite(&readout_index0,sizeof(quint32),1,outf);
		fwrite(&dwell_time0,sizeof(float),1,outf);
		for (int i=0; i<29; i++) {
			quint32 tmp=0;
			fwrite(&tmp,sizeof(quint32),1,outf);
		}
		for (int i=0; i<N0; i++) {
			float tmp_real=RO.real.value(i);
			float tmp_imag=RO.imag.value(i);
			fwrite(&tmp_real,sizeof(float),1,outf);
			fwrite(&tmp_imag,sizeof(float),1,outf);
		}
	}
	fclose(outf);
	
	printf("\n");
	
	return 0;
	
}
