
#include "vmrisimulator.h"
#include <stdio.h>
#include "qjson.h"
#include "textfile.h"

void usage() {
	printf("usage: vmriengine [isochromats.json] [simblocks.json] [output.raw]\n");
}

int main(int argc,char **argv) {
	
	if (argc-1<3) {
		usage();
		return -1;
	}
	
	QString isochromats_file=QString(argv[1]);
	QString simblocks_file=QString(argv[2]);
	QString output_file=QString(argv[3]);
	
	QMap<QString,QVariant> isochromats0=parseJSON(read_text_file(isochromats_file)).toMap();
	QMap<QString,QVariant> simblocks0=parseJSON(read_text_file(simblocks_file)).toMap();
	
	VmriSimulator VS;
	VS.setT1(isochromats0["T1"].toDouble());
	VS.setT2(isochromats0["T2"].toDouble());
	VS.setIsochromats(isochromats0);
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
	
	return 0;
	
}
