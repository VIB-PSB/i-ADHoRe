#include "Settings.h"
#include "AlignDataSet.h"

#include "debug/FileException.h"

#ifdef DEBUG
#include "debug/stopwatch.h"
#endif

#include <iostream>
#include <cstdlib>
#include <cmath>
using std::cerr;
using std::endl;


int main (int argc, char** argv) {
	#ifdef DEBUG
	Stopwatch timer;
	#endif

	cerr << endl;
	cerr << "\33[1;32mThis is i-ALIGN 0.1 (BETA)\33[0m" << endl;
	cerr << "Written by Sebastian Proost based on i-ADHoRE" << endl;
	cerr << "Algorithm designed by Koen Janssens, Cedric Simillion, Klaas Vandepoele, Yvan Saeys and Yves Van de Peer" << endl;
	cerr << "(c) 2008, Flanders Interuniversity Institute of Biotechnology, VIB" << endl << endl << endl;

	if (argc != 2) {
		cerr << "Usage: " << argv[0] << " [configuration file]" << endl << endl << endl;
	}
	else {
		try {
			Settings settings(argv[1]);

			AlignDataSet dataset(settings);

			dataset.mapGenes();

			dataset.remapTandems();

			dataset.align();

			dataset.output();

			cerr << endl << endl << "All Done!" << endl << endl << endl;
		}
		catch (const FileException& fe) {
			cerr << fe.what() << endl << endl;
		}
	}
	#ifdef DEBUG
	timer.stop();

	cerr << "time exceeded: " << timer.time() << endl << endl;
	#endif

	return 0;
}
