#include "parallel.h"
#include "Settings.h"
#include "DataSet.h"

#include "debug/FileException.h"

#ifdef DEBUG
#include "debug/stopwatch.h"
#endif

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "util.h"

using std::cout;
using std::cerr;
using std::endl;

int main (int argc, char** argv) {

    ParToolBox::init(&argc, &argv);

    // only show output from process 0
    ParToolBox::silenceOthers(0);

    // check that the number of program arguments is correct
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " [configuration file]" << endl;

        ParToolBox::destroy();
        exit(EXIT_FAILURE);
    }

    // seed the randomizer
    srand((unsigned)time(0));

    cout << endl;
    cout << "This is i-ADHoRe v3.0." << endl;
    cout << "Copyright (c) 2002-2010, Flanders Interuniversity "
        "Institute for Biotechnology, VIB." << endl;
    cout << "Algorithm designed by Klaas Vandepoele, Cedric Simillion, "
        "Jan Fostier, Dieter De Witte,\nKoen Janssens, Sebastian Proost, "
        "Yvan Saeys and Yves Van de Peer." << endl << endl;

    cout << "Process " << ParToolBox::getProcID()+1 << "/"
         << ParToolBox::getNumProcesses() << " is alive on "
         << ParToolBox::getProcName() << "." << endl << endl;

    try {
        Settings settings(argv[1]);
        DataSet dataset(settings);

        dataset.mapGenes();
        dataset.remapTandems();

        // create output directory and generate empty files
        if (ParToolBox::getProcID() == 0) {
            dataset.prepareOutput();
            dataset.outputGenes();
        }

        // create a permutation to order the genelists from big to small
        dataset.sortGeneLists();

        switch (settings.getClusterType()) {
            case Collinear:
                cout << "Collinear Search" << endl;
                break;
            case Cloud:
                cout << "Synthenic Cloud Search" << endl;
                break;
            case Hybrid:
                cout << "Hybrid Search (first Collinear, then Cloud Search on what is remaining in GHM)" << endl;
                break;
        }

        cout << "Level 2 multiplicon detection..."; cout.flush();
        Util::startChrono();
        dataset.parallelLevel2ADHoReDyn();
        cout << "\tdone. (time: " << Util::stopChrono() << "s)" << endl;

        cout << "Profile detection..."; cout.flush();
        if (!settings.verboseOutput())  // enable silenced mode
            ParToolBox::silenceAll();
        Util::startChrono();
        dataset.profileDetection();

        if (!settings.verboseOutput())
            ParToolBox::restoreOutput(0);
        if (settings.verboseOutput())
            cout << "Time for Higher Level Detection: "
                 << Util::stopChrono() << "s." << endl;
        else
            cout << "\t\tdone. (time: " << Util::stopChrono() << "s)" << endl;

        if (ParToolBox::getProcID() == 0)
            dataset.output();
        cout << endl << endl << "All Done!  Bye..." << endl << endl << endl;
    }
    catch (const FileException& fe) {
        cerr << fe.what() << endl << endl;
    }

    ParToolBox::destroy();
    return EXIT_SUCCESS;
}
