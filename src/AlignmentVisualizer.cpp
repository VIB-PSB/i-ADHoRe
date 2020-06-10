#include "Settings.h"
#include <cassert>
#include <string>
#include <map>
#include "DataSet.h"
#include "PostProcessor.h"
using namespace std;

int main (int argc, char** argv) {


    assert(argc==3);
    cout << "i-Visualize module started" << endl;
    Settings settings(argv[1]);
    int mID=atoi(argv[2]);


    DataSet dataset(settings);

    dataset.mapGenes();
    dataset.remapTandems();
    PostProcessor postprocessor(mID,settings.getOutputPath(), &dataset, settings.useFamily());
    cout << "Postprocessor built for multiplicon " << mID << " with tandemGap= " << settings.getTandemGap() << endl;
    postprocessor.visualizeMultiplicon(settings.getTandemGap());
    //postprocessor.printMultiplicon();
    cout << "leaving i-Visualize" << endl;

}
