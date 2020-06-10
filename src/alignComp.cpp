#include "alignComp.h"

void AlignScore::addScore(uint level, uint nAP, uint nAlAP, uint nHom,
                          uint nAlHom, uint lAl, double time)
{
    // if the size of the containers is too small: expand them
    if (level >= numAP.size()) {
        numAP.resize(level+1, 0u);
        numAlAP.resize(level+1, 0u);
        numHom.resize(level+1, 0u);
        numAlHom.resize(level+1, 0u);
        lengthAl.resize(level+1, 0u);
        numProfiles.resize(level+1, 0u);
        timeProfiles.resize(level+1, 0.0);
    }

    numAP[level] += nAP;
    numAlAP[level] += nAlAP;
    numHom[level] += nHom;
    numAlHom[level] += nAlHom;
    lengthAl[level] += lAl;
    numProfiles[level]++;
    timeProfiles[level] += time;
}

void AlignScore::reset()
{
    numAP.clear();
    numAlAP.clear();
    numHom.clear();
    numAlHom.clear();
    lengthAl.clear();
    numProfiles.clear();
    timeProfiles.clear();
}
