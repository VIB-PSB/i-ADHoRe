#ifdef HAVE_MPI
    #include <mpi.h>
#endif

#include "DataSet.h"

#include "AnchorPoint.h"
#include "BaseCluster.h"
#include "SynthenicCloud.h"
#include "GeneList.h"
#include "GenePairs.h"
#include "GeneFamily.h"
#include "Multiplicon.h"
#include "Profile.h"
#include "ListElement.h"
#include "Gene.h"
#include "GHMProfile.h"
#include "Settings.h"

#include <cassert>
#include "util.h"
#include "parallel.h"

#include <climits>

#include <pthread.h>

#include "AlignmentDrawer.h"

#include "alignment/GGAligner2.h"

void DataSet::level2ADHoRe(int firstItem, int nItems, int threadID,
                           vector<Multiplicon*> &mpl_output, vector<SynthenicCloud*>& scl_output) const
{
    bool col,cloud;

    switch (settings.getClusterType()) {
        case Collinear: cloud=false; col=true;  break;
        case Cloud:     cloud=true;  col=false; break;
        case Hybrid:    cloud=true;  col=true;  break;
    }

    for (int i = 0; i < workForThread[threadID].size(); i++) {

        uint x = workForThread[threadID][i].first;
        uint y = workForThread[threadID][i].second;

        lluint w = genelists[x]->getSize() * genelists[y]->getSize();

        //cerr << "I am proc " << ParToolBox::getProcID() << ", thread " <<
        //d threadID << ", doing: " << x << " / " << y << ", weight = " << w << endl;

        //necessary to use different vector since output contains also mps of other GHMs
        vector<Multiplicon*> multipliconsColSearch;

        if (col)
        {
            GHM ghm (*genelists[x], *genelists[y],false);
            ghm.buildMatrix(settings.useFamily());

            ghm.run(settings);
            ghm.getMultiplicons(mpl_output);

            if (settings.showGHM(x,y))
                ghm.visualizeGHM(settings.getOutputPath());

            if (cloud)
                ghm.getMultiplicons(multipliconsColSearch);
        }

        if (cloud)
        {
            GHM ghm (*genelists[x], *genelists[y],true);
            ghm.buildMatrix(settings.useFamily());

            //remove AP already in baseClusters from Col Search, in cloudmode nothing is removed
            ghm.removeAPFromMultiplicons(multipliconsColSearch);

            ghm.run(settings);
            ghm.getClouds(scl_output);

            if (settings.showGHM(x,y)) ghm.visualizeGHM(settings.getOutputPath());
        }

        multipliconsColSearch.clear();

        // update x and y
        y++;
        if (y == (int)genelists.size()) {
            x++;
            y = x;
        }
    }

    workForThread[threadID].clear();
}

uint DataSet::getProcForPackage(const vector<uint64_t> &weightPerProc)
{
    uint sp = 0;
    uint64_t sw = weightPerProc[sp];

    for (uint p = 0; p < weightPerProc.size(); p++) {
        if (weightPerProc[p] < sw) {
            sw = weightPerProc[p];
            sp = p;
        }
    }

    return sp;
}

bool DataSet::getSomeWorkL2(int &firstItem, int &nItems, int threadID,
                            bool augmentWIP, bool getSmallWP)
{
    // make sure no race conditions occur in this part
    pthread_mutex_lock (&queueMutex);

    // check whether we're ready
    if (cX >= genelists.size()) {
        pthread_mutex_unlock (&queueMutex);
        return false;
    }

    // determine a good nItems to process based on the weight of the work
    lluint targetWeight = (totalWeight - weightSoFar) * packagePerc /
        (double)(settings.getNumThreads() * ParToolBox::getNumProcesses());

    if (targetWeight < minWeight)
        targetWeight = minWeight;

    lluint thisWeight = 0;

    while (true) {

        uint lX = indexToList[cX];
        uint lY = indexToList[cY];

        GeneList &glX = *genelists[lX];
        GeneList &glY = *genelists[lY];
        lluint weight = glX.getSize() * glY.getSize();
        uint proc = getProcForPackage(weightPerProc);
        weightPerProc[proc] += weight;

        if (proc == ParToolBox::getProcID()) {
            workForThread[threadID].push_back(pair<uint, uint>(lX, lY));
            thisWeight += weight;
        }

        weightSoFar += weight;

        if ((cX - cY) < 2) {
            // check whether we're ready
            if ((cX == cY) && (cX == genelists.size()-1)) {
                cX = genelists.size();
                break;
            }
            int newD = cX+cY+1;
            cX = (newD < genelists.size()) ? newD : genelists.size() - 1;
            cY = (newD < genelists.size()) ? 0 : newD - genelists.size() + 1;
        } else {
            cX--;
            cY++;
        }

        if (thisWeight > targetWeight)
            break;
    }

    if (augmentWIP)
        workInProgress++;

    pthread_mutex_unlock (&queueMutex);
    return true;
}

void DataSet::parallelLevel2ADHoReDyn()
{
    // initialize load distribution parameters
    workFunction = &DataSet::level2ADHoRe;
    getSomeWorkFunction = &DataSet::getSomeWorkL2;
    cX = cY = 0;
    weightPerProc.clear();
    weightPerProc.resize(ParToolBox::getNumProcesses(), 0);
    workForThread = new std::vector<std::pair<uint, uint> >[settings.getNumThreads()];

    packagePerc = 0.05;
    minPackagePerc = 0.005;

    totalWeight = 0;
    for (uint x = 0; x < genelists.size(); x++) {
        lluint sx = genelists[x]->getRemappedElementsLength();
        for (unsigned int y = x; y < genelists.size(); y++) {
            lluint sy = genelists[y]->getRemappedElementsLength();
            totalWeight += sx * sy;
        }
    }

    minWeight = minPackagePerc * totalWeight /
        ((double)settings.getNumThreads() * ParToolBox::getNumProcesses()) + 1;
    weightSoFar = 0;

    createThreadPool();
    wakeThreads();

    // finish the work packet
    finishWorkPacket();
    finishWorkerThreads();

#ifdef HAVE_MPI
    // serialize the locally found multiplicons and clouds (M multiplicon, C clouds)
    int thisBuffSizeM = Multiplicon::getPackSize(localMultiplicons);
    char *bufferM = new char[thisBuffSizeM];
    Multiplicon::packMultiplicons(localMultiplicons, bufferM);

    int thisBuffSizeC= SynthenicCloud::getPackSize(localClouds);
    char *bufferC = new char[thisBuffSizeC];
    SynthenicCloud::packSynthenicClouds(localClouds,bufferC);

    int thisProc = ParToolBox::getProcID();
    int nProc = ParToolBox::getNumProcesses();

    // Communicate the buffer size to all the processes
    int *buffSizeM = new int [nProc];
    MPI_Allgather(&thisBuffSizeM, 1, MPI_INT,
                  buffSizeM, 1, MPI_INT, MPI_COMM_WORLD);
    int *buffSizeC = new int [nProc];
    MPI_Allgather(&thisBuffSizeC, 1, MPI_INT,
                  buffSizeC, 1, MPI_INT, MPI_COMM_WORLD);

    // Calculate the displacement
    int *displM = new int[nProc]; displM[0] = 0;
    int *displC = new int[nProc]; displC[0] = 0;
    for (int i = 1; i < nProc; i++) {
        displM[i] = displM[i-1]+buffSizeM[i-1];
        displC[i] = displC[i-1]+buffSizeC[i-1];
    }
    // Calculate the total buffer size
    int totalBuffSizeM = 0;
    int totalBuffSizeC = 0;
    for (int i = 0; i < nProc; i++) {
        totalBuffSizeM+= buffSizeM[i];
        totalBuffSizeC+= buffSizeC[i];
    }
    char *recvBufferM = new char[totalBuffSizeM];
    char *recvBufferC = new char[totalBuffSizeC];

    MPI_Allgatherv(bufferM, thisBuffSizeM, MPI_CHAR,
                   recvBufferM, buffSizeM, displM, MPI_CHAR, MPI_COMM_WORLD);
    MPI_Allgatherv(bufferC, thisBuffSizeC, MPI_CHAR,
                   recvBufferC, buffSizeC, displC, MPI_CHAR, MPI_COMM_WORLD);

    // use this construction to make sure that the order in which the
    // multiplicons are stored is the same for all processes
    for (int i = 0; i < nProc; i++) {
        if (i == thisProc) {
            addMultiplicons(localMultiplicons);
            addClouds(localClouds);
        } else {
            vector<Multiplicon*> mplicons;
            Multiplicon::unpackL2Multiplicons(recvBufferM + displM[i], mplicons,
                                              genelists, settings.useFamily());
            addMultiplicons(mplicons);

            vector<SynthenicCloud*> sclouds;
            SynthenicCloud::unpackSynthenicClouds(recvBufferC + displC[i], sclouds);
            addClouds(sclouds);
        }
    }

    delete [] recvBufferM;
    delete [] displM;
    delete [] buffSizeM;
    delete [] bufferM;

    delete [] recvBufferC;
    delete [] displC;
    delete [] buffSizeC;
    delete [] bufferC;
#else    // if we don't have MPI
    addMultiplicons(localMultiplicons);
    addClouds(localClouds);
#endif

    localMultiplicons.clear();
    localClouds.clear();
    sortByMultipliconSize(multiplicons);

    destroyThreadPool();

    delete [] workForThread;
}
