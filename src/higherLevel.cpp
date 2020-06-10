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

void DataSet::profileSearch(int firstItem, int nItems, int threadID,
                            vector<Multiplicon*>& target, vector<SynthenicCloud*>& dummyVar) const
{
    for (int y = firstItem; y < firstItem + nItems; y++) {
        double startTime = Util::getTime();

        int lst = indexToList[y];
        GeneList &gl = *genelists[lst];
        if (gl.getSize() - gl.getNumberOfMaskedElements() >=
            (unsigned int)settings.getAnchorPoints()) {
            GHMProfile ghm (*profileUC, gl);
            ghm.buildMatrix();
            ghm.run(settings);
            ghm.getMultiplicons(target);
        }
    }
}

bool DataSet::getSomeWorkHL(int &firstItem, int &nItems, int threadID,
                          bool augmentWIP, bool getSmallWP)
{
    // make sure no race conditions occur in this part
    pthread_mutex_lock (&queueMutex);

    // if there is no more work, unlock the mutex and exit
    if (currentIndex > finalIndex) {
        pthread_mutex_unlock (&queueMutex);
        return false;
    }

    firstItem = currentIndex;

    // determine a good nItems to process based on the weight of the work
    lluint targetWeight = (totalWeight - weightSoFar) *
        packagePerc / (double)settings.getNumThreads();

    if (targetWeight < minWeight)
        targetWeight = minWeight;

    // it == final element that will be processed
    map<lluint, uint>::iterator it = cumWeightToIndex.lower_bound(weightSoFar + targetWeight);
    if (it == cumWeightToIndex.end()) it--;

    nItems = it->second - firstItem + 1;
    if (firstItem + nItems - 1 >= finalIndex)
        nItems = finalIndex - firstItem + 1;

    weightSoFar = it->first;
    currentIndex += nItems;

    //cerr << "Proc " << ParToolBox::getProcID() << " from " << firstItem << " to " << currentIndex+nItems-1 << endl;

    if (augmentWIP)
        workInProgress++;

    pthread_mutex_unlock (&queueMutex);
    return true;
}

void DataSet::profileDetection()
{
    workFunction = &DataSet::profileSearch;
    getSomeWorkFunction = &DataSet::getSomeWorkHL;

    cout << endl;

    packagePerc = 0.1;
    minPackagePerc = 0.005;

    // create the weights and a mapping list -> pair
    vector<lluint> weights;
    weights.reserve(genelists.size());

    for (unsigned int l = 0; l < genelists.size(); l++)
        weights.push_back(genelists[l]->getRemappedElementsLength());

    vector<lluint> indexToWeight;
    ParToolBox::statPartWorklSort(weights, ParToolBox::getNumProcesses(),
                                  startPos, indexToWeight, indexToList,
                                  cumWeight, cumWeightToIndex);

    createThreadPool();

    for (unsigned int i = 0; i < multiplicons.size(); i++)
        multiplicons_to_evaluate.push_front(multiplicons[i]);

    unsigned int profile_id = 1;

    // =============================================

    double alignTime = 0.0, flushTime = 0.0;

    while (!multiplicons_to_evaluate.empty()) {
        Multiplicon* multiplicon = multiplicons_to_evaluate.front();

        cout << multiplicons_to_evaluate.size() << " multiplicons to evaluate "
             << "- evaluating level " << multiplicon->getLevel()
             << " multiplicon... ";

        // remove the multiplicon that has been evaluated
        multiplicons_to_evaluate.pop_front();

        try {
            Util::startChrono();
            GeneList &lY = *genelists[multiplicon->getYObjectID()];

            multiplicon->createProfile(profile_id);
            const Profile *profile = multiplicon->getProfile();

            if (settings.getCompareAligners())
                multiplicon->compareAligners(NW, GG, RA, RC, RAC, LL, LLBS, LS);

            if ((multiplicon->getLevel() == 2) && (!settings.level2Only())) {
                if (allMasked(lY, multiplicon->getBeginY(), multiplicon->getEndY())) {
                    multiplicon->align(settings.getAlignmentMethod(),
                                       settings.getMaxGapsInAlignment());
                    throw ProfileException("level-2 multiplicon is redundant");
                }
            }

            // masking check
            if (allMasked(lY, multiplicon->getBeginY(), multiplicon->getEndY()))
                throw ProfileException("all elements masked");

            // align the multiplicon in a profile
            multiplicon->align(settings.getAlignmentMethod(),
                               settings.getMaxGapsInAlignment());

            // check the quality of higher-level profiles
            if (multiplicon->getLevel() > 2)
                multiplicon->checkAlignment(settings.getAnchorPoints());

            // mask the X-genelist
            if (multiplicon->getLevel() == 2) {
                GeneList &lX = *genelists[multiplicon->getXObjectID()];
                lX.mask(multiplicon->getBeginX(), multiplicon->getEndX());
            }

            // mask the Y-genelist
            lY.mask(multiplicon->getBeginY(), multiplicon->getEndY());

            double alignmentTime = Util::stopChrono();
            alignTime += alignmentTime;

            if (!settings.level2Only()) {
                vector<Multiplicon*> new_multiplicons;
                parallelProfileSearch(profile, new_multiplicons);
                sortByMultipliconSize(new_multiplicons);
                cout << new_multiplicons.size() << " new multiplicons found.";

                //adding new found multiplicons to the 'to be evaluated' queue
                for (unsigned int i = 0; i < new_multiplicons.size(); i++) {
                    new_multiplicons[i]->setParentID(multipliconID);
                    multiplicons_to_evaluate.push_front(new_multiplicons[i]);
                }
            }

            cout << endl;

            //store the multiplicon as evaluated
            multiplicon->setIsRedundant( false );
            multiplicon->setId(multipliconID++);
            evaluated_multiplicons.push_back(multiplicon);

            profile_id++;

        } catch(const ProfileException& e) {
            alignTime += Util::stopChrono();
            cout << e.what() << endl;

            if ( multiplicon->getLevel() == 2 ) {
                multiplicon->setIsRedundant( true );
                multiplicon->setId(multipliconID++);
                evaluated_multiplicons.push_back(multiplicon);
            }
            else {
                delete multiplicon;
            }
        }

        double startTime = Util::getTime();
        if ((evaluated_multiplicons.size() >= settings.getFlushOutput())
            or (clouds.size()>=settings.getFlushOutput())) {
            if (ParToolBox::getProcID() == 0){
                flushOutput();
                if (settings.writeStatistics()) cout << "WARNING: statistics won't be correct when output flushing is used!!" << endl;
            }

            //delete multiplicons
            vector<Multiplicon*>::const_iterator itM = evaluated_multiplicons.begin();
            for ( ; itM != evaluated_multiplicons.end(); itM++)
                delete (*itM);
            evaluated_multiplicons.clear();

            //delete clouds
            vector<SynthenicCloud*>::const_iterator itC=clouds.begin();
            for ( ; itC != clouds.end(); itC++)
                delete (*itC);
            clouds.clear();
        }

        flushTime += Util::getTime() - startTime;
    }

    // uncomment to get alignment time and flushing times
    /*cerr << "Total alignment time " << alignTime << endl;
    cerr << "Total flushing time " << flushTime << endl;*/

    destroyThreadPool();

    // flush final output
    if (ParToolBox::getProcID() == 0)
    {
        statistics();
        flushOutput();
    }
    deleteEvaluatedMultiplicons();
    deleteEvaluatedClouds();
}

void DataSet::parallelProfileSearch(const Profile *profile,
                                    vector<Multiplicon*>& target)
{
    // the the pointer to the profile under consideration
    profileUC = profile;

    uint thisProc = ParToolBox::getProcID();
    uint nProc= ParToolBox::getNumProcesses();

    currentIndex = startPos[thisProc];
    finalIndex = (thisProc == nProc -1) ? genelists.size() - 1 :
        startPos[thisProc+1] - 1;
    totalWeight = cumWeight[finalIndex];
    weightSoFar = (thisProc == 0) ? 0 : cumWeight[currentIndex - 1];

    minWeight = minPackagePerc * cumWeight.back() /
        (double)settings.getNumThreads() + 1;

    wakeThreads();

    // finish the work packet
    finishWorkPacket();

    // synchronize worker thread to make sure ALL work is finished
    finishWorkerThreads();

#ifdef HAVE_MPI
    int thisBuffSize = Multiplicon::getPackSize(localMultiplicons);
    char *buffer = new char[thisBuffSize];
    Multiplicon::packMultiplicons(localMultiplicons, buffer);

    // Communicate the buffer size to all the processes
    int *buffSize = new int [nProc];
    MPI_Allgather(&thisBuffSize, 1, MPI_INT, buffSize,
                  1, MPI_INT, MPI_COMM_WORLD);

    // Calculate the displacement
    int *displ = new int[nProc]; displ[0] = 0;
    for (int i = 1; i < nProc; i++)
        displ[i] = displ[i-1]+buffSize[i-1];

    // Calculate the total buffer size
    int totalBuffSize = 0;
    for (int i = 0; i < nProc; i++)
        totalBuffSize += buffSize[i];

    char *recvBuffer = new char[totalBuffSize];
    MPI_Allgatherv(buffer, thisBuffSize, MPI_CHAR, recvBuffer,
                   buffSize, displ, MPI_CHAR, MPI_COMM_WORLD);

    // use this construction to make sure that the order in which the
    // multiplicons are stored is the same for all processes
    for (int i = 0; i < nProc; i++) {
        if (i == thisProc) {
            target.insert(target.end(), localMultiplicons.begin(),
                          localMultiplicons.end());
        } else {
            vector<Multiplicon*> mplicons;
            Multiplicon::unpackHLMultiplicons(recvBuffer + displ[i], mplicons,
                                              *profile, genelists, settings.useFamily());
            target.insert(target.end(), mplicons.begin(), mplicons.end());
        }
    }

    delete [] recvBuffer;
    delete [] displ;
    delete [] buffSize;
    delete [] buffer;
#else    // don't use MPI
    target.insert(target.end(), localMultiplicons.begin(),
                  localMultiplicons.end());
#endif

    localMultiplicons.clear();
}
