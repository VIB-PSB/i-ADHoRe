#include "parallel.h"
#include <cmath>
#include <cassert>
#include <iostream>

using namespace std;

char ParToolBox::procName[MAX_PROCESSOR_NAME] = "localhost";
int ParToolBox::nProc = 1;
int ParToolBox::thisProc = 0;
streambuf* ParToolBox::sb = NULL;

void ParToolBox::init(int *argc, char ***argv)
{
#ifdef HAVE_MPI
    // check whether MPI is already initialized
    int MPI_hasInit;
    MPI_Initialized(&MPI_hasInit);

    // init MPI if necessary
    if (!MPI_hasInit)
        MPI_Init(argc, argv);

    // get process info
    MPI_Comm_rank(MPI_COMM_WORLD, &thisProc);
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);

    int dummy;
    MPI_Get_processor_name(procName, &dummy);
#endif
}

void ParToolBox::destroy()
{
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}

void ParToolBox::statPartitionWorkload(const vector<double> &weight,
                                       vector<int> &startPos, int nProc)
{
    assert (nProc >= 1);

    // prepare the partition vector
    startPos.resize(nProc, weight.size());
    if (weight.empty()) return;

    startPos[0] = 0;

    // calculate the total workload
    double totalWeight = 0.0;
    vector<double>::const_iterator it;
    for (it = weight.begin(); it != weight.end(); it++)
        totalWeight += *it;

    double targetWeight = totalWeight / nProc;
    double currWeight = weight[0];
    int currProc = 0, cnt = 1;

    // partition the workload among the processes
    for (it = weight.begin() + 1 ; it != weight.end(); it++, cnt++) {

        double currDiff = abs(currWeight - targetWeight);
        double newDiff = abs(currWeight + *it - targetWeight);

        if (currDiff < newDiff) {
            currProc++;
            startPos[currProc] = cnt;

            currWeight = *it;

            // check if we are ready here
            if (currProc == (nProc-1)) return;
        } else {
            currWeight += *it;
        }
    }
}

uint getSmallestPackage(const vector<uint64_t> &weightPerPackage)
{
    uint sp = 0;
    uint64_t sw = weightPerPackage[sp];

    for (uint p = 0; p < weightPerPackage.size(); p++) {
        if (weightPerPackage[p] < sw) {
            sw = weightPerPackage[p];
            sp = p;
        }
    }

    return sp;
}

void ParToolBox::statPartWorklSort(const vector< uint64_t >& weights,
                                   uint nPackages,
                                   vector< uint >& startIndex,
                                   vector< uint64_t >& indexToWeight,
                                   vector< uint >& indexToList,
                                   vector< uint64_t>& cumWeight,
                                   map< uint64_t, uint>& cumWeightToIndex)
{
    // reset output variables
    startIndex.clear();
    startIndex.resize(nPackages);
    indexToWeight.clear();
    indexToWeight.resize(weights.size());
    indexToList.clear();
    indexToList.resize(weights.size());
    cumWeight.clear();
    cumWeight.reserve(weights.size());
    cumWeightToIndex.clear();

    // first, sort the weights
    multimap<uint64_t, uint > weightToList;
    for (uint l = 0; l < weights.size(); l++)
        weightToList.insert(pair<uint64_t, uint>(weights[l], l));

    vector<uint64_t> weightPerPackage(nPackages, 0);
    vector<uint> indicesPerPackage(nPackages, 0);

    // calculate the weight per package and the number of indices per package
    multimap<uint64_t, uint>::reverse_iterator it = weightToList.rbegin();
    for ( ; it != weightToList.rend(); it++) {
        uint sp = getSmallestPackage(weightPerPackage);
        weightPerPackage[sp] += it->first;
        indicesPerPackage[sp]++;
    }

    // calculate the start position of the packages
    for (uint p = 0, currStartPos = 0; p < nPackages; p++) {
        startIndex[p] = currStartPos;
        currStartPos += indicesPerPackage[p];
    }

    vector<uint> currPos = startIndex;
    fill(weightPerPackage.begin(), weightPerPackage.end(), 0);

    // finally, calculate the sorted indexToWeight and indexToList
    for (it = weightToList.rbegin() ; it != weightToList.rend(); it++) {
        uint sp = getSmallestPackage(weightPerPackage);
        weightPerPackage[sp] += it->first;
        indexToList[currPos[sp]] = it->second;
        indexToWeight[currPos[sp]] = it->first;
        currPos[sp]++;
    }

    // create the cumulative weight structures
    uint64_t cW = 0;
    for(uint i = 0; i < weights.size(); i++) {
        cW += indexToWeight[i];
        cumWeight.push_back(cW);
        cumWeightToIndex.insert(pair<uint64_t, uint>(cW, i));
    }
}

void ParToolBox::dynPartitionWorkload(const std::vector< double >& weight,
                                      std::vector< uint >& startPos,
                                      int nProc, double packagePerc,
                                      double minPackagePerc)
{
    startPos.clear();

    if (nProc == 1) {
        startPos.push_back(0);
        return;
    }

    // calculate the total workload
    double totalWeight = 0.0;
    vector<double>::const_iterator it;
    for (it = weight.begin(); it != weight.end(); it++)
        totalWeight += *it;

    double minWeight = minPackagePerc * totalWeight / nProc;

    double weightLeft = totalWeight;
    uint currentIndex = 0;

    while (currentIndex < weight.size()) {

        startPos.push_back(currentIndex);

        // determine a good nItems to process based on the weight of the work
        double targetWeight = weightLeft * packagePerc / nProc;

        if (targetWeight < minWeight)
            return;

        double thisWeight = 0.0;
        for ( ; currentIndex < weight.size(); currentIndex++) {
            if (thisWeight >= targetWeight) break;
            thisWeight += weight[currentIndex];
        }

        weightLeft -= thisWeight;
    }
}
