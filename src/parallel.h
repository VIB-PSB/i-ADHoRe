#ifndef PARALLEL_H
#define PARALLEL_H

#ifdef HAVE_MPI
#include <mpi.h>
#define MAX_PROCESSOR_NAME MPI_MAX_PROCESSOR_NAME
#else
#define MAX_PROCESSOR_NAME 10
#endif

#include <vector>
#include <map>
#include <stdint.h>
#include <cstdio>
#include <iostream>

typedef unsigned int uint;

class ParToolBox
{
public:
    /*
     * Initialize MPI
     * @param argc Pointer to the number of command line arguments
     * @param argv Pointer to the command line arguments
     */
    static void init(int *argc, char ***argv);

    /*
     * Destroy MPI
     */
    static void destroy();

    /*
     * Static partition tasks among a set of processes
     * @param weights A vector of weight for each task (input)
     * @param proc A vector which contains the starting position of each
     * task for each process.  The starting position for process 0 is
     * always 0.  If there are more processes than task, then
     * startPos[p] = #tasks for p >= #tasks.
     * (output, existing content will be overwritten)
     * @param nProc Number of processes to partition between (optional)
     */
    static void statPartitionWorkload(const std::vector<double> &weight,
                                      std::vector<int> &startPos,
                                      int nProc);

    /*
     * Dynamic partition of tasks among a set of processes
     */
    static void dynPartitionWorkload(const std::vector< double >& weight,
                                     std::vector<uint>& startPos,
                                     int nProc, double packagePerc,
                                     double minPackagePerc);

    static void statPartWorklSort(const std::vector< uint64_t >& weights,
                                  uint nPackages,
                                  std::vector< uint >& startIndex,
                                  std::vector< uint64_t >& indexToWeight,
                                  std::vector< uint >& indexToList,
                                  std::vector< uint64_t>& cumWeight,
                                  std::map< uint64_t, uint>& cumWeightToIndex);

    /*
     * Get the process name
     */
    static char* getProcName() {
        return procName;
    }

    /*
     * Get the total number of processes
     */
    static uint getNumProcesses() {
        return nProc;
    }

    /*
     * Get the identifier for the current process
     */
    static int getProcID() {
        return thisProc;
    }

    static void silenceOthers(int procID) {
        if (thisProc != procID)
            std::cout.rdbuf(NULL);
    }

    static void silenceAll() {
        sb = std::cout.rdbuf();
        std::cout.rdbuf(NULL);
    }

    static void restoreOutput(int procID) {
        if (thisProc == procID)
            std::cout.rdbuf(sb);
    }

private:

    static char procName[MAX_PROCESSOR_NAME];
    static int nProc;   // number of parallel processes
    static int thisProc;    // identifier for this process
    static std::streambuf *sb;
};

#endif
