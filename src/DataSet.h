#ifndef __DATASET_H
#define __DATASET_H

class Settings;
class GeneList;
class Multiplicon;
class SynthenicCloud;
class GenePairs;
class ListElement;
class Gene;
class Profile;
class GHMProfile;
class ListFile;

class PackingTest;
class GapsTest;

class DataSet;

#include "debug/FileException.h"

#include "headers.h"
#include "alignComp.h"

#include <stdint.h>

extern "C" void* startThread(void *args);

typedef uint64_t lluint;

typedef struct {
    DataSet *dataset;
    int threadID;
} ThreadArgs;

class DataSet {

public:
    ///////////////////////////////
    //CONSTRUCTORS AND DESTRUCTOR//
    ///////////////////////////////

    /*
    *constructs a dataset object and gets the properties out of the settings object.
    *creates genelist objects for every listfile
    */
    DataSet(const Settings& sett);

    /*
    *destructor
    */
    ~DataSet();


    //////////////////
    //PUBLIC METHODS//
    //////////////////

    /*
    *reads alle genepairs from the blast table file and removes the not
    *necessary pairs
    */
    void mapGenes();
    void getGenePairs();
    void getGeneFamilies();

    /*
    *remaps pairs and indirect pairs for every genelist on one gene
    */
    void remapTandems();

    /**
     * Run a portion of the level-2 ADHoRe in thread-safe manner
     * @param beginX First x-index of the genelist to consider
     * @param beginY First y-index of the genelist to consider
     * @param nItem Number of items to process
     * @param output Found multiplicons (output)
     */
    void level2ADHoRe(int firstItem, int nItems, int threadID,
                      vector<Multiplicon*> &output, vector<SynthenicCloud*>& scl_output) const;

    void wakeThreads();

    void finishWorkPacket();

    void finishWorkerThreads();

    /**
     * Start the dynamic load balanced MPI/multi-threaded version of the level-2 ADHoRe
     */
    void parallelLevel2ADHoReDyn();

    /**
     * Get some workload to process in a thread
     * @param firstItem First item to process (output)
     * @param nItem Number of items to process (output)
     * @return True if there is still work to perform, false otherwise
     */
    bool getSomeWorkHL(int &firstItem, int &nItems, int threadID,
                     bool augmentWIP, bool getSmallWP = false);

    bool getSomeWorkL2(int &firstItem, int &nItems, int threadID,
                       bool augmentWIP, bool getSmallWP = false);

    //bool getSomeWorkL2(int threadID, bool augmentWIP);

    /*
    *detects higher level multiplicons by creating profiles from each multiplicon
    */
    void profileDetection();

    /*
    *outputs the results into text files
    */
    void output();

    /**
     * Flushes the output so far into text files
     */
    void flushOutput();

    /**
     * Create output directory and empty output files
     */
    void prepareOutput();

    /**
     * Output the genes.txt output file
     */
    void outputGenes();

    /*
    *produces a log file with general statistics
    */
    void statistics();

    /**
     * Convert an index to an x and y position in a upper triangular matrix
     * @param i Index (input)
     * @param N Dimension of the matrix
     * @param x x-position (output)
     * @param y y-position (output)
     */
    static void indexToXY(int i, int N, int &x, int &y);

    /**
     * Find a certain gene given its numerical geneID
     * @param geneID Numerical and unique gene identifier
     * @return Reference to the corresponding gene
     */
    const Gene& getGene(unsigned int geneID);

    void sortGeneLists();

    //reads genomes from genelists --> TODO: fill this list when the ini file is read !
    void getGenomes();

    void (DataSet::*workFunction)(int firstItem, int nItems, int nThreads,
                                  vector<Multiplicon*> & mpl_output, vector<SynthenicCloud*>& scl_output) const;
    bool (DataSet::*getSomeWorkFunction)(int &firstItem, int &nItems, int threadID,
                                 bool augmentWIP, bool getSmallWP);

    GeneList* getGeneList(const string& listName, const string& genomeName) const;
    
private:
    ///////////////////
    //PRIVATE METHODS//
    ///////////////////

    /**
     * Check whether a portion of a gene list is completely masked
     * @param list Reference to the list under consideration
     * @param begin Offset to the first element in the gene list
     * @param end Offset to the final element in the gene list
     */
    bool allMasked(const GeneList &list, int begin, int end);

    /*
     * Unpack a buffer of packed multiplicons
     */
    void unpackMultiplicons(const char *packedMult);

    /*
    *adds the vector of multiplicons to the existing vector of multiplicons
    */
    void addMultiplicons(const vector<Multiplicon*>& multiplicons);

    void addMultiplicons2(const vector<Multiplicon*>& multiplicons);

    void addClouds(const vector<SynthenicCloud*>& clouds);
    /*
    *sorts the vector of multiplicons by lowest discrete pseudo distance first, and then by
    *the lowest multiplicon size
    */
    void sortByMultipliconSize(vector<Multiplicon*>& multiplicons) const;

    /**
     * Search a profile against a number of genelists
     * @param firstItem First genelist to consider
     * @param nItem Number of genelists to consider
     * @param target Vector to store the found multiplicons (output)
     * @param dummyVar this variable is just added to make be compatible with function pointer
     */
    void profileSearch (int firstItem, int nItems, int threadID,
                        vector<Multiplicon*>& target, vector<SynthenicCloud*>& dummyVar) const;

    void parallelProfileSearch(const Profile *profile,
                               vector<Multiplicon*>& target);

    //returns the size of the specified remapped genome
    int getRemappedGenomeSize(string genome);
    //returns the size of the specified gene list
    int getRemappedGeneListSize(string genome, string list);

    //sub-functions of the statistics
    void getDuplicatedPortions();

    void getCollinearPortions();

    /**
     * Delete the memory associated with evaluated multiplicons
     */
    void deleteEvaluatedMultiplicons();
    void deleteEvaluatedClouds();

    /**
    * Visualize Alignments with AlignmentDrawer class (SVG)
    */
    void visualizeAlignedProfiles();

    /**
    * Prints profiles via GeneID followed by links (txt format)
    */
    void printProfiles();

    uint getProcForPackage(const vector<uint64_t> &weightPerProc);

    void runWorkFunction(int firstItem, int nItems, int threadID,
                         vector<Multiplicon*> &multiplicons, vector<SynthenicCloud*>& clouds) {
        (this->*workFunction)(firstItem, nItems, threadID, multiplicons,clouds);
    }

    bool runGetSomeWork(int &firstItem, int &nItems, int threadID,
                        bool augmentWIP, bool getSmallWP = false) {
        return (this->*getSomeWorkFunction)(firstItem, nItems, threadID,
                                    augmentWIP, getSmallWP);
    }

    void createThreadPool();
    void destroyThreadPool();


    void flushCollinear();
    void flushClouds();

    int max(int a, int b) {
        return (a > b) ? a : b;
    }

    int min(int a, int b) {
        return (a > b) ? b : a;
    }



    //////////////
    //ATTRIBUTES//
    //////////////

    //settings object containing all information necessary for the algorithm
    const Settings& settings;

    void readDataSet();
    //vector containing all genomes, is filled when getGenomes is called
    list<string> genomes;
    //vector containing all the genelists
    vector<GeneList*> genelists;
    //vector containing all calculated multiplicons by the level 1 algorithm
    vector<Multiplicon*> multiplicons;
    //queue containing all multiplicons that need to be evaluated into a profile
    deque<Multiplicon*> multiplicons_to_evaluate;
    //vector containing all multiplicons that have been evaluated
    vector<Multiplicon*> evaluated_multiplicons;

    vector<Multiplicon *> masterTMultiplicons;
    vector<SynthenicCloud*> masterTClouds;

    GenePairs *genepairs;

    // threading information

    pthread_t *threads;
    ThreadArgs *threadArgs;
    pthread_cond_t workerCond, masterCond;
    pthread_mutex_t workerMutex, queueMutex, mpcMutex;

    // load balancing parameters
    double packagePerc;
    double minPackagePerc;

    lluint totalWeight;
    lluint weightSoFar;
    lluint minWeight;

    // level 2 load balancing variables
    uint cX, cY;
    std::vector<lluint> weightPerProc;
    std::vector<uint> indexToList;
    std::vector<std::pair<uint, uint> > *workForThread;

    // higher level load balancing variables
    int currentIndex, finalIndex;
    vector<Multiplicon*> localMultiplicons;
    vector<SynthenicCloud*> localClouds;

    int workInProgress;
    bool destroyTP;

    double flushTime;

    // higher level load balancing
    std::vector<uint> startPos;
    std::vector<lluint> cumWeight;
    std::map<lluint, uint> cumWeightToIndex;

    // vector holding the allocated profiles
    // vector<const Profile *> profiles;

    // aux map for mapping between geneID and the gene itself
    map<int, GeneList*> geneIDmap;

    const Profile *profileUC;

    // output files
    std::string geneFile;
    std::string multipliconsFile;
    std::string baseclustersFile;
    std::string anchorpointsFile;
    std::string mplpairsFile;
    std::string segmentsFile;
    std::string listElementsFile;
    std::string alignmentFile;
    std::string synthenicCloudsFile;
    std::string cloudAnchorPointsFile;

    uint multipliconID;
    uint baseclusterID;
    uint anchorpointID;
    uint pairID;
    uint segmentID;
    uint elementID;
    uint cloudID;

    AlignScore NW, GG, RA, RC, RAC, LL, LLBS, LS;

    friend class PackingTest;
    friend class GapsTest;

    friend void* startThread(void *args);

    int nThreads;   // number of spawned (i.e. extra) threads

    //vector containing all SynthenicClouds (level2)
    vector<SynthenicCloud*> clouds;
};

#endif
