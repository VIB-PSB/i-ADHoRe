#ifndef POSTPROCESSOR_H
#define POSTPROCESSOR_H


#include <map>
#include "DataSet.h"
#include "Multiplicon.h"

#define Debug(x) cout << #x << " " << x << endl;

using namespace std;

const string segmentFile="segments.txt";
const string listElementsFile="list_elements.txt";

/**
 * Structure used for storing segment info: identifier, genelistname, genomename,
 * the first gene and the last gene of the segment
 */
struct segInfo
{
    int id;
    string listName;
    string genomeName;
    string firstGene;
    string lastGene;
};

/**
 * Stores info of listelements, being the genename and the position in the alignment 
 */
struct listElInfo
{
    string gene;
    int pos;
};


class PostProcessor {

public:
    /**
     * Constructor
     */
    PostProcessor(int mID, const string& path, DataSet* data, bool useFamily);

    /**
     * Destructor
     */
    ~PostProcessor();

    /**
     * Create SVG file with multiplicon
     */
    void visualizeMultiplicon(int tandemG );

    /**
     * Print the multiplicon like it is stored in the ouput file (for debugging purposes only)
     */
    void printMultiplicon();

private:

    int multipliconID;
    string outputPath;
    DataSet* dataset;
    bool useFamilies;

    vector<GeneList*> segments;
    vector<segInfo> segmentInfo;
    set<Link> homologs;

//PRIVATE METHODS
    vector<segInfo> getSegmentInfo();

    vector<vector<listElInfo> > getListElements(const vector<segInfo>& sInfo);

    /**
     * Matches the adresses of the listElements (stored in GeneLists in DataSet) to the segment loaded
     * from the output file, this way a 2d table is generated containing all segments WITH gaps
     * This method takes into account the possible section inversion of the genelist
     * @param maxPos maximum position of a gene in the segment
     */
    void createSegment(segInfo& sI, vector<listElInfo>& lEs, int maxPos);

    /**
     * Creates a set of links between the listElements which are homologs
     */
    void createLinks();

};

#endif
