#ifndef __GHM_H
#define __GHM_H

#include "GeneList.h"
#include "bmp/bmp.h"
#include "bmp/grafix.h"


class Settings;
class AnchorPoint;
class BaseCluster;
class SynthenicCloud;
class Multiplicon;
class Cluster;




#define Print(x) cout << x << endl;

#ifdef DEBUG
#include "debug/pngwriter.h"
#endif

#include "headers.h"

enum Orientation { OPP_ORIENT = 0, SAME_ORIENT = 1, MIXED_ORIENT = 0};



class GHM
{

public:
    ///////////////////////////////
    //CONSTRUCTORS AND DESTRUCTOR//
    ///////////////////////////////

    /**
    * Constructs a GHM and fills in the matrix
    *
    * @param xObject The GeneList object corresponding to the X-axis
    * @param yObject The GeneList object corresponding to the Y-axis
    */
    GHM(const GeneList& xObject, const GeneList& yObject, bool cloudSearch=false);

    /**
    * Destructor
    */
    virtual ~GHM();

    //////////////////
    //PUBLIC METHODS//
    //////////////////

    /**
    * Runs the basic ADHoRe algorithm on the gene homology matrix.
    *
    * @param sett Settings object with the settings needed to run the algorithm
    */
    void run(const Settings& sett);

    /**
    * This creates every matching point in the matrix,
    *
    */
    virtual void buildMatrix(bool useFamiles=false);
    void buildMatrixFast();

    /**
    * Preprocesses the multiplicons and fills in the vector with all multiplicons
    *
    * @param useFamily True if we're using gene families
    * @param minHomologs Minimum number of homologs per segment
    */
protected:
    virtual void setMultiplicons(bool useFamily,
                                 int minHomologs);
public:
    void getMultiplicons(vector<Multiplicon*>& mps) const;

    /**
    * This fills in the vector with synthenic clouds + adds GeneInformation to APs
    *
    * @param scv Synthenic Cloud vector that gets filled in with every cloud (ptr)
    */
    virtual void getClouds(vector<SynthenicCloud*>& scv) const;

    /**
    * Returns the total number of homologous points in the matrix
    *
    * @return The total number of homologous points in the matrix
    */
    unsigned int getNumberOfPoints() const {
        return count_points[0] + count_points[1];
    }

    /**
    * Returns the total number of the y_list that have a homolog with the x_object
    *
    * @return The total number of genes of the y_object that have a homolog with the x_object
    */
    //unsigned int getNumberOfYHits(); //NOTE this isn't used anywhere!

    /**
    * Returns the length of the X-axis
    *
    * @return The length of the X-axis
    */
    virtual unsigned int getXObjectSize() const {
        return x_object.getSize();
    }

    /**
    * Returns the length of the Y-axis
    *
    * @return The length of the Y-axis
    */
    virtual unsigned int getYObjectSize() const {
        return y_object.getSize();
    }

    /**
    * Returns whether or not in any orientation class there is a value in the
    * cell at a given coordinate (x,y) in a GHM.
    *
    * @param x The x-coordinate
    * @param y The y-coordinate
    *
    * @return Whether or not the specified cell corresponds to a homologous point
    */
    bool identity(int x, int y);

    /**
    * Visualizes the GHM into a png image file, with the name of the lists in the filename
    *
    * @param gapSize The maximum number of space between two homologous points
    * @param outputPath The directory where the GHM should be plotted to
    */
//#ifdef DEBUG
//    void plotGHM(int gapSize, const string& outputPath);
//#endif

    bool isSynthenicCloudSearch() const
    {
      return isCloudSearch;
    }

    /**
    * Generates a bitmap with AP in white, AP in clusters in blue, and filtered AP in red
    */
    void visualizeGHM(const std::string& output_path) const;

    /**
    * In case of hybrid search first AP found in collinear clusters are removed from GHM before alg runs
    */
    void removeAPFromMultiplicons(const vector<Multiplicon*>& mps);

private:
    ///////////////////
    //PRIVATE METHODS//
    ///////////////////

    /**
    * Runs the i-Adhore algorithm for Collinear clusters
    */
    void runCollinear(const Settings& settings);

    /**
    *calculates a number of properties used in the statistical validation of clusters
    */
    void prepareForStatisticalValidation();

    /**
    * Creates initial baseclusters
    * @param gap Gap size (criterium for adding new points to the BC)
    * @param orientation Orientation of the BC
    * @param qValue Qualitiy assessment criterium of the created BC
    */
    void seedBaseClusters(int gap, bool orientation, double qValue);

    /**
     * Creates initial baseclusters starting from a singleton basecluster
     * @param basecluster Basecluster object (input / output)
     * @param gap Gap size (criterium for adding new points to the BC)
     */
    void seedBaseCluster(BaseCluster *basecluster, int gap) const;

    /**
     * Searches for AP to be inserted into a BC that are within gap range
     * @param gap Gap size (criterium for adding new points to the BC)
     * @param clusterOrientation Orientation of the clusters
     * @param ghmOrientation Orientation of the GHM
     * @param qValue Qualitiy assessment criterium of the enriched BC
     */
    void enrichClusters(int gap, bool clusterOrientation,
                        bool ghmOrientation, double qValue);

    /**
     * Merges pairs of baseclusters of the same orientation class
     * @param gap Maximum gap size allowed between the two clusters
     * @param orientation Specifier for the orientation class
     * @param qValue Qualitiy assessment criterium
     */
    void joinClusters(int gap, bool orientation, double qValue);

    /**
     * Searches for the closest clusters within specified orientation classes
     * @param gap Maximum gap size allowed between the two clusters
     * @param orientI Orientation class of first cluster
     * @param orientJ Orientation class of second cluster
     * @param closestI Iterator pointing to the found cluster (output)
     * @param closestJ Iterator pointing to the found cluster (output)
     * @param qValue Qualitiy assessment criterium
     * @return True of two clusters are found, false otherwise
     */
    bool getClosestClusters(int gap, bool orientI, bool orientJ,
                            vector<BaseCluster*>::iterator& closestI,
                            vector<BaseCluster*>::iterator& closestJ,
                            double qValue);

    /**
    *performs statistical filtering on all baseclusters.
    *clusters with a probability higher than the specified cutoff are removed.
    *the probability is stored in the remaining clusters
    */
    void filterBaseClusters(double probCutoff);
    void filterBaseClustersBF(double probCutoff);
    void filterBaseClustersFDR(double probCutoff);

    /**
    *performs metaclustering i.e. grouping clusters together.
    *this is done by first trying to merge clusters from one orientation class
    *with clusters from the other orientation class by twisting them.
    *clusters of the same class are grouped together if they are close enough
    *within the specified clustergap.
    *these groups are called metaclusters and are stored into multiplicon objects.
    */
    void clusterClusters(int clusterGap, double qValue, unsigned int cntAnchorpoints);

/*    /** //NOTE Just doublechecks, no longer useful!
    *the checks if all metaclusters of the individual baseclusters still allows to
    *detect a possible alignment path later on.
    *this is done by checking if there are two clusters of the same orientation
    *that have overlapping coordinates but are not located withing each others
    *interval.
    */
    //void checkMetaclusterConfiguration();

    /**
    *twists the clusters of the specified orientation (mirroring around their central y-value)
    *this changes the individual anchorpoints' y-values
    */
    void twistClusters(bool orientation);

    /**
    *adds a basecluster to the 2 dimensional vector contains baseclusters
    */
    void addCluster(Multiplicon& cluster);

    /**
    *for debugging purposes, this writes out every basecluster to the console
    */
    //void write(); //NOTE DEBUG

//*************
//CLOUD METHODS
//*************

    /**
    * Runs the i-Adhore algorithm for Synthenic clouds 
    */
    void runSyntheny(const Settings& settings);

    /**
    * Creates intial clouds
    * @param gap Gap size (criterium for adding new points to the SC) ->extension over which is the bounding box is resized
    * @param bf if bruteforce=true then the exact distances between the ap is calculated
    */
    void condenseClouds(uint gap, bool bf);

    /**
    * In contrast to the equivalent seedBaseCluster condenseCloud can also inflate nonunit size clouds
    * @param gap Gap size (criterium for adding new points to the SC) ->extension over which is the bounding box is resized
    * @param sCloud Cloud to be inflated
    * @param foundNewAP New AP found during condensing process
    * @param bf if bruteforce=true then the exact distances between the ap is calculated
    */
    void condenseCloud(SynthenicCloud& sCloud, uint gap, vector<AnchorPoint>& foundNewAP, bool bf);



    /**
    * Remove newly found AP from GHM
    * NOTE in case of condenseClouds the first AP will not be in the foundAP queue, it will be manually removed
    * to not disturb the iterators there!
    *@param foundNewAP vector with anchorPoints to be removed
    */
    void removeAddedAnchorPoints(vector<AnchorPoint>& foundNewAP);

    /**
    * Search for APs in rectangular area in GHM 
    *@param loX left boundary
    *@param hiX right boundary
    *@param loY lower boundary 
    *@param hiY upper boundary
    *@param sCloud found AP will be added to sCloud
    */
    void addAPFromSearchBox(int loX, int hiX, int loY, int hiY, SynthenicCloud& sCloud, vector<AnchorPoint>& foundExtraAP);

    //NOTE Because AP in searchbox might still be further away than the gap, we could for the foundAP calculate their
    // actual distance to the cloud! (this is implemented in this function)
    void addAPFromSearchBoxBF(int gap, int loX, int hiX, int loY, int hiY, SynthenicCloud& sCloud, vector<AnchorPoint>& foundExtraAP);


    /**
    * Tries to add new AP to the already existing Clouds
    * @param gap Gap size (criterium for adding new points to the SC)
    * @param bf if bruteforce=true then the exact distances between the ap is calculated
    * complexity is than quadratic but the avalanche effect is definitely absent
    */
    void inflateClouds(uint gap, bool bf);

    /**
    * Combines clouds which lie closely together
    * In inflate and condensecloud AP are added through inflation of the bounding box
    * In mergeClouds the actual distance will be calculated and will have to be smaller than gap
    * @param gap Gap size (criterium for adding new points to the SC)
    */
    void mergeClouds(uint clustergap, bool bf);

    /**
    * The difference with mergeClouds is that always the exact distance will be calculated between the clouds
    * which is alot more time intensive, and not really useful
    */
    void mergeCloudsBruteForce(uint gap, bool bf);

    /**
    * Calculate distances between the corners of the bounding boxes
    *@return minimum corner distance
    */
    uint estimateCloudKspd(SynthenicCloud& scloud1, SynthenicCloud& scloud2) const;

    /**
    * For both Clusters the minimal distance between AP is calculated
    * @param frameThickness Only AP which lie within a distance of frameThickness from the box edge are taken into consideration
    * this might reduce the number of kspds needed to calculate 
    * @Return The returned distance is irrelevant if it is larger than 2*frameThickness, otherwise it returns the minimal distance
    * between two cloud's AP
    * NOTE: don't use when clouds overlap, it will return wrong result
    */
    uint calculateMinimalKspdBetweenClouds(SynthenicCloud& scloud1, SynthenicCloud& scloud2, uint frameThickness) const;

    /**
    * Calculates the exact distance between two clusters by calculating the distance between all AP pairs
    * @return Returns exact distance O(AP^2) time!
    */
    uint calculateMinimalKspdBetweenCloudsBruteForce(SynthenicCloud& scloud1, SynthenicCloud& scloud2) const;

    /**
    * @return true when bounding boxes of two clouds overlap
    */
    bool boxOverlap(SynthenicCloud& sCloud1, SynthenicCloud& sCloud2) const;

    /**
    * Checks whether AP from one cloud reside in the bounding box of another cloud
    * @return true if AP are found in the cross section of both boxes
    *NOTE always first check for boxOverlap!
    */
    bool cloudsOverlap(SynthenicCloud& sCloud1, SynthenicCloud& sCloud2) const;

    /**
    * Remove AP form cloud2 and add them to cloud1, remove cloud2 from vector
    */
    void joinClouds(list<SynthenicCloud*>::iterator cloudIt1, list<SynthenicCloud*>::iterator cloudIt2, bool bf);

    /**
    * Searches in bounding box of 2 merged clouds if there are points in GHM which are not inserted into this new cloud
    * @param sCloud Union of two clouds (note that the difference between the bounding box of sCloud and the union of the bbs
    * of the joined clouds is not empty => this enrichment is necessary!
    */
    void enrichMergedClouds(SynthenicCloud& sCloud);

    /**
    * filtering clouds based on: Density, Binomial distribution or BinomialDistributionCorrected (own idea)
    * BF=Bonferoni multiple hypothesis correction, FDR false discovery rate mhc
    */
    void filterCloudsBinomialD(double probCutoff);
    void filterCloudsBinomialDCorr(double probCutoff);

    void filterCloudsBinomialDBF(double probCutoff);
    void filterCloudsBinomialDCorrBF(double probCutoff);

    void filterCloudsBinomialDFDR(double probCutoff);
    void filterCloudsBinomialDCorrFDR(double probCutoff);
    //void filterCloudsDensityCriterium(double density);

    /**
    * @return number of AP/number of size of GHM 
    */
    double calculateAPDensity() const;

//*************************
//GHM Visualization METHODS
//*************************

    void visualizeSynthenicClouds(const std::string& output_path) const;
    void visualizeBaseClusters(const std::string& output_path) const;   
#ifdef HAVE_PNG
    void visualizeSynthenicCloudsPNG(const std::string& output_path) const;
    void visualizeBaseClustersPNG(const std::string& output_path) const;
#endif

protected:

    ////////////////////////
    //PROTECTED ATTRIBUTES//
    ////////////////////////

    //the genelist object representing the x coordinates
    const GeneList& x_object;

    //the GeneList object representing the y coordinates
    const GeneList& y_object;

    //the current level
    int level;

    //vector of 2 homology matrices (one for every orientation)
    vector<map<int, set<int> > > matrix;

    //vector containing the multiplicons
    vector<Multiplicon*> multiplicons;

    //used for statistical validation of the clusters
    unsigned long long area;

    //the total number of homologous points in the matrix
    unsigned long long count_points[2];

    //2-dimensional vector containing the baseclusters for each orientation class
    vector<BaseCluster*> baseclusters[2];

    // baseclusters that were filtered for each orientation class
    vector<BaseCluster*> filteredBC[2];

    //a boolean meaning if the 2 genelists are identical or not
    bool identical;

    //vector containing a list of genelists, when multiple ghm's are merged
    vector<const GeneList*> x_object_merged;

    bool active;

    //list containing synthenic clouds found in GHM
    //due to cloudmerging it might be more efficient to use a list, since random acces is not needed
    list<SynthenicCloud*> sClouds; 

    vector<SynthenicCloud*> filteredSC;

    //true is cloudSearch, false if collinear search
    bool isCloudSearch;

};

#endif
