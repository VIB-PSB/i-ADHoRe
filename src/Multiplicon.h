#ifndef __MULTIPLICON_H
#define __MULTIPLICON_H

#include "Cluster.h"
#include "headers.h"
#include "GeneList.h"
#include "Settings.h"
#include "alignComp.h"
#include <cassert>
#include <set>

class GeneList;
class Profile;
class BaseCluster;
class Settings;
class AlignmentException;

// ============================================================================
// LINK CLASS
// ============================================================================

class Link
{

public:
    Link(int segmentX_, int segmentY_, int geneXID_,
         int geneYID_, bool isAP_ = false) :
            segmentX(segmentX_), segmentY(segmentY_), geneXID(geneXID_),
            geneYID(geneYID_), isAP(isAP_), isAligned(false) {
        assert(segmentX_<segmentY_);
    }

    int segmentX;
    int segmentY;
    int geneXID;
    int geneYID;
    bool isAP;
    bool isAligned;

    friend bool operator<(const Link &lhs, const Link &rhs);
};

inline bool operator<(const Link &lhs, const Link &rhs)
{
    if (lhs.segmentX != rhs.segmentX)
        return (lhs.segmentX < rhs.segmentX);
    if (lhs.segmentY != rhs.segmentY)
        return (lhs.segmentY < rhs.segmentY);
    if (lhs.geneXID != rhs.geneXID)
        return (lhs.geneXID < rhs.geneXID);
    return (lhs.geneYID < rhs.geneYID);
}

// ============================================================================
// MULTIPLICON CLASS
// ============================================================================

class Multiplicon: public Cluster {

public:
    ///////////////////////////////
    //CONSTRUCTORS AND DESTRUCTOR//
    ///////////////////////////////

    /**
    * Constructs a multiplicon object with the corresponding genelist IDs
    */
    Multiplicon(int xObjectID, int yObjectID, int level);

    /**
     * Constructs a multiplicon object from a packed buffer
     * @param buffer Buffer that contains serialized multiplicon object
     * @param genelists Reference to the gene lists
     * @param useFamily True if we're using gene families
     */
    Multiplicon(const char *buffer, const vector<GeneList* >& genelists,
                bool useFamily);

    /**
     * Constructs a multiplicon object from a packed buffer
     * @param buffer Buffer that contains serialized multiplicon object
     * @param xObject xObject profile used to create this multiplicon
     * @param genelists Reference to the gene lists
     * @param useFamily True if we're using gene families
     */
    Multiplicon(const char *buffer, const Profile &xObject,
                const vector<GeneList* >& genelists, bool useFamily);

    /**
     * Destructor
     */
    ~Multiplicon();


    //////////////////
    //PUBLIC METHODS//
    //////////////////

    /*
    * Returns a pointer to the vector containing baseclusters
    */
    const vector<BaseCluster*>& getBaseClusters() const {
        return baseclusters;
    }

    /*
     * Adds the baseclusters from another multiplicon to this multiplicon.
     * the other multiplicon is also made empty
     */
    void addBaseClusters(Multiplicon& multiplicon);

    /*
    * Adds one single basecluster to this multiplicon
    */
    void addBaseCluster(BaseCluster& cluster);

    /*
    * Clears the multiplicon from all baseclusters but does not free the memory!!
    */
    void clear() {
        baseclusters.clear();
    }

    /*
    *calculates and returns the total number of anchorpoints from all baseclusters
    *in this multiplicon
    */
    unsigned int getCountAnchorPoints() const;

    /*
    *returns the lowest x-value from all anchorpoints
    */
    unsigned int getLowestX() const;

    /*
    *returns the highest x-value from all anchorpoints
    */
    unsigned int getHighestX() const;

    /*
    *returns the lowest y-value from all anchorpoints
    */
    unsigned int getLowestY() const;

    /*
    *returns the highest y-value from all anchorpoints
    */
    unsigned int getHighestY() const;

    /*
    *returns the level from the multiplicon
    */
    int getLevel() const {
        return level;
    }

    /**
     * Create a profile from this multiplicon
     * @param profileID Unique identifier for the profile
     */
    void createProfile(int profileID);

    /**
     * Align the profile associated with this multiplicon
     * @param alignMethod Method to be used for aligning
     * @param maxGaps Maximum number of gaps in the alignment
     */
    void align(const AlignmentMethod &alignMethod, int maxGaps);

    /**
     * Check the alignment of the profile associated with this multiplicon
     * @param minHomologs The minimum number of aligned homologs per segment
     */
    void checkAlignment(int minHomologs);

    /**
     * Returns the profile created from this multiplicon
     * @return A pointer to the profile in this multiplicon
     */
    const Profile* getProfile() const {
        return profile;
    }

    /**
     * For debugging purposes, this writes out every basecluster to the console
     */
    void write () const;

    /**
     * Gives the multiplicon an id to be possible to identify it
     * @param ID Unique identifier for this multiplicon
     */
    void setId(unsigned int ID) {
        this->multipliconID = ID;
    }

    /**
     * Returns the ID of the multiplicon
     * @return The ID of the multiplicon
     */
    unsigned int getId() const {
        return multipliconID;
    }

    void setIsRedundant(bool isRedundant_) {
        isRedundant = isRedundant_;
    }

    bool getIsRedundant() const {
        return isRedundant;
    }

    /**
     * Calculate and return the package size (in byte)
     */
    int getPackSize() const;

    /**
     * Pack data in a char stream
     * @param buffer Pre-allocated buffer to store the data in
     * @return The number of chars in the buffer that were used
     */
    int pack(char *buffer) const;

    /**
     * Get the number of chars to pack a number of multiplicons
     * @param mplicons Multiplicons to be packed
     * @return Number of chars to pack the multiplicons
     */
    static int getPackSize(const vector<Multiplicon*> &mplicons);

    /**
     * Pack multiplicons in a buffer
     * @param mplicons Multiplicons to be packed
     * @param buffer Pre-allocated buffer (output)
     */
    static void packMultiplicons(const vector<Multiplicon*> &mplicons,
                                 char* buffer);

    /**
     * Unpack multiplicons in a buffer
     * @param buffer Buffer that contains packed multiplicons
     * @param mplicons Multiplicons to be unpacked (output)
     * @param genelists Reference to the genelists
     * @param useFamily True if we're using gene families
     */
    static void unpackL2Multiplicons(const char *buffer,
                                     vector<Multiplicon*> &mplicons,
                                     const vector<GeneList *>& genelists,
                                     bool useFamily);

    /**
     * Unpack multiplicons in a buffer
     * @param buffer Buffer that contains packed multiplicons
     * @param mplicons Multiplicons to be unpacked (output)
     * @param profile Reference to the profile used to create this object
     * @param genelists Reference to the genelists
     * @param useFamily True if we're using gene families
     */
    static void unpackHLMultiplicons(const char *buffer,
                                     vector<Multiplicon*> &mplicons,
                                     const Profile &profile,
                                     const vector<GeneList *>& genelists,
                                     bool useFamily);

    friend bool operator==(const Multiplicon &lhs, const Multiplicon &rhs);

    friend bool operator!=(const Multiplicon &lhs, const Multiplicon &rhs) {
        return !(lhs == rhs);
    }

    void setParentID(uint ID) {
        parentID = ID;
    }

    uint getParentID() const {
        return parentID;
    }

    /**
     * Get an iterator to the first homologous pair
     * @return An iterator to the first homologous pair
     */
    set<Link>::const_iterator getHomBegin() const {
        return homologs.begin();
    }

    /**
     * Get an iterator past the final homologous pair
     * @return An iterator past the final homologous pair
     */
    set<Link>::const_iterator getHomEnd() const {
        return homologs.end();
    }

    /**
     * Get the xSegments
     */
    const vector<GeneList*>& getXSegments() const {
        return xSegments;
    }

    /**
     * Get the ySegments
     */
    GeneList* getYSegment() const {
        return ySegment;
    }

    /**
     * Extract the relevant x-object portion of a genelist
     * @param xObject Genelist that was used as x-object in GHM L2 search
     */
    void extractXObject(const GeneList &xObject);

    /**
     * Extract the relevant x-object portion of a profile
     * @param xObject Profile that was used as x-object in GHM profile search
     */
    void extractXObject(const Profile &xObject);

    /**
     * Extract the relevant y-object portion of a genelist
     * @param yObject Genelist that was used in the GHM profile search
     */
    void extractYObject(const GeneList& yObject);

    /**
     * Create a list homologous links
     * @param useFamily True if we're using gene families
     * @param minHomologs The minimum number of homologs per segment
     * @return True if each segment contains at least minHomologs homologs
     */
    bool createHomologs(bool useFamily, int minHomologs);

    /**
     * Get a reference to the homologous genes
     */
    const std::set<Link>& getHomologs() const {
        return homologs;
    }

    /**
     * Align this profile using various alignment procedures and log results
     * @param NW Score for Needleman-Wunsch (all links equal weight)
     * @param GG Score for the original Greedy-Graph aligner
     * @param RA Score for the random aligner
     * @param RC Score for the random conflict aligner
     * @param RAC Score for the random active conflict aligner
     * @param LL Score for the longest link aligner
     * @param LLBS Score for the lowest lower bound score aligner
     * @param LS Score for the lowest score aligner
     */
    void compareAligners(AlignScore &NW, AlignScore &GG,
                         AlignScore &RA, AlignScore &RC,
                         AlignScore &RAC, AlignScore &LL,
                         AlignScore &LLBS, AlignScore &LS);

private:

   /**
    * Copy constructor
    */
    Multiplicon(const Multiplicon& multiplicon) {};

    /**
    * Assignment operator
    */
    void operator=(const Multiplicon& multiplicon) {};

    /**
     * Remove the homologous links that are not contained in the xSegment
     */
    void pruneLinks();

    /**
     * Add the homologous links due to the new ySegment
     */
    void addLinks(bool useFamily);

    //////////////
    //ATTRIBUTES//
    //////////////

    // the level of the multiplicon
    int level;

    // vector containing the baseclusters inside this multiplicon
    vector<BaseCluster*> baseclusters;

    // the segments contained in the multiplicon
    vector<GeneList*> xSegments;

    // the y-segment
    GeneList *ySegment;

    // the profile corresponding to this multiplicon
    Profile* profile;

    // set of homologous pairs within this multiplicon
    std::set<Link> homologs;

    // a boolean indicating whether or not a multiplicon is redundant
    bool isRedundant;

    uint multipliconID;
    uint parentID;
};

#endif
