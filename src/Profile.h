#ifndef __PROFILE_H
#define __PROFILE_H

class Multiplicon;

#include "GeneList.h"
#include "AnchorPoint.h"
#include "alignComp.h"
#include "Multiplicon.h"
#include "alignment/Node.h"

#include <cassert>

// ============================================================================
// PROFILE EXCEPTION CLASS
// ============================================================================

class ProfileException : public runtime_error
{
    public:
        ProfileException(const string& msg = "") : runtime_error(msg) {}
};

// ============================================================================
// PROFILE CLASS
// ============================================================================

class Profile {

public:

    /**
     * Construct a profile from a multiplicon
     */
    Profile(const Multiplicon& m, unsigned int _id);

    /**
     * Aligns the segments in the profile
     * @param alignMethod Method to be used for aligning
     * @param maxGaps Maximum number of gaps in the alignment
     */
    void align(const AlignmentMethod &alignMethod, int maxGaps);

    /**
     * Cut the aligned multiplicon
     */
    void cutAlignment();

    /**
     * Align this profile using various alignment procedures and log results
     * @param NW Score for Needleman-Wunsch (all links equal weight)
     * @param GG Score for the original Greedy-Graph aligner
     * @param RA Score for the random aligner
     * @param RC Score for the random conflict aligner
     * @param RAC Score for the random active conflict aligner
     * @param LL Score for the longest link aligner
     * @param LLBS Score for the lowest estimated score aligner
     * @param LS Score for the lowest score aligner
     */
    void compareAligners(AlignScore &NW, AlignScore &GG,
                         AlignScore &RA, AlignScore &RC,
                         AlignScore &RAC, AlignScore &LL,
                         AlignScore &LLBS, AlignScore &LS);

    /**
     * Returns the segments of the profile that have been aligned, gaps included
     * @return The segments of the profile that have been aligned
     */
    const vector<GeneList*>& getSegments() const {
        return segments;
    }

    /**
     * Returns the size of the aligned segments
     * @return The size of the aligned segments
     */
    unsigned int getSize() const {
        if (segments.size() > 0)
            return segments[0]->getSize();
        else
            return 0;
    }

    /**
     * Gives the profile a unique id
     * @param id Unique id for the profile
     */
    void setId(unsigned int id) {
        this->id = id;
    }

    /**
     * Get the id of the profile
     * @return The id of the profile
     */
    unsigned int getId() {
        return id;
    }

    /**
     * Get an iterator to the first homologous pair
     * @return An iterator to the first homologous pair
     */
    set<Link>::const_iterator getHomBegin() const {
        return multiplicon.getHomBegin();
    }

    /**
     * Get an iterator past the final homologous pair
     * @return An iterator past the final homologous pair
     */
    set<Link>::const_iterator getHomEnd() const {
        return multiplicon.getHomEnd();
    }

    /**
     * Print a number of aligned anchor points to the screen
     * @param profile Profile under consideration
     */
    static void printProfile(const std::vector<GeneList*>& profile);

    /**
     * Check whether an anchorpoint still exists in the profile
     * @param geneXID Identifier for the X gene
     * @param geneYID Identifier for the Y gene
     */
    bool findAnchorPoint(int geneXID, int geneYID) const;

    /**
     * Returns the level of the profile
     */
    int getLevel() const {
        return segments.size();
    }

    /**
     * Get the parent multiplicon this profile was created from
     * @return The parent multiplicon
     */
    const Multiplicon& getMultiplicon() const {
        return multiplicon;
    }

    /**
     * Checks the alignment of the profile
     * @param minHomologs The minimum number of aligned homologs per segment
     */
    void checkAlignment(int minHomologs);

private:

    /**
     * Delete the memory associated with some segments
     * @param segments Segments under consideration
     */
    static void clearSegments(vector<GeneList*>& segments);

    /**
     * Deep copy some segments
     * @param src Source segments to be copied
     * @param dst Destination where the segments should be copied to
     */
    static void deepCopySegments(const vector<GeneList*>& src,
                                 vector<GeneList*>& dst);

    /**
     * Permutates the unaligned_y_object of the profile according to the
     * orientation and was_twisted attribute of each basecluster of its
     * multiplicon
     */
    void applyYListPermutation();

    /**
     * Mark the homologs that have been aligned
     * @param segments Segments under consideration
     */
    void markAlignedHomologs(vector<GeneList*> &segments);

    /**
     * Checks if the last segment added wasn't completely masked.
     * if so, an exception will be thrown. if not, the entire segment
     * will be completely unmasked
     */
    void checkMasking(bool level2Only);

    /**
     * Calculates the alignment score based on the number of homologs aligned and gaps
     */
    double calculateAlignmentScore(vector<GeneList>& lists);

    void createNodes(const set<Link> &links,
                     const vector<GeneList*> &segments,
                     vector<vector<Node *> >& nodes,
                     bool onlyAP, bool priorityAP = false) const;

    void destroyNodes(vector<vector<Node *> >& nodes) const;

    //////////////
    //ATTRIBUTES//
    //////////////

    //the multiplicon where the profile is created from
    const Multiplicon& multiplicon;

    //profile is valid or not
    bool valid;

    // the segments that have been aligned
    vector<GeneList*> segments;

    //the y genelist where the profile is created from
    GeneList *unaligned_y_list;

    //id
    unsigned int id;

    // alignment statistics
    unsigned int numAlHom, numAP, numAlAP, numHom;

    // alignment statistics per segment
    vector<int> numAlHomV, numAlAPV;
};

#endif
