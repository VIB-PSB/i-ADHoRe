#ifndef __NWALIGNER_H
#define __NWALIGNER_H

#include "Aligner.h"

class NWAligner : public Aligner
{
public:
    /**
     * Default constructor
     * @param profile Reference to the profile
     * @param gapScore Score for a gap
     * @param misScore Score for aligning two non homologous pairs
     * @param APScore Score for aligning two anchorpoints
     * @param homScore Score for aligning two homologous gene (no AP)
     * @param maxGaps Maximum number of gaps inserted in alignment
     */
    NWAligner(const set<Link> &homologs, int gapScore, int misScore,
              int APScore, int homScore, int maxGaps);

    /**
     * Align a number of segments (the last segment is unaligned)
     * @param segments Reference to the segments to align (input / output)
     */
    void align(vector<GeneList*>& segments);

private:

    /**
     * Check whether two gene in the segments list are pairs
     * @param x x-Coordinate of the gene
     * @param y y-Coordinate of the gene
     * @param segments Segments (input)
     * @param aHit True if the two genes are pairs (output)
     * @return Score for the match / mismatch
     */
    int hit(int x, int y, const vector<GeneList*>& segments, bool &aHit);


    /**
     * Create a map of scores
     * @param segments Segments (input)
     */
    void createHitmap(const vector<GeneList*>& segments);

    const set<Link> &homologs;
    map<int, int> hitmap;

    // NW scoring system
    int gapScore;
    int misScore;
    int APScore;
    int homScore;

    // other parameters
    uint maxGaps;
};

#endif
