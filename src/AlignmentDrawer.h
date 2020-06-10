#ifndef ALIGNMENTDRAWER_H
#define ALIGNMENTDRAWER_H

#include "headers.h"
#include "Profile.h"
#include "DataSet.h"
#include <iostream>

#include "SvgWriter.h"

class Conflict
{
public:

    /**
    * Constructor
    */
    Conflict(int i, int j, int k, int l) : indexOfSeg1(i), posInSeg1(j), indexOfSeg2(k), posInSeg2(l){};

    /**
    * Copy Constructor
    */
    Conflict(const Conflict& c);

    /**
    * Assignment operator
    */
    Conflict& operator=(const Conflict& c);

    int indexOfSeg1;
    int posInSeg1;
    int indexOfSeg2;
    int posInSeg2;

private:
    Conflict(){};
};

class AlignmentDrawer
{
public:

    /**
    * Constructor
    */
    AlignmentDrawer(const Profile* prof);

    /**
     * Constructor used in postprocessor
     */
    AlignmentDrawer(const vector<GeneList*>& segs);

    /**
    * Destructor
    */
    ~AlignmentDrawer();

    bool buildColorMatrix(DataSet* data, int tandemG);

    bool buildColorMatrixPostP(DataSet* data, int tandemG, set<Link>& homologs);

    void generateAlignmentSVG(const std::string filename);

    bool visualizeConflicts;
private:

    /**
    * matrix containing an int >0 for homologs, 0 no homologs, -1 for a gap
    */
    vector<vector<int> > colorMatrix;

    vector<GeneList*> segments;
    const Profile* profPtr;
    vector<Color> colors;
    vector<Conflict> conflicts;

    SvgWriter* svg;

    bool badProfile;
    int numColorsNeeded;
    double border;
    double squareWidth;
    double squareHeight;
    const double horSquareDist;
    const double vertSquareDist;

/**PRIVATE METHODS*/

    AlignmentDrawer() : border(0.0), squareWidth(0.0), squareHeight(0.0),horSquareDist(0.0), vertSquareDist(0.0) {};
    AlignmentDrawer(const AlignmentDrawer& a) : border(0.0), squareWidth(0.0), squareHeight(0.0),horSquareDist(0.0), vertSquareDist(0.0) {};
    void operator=(const AlignmentDrawer& a){};

    /**
    * Resize colorMatrix and put all values to -1 (gap)
    */
    void initializeColorMatrix();

    int findGenePosInSegment(int indexSeg, const Gene& g);

    void generateColorVector();

    void drawConflicts();

    void coordCenterAPij(int i, int j, double& x, double& y) const;

    bool buildColorMatrixOverlappingCode(DataSet* data, int tandemG,
                                         set<Link>::const_iterator start, set<Link>::const_iterator stop);

    void setBorder();
    void setBoxSize();
};

#endif
