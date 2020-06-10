#ifndef __GGAligner_H
#define __GGAligner_H

class Node;
class GeneList;
#include "Aligner.h"
#include "../headers.h"

class GGAligner : public Aligner
{

public:
        /**
         * Creates an GGAligner object
         * @param gapSize The maximal space between two homologous elements
         */
	GGAligner (int gapSize);

        /**
         * Aligns the segments and fills in the aligned elements for every genelist into the vector
         * @param segments The vector that is filled in with the aligned elements for every genelist
         */
	void align(vector<GeneList*>& segments);

	/**
	 * Return the number of aligned anchor points after the aligning
	 * @return The number of aligned anchor points after the aligning
	 */
	int getNumAlignedPoints() const;

private:

        //segments that are getting aligned
	vector<vector<Node> > segments;

	//false is a gap, true is not a gap
	vector<vector<int> > alignings;

	//maximal distance between 2 equal nodes
        int gap;

	friend std::ostream& operator<<(std::ostream& os, const GGAligner& a);
};

#endif
