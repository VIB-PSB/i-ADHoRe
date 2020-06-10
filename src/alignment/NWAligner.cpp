#include "NWAligner.h"

#include "../GeneList.h"
#include "../ListElement.h"
#include "AlignmentException.h"

NWAligner::NWAligner(const set<Link> &homologs_, int gapScore_, int misScore_,
                     int APScore_, int homScore_, int maxGaps_) :
                     homologs(homologs_), gapScore(gapScore_),
                     misScore(misScore_), APScore(APScore_),
                     homScore(homScore_), maxGaps(maxGaps_)
{

}

int max(int a, int b)
{
    return (a > b) ? a : b;
}

int NWAligner::hit(int x, int y, const vector<GeneList*>& segments, bool &aHit)
{
    int numX = segments[0]->getSize();
    int numID = y*numX + x; // a unique number for each (x,y) pair

    map<int, int>::iterator it = hitmap.find(numID);
    if (it != hitmap.end()) {
        aHit = true;
        return it->second;
    }

    aHit = false;
    return misScore;
}

void NWAligner::createHitmap(const vector<GeneList*>& segments)
{
    hitmap.clear();
    int numX = segments[0]->getSize();

    // create a mapping between numID and x, y
    map<int, int> numIDToCoord;
    for (uint i = 0; i < segments.size(); i++) {
        for (uint j = 0; j < segments[i]->getSize(); j++) {
            ListElement *le = segments[i]->getRemappedElements()[j];
            if (le->isGap()) continue;

            int numID = le->getNumID();
            numIDToCoord[numID] = j;
        }
    }

    set<Link>::iterator it;
    for (it = homologs.begin(); it != homologs.end(); it++) {
        if (it->segmentY != (segments.size() - 1)) continue;

        int x = numIDToCoord[it->geneXID];
        int y = numIDToCoord[it->geneYID];
        int numID = y*numX + x; // a unique number for each (x,y) pair

        int score = (it->isAP) ? APScore : homScore;

        map<int, int>::iterator s = hitmap.find(numID);
        if (s != hitmap.end())
            s->second = max(score, s->second);
        else
            hitmap[numID] = score;
    }
}

void NWAligner::align(vector<GeneList*>& segments)
{
    // create a map of (x, y) scores
    createHitmap(segments);

    // size of the objects
    int Y = segments.size() - 1;
    int max_x = segments[0]->getSize();
    int max_y = segments[Y]->getSize();

    // scoring matrix
    vector< vector<int> > matrix;
    // initialize the matrix
    matrix.resize(max_x+1);
    for (int i = 0; i <= max_x; i++) {
        matrix[i].resize(max_y+1, 0);
    }

    // fill in the first colum and row
    for (int x = 0; x <= max_x; x++)
        matrix[x][0] = x * gapScore;
    for (int y = 0; y <= max_y; y++)
        matrix[0][y] = y * gapScore;

    // fill in the rest of the matrix
    for (int x = 1; x <= max_x; x++) {
        for (int y = 1; y <= max_y; y++) {
            // calculate gap scores
            int up_score = matrix[x][y-1] + gapScore;    // choice 3
            int left_score = matrix[x-1][y] + gapScore;  // choice 2

            // calculate match score                // choice 1
            bool aHit;
            int diagonal_score = matrix[x-1][y-1] + hit(x-1,y-1, segments, aHit);

            // find the best score
            if (diagonal_score >= up_score) {
                if (diagonal_score >= left_score) {     // choice 1
                    matrix[x][y] = diagonal_score;
                }
                else {                                  // choice 2
                    matrix[x][y] = left_score;
                }
            }
            else {
                if (up_score >= left_score) {           // choice 3
                    matrix[x][y] = up_score;
                }
                else {                                  // choice 2
                    matrix[x][y] = left_score;
                }
            }
        }
    }

    int x = max_x;
    int y = max_y;

    vector<int> pathX;
    vector<int> pathY;

    while (x > 0 && y > 0) {
        bool aHit;
        int S = hit(x-1, y-1, segments, aHit);
        if (aHit) {
            pathX.push_back(x-1);
            pathY.push_back(y-1);
        }
        int score = matrix[x][y];
        int scoreDiag = matrix[x-1][y-1];
        int scoreUp = matrix[x][y-1];
        int scoreLeft = matrix[x-1][y];

        if (score == scoreDiag + S) {
            x--;
            y--;
        } else if (score == scoreLeft + gapScore) {
            x--;
        } else if (score == scoreUp + gapScore) {
            y--;
        } else {
            // we should never be here!
            assert(false);
        }
    }

    reverse(pathX.begin(), pathX.end());
    reverse(pathY.begin(), pathY.end());

    int gapsX = 0, gapsY = 0;
    int firstX = 0, firstY = 0;

    for (uint i = 0; i < pathX.size(); i++) {
        // insert gaps in pathX
        if ((pathX[i] + gapsX) < (pathY[i] + gapsY)) {
            int numGaps = pathY[i] + gapsY - pathX[i] - gapsX;
            int lastX = pathX[i] + gapsX;
            for (unsigned int j = 0; j < segments.size() - 1; j++) {
                if (numGaps > maxGaps)
                    throw AlignmentException("alignment failed...too many gaps in profile");

                    segments[j]->introduceGaps(firstX, lastX, numGaps);
            }
            gapsX += numGaps;
        }

        // insert gaps in pathY
        if ((pathX[i] + gapsX) > (pathY[i] + gapsY)) {
            int numGaps = pathX[i] + gapsX - pathY[i] - gapsY;
            int lastY = pathY[i] + gapsY;

            if (numGaps > maxGaps)
                throw AlignmentException("alignment failed...too many gaps in profile");

            segments.back()->introduceGaps(firstY, lastY, numGaps);
            gapsY += numGaps;
        }
        firstX = pathX[i] + gapsX;
        firstY = pathY[i] + gapsY;
    }

    //make sure every list is equal in size -> insert gaps past the end
    unsigned int largest = 0;
    for (unsigned int i = 0; i < segments.size(); i++)
        largest = max(largest, (int)segments[i]->getSize());

    for (unsigned int i = 0; i < segments.size(); i++) {
        uint size = segments[i]->getSize();
        if (size == largest) continue;

        if ((largest - size) > maxGaps)
            throw AlignmentException("alignment failed...too many gaps in profile");

        segments[i]->introduceGaps(size, size, largest - size);
    }
}
