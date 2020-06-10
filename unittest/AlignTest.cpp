#include <gtest/gtest.h>
#include <vector>
#include "../src/alignment/GGAligner2.h"
#include "../src/alignment/Node.h"
#include "../src/GeneList.h"

using namespace std;

class AlignTest : public ::testing::Test
{
protected:
    virtual void SetUp();
    virtual void TearDown();

    GGAligner2 *aligner;
    vector<vector<Node*> > nodes;
    vector<GeneList*> segments;
};

void AlignTest::SetUp()
{
    int numSeg = 6;
    int segLength = 2;

    for (int i = 0; i < numSeg; i++) {
        segments.push_back(new GeneList(segLength));
        nodes.push_back(vector<Node*>(segLength));
        for (int j = 0; j < segLength; j++)
            nodes[i][j] = new Node(&segments[i]->getLe(j), i, j, j);
    }

    // 5 - 2

    // 1-0
    nodes[0][0]->addLink(*nodes[1][0], 1.0);
    nodes[1][0]->addLink(*nodes[0][0], 1.0);

    nodes[0][0]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[0][0], 1.0);

    nodes[0][0]->addLink(*nodes[3][0], 1.0);
    nodes[3][0]->addLink(*nodes[0][0], 1.0);

    nodes[0][0]->addLink(*nodes[4][0], 1.0);
    nodes[4][0]->addLink(*nodes[0][0], 1.0);

    // 0-1
    nodes[0][1]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[0][1], 1.0);

    nodes[0][1]->addLink(*nodes[2][0], 1.0);
    nodes[2][0]->addLink(*nodes[0][1], 1.0);

    nodes[0][1]->addLink(*nodes[3][1], 1.0);
    nodes[3][1]->addLink(*nodes[0][1], 1.0);

    nodes[0][1]->addLink(*nodes[4][1], 1.0);
    nodes[4][1]->addLink(*nodes[0][1], 1.0);

    // 1-0
    nodes[1][0]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[1][0], 1.0);

    nodes[1][0]->addLink(*nodes[3][0], 1.0);
    nodes[3][0]->addLink(*nodes[1][0], 1.0);

    nodes[1][0]->addLink(*nodes[4][0], 1.0);
    nodes[4][0]->addLink(*nodes[1][0], 1.0);

    // 1-1
    nodes[1][1]->addLink(*nodes[2][0], 1.0);
    nodes[2][0]->addLink(*nodes[1][1], 1.0);

    nodes[1][1]->addLink(*nodes[3][1], 1.0);
    nodes[3][1]->addLink(*nodes[1][1], 1.0);

    nodes[1][1]->addLink(*nodes[4][1], 1.0);
    nodes[4][1]->addLink(*nodes[1][1], 1.0);

    // 2-0
    nodes[2][0]->addLink(*nodes[3][1], 1.0);
    nodes[3][1]->addLink(*nodes[2][0], 1.0);

    nodes[2][0]->addLink(*nodes[4][1], 1.0);
    nodes[4][1]->addLink(*nodes[2][0], 1.0);

    // 2-1
    nodes[2][1]->addLink(*nodes[3][0], 1.0);
    nodes[3][0]->addLink(*nodes[2][1], 1.0);

    nodes[2][1]->addLink(*nodes[4][0], 1.0);
    nodes[4][0]->addLink(*nodes[2][1], 1.0);

    // 3-0
    nodes[3][0]->addLink(*nodes[4][0], 1.0);
    nodes[4][0]->addLink(*nodes[3][0], 1.0);

    // 3-1
    nodes[3][1]->addLink(*nodes[4][1], 1.0);
    nodes[4][1]->addLink(*nodes[3][1], 1.0);

    // 5-0
    nodes[5][0]->addLink(*nodes[0][1], 1.0);
    nodes[0][1]->addLink(*nodes[5][0], 1.0);

    nodes[5][0]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[5][0], 1.0);

    nodes[5][0]->addLink(*nodes[2][0], 1.0);
    nodes[2][0]->addLink(*nodes[5][0], 1.0);

    nodes[5][0]->addLink(*nodes[3][1], 1.0);
    nodes[3][1]->addLink(*nodes[5][0], 1.0);

    nodes[5][0]->addLink(*nodes[4][1], 1.0);
    nodes[4][1]->addLink(*nodes[5][0], 1.0);


    // 2 - 3 (2 direct conflicts vs. one)
   /* nodes[0][0]->addLink(*nodes[1][2], 1.0);
    nodes[1][2]->addLink(*nodes[0][0], 1.0);

    nodes[1][0]->addLink(*nodes[0][1], 1.0);
    nodes[0][1]->addLink(*nodes[1][0], 1.0);

    nodes[1][1]->addLink(*nodes[0][2], 1.0);
    nodes[0][2]->addLink(*nodes[1][1], 1.0);*/

    // 2 - 3 (2 direct conflicts vs. one)
    /*nodes[0][0]->addLink(*nodes[1][0], 1.0);
    nodes[1][0]->addLink(*nodes[0][0], 1.0);

    nodes[1][0]->addLink(*nodes[0][1], 1.0);
    nodes[0][1]->addLink(*nodes[1][0], 1.0);

    nodes[0][0]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[0][0], 1.0);*/

  /*  nodes[1][0]->addLink(*nodes[0][1], 1.0);
    nodes[0][1]->addLink(*nodes[1][0], 1.0);*/

  /*  nodes[2][0]->addLink(*nodes[0][2], 1.0);
    nodes[0][2]->addLink(*nodes[2][0], 1.0);*/

    /*nodes[0][0]->addLink(*nodes[1][0], 1.0);
    nodes[1][0]->addLink(*nodes[0][0], 1.0);

    nodes[0][0]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[0][0], 1.0);

    nodes[1][0]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[1][0], 1.0);

    nodes[0][1]->addLink(*nodes[2][0], 1.0);
    nodes[2][0]->addLink(*nodes[0][1], 1.0);*/

    //nodes[1][0]->addLink(*nodes[2][0], 1.0);
   // nodes[2][0]->addLink(*nodes[1][0], 1.0);

   // nodes[1][1]->addLink(*nodes[2][1], 1.0);
  //  nodes[2][1]->addLink(*nodes[1][1], 1.0);

    // 3 - 3 (2 direct conflicts vs. one)
    /*nodes[0][0]->addLink(*nodes[2][2], 1.0);
    nodes[2][2]->addLink(*nodes[0][0], 1.0);

    nodes[1][0]->addLink(*nodes[0][1], 1.0);
    nodes[0][1]->addLink(*nodes[1][0], 1.0);

    nodes[1][1]->addLink(*nodes[2][0], 1.0);
    nodes[2][0]->addLink(*nodes[1][1], 1.0);

    nodes[1][2]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[1][2], 1.0);*/

    // 3 - 3 (2 direct conflicts vs. one)
    /*nodes[0][0]->addLink(*nodes[2][2], 1.0);
    nodes[2][2]->addLink(*nodes[0][0], 1.0);

    nodes[1][0]->addLink(*nodes[0][1], 1.0);
    nodes[0][1]->addLink(*nodes[1][0], 1.0);

    nodes[1][1]->addLink(*nodes[2][0], 1.0);
    nodes[2][0]->addLink(*nodes[1][1], 1.0);

    nodes[1][2]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[1][2], 1.0);

    nodes[1][4]->addLink(*nodes[2][2], 0.5);
    nodes[2][2]->addLink(*nodes[1][4], 0.5);

    nodes[1][3]->addLink(*nodes[0][2], 1.0);
    nodes[0][2]->addLink(*nodes[1][3], 1.0);*/


    //
    /*nodes[0][0]->addLink(*nodes[1][0], 1.0);
    nodes[1][0]->addLink(*nodes[0][0], 1.0);

    nodes[0][0]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[0][0], 1.0);*/

/*    nodes[0][0]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[0][0], 1.0);

    nodes[0][1]->addLink(*nodes[2][2], 1.0);
    nodes[2][2]->addLink(*nodes[0][1], 1.0);

    nodes[1][0]->addLink(*nodes[0][2], 1.0);
    nodes[0][2]->addLink(*nodes[1][0], 1.0);

    nodes[2][0]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[2][0], 1.0);*/

    // 4 - 5
    /*  nodes[0][0]->addLink(*nodes[2][1], 1.0);
        nodes[2][1]->addLink(*nodes[0][0], 1.0);

        nodes[1][0]->addLink(*nodes[0][1], 1.0);
        nodes[0][1]->addLink(*nodes[1][0], 1.0);

        nodes[2][0]->addLink(*nodes[1][1], 1.0);
        nodes[1][1]->addLink(*nodes[2][0], 1.0);

        nodes[3][0]->addLink(*nodes[2][3], 1.0);
        nodes[2][3]->addLink(*nodes[3][0], 1.0);

        nodes[3][1]->addLink(*nodes[2][4], 1.0);
        nodes[2][4]->addLink(*nodes[3][1], 1.0);

        nodes[2][2]->addLink(*nodes[3][2], 1.0);
        nodes[3][2]->addLink(*nodes[2][2], 1.0);*/

    // 3 - 3
    /*nodes[0][0]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[0][0], 1.0);

    nodes[0][1]->addLink(*nodes[2][2], 1.0);
    nodes[2][2]->addLink(*nodes[0][1], 1.0);

    nodes[1][0]->addLink(*nodes[0][2], 1.0);
    nodes[0][2]->addLink(*nodes[1][0], 1.0);

    nodes[2][0]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[2][0], 1.0);*/

    // 3 - 3
    /*nodes[0][0]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[0][0], 1.0);

    nodes[0][0]->addLink(*nodes[2][2], 1.0);
    nodes[2][2]->addLink(*nodes[0][0], 1.0);

    nodes[1][0]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[1][0], 1.0);

    nodes[2][0]->addLink(*nodes[1][2], 1.0);
    nodes[1][2]->addLink(*nodes[2][0], 1.0);*/

    // 3 - 2
    /*nodes[0][0]->addLink(*nodes[1][0], 1.0);
    nodes[1][0]->addLink(*nodes[0][0], 1.0);

    nodes[0][0]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[0][0], 1.0);

    nodes[2][0]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[2][0], 1.0);

    nodes[1][0]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[1][0], 1.0);*/

    //nodes[2][0]->addLink(*nodes[1][2], 1.0);
    //nodes[1][2]->addLink(*nodes[2][0], 1.0);


    //nodes[1][1]->addLink(*nodes[0][2], 1.0);
    //nodes[0][2]->addLink(*nodes[1][1], 1.0);

    // CASE 1
    /*nodes[0][0]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[0][0], 1.0);

    nodes[1][0]->addLink(*nodes[2][1], 1.0);
    nodes[2][1]->addLink(*nodes[1][0], 1.0);

    nodes[2][0]->addLink(*nodes[0][1], 1.0);
    nodes[0][1]->addLink(*nodes[2][0], 1.0);*/

    // CASE 2
    /*nodes[0][0]->addLink(*nodes[2][0], 1.0);
    nodes[2][0]->addLink(*nodes[0][0], 1.0);

    nodes[0][0]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[0][0], 1.0);

    nodes[2][1]->addLink(*nodes[1][0], 1.0);
    nodes[1][0]->addLink(*nodes[2][1], 1.0);*/

    // Immediate direct conflict test
    /*nodes[0][0]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[0][0], 1.0);

    nodes[1][0]->addLink(*nodes[0][1], 1.0);
    nodes[0][1]->addLink(*nodes[1][0], 1.0);

    nodes[1][2]->addLink(*nodes[2][3], 1.0);
    nodes[2][3]->addLink(*nodes[1][2], 1.0);

    nodes[1][3]->addLink(*nodes[2][4], 1.0);
    nodes[2][4]->addLink(*nodes[1][3], 1.0);

    nodes[2][0]->addLink(*nodes[1][4], 1.0);
    nodes[1][4]->addLink(*nodes[2][0], 1.0);

    nodes[2][1]->addLink(*nodes[1][5], 1.0);
    nodes[1][5]->addLink(*nodes[2][1], 1.0);

    nodes[2][2]->addLink(*nodes[1][6], 1.0);
    nodes[1][6]->addLink(*nodes[2][2], 1.0);*/


    /*nodes[0][0]->addLink(*nodes[1][0], 1.0);
    nodes[1][0]->addLink(*nodes[0][0], 1.0);

    nodes[1][0]->addLink(*nodes[0][2], 1.0);
    nodes[0][2]->addLink(*nodes[1][0], 1.0);

    nodes[0][1]->addLink(*nodes[1][1], 1.0);
    nodes[1][1]->addLink(*nodes[0][1], 1.0);*/

    aligner = new GGAligner2(1000, nodes, CH_LS);
    aligner->align(segments);
}

void AlignTest::TearDown()
{

}

TEST_F(AlignTest, Conflict1) {

}