#include <gtest/gtest.h>
#include "../src/Settings.h"
#include "../src/DataSet.h"
#include "../src/GeneList.h"
#include "../src/Multiplicon.h"
#include "../src/BaseCluster.h"
#include "../src/AnchorPoint.h"

using namespace std;

class PackingTest : public ::testing::Test
{
protected:
	virtual void SetUp();

	virtual void TearDown();

	Settings *settings;
	DataSet *dataset;
	vector<Multiplicon *> mplicons, mplicons2;
	char *buffer;
};

void PackingTest::SetUp()
{
	settings = new Settings("testset.ini");
	dataset = new DataSet(*settings);
	dataset->mapGenes();
	dataset->remapTandems();

	GeneList &gl = *dataset->genelists[0];
	for (int i = 0; i < 5; i++) {
		Multiplicon *mpl = new Multiplicon(gl.getID(), gl.getID(), 2);
		mplicons.push_back(mpl);
		for (int j = 0; j < 5; j++) {
			BaseCluster *bc = new BaseCluster(true);
			bc->setMultiplicon(*mpl);
			mpl->addBaseCluster(*bc);
			for (int k = 0; k < 5; k++)
				bc->addAnchorPoint(j, k);
		}
	}

	int buffSize = Multiplicon::getPackSize(mplicons);
	buffer = new char[buffSize];
	//Multiplicon::packMultiplicons(mplicons, buffer);
	//Multiplicon::unpackMultiplicons(buffer, mplicons2);
}

void PackingTest::TearDown()
{
	for (int i = 0; i < mplicons.size(); i++) {
		delete mplicons[i];
		delete mplicons2[i];
	}

	delete [] buffer;
	delete dataset;
	delete settings;
}

TEST_F(PackingTest, MultipliconPackingTest) {
	EXPECT_EQ(mplicons.size(), mplicons2.size());
	for (int i = 0; i < mplicons.size(); i++) {
		EXPECT_EQ(*mplicons[i] == *mplicons2[i], true);
	}
}
