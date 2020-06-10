#include <gtest/gtest.h>
#include "../src/Settings.h"
#include "../src/DataSet.h"
#include "../src/GeneList.h"
#include "../src/ListElement.h"

using namespace std;

class GapsTest : public ::testing::Test
{
protected:
	virtual void SetUp();

	virtual void TearDown();

	Settings *settings;
	DataSet *dataset;
	vector<Multiplicon *> mplicons, mplicons2;
	char *buffer;
};

void GapsTest::SetUp()
{
	settings = new Settings("testset.ini");
	dataset = new DataSet(*settings);
	dataset->mapGenes();
	dataset->remapTandems();

	dataset->genelists[0]->introduceGaps(0, 10, 7);
	for (int i = 0; i < 18; i++)
		if (dataset->genelists[0]->getRemappedElements()[i]->isGap())
			cout << "_";
		else
			cout << "X";
	cout << endl;

	dataset->genelists[0]->removeGaps();
	dataset->genelists[0]->introduceGaps(0, 1, 5);
	for (int i = 0; i < 10; i++)
		if (dataset->genelists[0]->getRemappedElements()[i]->isGap())
			cout << "_";
		else
			cout << "X";
	cout << endl;

	dataset->genelists[0]->removeGaps();
	dataset->genelists[0]->introduceGaps(0, 0, 5);
	for (int i = 0; i < 10; i++)
		if (dataset->genelists[0]->getRemappedElements()[i]->isGap())
			cout << "_";
		else
			cout << "X";
	cout << endl;
}

void GapsTest::TearDown()
{
	delete dataset;
	delete settings;
}

TEST_F(GapsTest, IntroduceGapsTest) {

}
