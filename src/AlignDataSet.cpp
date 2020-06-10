#include "AlignDataSet.h"

#include "AnchorPoint.h"
#include "BaseCluster.h"
#include "GeneList.h"
#include "GenePairs.h"
#include "GeneFamily.h"
#include "Multiplicon.h"
#include "Profile.h"
#include "ListElement.h"
#include "Gene.h"
#include "GHMProfile.h"
#include "Settings.h"
#include "alignment/NWAligner.h"
#include "alignment/GGAligner.h"
#include "alignment/Aligner.h"
#include "alignment/AlignmentException.h"



AlignDataSet::AlignDataSet(Settings& sett) : settings(sett)  {

	cerr << "Creating dataset...";

	list<ListFile> listfiles = settings.getListFiles();

	list<ListFile>::iterator it;
	for (it = listfiles.begin(); it != listfiles.end(); it++) {
		genelists.push_back(new GeneList(it->getListName(), it->getGenomeName(), it->getFileName()));
	}

	cerr << "\t\tDone!" << endl;
}

AlignDataSet::~AlignDataSet() {
	vector<Multiplicon*>::iterator itm = evaluated_multiplicons.begin();
	while (itm != evaluated_multiplicons.end()) {
		//delete profile
		if ((*itm)->getProfile() != NULL) {
			delete (*itm)->getProfile();
		}

		//delete multiplicon
		delete (*itm);
		itm = evaluated_multiplicons.erase(itm);
	}
}

void AlignDataSet::mapGenes() {

	if (settings.useFamily() == true)
	{
		cerr << "Mapping gene families...";
		getGeneFamilies();

	} else {
		cerr << "Mapping gene pairs...";
		getGenePairs();
	}

	cerr << "\t\tDone!" << endl;
}

void AlignDataSet::getGenePairs() {

	try {
		GenePairs genepairs (settings.getBlastTable());

		unsigned int count = 0;
		for (unsigned int i = 0; i < genelists.size(); i++) {
			count += genelists[i]->getElementsLength();
		}

		hash_map<string,Gene,stringhash> genes ((unsigned int)(count*1.25));
		for (unsigned int i = 0; i < genelists.size(); i++) {
			vector<ListElement*>& list = genelists[i]->getElements();
			for (unsigned int j = 0; j < list.size(); j++) {
				//genes[string(list.at(j)->getGene().getID())] = list[j]->getGene();
				genes.insert(pair<string, Gene> (list.at(j)->getGene().getID(), list[j]->getGene()));
			}
		}
		for (unsigned int i = 0; i < genelists.size(); i++) {
			vector<ListElement*>& list = genelists[i]->getElements();
			for (unsigned int j = 0; j < list.size(); j++) {
				hash_set<string,stringhash>* pairs = genepairs.getPairsOf(string(list.at(j)->getGene().getID()));

				if (pairs != 0) {
					hash_set<string, stringhash>::iterator it;

					vector<hash_set<string, stringhash>::iterator> iterators;

					for (it = pairs->begin(); it != pairs->end(); it++) {
						//hash_map<string,Gene*,stringhash>::iterator iter = genes.find(*it);
						if (genes.find(*it) == genes.end()) {
							//this gives problems on the hagrid server, older gcc compiler?
							//pairs->erase(it);
							iterators.push_back(it);
						}
					}
					for (unsigned int k = 0; k < iterators.size(); k ++) {
						pairs->erase(iterators[k]);
					}
				}

				list.at(j)->getGene().setPairs(*pairs);
			}
		}

	}
	catch (const FileException& f) {
		throw f;
	}
}

void AlignDataSet::getGeneFamilies() {

	try {

		GeneFamily genefamily(settings.getBlastTable());

		unsigned int count = 0;
		for (unsigned int i = 0; i < genelists.size(); i++) {
			count += genelists[i]->getElementsLength();
		}

		for (unsigned int i = 0; i < genelists.size(); i++) {
			vector<ListElement*>& list = genelists[i]->getElements();
			for (unsigned int j = 0; j < list.size(); j++) {

				string ID(list.at(j)->getGene().getID());

				string fam = genefamily.getFamilyOf(ID);
				list.at(j)->getGene().setFamily(genefamily.getFamilyOf(ID));

				//cerr << "\33[1;32m" << ID << "\33[0m\t" << fam << "\n";

			}

		}

	}
	catch (const FileException& f) {
		throw f;
	}
}

void AlignDataSet::remapTandems() {
	cerr << "Remapping tandem duplicates";
	for (unsigned int i = 0; i < genelists.size(); i++) {
		//genelists[i].remapTandems(settings.getGapSize());
		genelists[i]->remapTandems(settings.getTandemGap(), settings.useFamily());
		cerr << '.';
	}
	cerr << "\tDone!" << endl;
}

void AlignDataSet::align() {
	cerr << "Aligning gene lists";

	AlignedLists.push_back(genelists[0]);

	for(unsigned int i = 1; i < genelists.size(); i++)
	{
		AlignedLists.push_back(genelists[i]);

		/*try {
		Aligner* aligner;


		if (settings.getAlignmentMethod() == settings.NeedlemanWunsch) {
       cerr << "NWAligner not supported yet...Crash and Burn..." <<
       endl;
		}
		else if (settings.getAlignmentMethod() == settings.GreedyGraphbased) {
       aligner = new GGAligner(settings.getGapSize());
		}

		aligner->align(AlignedLists);

		delete aligner;
		}
		catch (const AlignmentException& e) {
       throw e;
		}*/
		cerr << ".";
		}
		try {
			Aligner* aligner = NULL;


			if (settings.getAlignmentMethod() == NeedlemanWunsch) {
				cerr << "NWAligner not supported yet...Crash and Burn..." <<
				endl;
			}
			else if (settings.getAlignmentMethod() == GreedyGraphbased) {
				aligner = new GGAligner(settings.getGapSize());
			}

			aligner->align(AlignedLists);

			delete aligner;
		}
		catch (const AlignmentException& e) {
			throw e;
		}

		cerr << "\tDone!" << endl;
	}

	void AlignDataSet::addMultiplicons(const vector<Multiplicon*>& mplicons) {
		for (unsigned int i = 0; i < mplicons.size(); i++) {
			multiplicons.push_back(mplicons[i]);
		}
	}

	void AlignDataSet::sortByMultipliconSize(vector<Multiplicon*>& mplicons) const {
		//first sort all on dpd
		for (unsigned int i = 1; i < mplicons.size(); i++) {
			Multiplicon* tmp = mplicons[i];
			unsigned int j = i;
			while (j > 0 && tmp->dpd() > mplicons[j-1]->dpd()) {
				mplicons[j] = mplicons[j-1];
				j--;
			}
			mplicons[j] = tmp;
		}

		//then sort on anchorpoints per multiplicon
		for (unsigned int i = 1; i < mplicons.size(); i++) {
			Multiplicon* tmp = mplicons[i];
			unsigned int j = i;
			while (j > 0 && tmp->getCountAnchorPoints() < mplicons[j-1]->getCountAnchorPoints()) {
				mplicons[j] = mplicons[j-1];
				j--;
			}
			mplicons[j] = tmp;
		}
	}

	void AlignDataSet::profileSearch(Profile& profile, vector<Multiplicon*>& target) {
		for (unsigned int i = 0; i < genelists.size(); i++) {
			if (genelists[i]->getSize() - genelists[i]->getNumberOfMaskedElements() >= (unsigned int)settings.getAnchorPoints()) {
				GHMProfile ghm (profile, *genelists[i]);
				ghm.buildMatrix();
				ghm.run(settings);
				//ghm.plotGHM(settings.getGapSize(), "datasets/bugsy/output_koen/profiles");
				ghm.getMultiplicons(target);
			}
		}
	}

	void AlignDataSet::output() {

		GeneFamily families("monocots.fam");

		cerr << "Generating output files...";


		//cout << "<table>" << endl;
		for (unsigned int j = 0; j < AlignedLists.size(); j++)
		{
			//cout << "<tr>" << endl;
			for (unsigned int i = 0; i < AlignedLists[j]->getSize(); i++)
			{
				string part = "000000";
				string famname = "";
				string name = AlignedLists[j]->getGeneName(i);
				int t = 0;

				if (families.hasFamily(AlignedLists[j]->getGeneName(i)))
				{
					famname = families.getFamilyOf(AlignedLists[j]->getGeneName(i));
					part = famname.substr(3);
					t = atoi(part.c_str())*3;
				}

				//cout << "<td bgcolor ='#" << t <<"'>";
				//cout << "<td>";
				cout << j << ":" << i << ":" << famname << ":" << name << endl;
				//cout << "</td>";
			}
			//cout << endl << "</tr>" << endl;
		}
		//cout << "</table>" << endl;



		cerr << "Done!" << endl;
	}

