#include "GenePairs.h"

GenePairs::GenePairs(const string& blastTableFile) 
{
        // only read the gene pairs if they haven't been read before
	if (pairs.size() != 0) return;

        ifstream fin(blastTableFile.c_str());
        if (!fin) throw FileException ("Error opening BLAST table file: " +
						blastTableFile);

        // count the number of lines in the input file
        int c_lines = 0;
        for (string dummy; getline(fin, dummy); c_lines++);
	int initHashSize = (int)(c_lines*0.0005);

        // reset input file stream to begin position
        fin.clear();
        fin.seekg(0);

	string A, B;

	while(readLine(fin, A, B)) {

                if (A == B) continue;   //only insert when not the same gene

		if (pairs.find(A) != pairs.end()) {		//insert A-B
			pairs[A]->insert(B);
		} else {
			hash_set<string, stringhash>* tab =
				new hash_set<string,stringhash>(initHashSize);
			tab->insert(B);
			pairs[A] = tab;
		}

		if (pairs.find(B) != pairs.end()) {		//insert B-A
			pairs[B]->insert(A);
		} else {
			hash_set<string, stringhash>* tab =
				new hash_set<string,stringhash>(initHashSize);
			tab->insert(A);
			pairs[B] = tab;
		}
	}

        fin.close();
}

GenePairs::~GenePairs()
{
 	hash_map<string, hash_set<string,stringhash>*,stringhash>::const_iterator it = pairs.begin();
	for ( ; it != pairs.end(); it++)
		delete it->second;
}

bool GenePairs::readLine(ifstream& fin, string &geneA, string &geneB)
{
        getline(fin, geneA, '\t');
        if (!fin.good()) return false;
        getline(fin, geneB);

        // perhaps we are reading the final empty line(s) of the file
        if (geneA.empty() || geneB.empty()) return false;

        char charA = geneA[geneA.length()-1];  // final character in geneA
        if ((charA == '+') || (charA == '-')) geneA.erase(geneA.length()-1);

        char charB = geneB[geneB.length()-1];  // final character in geneB
        if ((charB == '+') || (charB == '-')) geneB.erase(geneB.length()-1);

        return true;
}

hash_set<string,stringhash>* GenePairs::getPairsOf(const string& geneName)
{
	if (pairs.find(geneName) != pairs.end()) {
		return pairs[geneName];
	}
	else {
		return NULL;
	}
}
