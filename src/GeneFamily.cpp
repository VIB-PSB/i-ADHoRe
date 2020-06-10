#include "GeneFamily.h"

GeneFamily::GeneFamily(const string& blastTableFile) 
{
	//only read the gene pairs if they haven't been read before
	if (pairs.size() != 0) return;

	ifstream fin(blastTableFile.c_str());
	if (!fin) throw FileException ("Error opening BLAST table file: " + 
						blastTableFile);

	string A, F;

	while(readLine(fin, A, F))
		pairs[A] = F;

	fin.close();
}

bool GeneFamily::readLine(ifstream& fin, string &geneA, string &geneB)
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

const string& GeneFamily::getFamilyOf(const string& geneName)
{
	if (pairs.find(geneName) != pairs.end()) {
		return pairs[geneName];
	}
	else {
		cerr << geneName << " not found in blast table. EXITING...." 
		     << endl << endl;
		exit(EXIT_FAILURE);
	}
}
