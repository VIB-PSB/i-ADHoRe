#ifndef __GENEFAMILY_H
#define __GENEFAMILY_H

#include "debug/FileException.h"
#include "headers.h"

class GeneFamily
{
public:
	////////////////
	//CONSTRUCTORS//
	////////////////

	/*
	 *constructor that reads the blasttable file and builds up the hashtable
	 */
	GeneFamily(const string& blastTableFile);

	//////////////////
	//PUBLIC METHODS//
	//////////////////

	/*
	 * returns a pointer to the hashtable containing the genes
	 * corresponding with the specified genome
	 */
	const string& getFamilyOf(const string& geneName);

	bool hasFamily(const string& geneName) {
		return (pairs.find(geneName) != pairs.end());
	}

private:
	///////////////////
	//PRIVATE METHODS//
	///////////////////
		
	/*
	 * reads a single line from the blast table file and puts the 2
	 * genes into the 2 strings
	 */
	bool readLine(ifstream& fin, string &geneA, string &geneB);

	//////////////
	//ATTRIBUTES//
	//////////////

	//the hash_map containing all pairs read from the blast table file
	hash_map<string, string, stringhash> pairs;

};


#endif
