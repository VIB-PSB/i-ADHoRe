#ifndef __GENEPAIRS_H
#define __GENEPAIRS_H

#include "debug/FileException.h"
#include "headers.h"

class GenePairs
{
public:
	////////////////
	//CONSTRUCTORS//
	////////////////

	/**
	 * Constructor that reads the blasttable file and builds up the hashtable
	 */
	GenePairs(const string& blastTableFile);

	/**
	 * Destructor
	 */
	~GenePairs();

	//////////////////
	//PUBLIC METHODS//
	//////////////////

	/**
	 * returns a pointer to the hashtable containing the genes 
	 * corresponding with the specified genome
	 */
	hash_set<string,stringhash>* getPairsOf(const string& geneName);

private:
	///////////////////
	//PRIVATE METHODS//
	///////////////////

	/*
	 * reads a single line from the blast table file and puts the 2
	 * genes into the 2 strings
	 */
	static bool readLine(ifstream& fin, string &geneA, string &geneB);

	//////////////
	//ATTRIBUTES//
	//////////////

	// the hash_map with all pairs read from the blast table file
	hash_map<string, hash_set<string,stringhash>*,stringhash> pairs;
};

#endif
