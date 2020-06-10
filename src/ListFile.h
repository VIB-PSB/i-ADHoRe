#ifndef __LISTFILE_H
#define __LISTFILE_H

#include "headers.h"

class ListFile
{

public:
	///////////////
	//CONSTRUCTOR//
	///////////////

	/**
	* Constructs a listfile object and fills in the strings
	*/
	ListFile(const string& genomeName, const string& listName,
		 const string& fileName) : genomename(genomeName),
		 listname(listName), filename(fileName) { }


	//////////////////
	//PUBLIC METHODS//
	//////////////////

	/**
	* Returns the genomename of the listfile
	*/
	const string& getGenomeName() const { return genomename; }

	/**
	* Returns the listname
	*/
	const string& getListName() const { return listname; }

	/**
	* Returns the filename of the listfile
	*/
	const string& getFileName() const { return filename; }

private:
	//////////////
	//ATTRIBUTES//
	//////////////

	const string genomename;
	const string listname;
	const string filename;
};

#endif
