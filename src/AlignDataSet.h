#ifndef __ALIGNDATASET_H
#define __ALIGNDATASET_H

class Settings;
class GeneList;
class Multiplicon;
class GenePairs;
class ListElement;
class Gene;
class Profile;
class GHMProfile;
class ListFile;

#include "debug/FileException.h"

#include "headers.h"



class AlignDataSet {
	public:
		///////////////////////////////
		//CONSTRUCTORS AND DESTRUCTOR//
		///////////////////////////////

		/*
		 *constructs a dataset object and gets the properties out of the settings object.
		 *creates genelist objects for every listfile
		 */
		AlignDataSet(Settings& sett);

		/*
		 *destructor
		 */
		~AlignDataSet();


		//////////////////
		//PUBLIC METHODS//
		//////////////////

		/*
		 *reads alle genepairs from the blast table file and removes the not
		 *necessary pairs
		 */
	    void mapGenes();
		void getGenePairs();
		void getGeneFamilies();

		/*
		 *remaps pairs and indirect pairs for every genelist on one gene
		 */
		void remapTandems();

		/*
		 *aligns loaded lists
		 */
		void align();

		void output();

	private:
		///////////////////
		//PRIVATE METHODS//
		///////////////////

		/*
		 *adds the vector of multiplicons to the existing vector of multiplicons
		 */
		void addMultiplicons(const vector<Multiplicon*>& multiplicons);

		/*
		 *sorts the vector of multiplicons by lowest discrete pseudo distance first, and then by
		 *the lowest multiplicon size
		 */
		void sortByMultipliconSize(vector<Multiplicon*>& multiplicons) const;

		/*
		 *Takes a profile object as input and searches all genelists for
         *matching segments. These segments are added and returned as new
         *multiplicons into the target vector
		 */
		void profileSearch(Profile& profile, vector<Multiplicon*>& target);

		//GeneList getGeneLists();
		//void storeMultiplicon();
		//Multiplicon* nextMultiplicon();
		//int getMultipliconsLeft();

		//Multiplicon* getEvaluatedMultiplicon();
		//AnchorPoint* getAnchorPoints();


		//////////////
		//ATTRIBUTES//
		//////////////

		//settings object containing all information necessary for the algorithm
		Settings& settings;

		void readDataSet();
		//vector containing all the genelists
		vector<GeneList*> genelists;
		vector<GeneList*> AlignedLists;
		//vector containing all calculated multiplicons by the level 1 algorithm
		vector<Multiplicon*> multiplicons;
		//queue containing all multiplicons that need to be evaluated into a profile
        deque<Multiplicon*> multiplicons_to_evaluate;
        //vector containing all multiplicons that have been evaluated
		vector<Multiplicon*> evaluated_multiplicons;

};



#endif
