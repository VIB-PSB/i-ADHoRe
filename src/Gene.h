#ifndef __GENE_H
#define __GENE_H

#include "headers.h"



class Gene {

public:
    ///////////////////////////////
    //CONSTRUCTORS AND DESTRUCTOR//
    ///////////////////////////////

    /*
    *constructs the gene and fills in the attributes
    */
    Gene(const string& ID, const string& genomeName, const int coordinate,
         const bool orientation);

    /*
    *default constructor (gap)
    */
    Gene();

    //operator the combination genomename && ID should be unique !!
    bool operator== (Gene other) const {
        return (genomename == other.genomename && ID == other.ID);
    }


    //////////////////
    //PUBLIC METHODS//
    //////////////////

    /*
    *returns the genome name
    */
    const string& getGenomeName() const {
        return genomename;
    }

    /*
    *returns the ID of this gene
    */
    const string& getID() const {
        return ID;
    }

    /*
    *sets the pairs corresponding this gene
    */
    void setPairs(const hash_set<string,stringhash>& table) {
        pairs = &table;
    }

    /*
    *returns true if this gene has pairs
    */
    bool hasPairs() const {
        return (pairs != NULL || has_gf);
    }

    /*
    *returns true if the specified gene is a pair with this gene
    */
    bool isPairWith(const Gene& gene) const;

    /*
    *returns true if the specified gene is an indirect pair with this gene
    */
    bool isIndirectPairWith(const Gene& gene) const;

    /*
    *remaps the gene onto the other gene
    */
    void remapTo(Gene& gene);

    /*
    *returns true if this gene was remapped
    */
    bool isRemapped() const {
        return remapped == true;
    }

    /*
    *returns the original coordinate
    */
    int getCoordinate() const {
        return coordinate;
    }

    /*
    *returns the remapped coordinate
    */
    int getRemappedCoordinate() const {
        return remapped_coordinate;
    }

    /*
    *substracts the remapped coordinate with the specified argument
    */
    void setRemappedCoordinate(const int substract) {
        remapped_coordinate = coordinate - substract;
    }

    /*
    *returns true if this gene is the tandem representative if a tandem was formed
    */
    bool isTandemRepresentative() const {
        return is_tandem_representative;
    }

    /*
    *return true if the gene forms a tandem
    */
    bool isTandem() const {
        return is_tandem;
    }

    /*
    *returns the gene that is tandem representative if this gene was remapped
    */
    const Gene& tandemRepresentative() const {
        return *tandem_representative;
    }

    /*
    *returns the orientation class
    */
    bool getOrientation() const {
        return orientation;
    }


    /*
    *returns the gene family ID of this gene
    */
    const string& getFamily() const {
        return gf_id;
    }

    void setFamily(string fam)
    {
        gf_id = fam;
        has_gf = true;
    }


private:
    //////////////
    //ATTRIBUTES//
    //////////////

    string ID;
    string genomename;
    int coordinate;
    bool orientation;
    int remapped_coordinate;
    bool is_tandem;
    bool is_tandem_representative;
    bool remapped;
    string gf_id;
    bool has_gf;

    const Gene* tandem_representative; // pointer to the representative
    const hash_set<string, stringhash>* pairs; // pointer to the pairs
};

#endif
