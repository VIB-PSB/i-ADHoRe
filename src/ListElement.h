#ifndef __LISTELEMENT_H
#define __LISTELEMENT_H

class Gene;
#include "Gene.h"
#include "headers.h"

class ListElement {

public:
    ///////////////////////////////
    //CONSTRUCTORS AND DESTRUCTOR//
    ///////////////////////////////

    /*
    *default constructor creates the object as a gap with no meaningfull attributes
    */
    ListElement();

    /*
    *constructs the list element
    */
    ListElement(Gene& gene, bool orientation, bool masked);

    //////////////////
    //PUBLIC METHODS//
    //////////////////

    /*
    *constructs the list element
    */
    const Gene& getGene() const;

    /*
    *constructs the list element
    */
    Gene& getGene();

    /*
    *fills in the queue with integers corresponding to pairs with the listelement
    */
    void matchingPositions(const vector<ListElement*>& list,
                           queue<int>& q) const;

    /*
    *returns true if the list element is masked
    */
    bool isMasked() const;

    /*
    *sets if the list element should be masked
    */
    void setMasked(bool m) {
        masked = m;
    }

    /*
    *returns the orientation class of the list element
    */
    bool getOrientation() const;

    /*
    *inverts the orientation
    */
    void invertOrientation() {
        orientation = (!orientation);
    }

    /*
    *returns true if the listelement represents a gap
    */
    bool isGap() const {
        return gap;
    }

    void setNumID(int numID_) {
        numID = numID_;
    }

    int getNumID() const {
        return numID;
    }

    bool hasHomolog;    // element has homologs
    bool hasAP;         // element has anchor point
    bool hasAlHomolog;  // element has aligned homologs
    bool hasAlAP;       // element has aligned anchor point

private:
    //////////////
    //ATTRIBUTES//
    //////////////

    Gene gene;
    bool orientation;
    bool masked;
    bool gap;

    int numID;
};

#endif
