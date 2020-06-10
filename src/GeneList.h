#ifndef __GENELIST_H
#define __GENELIST_H

class ListElement;
class Gene;

#include "headers.h"


class GeneList {

public:
    ////////////////
    //CONSTRUCTORS//
    ////////////////

    /**
    * Constructs the genelist by parsing the listfile specified by the
    * filename and creating corresponding listelement objects
    * and gene objects
    */
    GeneList(const string& listName, const string& genomeName,
             const string& fileName);


    /**
    * Constructs a genelist from a segment of another genelist
    */
    GeneList(const GeneList& genelist, int begin, int end);

    /**
     * Default constructor
     */
    GeneList(int size);


    /**
     * Constructor to self assemble a genelist from listel pointer -> note that these will be deep copyed to avoid
     * memory leaks!!
     */
    GeneList(const string& listName, const string& genomeName, const vector<ListElement*>& segmentFromFile);

    /**
    * Destructor
    */
    ~GeneList();

    //////////////////
    //PUBLIC METHODS//
    //////////////////

    /*
    *performs remapping of tandem repeated genes in a gene list object.
    *Genes in a tandem are remapped onto the first gene which is marked as tandem_representative
    */
    void remapTandems(int gapSize, bool useFamily);

    /*
    *returns the total number of listelements that have not been remapped
    */
    int getElementsLength() const {
        return elements.size();
    }

    /*
    *returns the total number of listelements that have been remapped
    */
    int getRemappedElementsLength() const {
        return remapped_elements.size();
    }

    /*
    *returns the total number of listelements that have been remapped
    */
    unsigned int getSize() const {
        return remapped_elements.size();
    }

    ListElement& getLe(int id) const {
        return *remapped_elements.at(id);
    }

    vector<ListElement *>::const_iterator getLEBegin() const {
        return remapped_elements.begin();
    }

    vector<ListElement *>::const_iterator getLEEnd() const {
        return remapped_elements.end();
    }

    const string getGeneName(int i);

    /*
    *returns the list name
    */
    const string& getListName() const {
        return listname;
    }

    /*
    *returns the genome name
    */
    const string& getGenomeName() const {
        return genomename;
    }

    /*
    *returns a vector of all remapped list elements
    */
    const vector<ListElement*>& getRemappedElements() const {
        return remapped_elements;
    }

    /*
    *returns a vector of not remapped list elements
    */
    vector<ListElement*>& getElements() {
        return elements;
    }

    /*
    *returns the number of masked elements
    */
    int getNumberOfMaskedElements() const;

    /*
    *introduces between two positions in a segment a number of gaps
    */
    void introduceGaps(int begin, int end, int number_of_gaps);

    /*
    *introduces between two positions in a segment 1 gap
    */
    void introduceGap(int position);

    void pasteGap();

    /*
    *inverts the order of a section of this segment between the 2 integer arguments
    */
    void invertSection(int begin, int end);

    /*
    *masks a section of the genelist, from a begin to an end position
    */
    void mask(int begin, int end);

    /*
    *returns true if the genelist was contructed from a segment of another genelist
    */
    bool isSegment() const {
        return is_segment;
    }

    /*
    *for debugging purposes
    */
    void write();
    /*
    *for debugging purposes
    */
    void write2();

    /*
     * Set the ID of the genelist
     */
    void setID(int id_) {
        id = id_;
    };

    /*
     * Get the ID of the genelist
     */
    int getID() const {
        return id;
    }

    /**
     * Remove all elements with an index greater than the "last" index
     * @param last Final index to retain
     */
    void cutAfter(int last) {
        remapped_elements.erase(remapped_elements.begin() + last + 1,
                                remapped_elements.end());
    }

    /**
     * Remove all elements with an index smaller than the "first" index
     * @param last Final index to retain
     */
    void cutBefore(int first) {
        remapped_elements.erase(remapped_elements.begin(),
                                remapped_elements.begin() + first);
    }

    /**
     * Remove the gaps for the remapped elements
     */
    void removeGaps();

private:
    ///////////////////
    //PRIVATE METHODS//
    ///////////////////

    void remapTandemsPairs(int gapSize);
    void remapTandemsFamily(int gapSize);

    /**
    * Copy constructor
    */
    GeneList(const GeneList& genelist) {};

    /**
    * Assignment operator
    */
    void operator=(const GeneList& genelist) {};

    //////////////
    //ATTRIBUTES//
    //////////////

    bool is_segment;

    string listname;
    string genomename;

    vector<ListElement*> elements;
    vector<ListElement*> remapped_elements;

    int id;
};

#endif
