#ifndef __GHMPROFILE_H
#define __GHMPROFILE_H

#include "GHM.h"

class Profile;
class ListElement;
class Gene;


class GHMProfile : public GHM {

public:
    ///////////////////////////////
    //CONSTRUCTORS AND DESTRUCTOR//
    ///////////////////////////////

    /**
    * constructs a ghm from a profile and a genelist
    *
    * @param xObject The Profile object corresponding to the X-axis
    * @param yObject The GeneList object corresponding to the Y-axis
    */
    GHMProfile(const Profile& xObject, const GeneList& yObject);

    //////////////////
    //PUBLIC METHODS//
    //////////////////

    /**
    * Creates every matching point in the matrix,
    * the int argument indicates how many segments there are in the x_object,
    * for single genelists this is simply 1.
    */
    void buildMatrix();

    /**
     * Fills in the vector with every multiplicon
     * but first some additional data needs to be filled in
     * @param mps The vector that needs to be filled in with multiplicons
     * @param useFamily True if we're using gene families
     * @param minHomologs The minimum number of homologous genes
     */
    void setMultiplicons(bool useFamily,
                          int minHomologs);

    void getMultiplicons(vector<Multiplicon*>& mps) const;
#ifdef HAVE_PNG
    void visualizeGHM(const std::string& output_path) const;
#endif
private:

    //////////////
    //ATTRIBUTES//
    //////////////

    const Profile& x_object;
};

#endif
