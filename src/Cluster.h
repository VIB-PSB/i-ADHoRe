#ifndef __CLUSTER_H
#define __CLUSTER_H

class Cluster {

public:
    //////////////////
    //PUBLIC METHODS//
    //////////////////

    /**
     * Default constructor
     */
    Cluster() : x_objectID(-1), y_objectID(-1), begin_x(-1), end_x(-1),
            begin_y(-1), end_y(-1) {};

    /**
     * Constructor with x- and y-genelist IDs
     */
    Cluster(int x_objectID_, int y_objectID_) : x_objectID(x_objectID_),
            y_objectID(y_objectID_), begin_x(-1), end_x(-1), begin_y(-1),
            end_y(-1) {};

    /**
     * Fills in the begin and end values of the cluster
     */
    void setBounds() 
    {
        begin_x = getLowestX();
        begin_y = getLowestY();
        end_x = getHighestX();
        end_y = getHighestY();
    }

    /**
     * Returns the beginning x-coordinate
     */
    int getBeginX() const {
        return begin_x;
    }

    /**
     * Returns the ending x-coordinate
     */
    int getEndX() const {
        return end_x;
    }

    /**
     * Returns the beginning y-coordinate
     */
    int getBeginY() const {
        return begin_y;
    }

    /**
     * Returns the ending y-coordinate
     */
    int getEndY() const {
        return end_y;
    }

    /**
     * Returns for all anchorpoints the lowest x coordinate
     */
    virtual unsigned int getLowestX() const = 0;

    /**
     * Returns for all anchorpoints the highest x coordinate
     */
    virtual unsigned int getHighestX() const = 0;

    /**
     * Returns for all anchorpoints the lowest y coordinate
     */
    virtual unsigned int getLowestY() const = 0;

    /**
     * Returns for all anchorpoints the highest y coordinate
     */
    virtual unsigned int getHighestY() const = 0;

    /**
     * Return the ID of the x-genelist
     */
    int getXObjectID() const {
        return x_objectID;
    }

    /**
     * Return the ID of the y-genelist
     */
    int getYObjectID() const {
        return y_objectID;
    }

    /**
     * Calculates the discrete pseudo distance between 2 points
     */
    double dpd(double ref_x, double ref_y, double x, double y) const;

    /**
     * Calculates the discrete pseudo distance between to points
     */
    static int dpd(int refX, int refY, int x, int y);

    /**
     * Calculates the discrete pseudo distance for the begin and
     * end coordinates
     */
    double dpd() const;

    /**
     * Calculates the king step pseudo distance (minimum number of moves a king has to make on a chess board
     * to move from begin to end coordinate) between two coordinates
     */
    static int kspd(const int& x1, const int& y1, const int &x2, const int& y2);

    void setXObjectID(int x){x_objectID=x;};
    void setYObjectID(int y){y_objectID=y;};

protected:
    ///////////////////////////////
    //CONSTRUCTORS AND DESTRUCTOR//
    ///////////////////////////////

    /**
    * Destructor
    */
    virtual ~Cluster() {}

    //////////////
    //ATTRIBUTES//
    //////////////

    int x_objectID;
    int y_objectID;

    int begin_x, end_x;
    int begin_y, end_y;

    /////////////////////
    //PROTECTED METHODS//
    /////////////////////
};

#endif
