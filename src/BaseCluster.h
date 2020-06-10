#ifndef __BASECLUSTER_H
#define __BASECLUSTER_H

#include "headers.h"
#include "Cluster.h"
#include "AnchorPoint.h"

class Multiplicon;

class BaseCluster: public Cluster {

public:
    ///////////////////////////////
    //CONSTRUCTORS AND DESTRUCTOR//
    ///////////////////////////////

    /**
    * Constructs a basecluster object with the specified orientation class
    */
    BaseCluster(const bool orient);

    /**
    * Constructs a basecluster object from a packed buffer
    */
    BaseCluster(const char *buffer);

    /**
    * Destructor
    */
    ~BaseCluster() {}

    //////////////////
    //PUBLIC METHODS//
    //////////////////

    /**
    * Creates an anchorpoint from the specified x and y values and is placed
    * into the vector of anchorpoints
    */
    void addAnchorPoint(int x, int y);

    void filterDuplicates();

    void addBackBone(int X, int Y);

    /**
    * Adds an already existing anchorpoint to the vector of anchorpoints
    */
    void addAnchorPoint(AnchorPoint& point);

    void addBackBone(AnchorPoint& point);

    /**
    * Returns an iterator pointing to the first anchorpoint
    */
    multiset<AnchorPoint>::const_iterator getAPBegin() {
        return anchorpoints.begin();
    }

    /**
    * Returns an iterator pointing past the final anchorpoint
    */
    multiset<AnchorPoint>::const_iterator getAPEnd() {
        return anchorpoints.end();
    }

    /**
    * Returns the total number of anchorpoints
    */
    unsigned int getCountAnchorPoints() const;

    /**
    * Returns for all anchorpoints the lowest x coordinate
    */
    unsigned int getLowestX() const;

    /**
    * Returns for all anchorpoints the highest x coordinate
    */
    unsigned int getHighestX() const;

    /**
    * Returns for all anchorpoints the lowest y coordinate
    */
    unsigned int getLowestY() const;

    /**
    * Returns for all anchorpoints the highest y coordinate
    */
    unsigned int getHighestY() const;

    /**
    * Updates all parameters describing the confidence interval of the basecluster
    */
    void updateStatistics();

    /**
    * Calculates the squared Pearson value for all x and y coordinates of all
    * anchorpoints, including the specified values.
    */
    double r_squared(int X = 0, int Y = 0) const;

    /**
    * Calculates the squared Pearson value for all x and y coordinates of all
    * anchorpoints, including all anchorpoints from the specified cluster
    */
    double r_squared(BaseCluster& cluster) const;

    /**
    * Calculates the average dpd between each following anchorpoint
    */
    double averageDPD() const;

    /**
    * Calculates the distance from the basecluster to the specified x and y coordinate
    */
    double distanceToPoint(int x, int y) const;

    /**
    * Returns a boolean saying the specified x and y coordinate is located inside
    * the confidence interval of the basecluster
    */
    bool inInterval(int x, int y);

    /**
    * Returns for a given x coordinate the upper and lower y coordinates of the confidence interval
    */
    void intervalBounds(int x, double& up, double& down);

    /**
    * Returns 1 if either the x or y coordinates overlap the specified cluster
    * Returns 2 if both the x and y coordinates overlap the specified cluster
    * Returns 0 if nothing overlaps
    */
    int overlappingCoordinates(BaseCluster& cluster) const;

    /**
    * Returns the orientation class of the basecluster
    */
    bool getOrientation() const {
        return orientation;
    }

    /**
    * Merges the basecluster with the specified basecluster which anchorpoints
    * are inserted into the basecluster.
    * the basecluster is deleted afterwards.
    */
    void mergeWith(BaseCluster& cluster);

    /**
    * Calculates the distance from this basecluster to the specified basecluster
    */
    double distanceToCluster(BaseCluster& cluster);


    /**
    * @return true if the APs that are in the bounding box of cluster1 are in the confidence interval of the second and
    * AND vice versa
    */
    bool partialOverlappingInterval(BaseCluster& cluster);

    /**
    * Returns a pointer to the multiplicon which the basecluster is part of
    */
    Multiplicon* getMultiplicon() const {
        return multiplicon;
    }

    /**
    * Sets the corresponding multiplicon
    */
    void setMultiplicon(Multiplicon& multiplicon);

    /**
    * Removes a corresponding multiplicon with this basecluster
    */
    void resetMultiplicon() {
        multiplicon = NULL;
    }

    /**
    * Calculates the probability to be generated by chance.
    * the area and number of points in the GHM needs to be passed as arguments
    */
    double calculateProbability(double area, double count_points,
                                int level);

    /**
    * Calculates the probability to be generated by chance.
    * the area and number of points in the GHM needs to be passed as arguments
    */
    //NOTE apparently this code is not used
    //double calculateProbability2(int area, int count_points, int level);

    /**
     * Uses Probability calculation designed for Cloud methods (explanation in synthenicCloud.h)
     */
    double calculateProbabilityBinomialD(double area, double count_points, int level) const;
    double calculateProbabilityBinomialDCorr(double area, double count_points, int level) const;

    /**
    * Sets the specified probability
    */
    void setRandomProbability(double probability) {
        random_probability = probability;
    }

    /**
    * Returns the probability
    */
    double getRandomProbability() const {
        return random_probability;
    }

    /**
    * Returns a boolean saying if the basecluster was twisted or not
    */
    bool wasTwisted() const {
        return was_twisted;
    }

    /**
    * Fills in the twisted boolean with the specified argument
    */
    void wasTwisted(bool twisted) {
        was_twisted = twisted;
    }

    /**
    * Mirrors the y-value of each anchorpoint around the central y-value of the basecluster
    */
    void twistCluster();

    double getXEnd1() const {
        return x_end1;
    }
    double getXEnd2() const {
        return x_end2;
    }
    double getYEnd1() const {
        return y_end1;
    }
    double getYEnd2() const {
        return y_end2;
    }

    /**
    * For debuging purposes, this writes information about the basecluster
    * and each anchorpoint to the console.
    */
    void write() const;

    /*
     * Calculate and return the package size (in byte)
     */
    int getPackSize() const;

    /*
     * Pack data in a char stream
     * @param buffer Pre-allocated buffer to store the data in
     * @return The number of chars in the buffer that were used
     */
    int pack(char *buffer) const;

    friend bool operator==(const BaseCluster &lhs, const BaseCluster &rhs);

    friend bool operator!=(const BaseCluster &lhs, const BaseCluster &rhs) {
        return !(lhs == rhs);
    }

    /**
    * Return the ID of the x-genelist
    */
    int getXObjectID() const {
        cerr << "No x-object used in BaseCluster" << endl;
        exit(EXIT_FAILURE);
    }

    /**
    * Return the ID of the y-genelist
    */
    int getYObjectID() const {
        cerr << "No x-object used in BaseCluster" << endl;
        exit(EXIT_FAILURE);
    }

    void setID(int ID) {
        clusterID = ID;
    }

    uint getID() {
        return clusterID;
    }

private:
    ///////////////////
    //PRIVATE METHODS//
    ///////////////////

    /**
    * Calculates the best fit line through all anchorpoints of the basecluster
    */
    void regression();

    /**
     * distance of P3 to segment P1P2
     */
    double pointLineDpd(double x1, double y1, double x2,
                        double y2, double x3, double y3);

    //////////////
    //ATTRIBUTES//
    //////////////

    multiset<AnchorPoint> anchorpoints;
    set<AnchorPoint> backBone;
    Multiplicon* multiplicon;

    double random_probability;
    bool orientation;

    double a, b; //coefficients of the best fit line a+bx
    double avg_x; 
    double var_x;
    double mrss; //mean residual sum of squares
    double x_end1, x_end2, y_end1, y_end2; //coordinates of outer points of regression line

    bool was_twisted;

    uint clusterID;    // unique identifier for the basecluster
};

#endif
