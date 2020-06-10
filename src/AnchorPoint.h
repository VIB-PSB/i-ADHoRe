#ifndef __ANCHORPOINT_H
#define __ANCHORPOINT_H

#include <cstring>
#include "Cluster.h"

class AnchorPoint
{
public:
	///////////////////////////////
	//CONSTRUCTORS AND DESTRUCTOR//
	///////////////////////////////

	/**
	* Constructs an anchorpoint without known geneIDs
	*/
	AnchorPoint(int X, int Y, bool isReal) : geneXID(-1), geneYID(-1),
		x(X), y(Y), is_real(isReal) {}

	/**
	* Constructs an anchorpoint with known geneIDs
	*/
	AnchorPoint(int geneXID_, int geneYID_, int X, int Y, bool isReal) :
		geneXID(geneXID_), geneYID(geneYID_), x(X), y(Y),
		is_real(isReal) {}

	/**
	* Constructs an anchorpoint object from a packed buffer
	*/
	AnchorPoint(const char *buffer)
	{
		memcpy(&geneXID, buffer, sizeof(geneXID));
		buffer += sizeof(geneXID);
		memcpy(&geneYID, buffer, sizeof(geneYID));
		buffer += sizeof(geneYID);
		memcpy(&x, buffer, sizeof(x));
		buffer += sizeof(x);
		memcpy(&y, buffer, sizeof(y));
		buffer += sizeof(y);
		memcpy(&is_real, buffer, sizeof(is_real));
		is_real = is_real;
	}

	//////////////////
	//PUBLIC METHODS//
	//////////////////

	/**
	* Returns the x-gene ID
	*/
	int getGeneXID() const { return geneXID; }

	/**
	* Returns the y-gene ID
	*/
	int getGeneYID() const { return geneYID; }

	/**
	* Returns the x coordinate
	*/
	int getX() const { return x; }

	/**
	* Returns the y coordinate
	*/
	int getY() const { return y; }

	/**
	* Returns true if the anchorpoint is a real anchorpoint
	*/
	bool isRealAnchorPoint() const { return is_real; }

	/**
	* Mirrors the y-value around the center
	*/
	void twistY(int max_y, int min_y) { y = max_y + min_y - y; }

	/**
	* Sets the corresponding gene IDs of this anchorpoint
	*/
	void setGeneIDs(int geneXID_, int geneYID_) {
		geneXID = geneXID_;
		geneYID = geneYID_;
	}

	/**
	 * Return the size of the packed anchorpoint
	 */
	static int getPackSize() {
		return sizeof(int) + sizeof(int) + sizeof(int) +
			sizeof(int) + sizeof(bool);
	}

	double dpdToAP(const AnchorPoint& a) const {
        return Cluster::dpd(x,y,a.getX(),a.getY());
    }

	/**
	 * Pack data in a char stream
	 * @param buffer Pre-allocated buffer to store the data in
	 * @return The number of chars in the buffer that were used
	 */
	int pack(char *buffer) const {
		const char *bufferOrig = buffer;

		memcpy(buffer, &geneXID, sizeof(geneXID));
		buffer += sizeof(geneXID);
		memcpy(buffer, &geneYID, sizeof(geneYID));
		buffer += sizeof(geneYID);
		memcpy(buffer, &x, sizeof(x));
		buffer += sizeof(x);
		memcpy(buffer, &y, sizeof(y));
		buffer += sizeof(y);
		memcpy(buffer, &is_real, sizeof(is_real));
		buffer += sizeof(is_real);

		return buffer - bufferOrig;
	}

	friend bool operator==(const AnchorPoint &lhs, const AnchorPoint &rhs);
	friend bool operator!=(const AnchorPoint &lhs, const AnchorPoint &rhs) {
		return !(lhs == rhs);
	}
        friend bool operator<(const AnchorPoint &lhs, const AnchorPoint &rhs);
	friend bool operator>(const AnchorPoint &lhs, const AnchorPoint &rhs);

private:
	//////////////
	//ATTRIBUTES//
	//////////////

	int geneXID;	// unique identifier for the x-gene
	int geneYID;	// unique identifier for the y-gene
	int x, y;	// x and y position of the gene
	bool is_real;	// is it a real anchorpoint or not
};

inline bool operator==(const AnchorPoint &lhs, const AnchorPoint &rhs)
{
	return ((lhs.x == rhs.x) && (lhs.y == rhs.y) &&
		(lhs.is_real == rhs.is_real));
}

inline bool operator<(const AnchorPoint &lhs, const AnchorPoint &rhs)
{
    /*
	// a real anchorpoint is always smaller than a false one
	if (lhs.is_real != rhs.is_real)
		return (lhs.is_real);
	// if both are false, there is no ordening
	if (!lhs.is_real)
		return false;
    */
	// if both are real, sort them on x-value
	return (lhs.x < rhs.x);
}

inline bool operator>(const AnchorPoint &lhs, const AnchorPoint &rhs)
{
    /* //NOTE this ordening is no longer relevant, it also disturbs the bounding box calculation!
	// a real anchorpoint is always smaller than a false one
	if (lhs.is_real != rhs.is_real)
		return (!lhs.is_real);
	// if both are false, there is no ordening
	if (!lhs.is_real)
		return false;
    */
	// if both are real, sort them on x-value
	return (lhs.x > rhs.x);
}

#endif
