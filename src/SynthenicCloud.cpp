#include "SynthenicCloud.h"
#include "hpmath.h"
#include <cassert>

SynthenicCloud::SynthenicCloud()
{

}

SynthenicCloud::SynthenicCloud(const SynthenicCloud& sCloud)  : Cluster(sCloud.getXObjectID(),sCloud.getYObjectID())
, random_probability(-1.0)
{
    vector<AnchorPoint>::const_iterator it=sCloud.getAPBegin();
    for (; it!=sCloud.getAPEnd(); it++) {
        addAnchorPoint(*it);
    }
}

SynthenicCloud& SynthenicCloud::operator=(const SynthenicCloud& sCloud)
{
    //check for self-assignment
    if (&sCloud==this)
        return (*this);

    //destruct current object
    anchorPoints.clear();

    x_objectID=sCloud.getXObjectID();
    y_objectID=sCloud.getYObjectID();

    vector<AnchorPoint>::const_iterator it=sCloud.getAPBegin();
    for (; it!=sCloud.getAPEnd(); it++) {
        addAnchorPoint(*it);
    }
    return (*this);
}

SynthenicCloud::~SynthenicCloud()
{
    anchorPoints.clear();
}

void SynthenicCloud::addAnchorPoint (int x,int y)
{
    modifyBoundingBox(x,y);
    anchorPoints.push_back(AnchorPoint(x,y,true));
}

void SynthenicCloud::addAnchorPoint(const AnchorPoint& point)
{
    modifyBoundingBox(point.getX(),point.getY());
    anchorPoints.push_back(point);
}

uint SynthenicCloud::getCountAnchorPoints() const
{
    return anchorPoints.size();
}

uint SynthenicCloud::getLowestX() const
{
    if (anchorPoints.size()==0)
        return 0;

    uint lowestX=anchorPoints[0].getX();

    for (int i=1; i<anchorPoints.size(); i++) {
        const AnchorPoint& AP=anchorPoints[i];
        if (AP.getX()<lowestX)
            lowestX=AP.getX();
    }
    return lowestX;
}

uint SynthenicCloud::getLowestY() const
{
    if (anchorPoints.size()==0)
        return 0;

    uint lowestY=anchorPoints[0].getY();

    for (int i=1; i<anchorPoints.size(); i++) {
        const AnchorPoint& AP=anchorPoints[i];
        if (AP.getY()<lowestY) lowestY=AP.getY();
    }
    return lowestY;
}

uint SynthenicCloud::getHighestX() const
{
    if (anchorPoints.size()==0)
        return 0;
    uint highestX=anchorPoints[0].getX();

    for (int i=1; i<anchorPoints.size(); i++) {
        const AnchorPoint& AP=anchorPoints[i];
        if (AP.getX()>highestX) highestX=AP.getX();
    }
    return highestX;
}

uint SynthenicCloud::getHighestY() const
{
    if (anchorPoints.size()==0)
        return 0;

    uint highestY=anchorPoints[0].getY();

    for (int i=1; i<anchorPoints.size(); i++) {
        const AnchorPoint& AP=anchorPoints[i];
        if (AP.getY()>highestY) highestY=AP.getY();
    }
    return highestY;
}

double SynthenicCloud::calculateProbabilityBinomialD(double APDensity) const
{
    //parameters binomial distribution
    double p=APDensity;
    int n=calculateBoxHeight()*calculateBoxWidth();
    int x=getCountAnchorPoints();

    double poomp=p/(1.0-p);
    double Pimo=omxn(p,n);

    double F=Pimo; //cumulative distrib (first term)

    //calc F(x-1) recursively
    for (int i=1; i<x; i++){
        double Pi=Pimo*poomp*(n-i+1.0)/(1.0*i);
        F+=Pi;
        Pimo=Pi;
    }
    return 1.0-F; //P(x>=X)
}

double SynthenicCloud::calculateProbabilityBinomialDCorr(double APDensity) const
{
    //parameters corrected binomial distribution
    double p=APDensity;

    int a=calculateBoxWidth();
    int b=calculateBoxHeight();

    int minAB=std::min(a,b);
    int n=a*b;
    int x=getCountAnchorPoints();

    double numerator;
    double denominator; //normalization ->due to rescaling of probe universe

    double poomp=p/(1-p);
    double Pimo=omxn(p,n);

    double F=Pimo; //cumulative Binomial Distr!! (will be corrected later on)

    /*if (x>minAB)
        cout << "calculateProbabilityBinomialDCorr fails x>min(a,b)!!" << endl;*/
    //NOTE No warning necessary since this means the cloud is definitely dense enough
    //it means that there must somewhere be a tandem in the cloud, the cumulative probability will be
    //equal to one

    //calc numerator F(x-1) recursively
    for (int i=1; i<x; i++){
        double Pi=Pimo*poomp*(a-i+1.0)*(b-i+1.0)/(1.0*i);
        F+=Pi;
        Pimo=Pi;
    }

    numerator=F;

    for (int i=x; i<minAB; i++){
        double Pi=Pimo*poomp*(a-i+1.0)*(b-i+1.0)/(1.0*i);
        F+=Pi;
        Pimo=Pi;
    }
    denominator=F;
    return (1.0-numerator/denominator);
}

double SynthenicCloud::calculateCloudDensity() const
{
    int dimx,dimy;
    dimx=calculateBoxWidth();
    dimy=calculateBoxHeight();
    int maxAP= ( dimx>dimy ) ? dimy : dimx;

    return ( getCountAnchorPoints() *1.0 ) / ( maxAP*1.0 );
}

void SynthenicCloud::modifyBoundingBox(int x,int y)
{
    //first AP
    if (anchorPoints.size()==0) {
        begin_x=x;
        end_x=x;
        begin_y=y;
        end_y=y;
    } else {

        if (x>end_x) end_x=x;
        else if (x<begin_x) begin_x=x;

        if (y>end_y) end_y=y;
        else if (y<begin_y) begin_y=y;
    }
}

vector<AnchorPoint> SynthenicCloud::calcAPInOuterFrame(uint frameThickness) const
{
    vector<AnchorPoint> APInFrame;
    AnchorPoint AP(0,0,true);

    for (int i=0; i<anchorPoints.size(); i++) {

        AP=anchorPoints[i];
        if      (AP.getX() < begin_x+frameThickness) APInFrame.push_back(AP);
        else if (AP.getX() > end_x-frameThickness) APInFrame.push_back(AP);
        else if (AP.getY() < begin_y+frameThickness) APInFrame.push_back(AP);
        else if (AP.getY() > end_y-frameThickness) APInFrame.push_back(AP);
    }

    return APInFrame;
}

uint SynthenicCloud::distanceToCloud(int x, int y) const
{
    assert(anchorPoints.size()!=0);

    int distance=Cluster::kspd(x,y,anchorPoints[0].getX(),anchorPoints[0].getY());
    int testDist;
    for (int i=1; i<anchorPoints.size(); i++) {
        const AnchorPoint& AP=anchorPoints[i];
        testDist=Cluster::kspd(x,y,AP.getX(),AP.getY());
        if (distance>testDist) distance=testDist;
    }
    return distance;
}

bool SynthenicCloud::coordInCloudBox(int x, int y) const
{
    if ((x>=begin_x) and (x<=end_x)){
        if ((y>=begin_y) and (y<=end_y))
            return true;
        else
            return false;
    } else
        return false;
}

/**MPI functions**/

int SynthenicCloud::getPackSize() const
{

    int packSize = 0;

    //clusterData
    packSize += sizeof(begin_x) + sizeof(end_x) +
                sizeof(begin_y) + sizeof(end_y) +
                sizeof(x_objectID)+sizeof(y_objectID);

    packSize += sizeof(random_probability);


    packSize += sizeof(int); // for storing the number of anchorpoints

    //anchorPoints in cloud
    packSize += anchorPoints.size()*AnchorPoint::getPackSize();
    return packSize;
}


int SynthenicCloud::pack(char *buffer) const
{
    const char *bufferOrig = buffer;

    memcpy(buffer, &x_objectID, sizeof(x_objectID));
    buffer += sizeof(x_objectID);
    memcpy(buffer, &y_objectID, sizeof(y_objectID));
    buffer += sizeof(y_objectID);

    memcpy(buffer, &begin_x, sizeof(begin_x));
    buffer += sizeof(begin_x);
    memcpy(buffer, &end_x, sizeof(end_x));
    buffer += sizeof(end_x);
    memcpy(buffer, &begin_y, sizeof(begin_y));
    buffer += sizeof(begin_y);
    memcpy(buffer, &end_y, sizeof(end_y));
    buffer += sizeof(end_y);

    memcpy(buffer, &random_probability, sizeof(random_probability));
    buffer += sizeof(random_probability);

    int size = anchorPoints.size();
    memcpy(buffer, &size, sizeof(size));
    buffer += sizeof(size);

    for (int i=0; i<anchorPoints.size(); i++)
        buffer += anchorPoints[i].pack(buffer);

    int bufSize=buffer-bufferOrig;
    if (bufSize!=getPackSize()){
        cout << "bufSize= " << bufSize << endl;
        cout << "getPackSize()= " << getPackSize() << endl;
        cout << endl;
        cout << "anchorPoints.size= " << size << endl;
        cout << endl;
        cout << "sizeof(x_objectID)=" << sizeof(x_objectID) << endl;
        cout << "sizeof(y_objectID)= " << sizeof(y_objectID) << endl;
        cout << "sizeof(begin_x)= " << sizeof(begin_x) << endl;
        cout << "sizeof(end_x)= " << sizeof(end_x) << endl;
        cout << "sizeof(begin_y)= " << sizeof(begin_y) << endl;
        cout << "sizeof(end_y)= " << sizeof(end_y) << endl;
        cout << "sizeof(random_probability)= " << sizeof(random_probability) << endl;
        cout << endl;
        cout << "sizeof(size)= " << sizeof(size) << endl;
        cout << "sizeof(int)= " << sizeof(int) << endl;
        cout << endl;
        cout << "size of pack(buffer) statement = " << anchorPoints[0].pack(buffer) << endl;
        assert(bufSize==getPackSize());
    }

    return bufSize;

}


int SynthenicCloud::getPackSize(const vector<SynthenicCloud*> &clouds)
{
    int packSize = 0;
    packSize += sizeof(int);    // for storing size
    for (int i=0; i<clouds.size(); i++)
        packSize += clouds[i]->getPackSize();

    return packSize;
}


void SynthenicCloud::packSynthenicClouds(const vector<SynthenicCloud*> &clouds,
                                char* buffer)
{
    int nClouds = clouds.size();
    memcpy(buffer, &nClouds, sizeof(nClouds));
    buffer += sizeof(nClouds);

    for (int i=0; i<clouds.size(); i++)
        buffer += clouds[i]->pack(buffer);

}

void SynthenicCloud::unpackSynthenicClouds(const char *buffer,
                                    vector<SynthenicCloud*>& clouds)
{
    int nClouds = 0;
    memcpy(&nClouds, buffer, sizeof(nClouds));
    buffer += sizeof(nClouds);

    clouds.reserve(clouds.size() + nClouds);
    for (int i = 0; i < nClouds; i++) {
        clouds.push_back(new SynthenicCloud(buffer));
        buffer += clouds.back()->getPackSize();
    }


}

SynthenicCloud::SynthenicCloud(const char *buffer) : Cluster()
{

    memcpy(&x_objectID, buffer, sizeof(x_objectID));
    buffer += sizeof(x_objectID);
    memcpy(&y_objectID, buffer, sizeof(y_objectID));
    buffer += sizeof(y_objectID);

    memcpy(&begin_x, buffer, sizeof(begin_x));
    buffer += sizeof(begin_x);
    memcpy(&end_x, buffer, sizeof(end_x));
    buffer += sizeof(end_x);
    memcpy(&begin_y, buffer, sizeof(begin_y));
    buffer += sizeof(begin_y);
    memcpy(&end_y, buffer, sizeof(end_y));
    buffer += sizeof(end_y);

    memcpy(&random_probability, buffer, sizeof(random_probability));
    buffer += sizeof(random_probability);

    int nAP=0;
    memcpy(&nAP, buffer, sizeof(nAP));
    buffer += sizeof(nAP);
    for (int i=0; i<nAP; i++) {
       anchorPoints.push_back(AnchorPoint(buffer));
       buffer += AnchorPoint::getPackSize();
   }
}

double SynthenicCloud::getRandomProbability() const
{
   if (!probabilitySet()) std::cerr << "Cloud probability was not calculated!" << endl;
   return random_probability;
}

void SynthenicCloud::setRandomProbabilityBinomialD(double APDensity)
{
    random_probability=calculateProbabilityBinomialD(APDensity);
}

void SynthenicCloud::setRandomProbabilityBinomialDCorr(double APDensity)
{
    random_probability=calculateProbabilityBinomialDCorr(APDensity);
}
