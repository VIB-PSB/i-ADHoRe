#include "BaseCluster.h"

#include "AnchorPoint.h"
#include "hpmath.h"
#include <algorithm>

BaseCluster::BaseCluster(const bool orient) : multiplicon(NULL),
        orientation(orient), b(0.0), was_twisted(false)
{

}

BaseCluster::BaseCluster(const char *buffer) : multiplicon(NULL)
{
    memcpy(&random_probability, buffer, sizeof(random_probability));
    buffer += sizeof(random_probability);
    memcpy(&orientation, buffer, sizeof(orientation));
    buffer += sizeof(orientation);
    memcpy(&a, buffer, sizeof(a));
    buffer += sizeof(a);
    memcpy(&b, buffer, sizeof(b));
    buffer += sizeof(b);
    memcpy(&avg_x, buffer, sizeof(avg_x));
    buffer += sizeof(avg_x);
    memcpy(&var_x, buffer, sizeof(var_x));
    buffer += sizeof(var_x);
    memcpy(&mrss, buffer, sizeof(mrss));
    buffer += sizeof(mrss);
    memcpy(&x_end1, buffer, sizeof(x_end1));
    buffer += sizeof(x_end1);
    memcpy(&x_end2, buffer, sizeof(x_end2));
    buffer += sizeof(x_end2);
    memcpy(&y_end1, buffer, sizeof(y_end1));
    buffer += sizeof(y_end1);
    memcpy(&y_end2, buffer, sizeof(y_end2));
    buffer += sizeof(y_end2);
    memcpy(&was_twisted, buffer, sizeof(was_twisted));
    buffer += sizeof(was_twisted);
    memcpy(&begin_x, buffer, sizeof(begin_x));
    buffer += sizeof(begin_x);
    memcpy(&end_x, buffer, sizeof(end_x));
    buffer += sizeof(end_x);
    memcpy(&begin_y, buffer, sizeof(begin_y));
    buffer += sizeof(begin_y);
    memcpy(&end_y, buffer, sizeof(end_y));
    buffer += sizeof(end_y);

    int size = anchorpoints.size();
    memcpy(&size, buffer, sizeof(size));
    buffer += sizeof(size);

    for (int i = 0; i < size; i++) {
        anchorpoints.insert(AnchorPoint(buffer));
        buffer += AnchorPoint::getPackSize();
    }
}

void BaseCluster::addAnchorPoint(int X, int Y) {
    anchorpoints.insert(AnchorPoint(X, Y, true));
    b = 0.0; //indicates the statistics need to be updated
}

struct CompareAP {
    bool operator()(const AnchorPoint& left, const AnchorPoint& right) const
    {
        if (left.getX() != right.getX())
                return left.getX() < right.getX();
        return left.getY() < right.getY();
    }
};

void BaseCluster::filterDuplicates() {
    set<AnchorPoint, CompareAP> unique(anchorpoints.begin(), anchorpoints.end());
    anchorpoints = multiset<AnchorPoint>(unique.begin(), unique.end());
}

void BaseCluster::addBackBone(int X, int Y) {
    backBone.insert(AnchorPoint(X, Y, true));
    //b = 0.0; //indicates the statistics need to be updated
}

void BaseCluster::addAnchorPoint(AnchorPoint& point) {
    anchorpoints.insert(point);
    b = 0.0; //indicates the statistics need to be updated
}

unsigned int BaseCluster::getCountAnchorPoints() const
{
    unsigned int c = 0;

    multiset<AnchorPoint>::const_iterator e = anchorpoints.begin();
    for ( ; e != anchorpoints.end(); e++)
        if (e->isRealAnchorPoint())
            c++;

    return c;
}

unsigned int BaseCluster::getLowestX() const
{
    if (anchorpoints.size() == 0) return 0;

    return anchorpoints.begin()->getX();
}

unsigned int BaseCluster::getLowestY() const
{
    if (anchorpoints.size() == 0) return 0;

    multiset<AnchorPoint>::const_iterator e = anchorpoints.begin();
    int lowest_y = e->getY();

    for (e++; e != anchorpoints.end(); e++) {
        int y = e->getY();
        if (y < lowest_y)
            lowest_y = y;
    }
    return lowest_y;
}

unsigned int BaseCluster::getHighestX() const
{
    if (anchorpoints.size() == 0) return 0;

    return (--anchorpoints.end())->getX();
}

unsigned int BaseCluster::getHighestY() const {
    if (anchorpoints.size() == 0) return 0;

    multiset<AnchorPoint>::const_iterator e = anchorpoints.begin();
    int highest_y = e->getY();

    for (e++; e != anchorpoints.end(); e++) {
        int y = e->getY();
        if (y > highest_y)
            highest_y = y;
    }
    return highest_y;
}

void BaseCluster::setMultiplicon(Multiplicon& multiplicon) {
    this->multiplicon = &multiplicon;
}

double BaseCluster::r_squared(int X, int Y) const {
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
    unsigned int n = 0;

    multiset<AnchorPoint>::const_iterator e = anchorpoints.begin();
    for ( ; e != anchorpoints.end(); e++) {
        int x = e->getX();
        int y = e->getY();
        sum_x += x;
        sum_y += y;
        sum_xy += x*y;
        sum_x2 += x*x;
        sum_y2 += y*y;
        n++;
    }
    //if X and Y are given as arguments, count them too
    if (X != 0 || Y != 0) {
        sum_x += X;
        sum_y += Y;
        sum_xy += X*Y;
        sum_x2 += X*X;
        sum_y2 += Y*Y;
        n++;
    }

    double value = ((sum_x2 - sum_x * sum_x / n) * (sum_y2 - sum_y * sum_y / n));

    if (value == 0) {
        return -1;
    }
    else {
        double r = (sum_xy - (sum_x * sum_y / n)) / sqrt(value);
        return (r * r);
    }
}

double BaseCluster::r_squared(BaseCluster& cluster) const {
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
    unsigned int n = 0;

    multiset<AnchorPoint>::const_iterator e = anchorpoints.begin();
    for ( ; e != anchorpoints.end(); e++) {
        int x = e->getX();
        int y = e->getY();
        sum_x += x;
        sum_y += y;
        sum_xy += x*y;
        sum_x2 += x*x;
        sum_y2 += y*y;
        n++;
    }

    //count the cluster too
    e = cluster.anchorpoints.begin();
    for (; e != cluster.anchorpoints.end(); e++) {
        int x = e->getX();
        int y = e->getY();
        sum_x += x;
        sum_y += y;
        sum_xy += x*y;
        sum_x2 += x*x;
        sum_y2 += y*y;
        n++;
    }

    double value = ((sum_x2 - sum_x * sum_x / n) *
                    (sum_y2 - sum_y * sum_y / n));

    if (value == 0) {
        return -1;
    }
    else {
        double r = (sum_xy - (sum_x * sum_y / n)) / sqrt(value);
        return (r * r);
    }
}

double BaseCluster::averageDPD() const {

    double totaldpd = 0;

    multiset<AnchorPoint>::const_iterator prev = anchorpoints.begin();
    multiset<AnchorPoint>::const_iterator curr = anchorpoints.begin();
    for (curr++; curr != anchorpoints.end(); prev++, curr++) {
        totaldpd += dpd(prev->getX(), prev->getY(),
                        curr->getX(), curr->getY());
    }

    return totaldpd / anchorpoints.size();
}


void BaseCluster::updateStatistics() {

    // if it's already calculated, get out of here
    if (b != 0.0) return;

    int N=anchorpoints.size();
    regression();

    var_x= 0.0;
    mrss = 0.0;

    multiset<AnchorPoint>::const_iterator e = anchorpoints.begin();
    for ( ; e != anchorpoints.end(); e++) {
        int x = e->getX();
        int y = e->getY();
        double y_est = b * x + a;
        mrss += (y - y_est) * (y - y_est);
        var_x += (x - avg_x) * (x - avg_x);
    }

    mrss   = mrss  / (N-2);
    var_x  = var_x / (N-1);
    x_end1 = getLowestX();
    x_end2 = getHighestX();
    y_end1 = x_end1 * b + a;
    y_end2 = x_end2 * b + a;
}

void BaseCluster::regression() {

    double N=anchorpoints.size();
    if (N < 3.0) return;

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;

    multiset<AnchorPoint>::const_iterator e = anchorpoints.begin();
    for ( ; e != anchorpoints.end(); e++) {
        int x = e->getX();
        int y = e->getY();
        sum_x += x;
        sum_y += y;
        sum_xy += x*y;
        sum_x2 += x*x;
    }

    avg_x = sum_x/N;
    double avg_y = sum_y/N;

    var_x  = (sum_x2 - N*avg_x*avg_x) /(N-1.0);
    double cov_xy = (sum_xy - N*avg_x*avg_y) /(N-1.0);

    b = cov_xy / var_x;
    a = avg_y - b * avg_x;
}

double BaseCluster::distanceToPoint(int x, int y) const
{
    if (x >= x_end1 && x <= x_end2) {
        return 0;
    }
    else if ((double)x < x_end1) {
        return dpd((double)x, (double)y, x_end1, y_end1);
    }
    else {
        return dpd((double)x, (double)y, x_end2, y_end2);
    }
}

bool BaseCluster::inInterval(int x, int y) {
    double up, down;
    intervalBounds(x, up, down);
    if (y >= down && y <= up) {
        return true;
    }
    else {
        return false;
    }
}

void BaseCluster::intervalBounds(int x, double& up, double& down) {
    updateStatistics();
    double tnm2_95percent = 1.96; 
    //NOTE tdistribution table used from: http://www.sociology.ohio-state.edu/people/ptv/publications/p%20values/p_value_tables.html
    //99.9% = 3.29
   //99% = 2.58
   //95% = 1.96
    double n=anchorpoints.size();
    double predInterval =tnm2_95percent
                            *sqrt(mrss)
                            *sqrt(
                                1.0 + 1.0/n + (x-avg_x)*(x-avg_x)/(n-1)/var_x
                             );
    up = (a + b * x)+predInterval;
    down = (a + b * x)-predInterval;

//here it seems better to use the prediction interval: that is an estimate of y boundaries given a certain x
// the alternative is the use of a confidenceinterval which estimates puts boundaries for the mean value of y
//given x..
}

void BaseCluster::mergeWith(BaseCluster& cluster)
{
    multiset<AnchorPoint>::const_iterator e = cluster.anchorpoints.begin();
    for ( ; e != cluster.anchorpoints.end(); e++) {
        addAnchorPoint(e->getX(), e->getY());
    }

    for (e = cluster.backBone.begin(); e != cluster.backBone.end(); e++) {
        addBackBone(e->getX(), e->getY());
    }

    b=0.0;
    updateStatistics();
}

bool intersect(double x1, double y1, double x2, double y2,
               double x3, double y3, double x4, double y4)
{
    double ud = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
    if (ud != 0) {
        double ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / ud;
        double ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / ud;
        if (ua >= 0 && ua <= 1 && ub >= 0 && ub <= 1)
            return true;
    }
    return false;
}

double min(double a, double b)
{
    return (a < b) ? a : b;
}

double max(double a, double b)
{
    return (a > b) ? a : b;
}

double BaseCluster::pointLineDpd(double x1, double y1, double x2,
                                 double y2, double x3, double y3)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double x0, y0;
    if ((dx == 0) && (dy == 0)) {
        x0 = x1;
        y0 = y1;
    } else {

        //t is the coordinate of the orthogonal projection P0 of P3 on the line P1P2
        // P=P1+|P1P2|*t
        double t = ((x3 - x1) * dx + (y3 - y1) * dy) / (dx * dx + dy * dy);
        t = min(max(0,t),1); //t must be in [0..1]
        x0 = x1 + t * dx;
        y0 = y1 + t * dy;
    }

    return dpd(x0, y0, x3, y3);
}

double BaseCluster::distanceToCluster(BaseCluster& cluster)
{
    updateStatistics();
    cluster.updateStatistics();

    double lowest_distance = 0;

    if (intersect(x_end1, y_end1, x_end2, y_end2, cluster.x_end1,
                  cluster.y_end1, cluster.x_end2, cluster.y_end2))
        return 0;

    lowest_distance = pointLineDpd(x_end1, y_end1, x_end2, y_end2,
                                   cluster.x_end1, cluster.y_end1);

    double d = pointLineDpd(x_end1, y_end1, x_end2, y_end2,
                            cluster.x_end2, cluster.y_end2);
    lowest_distance = min(lowest_distance, d);

    d = pointLineDpd(cluster.x_end1, cluster.y_end1, cluster.x_end2,
                     cluster.y_end2, x_end1, y_end1);
    lowest_distance = min(lowest_distance, d);

    d = pointLineDpd(cluster.x_end1, cluster.y_end1, cluster.x_end2,
                     cluster.y_end2, x_end2, y_end2);
    lowest_distance = min(lowest_distance, d);

    //due to bascule instability we also calculate dpds between outerAPs

    const AnchorPoint& start1=*anchorpoints.begin();
    const AnchorPoint& stop1=*(--anchorpoints.end());
    const AnchorPoint& start2=*cluster.anchorpoints.begin();
    const AnchorPoint& stop2=*(--cluster.anchorpoints.end());

    d=start1.dpdToAP(start2);
    lowest_distance = min(lowest_distance, d);

    d=start1.dpdToAP(stop2);
    lowest_distance = min(lowest_distance, d);

    d=stop1.dpdToAP(start2);
    lowest_distance = min(lowest_distance, d);

    d=stop1.dpdToAP(stop2);
    lowest_distance = min(lowest_distance, d);

    //the dpds between regression line of 1 and outerAPs of 2

    d = pointLineDpd(x_end1, y_end1, x_end2, y_end2,
                                   start2.getX(),start2.getY() );
    lowest_distance = min(lowest_distance, d);

    d = pointLineDpd(x_end1, y_end1, x_end2, y_end2,
                                   stop2.getX(),stop2.getY() );
    lowest_distance = min(lowest_distance, d);

    //and vice versa

    d = pointLineDpd(cluster.x_end1, cluster.y_end1, cluster.x_end2, cluster.y_end2,
                                   start1.getX(),start1.getY() );
    lowest_distance = min(lowest_distance, d);

    d = pointLineDpd(cluster.x_end1, cluster.y_end1, cluster.x_end2, cluster.y_end2,
                                   stop1.getX(),stop1.getY() );
    lowest_distance = min(lowest_distance, d);


    return lowest_distance;
}

int BaseCluster::overlappingCoordinates(BaseCluster& cluster) const
{
    int overlap_status = 0;

    double ax1 = x_end1, ay1 = y_end1, ax2 = x_end2, ay2 = y_end2;
    double bx1 = cluster.x_end1, by1 = cluster.y_end1, bx2 = cluster.x_end2, by2 = cluster.y_end2;

    if (!orientation) {
        ay1 = y_end2;
        ay2 = y_end1;
    }
    if (!cluster.orientation) {
        by1 = cluster.y_end2;
        by2 = cluster.y_end1;
    }

    if (ax1 < bx2 && ax2 > bx1)
        overlap_status++;
    if (ay1 < by2 && ay2 > by1)
        overlap_status++;

    return overlap_status;
}

bool BaseCluster::partialOverlappingInterval(BaseCluster& c) {

    bool overlapping=true;
    updateStatistics();
    c.updateStatistics();

    multiset<AnchorPoint>::const_iterator e;

    //check if AP of cluster c that are in the bounding box of THIS cluster are within the confidenceinterval of THIS cluster

    int x1min, x1max;
    x1min=getLowestX();
    x1max=getHighestX();

    e=c.anchorpoints.begin();

    for ( ; e != c.anchorpoints.end(); e++) {

        int x=e->getX();
        int y=e->getY();
        if ( (x>=x1min) and (x<=x1max)) {
            if (!inInterval(x, y)) {
                overlapping = false;
                break;
            }

        } else
            ;
    }

    if (!overlapping) return overlapping;

    //check if AP of THIS cluster that are in the bounding box of cluster c are within the confidenceinterval of cluster c

    int x2min, x2max;
    x2min=c.getLowestX();
    x2max=c.getHighestX();

    e=anchorpoints.begin();

    for ( ; e != anchorpoints.end(); e++) {

        int x=e->getX();
        int y=e->getY();
        if ( (x>=x2min) and (x<=x2max)) {
            if (!c.inInterval(x,y)) {
                overlapping = false;
                break;
            }

        } else
            ;
    }

    return overlapping;
}

double BaseCluster::calculateProbability(double area, double points, int level)
{
    double density = (double)points / (double)area;
    double probability = 1.0;

    multiset<AnchorPoint>::const_iterator it1 = backBone.begin();
    it1++;
    multiset<AnchorPoint>::const_iterator it2 = backBone.begin();

    while (it1 != backBone.end())
    {
        double distance = dpd((it2)->getX(), (it2)->getY(),
                              (it1)->getX(), (it1)->getY());

        int c = (int)ceil(0.5 * distance * distance);

        // calculate p_i = 1.0 - pow(1.0 - density, c), this is
        // the change to find *at least one* homolog within
        // the DPD distance.
        double p_i = 1.0 - pow(1.0 - density, c);

        probability *= p_i;
        it1++;
        it2++;
    }

    it1 = anchorpoints.begin();
    it1++;
    it2 = anchorpoints.begin();

    while (it1 != anchorpoints.end()) {
        double distance = dpd((it2)->getX(), (it2)->getY(),
                              (it1)->getX(), (it1)->getY());

        if (distance == 0) {
            anchorpoints.erase(it1++);
        } else {
            it1++;
        }
    }

    return probability;
}

/*double BaseCluster::calculateProbability2(int area, int points, int level) {
    points = points - anchorpoints.size();

    double density = (double)points / area;
    double probability = 1.0;

    multiset<AnchorPoint>::const_iterator it1 = anchorpoints.begin();
    it1++;
    multiset<AnchorPoint>::const_iterator it2 = anchorpoints.begin();

    while (it1 != anchorpoints.end()) {

        double distance = dpd((it2)->getX(), (it2)->getY(),
                              (it1)->getX(), (it1)->getY());

        if (distance == 0) {
            anchorpoints.erase(it1++);
        }
        else {
            int c = (int)ceil(0.5 * distance * distance);

            // calculate p_i = 1.0 - pow(1.0 - density, c);
            //  double p_i = ompowopxn(-density, c);

            double p_i = c * density * pow(1.0 - density, c - 1);
            // calculate p_i = 1 - (1 - p_i)^(n) in high precision
            p_i = ompowopxn(-p_i, level - 1);

            probability *= p_i;

            it1++;
            it2++;
        }
    }

    // evaluate logNotPglobal = point*log(1 - probability) in high precision
    double logNotPglobal = points*logopx(-probability);

    // evaluate pGlobal = 1-exp(logNotPglobal) in high precision
    return omexpx(logNotPglobal);
}*/


double BaseCluster::calculateProbabilityBinomialD(double area, double count_points, int level) const
{
    //parameters binomial distribution
    double APDensity=count_points/area;
    double p=APDensity;
    int boxHeight=getHighestY()-getLowestY();
    int boxWidth=getHighestX()-getLowestX();
    int n=boxHeight*boxWidth;
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

double BaseCluster::calculateProbabilityBinomialDCorr(double area, double count_points, int level) const
{
    //parameters corrected binomial distribution
    double APDensity=count_points/area;
    double p=APDensity;

    int a=getHighestX()-getLowestX();
    int b=getHighestY()-getLowestY();

    int minAB=std::min(a,b);
    int n=a*b;
    int x=getCountAnchorPoints();

    double numerator;
    double denominator; //normalization ->due to rescaling of probe universe

    double poomp=p/(1-p);
    double Pimo=omxn(p,n);

    double F=Pimo; //cumulative Binomial Distr!! (will be corrected later on)



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

void BaseCluster::twistCluster() {
    int max_y = getHighestY();
    int min_y = getLowestY();

    // swapping the y-values might change the ordening and breaking
    // the iterator, hence, we must use a temporary object
    multiset<AnchorPoint> twistedAnchorpoints;

    multiset<AnchorPoint>::iterator e = anchorpoints.begin();
    for ( ; e != anchorpoints.end(); e++) {
        AnchorPoint ap = *e;
        ap.twistY(max_y, min_y);
        twistedAnchorpoints.insert(ap);
    }

    anchorpoints = twistedAnchorpoints;

    b = 0.0;
    updateStatistics();
}

void BaseCluster::write () const {
    //cout << "basecluster: x_end1: "<< x_end1 << "\tx_end2: "<< x_end2
    //   << "\ty_end1: " << y_end1 << "\ty_end2: "<< y_end2 << endl;
    cout << "basecluster: begin_x: "<< begin_x << "\tbegin_y: "<< begin_y
         << "\tend_x: " << end_x << "\tend_y: "<< end_y << endl;
    cout << "anchorpoints: " << endl;

    multiset<AnchorPoint>::const_iterator e = anchorpoints.begin();
    for ( ; e != anchorpoints.end(); e++) {
        cout << "------" << &(*e) << '\t' <<
             e->getX() << ' ' << e->getY() << endl;
    }
}

int BaseCluster::getPackSize() const
{
    int packSize = 0;

    packSize += sizeof(random_probability) + sizeof(orientation) +
                sizeof(a) + sizeof(b) + sizeof (avg_x) + sizeof(var_x) +
                sizeof(mrss) + sizeof(x_end1) + sizeof(x_end2) +
                sizeof(y_end1) + sizeof(y_end2) + sizeof(was_twisted) +
                sizeof(begin_x) + sizeof(end_x) +
                sizeof(begin_y) + sizeof(end_y);
    packSize += sizeof(int); // for storing the number of anchorpoints
    packSize += anchorpoints.size()*AnchorPoint::getPackSize();

    return packSize;
}

int BaseCluster::pack(char *buffer) const
{
    const char *bufferOrig = buffer;

    memcpy(buffer, &random_probability, sizeof(random_probability));
    buffer += sizeof(random_probability);
    memcpy(buffer, &orientation, sizeof(orientation));
    buffer += sizeof(orientation);
    memcpy(buffer, &a, sizeof(a));
    buffer += sizeof(a);
    memcpy(buffer, &b, sizeof(b));
    buffer += sizeof(b);
    memcpy(buffer, &avg_x, sizeof(avg_x));
    buffer += sizeof(avg_x);
    memcpy(buffer, &var_x, sizeof(var_x));
    buffer += sizeof(var_x);
    memcpy(buffer, &mrss, sizeof(mrss));
    buffer += sizeof(mrss);
    memcpy(buffer, &x_end1, sizeof(x_end1));
    buffer += sizeof(x_end1);
    memcpy(buffer, &x_end2, sizeof(x_end2));
    buffer += sizeof(x_end2);
    memcpy(buffer, &y_end1, sizeof(y_end1));
    buffer += sizeof(y_end1);
    memcpy(buffer, &y_end2, sizeof(y_end2));
    buffer += sizeof(y_end2);
    memcpy(buffer, &was_twisted, sizeof(was_twisted));
    buffer += sizeof(was_twisted);
    memcpy(buffer, &begin_x, sizeof(begin_x));
    buffer += sizeof(begin_x);
    memcpy(buffer, &end_x, sizeof(end_x));
    buffer += sizeof(end_x);
    memcpy(buffer, &begin_y, sizeof(begin_y));
    buffer += sizeof(begin_y);
    memcpy(buffer, &end_y, sizeof(end_y));
    buffer += sizeof(end_y);

    int size = anchorpoints.size();
    memcpy(buffer, &size, sizeof(size));
    buffer += sizeof(size);

    multiset<AnchorPoint>::const_iterator it = anchorpoints.begin();
    for ( ; it != anchorpoints.end(); it++)
        buffer += (*it).pack(buffer);

    return buffer - bufferOrig;
}

bool operator==(const BaseCluster &lhs, const BaseCluster &rhs)
{
    if (lhs.random_probability != rhs.random_probability) return false;
    if (lhs.orientation != rhs.orientation) return false;
    if (lhs.a != rhs.a) return false;
    if (lhs.b != rhs.b) return false;
    if (lhs.avg_x != rhs.avg_x) return false;
    if (lhs.var_x != rhs.var_x) return false;
    if (lhs.mrss != rhs.mrss) return false;
    if (lhs.x_end1 != rhs.x_end1) return false;
    if (lhs.x_end2 != rhs.x_end2) return false;
    if (lhs.y_end1 != rhs.y_end1) return false;
    if (lhs.y_end2 != rhs.y_end2) return false;
    if (lhs.was_twisted != rhs.was_twisted) return false;
    // base class comparison
    if (lhs.begin_x != rhs.begin_x) return false;
    if (lhs.end_x != rhs.end_x) return false;
    if (lhs.begin_y != rhs.begin_y) return false;
    if (lhs.end_y != rhs.end_y) return false;

    if (lhs.anchorpoints.size() != rhs.anchorpoints.size()) return false;
    multiset<AnchorPoint>::const_iterator it1 = lhs.anchorpoints.begin();
    multiset<AnchorPoint>::const_iterator it2 = rhs.anchorpoints.begin();

    for ( ; it1 != lhs.anchorpoints.end(); it1++, it2++)
        if (*it1 != *it2) return false;

    return true;
}
