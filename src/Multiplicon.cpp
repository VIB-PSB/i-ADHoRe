#include "AnchorPoint.h"
#include "Multiplicon.h"

#include "headers.h"
#include "GeneList.h"
#include "BaseCluster.h"
#include "ListElement.h"
#include "Profile.h"

typedef vector<ListElement* >::const_iterator VecListElementCIt;

Multiplicon::Multiplicon(int xObjectID, int yObjectID, int _level)
    : Cluster(xObjectID, yObjectID), level(_level), profile(NULL),
      multipliconID(0), isRedundant(false), parentID(0), ySegment(NULL)
{

}

Multiplicon::Multiplicon(const char *buffer,
                         const vector<GeneList* >& genelists,
                         bool useFamily) :
    profile(NULL), multipliconID(0), parentID(0), ySegment(NULL)
{
    memcpy(&x_objectID, buffer, sizeof(x_objectID));
    buffer += sizeof(x_objectID);
    memcpy(&y_objectID, buffer, sizeof(y_objectID));
    buffer += sizeof(y_objectID);
    memcpy(&level, buffer, sizeof(level));
    buffer += sizeof(level);
    memcpy(&isRedundant, buffer, sizeof(isRedundant));
    buffer += sizeof(isRedundant);
    memcpy(&begin_x, buffer, sizeof(begin_x));
    buffer += sizeof(begin_x);
    memcpy(&end_x, buffer, sizeof(end_x));
    buffer += sizeof(end_x);
    memcpy(&begin_y, buffer, sizeof(begin_y));
    buffer += sizeof(begin_y);
    memcpy(&end_y, buffer, sizeof(end_y));
    buffer += sizeof(end_y);

    int size = 0;
    memcpy(&size, buffer, sizeof(size));
    buffer += sizeof(size);

    baseclusters.reserve(size);
    for (int i = 0; i < size; i++) {
        BaseCluster *cluster = new BaseCluster(buffer);
        addBaseCluster(*cluster);
        buffer += cluster->getPackSize();
    }

    extractXObject(*genelists[x_objectID]);
    extractYObject(*genelists[y_objectID]);

    // ignore the return value of createHomologs
    createHomologs(useFamily, 0);
}

Multiplicon::Multiplicon(const char *buffer, const Profile &xObject,
                         const vector<GeneList* >& genelists, bool useFamily) :
    profile(NULL), multipliconID(0), parentID(0), ySegment(NULL)
{
    memcpy(&x_objectID, buffer, sizeof(x_objectID));
    buffer += sizeof(x_objectID);
    memcpy(&y_objectID, buffer, sizeof(y_objectID));
    buffer += sizeof(y_objectID);
    memcpy(&level, buffer, sizeof(level));
    buffer += sizeof(level);
    memcpy(&isRedundant, buffer, sizeof(isRedundant));
    buffer += sizeof(isRedundant);
    memcpy(&begin_x, buffer, sizeof(begin_x));
    buffer += sizeof(begin_x);
    memcpy(&end_x, buffer, sizeof(end_x));
    buffer += sizeof(end_x);
    memcpy(&begin_y, buffer, sizeof(begin_y));
    buffer += sizeof(begin_y);
    memcpy(&end_y, buffer, sizeof(end_y));
    buffer += sizeof(end_y);

    int size = 0;
    memcpy(&size, buffer, sizeof(size));
    buffer += sizeof(size);

    baseclusters.reserve(size);
    for (int i = 0; i < size; i++) {
        BaseCluster *cluster = new BaseCluster(buffer);
        addBaseCluster(*cluster);
        buffer += cluster->getPackSize();
    }

    extractXObject(xObject);
    extractYObject(*genelists[y_objectID]);

    // ignore the return value of createHomologs
    createHomologs(useFamily, 0);
}

Multiplicon::~Multiplicon()
{
    for (int i = 0; i < xSegments.size(); i++)
        delete xSegments[i];

    delete ySegment;

    vector<BaseCluster*>::iterator it = baseclusters.begin();
    for ( ; it != baseclusters.end(); it++)
        delete (*it);

    delete profile;
}

void Multiplicon::addBaseClusters(Multiplicon& multiplicon) {
    const vector<BaseCluster*>& clusters = multiplicon.getBaseClusters();
    for (unsigned int i = 0; i < clusters.size(); i++) {
        clusters[i]->setMultiplicon(*this);
        baseclusters.push_back(clusters[i]);
    }
    multiplicon.clear();
}

void Multiplicon::addBaseCluster(BaseCluster& cluster) {
    cluster.setMultiplicon(*this);
    baseclusters.push_back(&cluster);
}

unsigned int Multiplicon::getCountAnchorPoints() const {
    unsigned int size = 0;
    for (unsigned int i = 0; i < baseclusters.size(); i++) {
        size += baseclusters[i]->getCountAnchorPoints();
    }

    return size;
}

unsigned int Multiplicon::getLowestX() const {
    unsigned int lowest_x = baseclusters[0]->getLowestX();
    for (unsigned int j = 1; j < baseclusters.size(); j++) {
        if (baseclusters[j]->getCountAnchorPoints() > 0) {
            unsigned int x = baseclusters[j]->getLowestX();
            if (x < lowest_x)
                lowest_x = x;
        }
    }
    return lowest_x;
}

unsigned int Multiplicon::getLowestY() const {
    unsigned int lowest_y = baseclusters[0]->getLowestY();
    for (unsigned int j = 1; j < baseclusters.size(); j++) {
        if (baseclusters[j]->getCountAnchorPoints() > 0) {
            unsigned int y = baseclusters[j]->getLowestY();
            if (y < lowest_y)
                lowest_y = y;
        }
    }
    return lowest_y;
}

unsigned int Multiplicon::getHighestX() const {
    unsigned int highest_x = baseclusters[0]->getHighestX();
    for (unsigned int j = 1; j < baseclusters.size(); j++) {
        if (baseclusters[j]->getCountAnchorPoints() > 0) {
            unsigned int x = baseclusters[j]->getHighestX();
            if (x > highest_x)
                highest_x = x;
        }
    }
    return highest_x;
}

unsigned int Multiplicon::getHighestY() const {
    unsigned int highest_y = baseclusters[0]->getHighestY();
    for (unsigned int j = 1; j < baseclusters.size(); j++) {
        if (baseclusters[j]->getCountAnchorPoints() > 0) {
            unsigned int y = baseclusters[j]->getHighestY();
            if (y > highest_y)
                highest_y = y;
        }
    }
    return highest_y;
}

void Multiplicon::write () const {
    cout << "multiplicon (" << getLowestX() << ", " << getLowestY()
         << ") to (" << getHighestX() << ", " << getHighestY() << ")"
         << endl;
    cout << "begin_x " << begin_x << "\tbegin_y" <<  begin_y << "\tend_x "
         << end_x << "\tend_y" <<  end_y << endl;
    cout << "multiplicon baseclusters: " << endl;
    for (unsigned int k = 0; k < baseclusters.size(); k++) {
        baseclusters[k]->write();
    }
}

int Multiplicon::getPackSize() const
{
    int packSize = 0;

    packSize += sizeof(x_objectID) + sizeof(y_objectID) + sizeof(level) +
                sizeof(isRedundant) + sizeof(begin_x) +
                sizeof(end_x) + sizeof(begin_y) + sizeof(end_y);
    packSize += sizeof(int); // for storing the number of baseclusters
    vector<BaseCluster*>::const_iterator it = baseclusters.begin();
    for ( ; it != baseclusters.end(); it++)
        packSize += (*it)->getPackSize();

    return packSize;
}

int Multiplicon::pack(char *buffer) const
{
    const char *bufferOrig = buffer;

    memcpy(buffer, &x_objectID, sizeof(x_objectID));
    buffer += sizeof(x_objectID);
    memcpy(buffer, &y_objectID, sizeof(y_objectID));
    buffer += sizeof(y_objectID);
    memcpy(buffer, &level, sizeof(level));
    buffer += sizeof(level);
    memcpy(buffer, &isRedundant, sizeof(isRedundant));
    buffer += sizeof(isRedundant);
    memcpy(buffer, &begin_x, sizeof(begin_x));
    buffer += sizeof(begin_x);
    memcpy(buffer, &end_x, sizeof(end_x));
    buffer += sizeof(end_x);
    memcpy(buffer, &begin_y, sizeof(begin_y));
    buffer += sizeof(begin_y);
    memcpy(buffer, &end_y, sizeof(end_y));
    buffer += sizeof(end_y);

    int size = baseclusters.size();
    memcpy(buffer, &size, sizeof(size));
    buffer += sizeof(size);

    vector<BaseCluster*>::const_iterator it = baseclusters.begin();
    for ( ; it != baseclusters.end(); it++)
        buffer += (*it)->pack(buffer);

    return buffer - bufferOrig;
}

int Multiplicon::getPackSize(const vector<Multiplicon*> &mplicons)
{
    int packSize = 0;

    packSize += sizeof(int);    // for storing mplicons.size()
    vector<Multiplicon*>::const_iterator it = mplicons.begin();
    for ( ; it != mplicons.end(); it++)
        packSize += (*it)->getPackSize();

    return packSize;
}

void Multiplicon::packMultiplicons(const vector<Multiplicon*> &mplicons,
                                   char* buffer)
{
    int nMultiplicons = mplicons.size();
    memcpy(buffer, &nMultiplicons, sizeof(nMultiplicons));
    buffer += sizeof(nMultiplicons);

    vector<Multiplicon*>::const_iterator it = mplicons.begin();
    for ( ; it != mplicons.end(); it++)
        buffer += (*it)->pack(buffer);
}

void Multiplicon::unpackL2Multiplicons(const char *buffer,
                                       vector<Multiplicon*> &mplicons,
                                       const vector<GeneList *>& genelists,
                                       bool useFamily)
{
    int nMultiplicons = 0;
    memcpy(&nMultiplicons, buffer, sizeof(nMultiplicons));
    buffer += sizeof(nMultiplicons);

    mplicons.reserve(mplicons.size() + nMultiplicons);
    for (int i = 0; i < nMultiplicons; i++) {
        mplicons.push_back(new Multiplicon(buffer, genelists, useFamily));
        buffer += mplicons.back()->getPackSize();
    }
}

void Multiplicon::unpackHLMultiplicons(const char *buffer,
                                       vector<Multiplicon*> &mplicons,
                                       const Profile &xObject,
                                       const vector<GeneList *>& genelists,
                                       bool useFamily)
{
    int nMultiplicons = 0;
    memcpy(&nMultiplicons, buffer, sizeof(nMultiplicons));
    buffer += sizeof(nMultiplicons);

    mplicons.reserve(mplicons.size() + nMultiplicons);
    for (int i = 0; i < nMultiplicons; i++) {
        mplicons.push_back(new Multiplicon(buffer, xObject, genelists, useFamily));
        buffer += mplicons.back()->getPackSize();
    }
}

bool operator==(const Multiplicon &lhs, const Multiplicon &rhs)
{
    if (lhs.level != rhs.level) return false;
    if (lhs.isRedundant != rhs.isRedundant) return false;
    if (lhs.profile != rhs.profile) return false;
    // base class comparison
    if (lhs.begin_x != rhs.begin_x) return false;
    if (lhs.end_x != rhs.end_x) return false;
    if (lhs.begin_y != rhs.begin_y) return false;
    if (lhs.end_y != rhs.end_y) return false;
    if (lhs.x_objectID != rhs.x_objectID) return false;
    if (lhs.y_objectID != rhs.y_objectID) return false;

    if (lhs.baseclusters.size() != rhs.baseclusters.size()) return false;
    for (unsigned int i = 0; i < lhs.baseclusters.size(); i++)
        if (*lhs.baseclusters[i] != *rhs.baseclusters[i]) return false;

    return true;
}

void Multiplicon::extractXObject(const GeneList &xObject)
{
    assert (level == 2);

    // reserve one extra, because unaligned_y_list will be added
    // when aligning the structure
    xSegments.reserve(level);

    xSegments.push_back(new GeneList(xObject, getBeginX(), getEndX()));
}

void Multiplicon::extractXObject(const Profile &xObject)
{
    const vector<GeneList*>& profSeg = xObject.getSegments();
    assert (level == profSeg.size() + 1);

    // reserve one extra, because unaligned_y_list will be added
    // when aligning the structure
    xSegments.reserve(level);

    // extract the part that is relevant
    for (unsigned int i = 0; i < profSeg.size(); i++)
        xSegments.push_back(new GeneList(*profSeg[i], getBeginX(), getEndX()));

    // copy the homologs already stored in the profile's multiplicon
    homologs = xObject.getMultiplicon().homologs;
}

void Multiplicon::extractYObject(const GeneList &yObject)
{
    ySegment = new GeneList(yObject, getBeginY(), getEndY());
}

bool Multiplicon::createHomologs(bool useFamily, int minHomologs)
{
    // remove the homologs that are no longer valid (shortening of segment)
    pruneLinks();

    // add the extra homologs due to the y-segment
    addLinks(useFamily);

    // count the number of homologs per segment
    for (int i = 0; i < xSegments.size(); i++) {
        int count = 0;
        for (int j = 0; j < xSegments[i]->getSize(); j++) {
            if (xSegments[i]->getRemappedElements()[j]->hasHomolog)
                count++;
            if (count >= minHomologs)
                break;
        }

        if (count < minHomologs)
            return false;
    }

    for (int count = 0, j = 0; j < ySegment->getSize(); j++) {
        if (ySegment->getRemappedElements()[j]->hasHomolog)
            count++;
        if (count >= minHomologs)
            return true;
    }

    return false;
}

void Multiplicon::pruneLinks()
{
    // get out early (in case of level 2)
    if (homologs.empty())
        return;

    // create the listElements map (mapping between numID and ListElement)
    vector<map<int, ListElement*> > leMap(xSegments.size());
    for (int i = 0; i < xSegments.size(); i++) {
        VecListElementCIt e = xSegments[i]->getLEBegin();
        for ( ; e != xSegments[i]->getLEEnd(); e++) {
            ListElement &le = *(*e);
            le.hasHomolog = false;
            le.hasAP = false;
            leMap[i][le.getNumID()] = &le;
        }
    }

    // remove the homologs that are no longer applicable
    for (set<Link>::iterator e = homologs.begin(); e != homologs.end(); ) {
        map<int, ListElement*>::iterator itX, itY;
        int segX = e->segmentX;
        int segY = e->segmentY;
        itX = leMap[segX].find(e->geneXID);
        itY = leMap[segY].find(e->geneYID);
        if ((itX == leMap[segX].end()) || (itY == leMap[segY].end())) {
            homologs.erase(e++);
        } else {
            itX->second->hasHomolog = true;
            itY->second->hasHomolog = true;
            if (e->isAP) {
                itX->second->hasAP = true;
                itY->second->hasAP = true;
            }
            e++;
        }
    }
}

void Multiplicon::addLinks(bool useFamily)
{
    if (!useFamily) {
        // add the homologs between the xSegments and the ySegment
        VecListElementCIt itY = ySegment->getLEBegin();
        for ( ; itY < ySegment->getLEEnd(); itY++) {
            ListElement &eY = *(*itY);
            if (eY.isGap()) continue;
            if (!eY.getGene().hasPairs()) continue;

            for (int i = 0; i < xSegments.size(); i++) {
                VecListElementCIt itX = xSegments[i]->getLEBegin();
                for ( ; itX != xSegments[i]->getLEEnd(); itX++) {
                    ListElement &eX = *(*itX);
                    if (eX.isGap()) continue;
                    if (!eX.getGene().isPairWith(eY.getGene()))
                        continue;

                    // store the homolog gene IDs
                    homologs.insert(Link(i, xSegments.size(),
                                         eX.getNumID(), eY.getNumID()));
                    eX.hasHomolog = true;
                    eY.hasHomolog = true;
                }
            }
        }
    } else {    // we're using gene families
        map<string, set<pair<ListElement*, int> > > pairMap;

        // insert the elements of the profile
        for (int i = 0; i < xSegments.size(); i++) {
            for (int j = 0; j < xSegments[i]->getSize(); j++) {
                ListElement *le = xSegments[i]->getRemappedElements()[j];
                if (le->isGap()) continue;

                pairMap[le->getGene().getFamily()].insert(pair<ListElement*, int>(le, i));
            }
        }

        // insert the elements of the y-list
        for (int j = 0; j < ySegment->getSize(); j++) {
            ListElement *le = ySegment->getRemappedElements()[j];
            if (le->isGap()) continue;

            pairMap[le->getGene().getFamily()].insert(pair<ListElement*, int>(le, xSegments.size()));
        }

        // now extract the homologs
        for (int j = 0; j < ySegment->getSize(); j++) {
            ListElement *eY = ySegment->getRemappedElements()[j];
            if (eY->isGap()) continue;

            map<string, set<pair<ListElement*, int> > >::iterator it =
                pairMap.find(eY->getGene().getFamily());
            set<pair<ListElement*, int> > &group = it->second;

            set<pair<ListElement*, int> >::iterator e;
            for (e = group.begin(); e != group.end(); e++) {
                ListElement *eX = e->first;
                int segX = e->second;

                if (segX == xSegments.size()) continue;
                // to have compatibility with the pairwise case:
                if (eX->getNumID() == eY->getNumID()) continue;
                homologs.insert(Link(segX, xSegments.size(),
                                  eX->getNumID(), eY->getNumID()));

                eX->hasHomolog = true;
                eY->hasHomolog = true;
            }
        }
    }

    // mark the homologs that are AP
    int offX = getBeginX();
    int offY = getBeginY();
    vector<BaseCluster*>::const_iterator bc = getBaseClusters().begin();
    for ( ; bc != getBaseClusters().end(); bc++) {
        multiset<AnchorPoint>::const_iterator e = (*bc)->getAPBegin();
        for ( ; e != (*bc)->getAPEnd(); e++) {
            // get the y-element
            ListElement &eY = ySegment->getLe(e->getY()-offY);
            assert(!eY.isGap());
            assert (eY.getNumID() == e->getGeneYID());

            // find the corresponding x-element
            bool found = false;
            for (int i = 0; i < xSegments.size(); i++) {
                ListElement &eX = xSegments[i]->getLe(e->getX()-offX);
                if (eX.isGap()) continue;
                if (eX.getNumID() != e->getGeneXID()) continue;
                assert (eX.getGene().isPairWith(eY.getGene()));

                // mark the elements
                eX.hasAP = true;
                eY.hasAP = true;

                // store the AP gene IDs
                Link target(i, xSegments.size(),
                            e->getGeneXID(), e->getGeneYID());
                set<Link>::iterator lnk = homologs.find(target);
                assert(lnk != homologs.end());
                const_cast<Link&>(*lnk).isAP = true;

                found = true;
                break;
            }
            assert(found);
        }
    }
}

void Multiplicon::createProfile(int profileID)
{
    profile = new Profile(*this, profileID);
}

void Multiplicon::align(const AlignmentMethod &alignMethod, int maxGaps)
{
    assert(profile != NULL);
    profile->align(alignMethod, maxGaps);
}

void Multiplicon::compareAligners(AlignScore &NW, AlignScore &GG,
                                  AlignScore &RA, AlignScore &RC,
                                  AlignScore &RAC, AlignScore &LL,
                                  AlignScore &LLBS, AlignScore &LS)
{
    assert(profile != NULL);
    profile->compareAligners(NW, GG, RA, RC, RAC, LL, LLBS, LS);
}

void Multiplicon::checkAlignment(int minHomologs)
{
    assert(profile != NULL);
    profile->checkAlignment(minHomologs);
}
