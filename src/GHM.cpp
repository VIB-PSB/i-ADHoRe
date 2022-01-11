#include "GHM.h"

#include "AnchorPoint.h"
#include "BaseCluster.h"
#include "SynthenicCloud.h"
#include "Multiplicon.h"
#include "ListElement.h"

#include "util.h"
#include <cassert>

using namespace std;

GHM::GHM (const GeneList& xObject, const GeneList& yObject, bool cloudsearch) :
        x_object(xObject), y_object(yObject), level(2), area(0), isCloudSearch(cloudsearch)
{
    identical = false;

    if (!x_object.isSegment() &&
            x_object.getListName() == y_object.getListName() &&
            x_object.getGenomeName() == y_object.getGenomeName())
        identical = true;
}

GHM::~GHM()
{
    for (int i = 0; i < 2; i++) {
        vector<BaseCluster*>::iterator it = filteredBC[i].begin();
        for ( ; it != filteredBC[i].end(); it++)
            delete *it;
    }
}

void GHM::buildMatrixFast()
{
    const vector<ListElement*>& xList = x_object.getRemappedElements();
    const vector<ListElement*>& yList = y_object.getRemappedElements();

    //build family tree
    map<string,vector<pair<ListElement*,uint> > > geneFamilyTree;

    for (uint y = 0; y < yList.size(); y++) {
        ListElement &yElement = *yList[y];
        if (yElement.isGap()) continue;

        Gene &geneY = yElement.getGene();
        if (!geneY.hasPairs()) continue;

        geneFamilyTree[geneY.getFamily()].push_back(make_pair<>(&yElement,y));
    }

    for (uint x = 0; x < xList.size(); x++) {
        const ListElement &xElement = *xList[x];
        if (xElement.isGap()) continue;

        const Gene &geneX = xElement.getGene();
        if (!geneX.hasPairs()) continue;

        vector<pair<ListElement*,uint> > &geneYVector=geneFamilyTree[geneX.getFamily()];

        for (uint i=0; i<geneFamilyTree[geneX.getFamily()].size(); i++) {
            int y = geneYVector[i].second;

            if (identical && (y > x)) continue;

            const ListElement &yElement = *(geneYVector[i].first);

            if (geneX.getID() == yElement.getGene().getID())
                continue; //identical genes are not considered valid pairs

            if (isCloudSearch) {
                matrix[MIXED_ORIENT][x].insert(y);
                count_points[MIXED_ORIENT]++;
            } else {
                if (xElement.getOrientation() == yElement.getOrientation()) {
                    matrix[SAME_ORIENT][x].insert(y);
                    count_points[SAME_ORIENT]++;
                } else {
                    matrix[OPP_ORIENT][x].insert(y);
                    count_points[OPP_ORIENT]++;
                }
            }
        }
    }
}

void GHM::buildMatrix (bool useFamilies)
{
    // reset the number of points in the GHM
    count_points[OPP_ORIENT] = 0;
    count_points[SAME_ORIENT] = 0;
    count_points[MIXED_ORIENT]= 0;

    int numMatrices= isCloudSearch ? 1 : 2;

    // reset the GHM datastructure
    matrix.clear();
    matrix.resize(numMatrices);

    if (useFamilies) {
        buildMatrixFast();
        return;
    }

    const vector<ListElement*>& xList = x_object.getRemappedElements();
    const vector<ListElement*>& yList = y_object.getRemappedElements();

    for (unsigned int x = 0; x < xList.size(); x++) {
        const ListElement &xElement = *xList[x];

        if (xElement.isGap()) continue;
        if (!xElement.getGene().hasPairs()) continue;

        std::queue<int> queue;
        xElement.matchingPositions(yList, queue);

        while (!queue.empty()) {
            // get and remove first element
            unsigned int y = queue.front();
            queue.pop();

            // in case the xList == yList,
            // only store the upper triangular part
            if (identical && (y > x)) continue;

            const ListElement &yElement = *yList[y];

            if (isCloudSearch)
            {
                matrix[MIXED_ORIENT][x].insert(y);
                count_points[MIXED_ORIENT]++;
            }
            else {
                if (xElement.getOrientation() == yElement.getOrientation())
                {
                    matrix[SAME_ORIENT][x].insert(y);
                    count_points[SAME_ORIENT]++;
                }
                else
                {
                    matrix[OPP_ORIENT][x].insert(y);
                    count_points[OPP_ORIENT]++;
                }
            }

        }
    }
}

void GHM::run(const Settings& settings)
{
    if (isCloudSearch) runSyntheny(settings);
    else runCollinear(settings);
}


void GHM::runCollinear(const Settings& settings)
{
    if ((count_points[0]+count_points[1]) < (unsigned int)settings.getAnchorPoints()) return;

    prepareForStatisticalValidation();
    int gapsizes [10];
    settings.getGapSizes(&gapsizes[0]);

    for (int i = 0; i < 10 && gapsizes[i] > 0; i++) {

        // seed, enrich and join for positive orientation
        seedBaseClusters(gapsizes[i], 0, settings.getQValue());
        enrichClusters(gapsizes[i], 0, 0, settings.getQValue());
        joinClusters(gapsizes[i], 0, settings.getQValue());

        // seed, enrich and join for negative orientation
        seedBaseClusters(gapsizes[i], 1, settings.getQValue());
        enrichClusters(gapsizes[i], 1, 1, settings.getQValue());
        joinClusters(gapsizes[i], 1, settings.getQValue());

        /*NOTE debug code: visulalize GHM in every iteration
        string output_path=settings.getOutputPath();
        char temp[3]; //
        sprintf(temp,"%d",i);
        string path=string(output_path)+string(temp);
        visualizeBaseClustersPNG(path.c_str());
        */
    }

    enrichClusters(settings.getGapSize(), 0, 0, settings.getQValue());
    enrichClusters(settings.getGapSize(), 1, 1, settings.getQValue());
    enrichClusters(settings.getGapSize(), 1, 0, settings.getQValue());
    enrichClusters(settings.getGapSize(), 0, 1, settings.getQValue());

    switch (settings.getMultHypCorMethod()) {
    case None:
        filterBaseClusters(settings.getProbCutoff());
        break;
    case Bonferroni :
        filterBaseClustersBF(settings.getProbCutoff());
        break;
    case FDR:
        filterBaseClustersFDR(settings.getProbCutoff());
        break;
    }

    //NOTE filteredBC will be used to plot GHM => Don't erase!
    /*for (int o = 0; o < 2; o++) {
        for (int i = 0; i < filteredBC[o].size(); i++)
            delete filteredBC[o][i];
        filteredBC[o].clear();
    }*/

    enrichClusters(settings.getGapSize(), 0, 0, settings.getQValue());
    enrichClusters(settings.getGapSize(), 1, 1, settings.getQValue());
    enrichClusters(settings.getGapSize(), 1, 0, settings.getQValue());
    enrichClusters(settings.getGapSize(), 0, 1, settings.getQValue());

    clusterClusters(settings.getClusterGap(), settings.getQValue(),
                    settings.getAnchorPoints());

    setMultiplicons(settings.useFamily(), settings.getAnchorPoints());
}

void GHM::prepareForStatisticalValidation() {
    //calculate the area of the ghm
    if (identical) {
        area = x_object.getSize() * (x_object.getSize() - 1) / 2;
    } else {
        //get the number of masked rows
        int masked_rows;
        if (level > 2) {
            masked_rows = y_object.getNumberOfMaskedElements();
        }
        else {
            masked_rows = 0;
        }
        area = x_object.getSize() * (y_object.getSize() - masked_rows);
    }
}

void GHM::seedBaseCluster(BaseCluster *basecluster, int gap) const
{
    assert(basecluster->getCountAnchorPoints() == 1);
    bool orientation = basecluster->getOrientation();
    const map<int, set<int> > &mat = matrix[orientation];

    int refX = basecluster->getAPBegin()->getX();
    int refY = basecluster->getAPBegin()->getY();

    while (true) {
        int loX = refX + 1;
        int hiX = refX + gap;
        int loY = (orientation) ? refY + 1 : refY - gap;
        int hiY = (orientation) ? refY + gap : refY - 1;

        map<int, set<int> >::const_iterator itX = mat.lower_bound(loX);
        map<int, set<int> >::const_iterator endX = mat.upper_bound(hiX);

        int closestX = 0, closestY = 0;
        int closestDpd = gap + 1;

        for ( ; itX != endX; itX++) {
            int x = itX->first;
            const set<int>& setY = itX->second;

            set<int>::const_iterator itY = setY.lower_bound(loY);
            set<int>::const_iterator endY = setY.upper_bound(hiY);

            for ( ; itY != endY; itY++) {
                int y = *itY;

                int dpd = BaseCluster::dpd(refX, refY, x, y);
                if (dpd < closestDpd) {
                    closestDpd = dpd;
                    closestX = x;
                    closestY = y;
                }
            }
        }

        if (closestDpd == gap + 1) return;

        basecluster->addAnchorPoint(closestX, closestY);
        basecluster->addBackBone(closestX, closestY);
        refX = closestX;
        refY = closestY;
    }
}

void GHM::seedBaseClusters(int gap, bool orientation, double qValue)
{
    map<int, set<int> > &mat = matrix[orientation];
    map<int, set<int> >::iterator itX;
    set<int>::const_iterator itY;
    multiset<AnchorPoint>::const_iterator AP;

    for (itX = mat.begin(); itX != mat.end(); ) {
        for (itY = itX->second.begin(); itY != itX->second.end(); ) {
            BaseCluster* basecluster = new BaseCluster(orientation);
            basecluster->addAnchorPoint(itX->first, *itY);
            basecluster->addBackBone(itX->first, *itY);
            seedBaseCluster(basecluster, gap);

            if (basecluster->getCountAnchorPoints() > 2 &&
                    basecluster->r_squared() >= qValue) {
                baseclusters[orientation].push_back(basecluster);
                basecluster->updateStatistics();

                // delete found APs from matrix except the first
                AP = basecluster->getAPBegin();
                for (AP++; AP != basecluster->getAPEnd(); AP++)
                    mat[AP->getX()].erase(AP->getY());
                // delete the first AP
                itX->second.erase(itY++);
            } else {
                delete basecluster;
                itY++;
            }
        }

        // it is possible that all AP for this x were deleted
        if (itX->second.empty())
            mat.erase(itX++);
        else
            itX++;
    }
}

void GHM::addCluster(Multiplicon& cluster)
{
    multiplicons.push_back(&cluster);
}

void GHM::enrichClusters(int gap, bool clusterOrientation,
                         bool ghmOrientation, double qValue)
{
    //make copy of the baseclusters vector
    vector<BaseCluster*> clusters (baseclusters[clusterOrientation]);
    map<int, set<int> > &mat = matrix[ghmOrientation];

    map<int, set<int> >::iterator itX;
    set<int>::iterator itY;

    while (!clusters.empty()) {
        vector<BaseCluster*> clustersNextIteration;
        clustersNextIteration.reserve(clusters.size());
        vector<bool> changedClusters (clusters.size(), false);

        for (itX = mat.begin(); itX != mat.end(); ) {
            for (itY = itX->second.begin(); itY != itX->second.end(); ) {
                double closestDistance = gap + 1;

                int closestClusterIndex = -1;
                vector<BaseCluster*>::iterator it = clusters.begin();
                for (int index = 0; it != clusters.end(); it++, index++) {
                    (*it)->updateStatistics();

                    double distance = (*it)->distanceToPoint(itX->first,*itY);
                    if (distance >= closestDistance) continue;

                    if ((*it)->r_squared(itX->first,*itY) >= qValue &&
                            (*it)->inInterval(itX->first,*itY))
                        // && distance <= 2*(*it)->averageDPD())
                    {
                        closestDistance = distance;
                        closestClusterIndex = index;
                    }
                }

                if (closestClusterIndex != -1) {
                    BaseCluster *closestCluster = clusters[closestClusterIndex];
                    closestCluster->addAnchorPoint(itX->first, *itY);
                    closestCluster->updateStatistics();

                    if (!changedClusters[closestClusterIndex]) {
                        clustersNextIteration.push_back(closestCluster);
                        changedClusters[closestClusterIndex] = true;
                    }

                    itX->second.erase(itY++);
                } else {
                    itY++;
                }
            }

            if (itX->second.empty())
                mat.erase(itX++);
            else
                itX++;
        }
        clusters = clustersNextIteration;
    }
}

void GHM::joinClusters(int gap, bool orient, double qValue)
{
    vector<BaseCluster*>::iterator i, j;
    while (getClosestClusters(gap, orient, orient, i, j, qValue)) {
        // merge the two clusters
        (*i)->mergeWith(**j);

        // remove the old basecluster
        delete *j;
        baseclusters[orient].erase(j);

    }
}

bool GHM::getClosestClusters(int gap, bool orientI, bool orientJ,
                             vector<BaseCluster*>::iterator& closestI,
                             vector<BaseCluster*>::iterator& closestJ,
                             double qValue)
{
    vector<BaseCluster*> &bcI = baseclusters[orientI];
    vector<BaseCluster*> &bcJ = baseclusters[orientJ];
    vector<BaseCluster*>::iterator itI, itJ;

    double closestDistance = gap + 1;

    for (itI = bcI.begin(); itI != bcI.end(); itI++) {
        itJ = (orientI == orientJ) ? itI + 1 : bcJ.begin();
        for ( ; itJ != bcJ.end(); itJ++) {
            BaseCluster* cI = *itI;
            BaseCluster* cJ = *itJ;

            if (cI->getMultiplicon() != NULL &&
                    cJ->getMultiplicon() != NULL &&
                    cI->getMultiplicon() == cJ->getMultiplicon())
                continue;

            double distance = cI->distanceToCluster(*cJ);

            if (distance < closestDistance) { //check distance

                if (cI->partialOverlappingInterval(*cJ)) { //check if one cluster resides in the overlappinginterval of the other
                    if (cI->r_squared(*cJ) >= qValue) { //check quality > qvalue for combined cluster
                        closestI = itI;
                        closestJ = itJ;
                        closestDistance = distance;
                    }
                }
            }
        }
    }

    return (closestDistance < gap + 1);
}

void GHM::filterBaseClusters(double probCutoff)
{
    for (int orient = 0; orient < 2; orient++) {
        vector<BaseCluster*>::iterator it = baseclusters[orient].begin();
        while (it != baseclusters[orient].end()) {
            double pGlobal = (*it)->calculateProbability(area, count_points[orient], level);
            (*it)->setRandomProbability(pGlobal);

            // if the cluster was generated by chance
            if (pGlobal > probCutoff) {
                multiset<AnchorPoint>::iterator AP;
                AP = (*it)->getAPBegin();
                for ( ; AP != (*it)->getAPEnd(); AP++)
                    matrix[orient][AP->getX()].insert(AP->getY());

                filteredBC[orient].push_back(*it);
                it = baseclusters[orient].erase(it);

            } else {
                it++;
            }
        }
    }
}

void GHM::filterCloudsBinomialDBF(double probCutoff)
{
    filterCloudsBinomialD(probCutoff/count_points[MIXED_ORIENT]);
}

void GHM::filterCloudsBinomialDCorrBF(double probCutoff)
{
    filterCloudsBinomialDCorr(probCutoff/count_points[MIXED_ORIENT]);
}

void GHM::filterBaseClustersBF(double probCutoff)
{
    for (int orient = 0; orient < 2; orient++) {
        vector<BaseCluster*>::iterator it = baseclusters[orient].begin();
        while (it != baseclusters[orient].end()) {
            double pGlobal = (*it)->calculateProbability(area, count_points[orient], level);
            (*it)->setRandomProbability(pGlobal);

            // if the cluster was generated by chance
            if (pGlobal > (probCutoff / count_points[orient])) {
                multiset<AnchorPoint>::iterator AP;
                AP = (*it)->getAPBegin();
                for ( ; AP != (*it)->getAPEnd(); AP++)
                    matrix[orient][AP->getX()].insert(AP->getY());

                filteredBC[orient].push_back(*it);
                it = baseclusters[orient].erase(it);
            } else {
                it++;
            }
        }
    }
}

bool pValueSort(const BaseCluster *left, const BaseCluster *right)
{
    return left->getRandomProbability() < right->getRandomProbability();
}

bool pValueSortClouds(const SynthenicCloud* left, const SynthenicCloud* right) {


    return left->getRandomProbability() < right->getRandomProbability();
}

void GHM::filterBaseClustersFDR(double probCutoff)
{
    for (int orient = 0; orient < 2; orient++) {
        vector<BaseCluster*>::iterator it = baseclusters[orient].begin();
        while (it != baseclusters[orient].end()) {
            double pGlobal = (*it)->calculateProbability(area, count_points[orient], level);
            (*it)->setRandomProbability(pGlobal);
            it++;
        }

        // sort the basecluster according to their pValue
        sort(baseclusters[orient].begin(), baseclusters[orient].end(), pValueSort);

        int k;
        for (k = 0; k < baseclusters[orient].size(); k++) {
            if (baseclusters[orient][k]->getRandomProbability() > (k+1)*probCutoff/count_points[orient])
                break;
        }

        // move the filtered baseclusters to filteredBC
        filteredBC[orient] = vector<BaseCluster*>(baseclusters[orient].begin() + k, baseclusters[orient].end());
        // and erase them from baseclusters
        baseclusters[orient].erase(baseclusters[orient].begin() + k, baseclusters[orient].end());
    }
}

void GHM::clusterClusters(int clusterGap, double qValue,
                          unsigned int cntAnchorpoints)
{
    bool join_flag = true;
    while (join_flag) {
        join_flag = false;
        for (int k = 0; k < 4; k++) {
            bool or_x = (k < 2) ? true : false;
            bool or_y = (k % 2 == 0) ? true : false;

            //twist the clusters
            if (or_x != or_y)
                twistClusters(or_y);
            bool was_twisted = (or_x != or_y);

            vector<BaseCluster*>::iterator itI, it_j;
            bool found = getClosestClusters(clusterGap, or_x, or_y, itI, it_j, qValue);

            //untwist the clusters
            if (or_x != or_y)
                twistClusters(or_y);

            if (found && (!was_twisted || !(*itI)->wasTwisted())) {
                BaseCluster* cI = *itI;
                BaseCluster* cJ = *it_j;
                cJ->wasTwisted(was_twisted);
                join_flag = true;

                if (cI->getMultiplicon() != NULL &&
                        cJ->getMultiplicon() != NULL) {
                    Multiplicon* tmp = cJ->getMultiplicon();
                    cI->getMultiplicon()->addBaseClusters(*cJ->getMultiplicon());
                    //cJ->resetMultiplicon();
                    tmp->clear();
                }
                else if (cI->getMultiplicon() != NULL ||
                         cJ->getMultiplicon() != NULL) {
                    if (cI->getMultiplicon() == NULL) {
                        BaseCluster* tmp = cI;
                        cI = cJ;
                        cJ = tmp;
                    }
                    cI->getMultiplicon()->addBaseCluster(*cJ);
                }
                else {
                    Multiplicon* multiplicon =
                        new Multiplicon(x_object.getID(), y_object.getID(), level);
                    multiplicon->addBaseCluster(*cI);
                    multiplicon->addBaseCluster(*cJ);
                    addCluster(*multiplicon);
                }
            }
        }
    }

    //Check the configuration of any obtained metaclusters
    //checkMetaclusterConfiguration();

    //Assign all baseclusters that were not assigned to a metacluster to their own multiplicon
    for (unsigned int k = 0; k < 2; k++) {
        for (unsigned int i = 0; i < baseclusters[k].size(); i++) {
            if (baseclusters[k][i]->getMultiplicon() == NULL) {
                Multiplicon* multiplicon =
                    new Multiplicon(x_object.getID(), y_object.getID(), level);
                baseclusters[k][i]->wasTwisted(false);
                multiplicon->addBaseCluster(*baseclusters[k][i]);
                addCluster(*multiplicon);
            }
        }
    }

    //Remove empty multiplicons from the ghm
    vector<Multiplicon*>::iterator it = multiplicons.begin();
    while (it != multiplicons.end()) {
        if ((*it)->getBaseClusters().size() == 0 ||
                (*it)->getCountAnchorPoints() < cntAnchorpoints) {
            delete *it;
            it = multiplicons.erase(it);
        }
        else {
            it++;
        }
    }
}

void GHM::twistClusters(bool orientation)
{
    vector<BaseCluster*> &bc = baseclusters[orientation];
    for (vector<BaseCluster*>::iterator i = bc.begin(); i != bc.end(); i++)
        (*i)->twistCluster();
}

/*void GHM::checkMetaclusterConfiguration()
{
    vector<Multiplicon*>::iterator it = multiplicons.begin();
    while (it != multiplicons.end()) {
        bool configuration_ok = true;

        const vector<BaseCluster*>& baseclusters = (*it)->getBaseClusters();
        for (int i = 0; i < (int)baseclusters.size() - 1; i++) {
            for (int j = i + 1; j < (int)baseclusters.size(); j++) {
                BaseCluster &bi = *baseclusters[i];
                BaseCluster &bj = *baseclusters[j];
                if (bi.getOrientation() == bj.getOrientation()) {
                    if (bi.overlappingCoordinates(bj) != 0 &&
                            !bi.partialOverlappingInterval(bj)) {
                        configuration_ok = false;
                        cout << "configuration not ok!" << endl;
                    }
                }
                if (!configuration_ok) break;
            }
            if (!configuration_ok) break;
        }
        if (!configuration_ok) {
            for (unsigned int i = 0; i < baseclusters.size(); i++)
                baseclusters[i]->resetMultiplicon();
            // we do not want the baseclusters in
            // the multiplicon to be deleted
            (*it)->clear();
            delete *it;
            it = multiplicons.erase(it);
        }
        else {
            it++;
        }
    }
}*/

void GHM::setMultiplicons(bool useFamily,
                          int minHomologs)
{
    const vector<ListElement*>& x_elements = x_object.getRemappedElements();
    const vector<ListElement*>& y_elements = y_object.getRemappedElements();

    for (unsigned int i = 0; i < multiplicons.size(); i++) {
        multiplicons[i]->setBounds();

        const vector<BaseCluster*>& baseclusters = multiplicons[i]->getBaseClusters();
        for (unsigned int j = 0; j < baseclusters.size(); j++) {
            baseclusters[j]->setBounds();

            multiset<AnchorPoint>::const_iterator e = baseclusters[j]->getAPBegin();
            for ( ; e != baseclusters[j]->getAPEnd(); e++) {
                int x = e->getX();
                int y = e->getY();

                // this const_cast is allowed, because setting the gene pointers
                // does not invalidate the multiset iterator
                const_cast<AnchorPoint&>(*e).setGeneIDs(x_elements[x]->getNumID(),
                                                        y_elements[y]->getNumID());
            }
        }
        multiplicons[i]->extractXObject(x_object);
        multiplicons[i]->extractYObject(y_object);

        // for level-2 multiplicons, we ignore the return value of createHomologs
        multiplicons[i]->createHomologs(useFamily, minHomologs);
    }
}

void GHM::getMultiplicons(vector<Multiplicon*>& mps) const
{
    for (unsigned int i = 0; i < multiplicons.size(); i++) {
        mps.push_back(multiplicons[i]);
    }
}

bool GHM::identity(int x, int y)
{
    bool found = false;
    int orient = 0;
    while (orient < 2 && !found) {
        if (matrix[orient].find(x) != matrix[orient].end()) {
            if (matrix[orient][x].find(y) != matrix[orient][x].end()) {
                found = true;
            }
        }
        orient++;
    }

    return found;
}

/*unsigned int GHM::getNumberOfYHits()
{
    unsigned int result = 0;

    for (int orient = 0; orient < 2; orient++) {
        map<int, set<int> >::iterator itX;
        for (itX = matrix[orient].begin(); itX != matrix[orient].end(); itX++) {
            result += itX->second.size();
        }
    }

    return result;
}*/

/*void GHM::write()
{
    map<int, set<int> >::iterator itX;
    for (itX = matrix[0].begin(); itX != matrix[0].end(); itX++) {
        set<int>::iterator itY;
        for (itY = itX->second.begin(); itY != itX->second.end(); itY++) {
            cerr << "or: " << false << "\tx: " << itX->first << "\ty: " << *itY <<endl;
        }
    }
}*/

/*****************************************************************************/
/*SYNTHENIC CLOUD FUNCTIONS                                                  */
/*****************************************************************************/

void GHM::runSyntheny(const Settings& settings)
{
    if (count_points[MIXED_ORIENT] < (unsigned int)settings.getAnchorPoints()) return;

    int gapsizes [10];
    settings.getCloudGapSizes(&gapsizes[0]);
    int clusterGap=settings.getCloudClusterGap();

    for (int i = 0; i < 10 && gapsizes[i] > 0; i++) {

        condenseClouds(gapsizes[i],settings.isBruteforce());
        inflateClouds(gapsizes[i],settings.isBruteforce());
        mergeClouds(clusterGap,settings.isBruteforce());
    }

    switch (settings.getCloudFilterMethod()) {

    case Binomial:

        switch (settings.getMultHypCorMethod()) {
        case None:
            filterCloudsBinomialD(settings.getProbCutoff());
            break;
        case Bonferroni :
            filterCloudsBinomialDBF(settings.getProbCutoff());
            break;
        case FDR:
            filterCloudsBinomialDFDR(settings.getProbCutoff());
            break;
        }
        break;

    case BinomialCorr:

        switch (settings.getMultHypCorMethod()) {
        case None:
            filterCloudsBinomialDCorr(settings.getProbCutoff());
            break;
        case Bonferroni :
            filterCloudsBinomialDCorrBF(settings.getProbCutoff());
            break;
        case FDR:
            filterCloudsBinomialDCorrFDR(settings.getProbCutoff());
            break;
        }
        break;
    }

}

void GHM::condenseClouds(uint gap, bool bf)
{
    map<int, set<int> > &mat = matrix[MIXED_ORIENT];
    map<int, set<int> >::iterator itX;
    set<int>::const_iterator itY;

    vector<AnchorPoint> APRecycleBin; //APs to be removed from GHM

    for (itX = mat.begin(); itX != mat.end(); ) {
        for (itY = itX->second.begin(); itY != itX->second.end(); ) {

            SynthenicCloud* sCloud = new SynthenicCloud();
            sCloud->addAnchorPoint(itX->first, *itY);

            condenseCloud(*sCloud, gap, APRecycleBin,bf); //note that first AP will not be added to APRecycleBin!


            if (sCloud->getCountAnchorPoints() > 2) {
                sClouds.push_back(sCloud);

                removeAddedAnchorPoints(APRecycleBin);

                // delete the first AP
                itX->second.erase(itY++);

                if (bf) { //scan bounding box for extra dots
                    addAPFromSearchBoxBF(gap,sCloud->getBeginX(),sCloud->getEndX()
                        ,sCloud->getBeginY(),sCloud->getEndY(),*sCloud,APRecycleBin);
                    removeAddedAnchorPoints(APRecycleBin);
                }

            } else {

                delete sCloud;
                APRecycleBin.clear();
                itY++;
            }
        }
        // it is possible that all AP for this x were deleted
        if (itX->second.empty()) mat.erase(itX++);
        else itX++;
    }
}

void GHM::condenseCloud(SynthenicCloud& sCloud, uint gap, vector<AnchorPoint>& foundNewAP, bool bf)
{
    map<int, set<int> > &mat = matrix[MIXED_ORIENT];

    int numberOfAP;

    //search for new AP in a frame excluding the bounding box of the present APs
    do {
        numberOfAP=sCloud.getCountAnchorPoints();

        //determine boundaries of current bounding box
        int xmin=sCloud.getBeginX();
        int xmax=sCloud.getEndX();
        int ymin=sCloud.getBeginY();
        int ymax=sCloud.getEndY();

        if (!bf) {
            //search above
            addAPFromSearchBox(xmin-gap,xmax+gap,ymax+1,ymax+gap,sCloud,foundNewAP);
            //search below
            addAPFromSearchBox(xmin-gap,xmax+gap,ymin-gap,ymin-1,sCloud,foundNewAP);
            //search left
            addAPFromSearchBox(xmin-gap,xmin-1,ymin,ymax,sCloud,foundNewAP);
            //search right
            addAPFromSearchBox(xmax+1,xmax+gap,ymin,ymax,sCloud,foundNewAP);
        } else {
            //search above
            addAPFromSearchBoxBF(gap,xmin-gap,xmax+gap,ymax+1,ymax+gap,sCloud,foundNewAP);
            //search below
            addAPFromSearchBoxBF(gap,xmin-gap,xmax+gap,ymin-gap,ymin-1,sCloud,foundNewAP);
            //search left
            addAPFromSearchBoxBF(gap,xmin-gap,xmin-1,ymin,ymax,sCloud,foundNewAP);
            //search right
            addAPFromSearchBoxBF(gap,xmax+1,xmax+gap,ymin,ymax,sCloud,foundNewAP);
        }

    } while (sCloud.getCountAnchorPoints()-numberOfAP > 0); //while new points are found
}

void GHM::removeAddedAnchorPoints(vector<AnchorPoint>& APRecycleBin)
{
    map<int, set<int> > &mat = matrix[MIXED_ORIENT];

    for (int i=0; i<APRecycleBin.size(); i++) {
        AnchorPoint& AP=APRecycleBin.at(i);
        mat[AP.getX()].erase(AP.getY()); /*NOTE (is actually =="find and remove", so no risk at segmentation faults due to
        multiple removements!*/
    }
    APRecycleBin.clear();
}

void  GHM::addAPFromSearchBox(int loX, int hiX, int loY, int hiY, SynthenicCloud& sCloud, vector<AnchorPoint>& foundAP)

{
    map<int, set<int> > &mat = matrix[MIXED_ORIENT];
    map<int, set<int> >::const_iterator itX = mat.lower_bound(loX);

    int coX, coY; //coordinates of AP found in search box

    for ( ; itX != mat.upper_bound(hiX); itX++) {
        const set<int>& setY = itX->second;
        set<int>::const_iterator itY = setY.lower_bound(loY);

        for ( ; itY !=setY.upper_bound(hiY); itY++) {
            coX = itX->first;
            coY = *itY;

            sCloud.addAnchorPoint(coX,coY);
            foundAP.push_back(AnchorPoint(coX,coY,true));
        }
    }
}

void GHM::addAPFromSearchBoxBF(int gap, int loX, int hiX, int loY, int hiY, SynthenicCloud& sCloud, vector< AnchorPoint >& foundAP)
{
    map<int, set<int> > &mat = matrix[MIXED_ORIENT];
    map<int, set<int> >::const_iterator itX = mat.lower_bound(loX);

    int coX, coY; //coordinates of AP found in search box

    for ( ; itX != mat.upper_bound(hiX); itX++) {
        const set<int>& setY = itX->second;
        set<int>::const_iterator itY = setY.lower_bound(loY);

        for ( ; itY !=setY.upper_bound(hiY); itY++) {
            coX = itX->first;
            coY = *itY;

            //calculate distance to cloudIt
            int dist=sCloud.distanceToCloud(coX,coY);

            if (dist<=gap) {
                sCloud.addAnchorPoint(coX,coY);
                foundAP.push_back(AnchorPoint(coX,coY,true));
            }
        }
    }
}


void GHM::inflateClouds(uint gap, bool bf)
{
    vector<AnchorPoint> APRecycleBin; //APs to be removed from GHM
    list<SynthenicCloud*>::const_iterator it=sClouds.begin();

    for (; it!=sClouds.end(); it++) {

        condenseCloud(*(*it), gap,APRecycleBin,bf);
        removeAddedAnchorPoints(APRecycleBin);
    }
}

void GHM::mergeClouds(uint clustergap, bool bf)
{
    int distEst; //estimate of distance between clouds
    int sV1,sH1; //box side horizontal and vertical of box 1 and 2
    int sV2,sH2;
    int max1;  //maximum side of box 1
    int max2;  //maximum side of box 2
    int max3;  //maximum side of the two boxes

    int dist;

    list<SynthenicCloud*>::iterator it1=sClouds.begin();

    for (; it1!=sClouds.end(); it1++) {
        SynthenicCloud& sCloud1=*(*it1);

        sH1=sCloud1.calculateBoxWidth();
        sV1=sCloud1.calculateBoxHeight();

        max1=max(sH1,sV1);

        list<SynthenicCloud*>::iterator it2=sClouds.begin();
        for (; it2!=sClouds.end(); ) {

            if (it1==it2) { //don't compare cloud with itself
                it2++;
                continue;
            }
            SynthenicCloud& sCloud2=*(*it2);

            if (boxOverlap(sCloud1,sCloud2)) {

                if (cloudsOverlap(sCloud1,sCloud2)) {
                    joinClouds(it1,it2,bf);
                    it2=sClouds.begin(); //let it2 start over to allow further merging
                } else {

                    //FIXME this should be calcBruteForce -> since the boxes overlap closest points
                    //are not necessarily in the outer frames!!!!
                    dist=calculateMinimalKspdBetweenCloudsBruteForce(sCloud1,sCloud2/*,clustergap*/);

                    if (dist<clustergap)
                    {
                        joinClouds(it1,it2,bf);
                        it2=sClouds.begin(); //let it2 start over to allow further merging
                    } else
                        it2++;
                }
            } else {
                distEst=estimateCloudKspd(sCloud1,sCloud2);

                sH2=sCloud2.calculateBoxWidth();
                sV2=sCloud2.calculateBoxHeight();

                max2=max(sH2,sV2);
                max3=max(max1,max2);

                //check if clouds could be close
                if (distEst< max((int)clustergap,max3/2)) {
                    dist=calculateMinimalKspdBetweenClouds(sCloud1,sCloud2,clustergap);
                    if (dist<clustergap) {

                        joinClouds(it1,it2,bf);
                        it2=sClouds.begin(); //let j start over to allow further merging
                    } else
                        it2++;
                } else
                    it2++;
            }
        }
    }
}

void GHM::mergeCloudsBruteForce(uint gap, bool bf)
{
    int dist;

    list<SynthenicCloud*>::iterator it1=sClouds.begin();

    for (; it1!=sClouds.end(); it1++) {

        list<SynthenicCloud*>::iterator it2=sClouds.begin();
        SynthenicCloud& sCloud1=*(*it1);

        for (; it2!=sClouds.end(); ) {

            if (it1==it2) { //don't compare cloud with itself
                it2++;
                continue;
            }

            SynthenicCloud& sCloud2=*(*it2);

            dist=calculateMinimalKspdBetweenCloudsBruteForce(sCloud1, sCloud2);

            if (dist<gap) {
                joinClouds(it1,it2,bf);
                it2=sClouds.begin(); //restart
            } else
                it2++;

        }
    }
}

uint GHM::estimateCloudKspd(SynthenicCloud& sCloud1, SynthenicCloud& sCloud2) const
{
    assert(sCloud1.getCountAnchorPoints()>0);
    assert(sCloud2.getCountAnchorPoints()>0);

    int cornerDistance[16];

    int x1min=sCloud1.getBeginX();
    int x1max=sCloud1.getEndX();
    int y1min=sCloud1.getBeginY();
    int y1max=sCloud1.getEndY();

    int x2min=sCloud2.getBeginX();
    int x2max=sCloud2.getEndX();
    int y2min=sCloud2.getBeginY();
    int y2max=sCloud2.getEndY();

    cornerDistance[0] =Cluster::kspd(x1min,y1min,x2min,y2min);
    cornerDistance[1] =Cluster::kspd(x1min,y1min,x2min,y2max);
    cornerDistance[2] =Cluster::kspd(x1min,y1min,x2max,y2min);
    cornerDistance[3] =Cluster::kspd(x1min,y1min,x2max,y2max);

    cornerDistance[4] =Cluster::kspd(x1min,y1max,x2min,y2min);
    cornerDistance[5] =Cluster::kspd(x1min,y1max,x2min,y2max);
    cornerDistance[6] =Cluster::kspd(x1min,y1max,x2max,y2min);
    cornerDistance[7] =Cluster::kspd(x1min,y1max,x2max,y2max);

    cornerDistance[8] =Cluster::kspd(x1max,y1min,x2min,y2min);
    cornerDistance[9] =Cluster::kspd(x1max,y1min,x2min,y2max);
    cornerDistance[10]=Cluster::kspd(x1max,y1min,x2max,y2min);
    cornerDistance[11]=Cluster::kspd(x1max,y1min,x2max,y2max);

    cornerDistance[12]=Cluster::kspd(x1max,y1max,x2min,y2min);
    cornerDistance[13]=Cluster::kspd(x1max,y1max,x2min,y2max);
    cornerDistance[14]=Cluster::kspd(x1max,y1max,x2max,y2min);
    cornerDistance[15]=Cluster::kspd(x1max,y1max,x2max,y2max);

    int minDistance=cornerDistance[0];

    for (int i=1; i<16; i++) {

        if (minDistance>cornerDistance[i])
            minDistance=cornerDistance[i];
    }
    return minDistance;
}

uint GHM::calculateMinimalKspdBetweenClouds(SynthenicCloud& sCloud1, SynthenicCloud& sCloud2, uint frameThickness) const
{
    int minKspd=2*frameThickness+1;

    vector<AnchorPoint> outerAPCloud1;
    vector<AnchorPoint> outerAPCloud2;

    outerAPCloud1=sCloud1.calcAPInOuterFrame(frameThickness);
    outerAPCloud2=sCloud2.calcAPInOuterFrame(frameThickness);

    int xAP1,yAP1;
    int xAP2,yAP2;

    int calcKspd;

    for (int i=0; i<outerAPCloud1.size(); i++) {
        xAP1=outerAPCloud1[i].getX();
        yAP1=outerAPCloud1[i].getY();

        for (int j=0; j<outerAPCloud2.size(); j++) {
            xAP2=outerAPCloud2[j].getX();
            yAP2=outerAPCloud2[j].getY();

            calcKspd=Cluster::kspd(xAP1,yAP1,xAP2,yAP2);

            if (calcKspd<minKspd)
                minKspd=calcKspd;

        }

    }
    return minKspd;
}

uint GHM::calculateMinimalKspdBetweenCloudsBruteForce(SynthenicCloud& sCloud1, SynthenicCloud& sCloud2) const
{
    vector<AnchorPoint>::const_iterator it1=sCloud1.getAPBegin();
    vector<AnchorPoint>::const_iterator it2=sCloud2.getAPBegin();

    vector<AnchorPoint>::const_iterator it1End=sCloud1.getAPEnd();
    vector<AnchorPoint>::const_iterator it2End=sCloud2.getAPEnd();

    unsigned int x1,y1,x2,y2; //AP coordinates in both clouds

    x1=it1->getX();
    y1=it1->getY();
    x2=it2->getX();
    y2=it2->getY();

    uint minDist = Cluster::kspd(x1,y1,x2,y2);
    uint dist;

    for (; it1!=it1End; it1++) {
        for (; it2!=it2End; it2++) {

            x1=it1->getX();
            y1=it1->getY();
            x2=it2->getX();
            y2=it2->getY();

            dist=Cluster::kspd(x1,y1,x2,y2);
            if (minDist > dist)
                minDist=dist;
        }
    }
    return minDist;
}

bool GHM::boxOverlap(SynthenicCloud& sCloud1, SynthenicCloud& sCloud2) const
{
    int x1min=sCloud1.getBeginX();
    int x1max=sCloud1.getEndX();
    int y1min=sCloud1.getBeginY();
    int y1max=sCloud1.getEndY();

    int x2min=sCloud2.getBeginX();
    int x2max=sCloud2.getEndX();
    int y2min=sCloud2.getBeginY();
    int y2max=sCloud2.getEndY();

    //check if corners of box 1 are inside of box2
    int box1CornersX[4]={x1min,x1min,x1max,x1max};
    int box1CornersY[4]={y1min,y1max,y1min,y1max};

    for (int i=0; i<4; i++) {
        if (sCloud2.coordInCloudBox(box1CornersX[i],box1CornersY[i]))
            return true;
    }

    //check vice versa
    int box2CornersX[4]={x2min,x2min,x2max,x2max};
    int box2CornersY[4]={y2min,y2max,y2min,y2max};

    for (int i=0; i<4; i++) {
        if (sCloud1.coordInCloudBox(box2CornersX[i],box2CornersY[i]))
            return true;
    }

    return false;
}

bool GHM::cloudsOverlap(SynthenicCloud& sCloud1, SynthenicCloud& sCloud2) const
{
    //search in cross section of two boxes for APs

    vector<AnchorPoint>::const_iterator itCloud1=sCloud1.getAPBegin();
    vector<AnchorPoint>::const_iterator itCloud2=sCloud2.getAPBegin();

    int x,y;

    //check if AP of cloud 1 are in cross section
    for (; itCloud1!=sCloud1.getAPEnd(); itCloud1++ ) {

        x=itCloud1->getX();
        y=itCloud1->getY();

        if ((x >= sCloud2.getBeginX()) and (x <= sCloud2.getEndX())) {

            if ((y >= sCloud2.getBeginY()) and (y <= sCloud2.getEndY()))
                return true;
        }

    }

    //check if AP of cloud 2 are in cross section
    for (; itCloud2!=sCloud2.getAPEnd(); itCloud2++ ) {

        x=itCloud2->getX();
        y=itCloud2->getY();

        if ((x >= sCloud1.getBeginX()) and (x <= sCloud1.getEndX())) {

            if ((y >= sCloud1.getBeginY()) and (y <= sCloud1.getEndY()))
                return true;
        }

    }
    return false;
}



void GHM::joinClouds(list<SynthenicCloud*>::iterator cloudIt1, list<SynthenicCloud*>::iterator cloudIt2, bool bf)
{
    vector<AnchorPoint>::const_iterator it2=(*cloudIt2)->getAPBegin();
    vector<AnchorPoint>::const_iterator it2End=(*cloudIt2)->getAPEnd();

    SynthenicCloud& totalCloud=*(*cloudIt1);

    for (; it2!=it2End; it2++)
        totalCloud.addAnchorPoint(*it2);

    delete (*cloudIt2); //free memory
    sClouds.erase(cloudIt2);
    if (!bf)
        enrichMergedClouds(totalCloud);

}

void GHM::enrichMergedClouds(SynthenicCloud& sCloud)
{
    vector<AnchorPoint> foundAP;
    addAPFromSearchBox(sCloud.getBeginX(), sCloud.getEndX(), sCloud.getBeginY(), sCloud.getEndY(), sCloud, foundAP);
    removeAddedAnchorPoints(foundAP);
}

void GHM::filterCloudsBinomialD(double probCutoff)
{
    double APDensity=calculateAPDensity();

    list<SynthenicCloud*>::iterator it=sClouds.begin();
    for (; it!=sClouds.end();) {
        double prob=(*it)->calculateProbabilityBinomialD(APDensity);

        if (prob>probCutoff) {
            filteredSC.push_back(*it);
            sClouds.erase(it++);
        }
        else it++;
    }
}

void GHM::filterCloudsBinomialDCorr(double probCutoff)
{
    double APDensity=calculateAPDensity();

    list<SynthenicCloud*>::iterator it=sClouds.begin();
    for (; it!=sClouds.end(); ) {
        double prob=(*it)->calculateProbabilityBinomialDCorr(APDensity);

        if (prob>probCutoff) {
            filteredSC.push_back(*it);
            sClouds.erase(it++);
        }
        else it++;
    }
}

/*void GHM::filterCloudsDensityCriterium(double minDensity)
{
    list<SynthenicCloud*>::iterator it=sClouds.begin();
    for (; it!=sClouds.end(); ){

        if ((*it)->calculateCloudDensity()<minDensity) {
            filteredSC.push_back(*it);
            sClouds.erase(it++);
        }
        else it++;
    }
}*/


void GHM::getClouds(vector<SynthenicCloud*>& scv) const
{
    if (!isCloudSearch) return;

    const vector<ListElement*>& x_elements = x_object.getRemappedElements();
    const vector<ListElement*>& y_elements = y_object.getRemappedElements();

    list<SynthenicCloud*>::const_iterator it=sClouds.begin();
    for (; it!=sClouds.end(); it++) {
        SynthenicCloud& sCloud=*(*it);
        vector<AnchorPoint>::const_iterator e=sCloud.getAPBegin();

        sCloud.setXObjectID(x_object.getID());
        sCloud.setYObjectID(y_object.getID());

        //set geneInfo
        for (; e!=sCloud.getAPEnd(); e++)
        {
            int x = e->getX();
            int y = e->getY();
            const_cast<AnchorPoint&>(*e).setGeneIDs(x_elements[x]->getNumID(), y_elements[y]->getNumID());
        }
        scv.push_back(*it);
    }
}

double GHM::calculateAPDensity() const
{
    int totAP=count_points[SAME_ORIENT]+count_points[OPP_ORIENT];
    const vector<ListElement*>& xList = x_object.getRemappedElements();
    const vector<ListElement*>& yList = y_object.getRemappedElements();
    int areaGHM=xList.size()*yList.size();
    return (1.0*totAP)/(1.0*areaGHM);
}

//*************************
//GHM Visualization METHODS
//*************************
void GHM::visualizeGHM(const std::string& output_path) const
{

#ifdef HAVE_PNG
    if (isCloudSearch) visualizeSynthenicCloudsPNG(output_path);
    else visualizeBaseClustersPNG(output_path);
    return;
#endif
    if (isCloudSearch) visualizeSynthenicClouds(output_path);
    else visualizeBaseClusters(output_path);
}

#ifdef HAVE_PNG
void GHM::visualizeSynthenicCloudsPNG(const std::string& output_path) const
{
    string filename = "SYNTHGHM_"+x_object.getListName()+"_"
                      +y_object.getListName()+".png";
    cout << "Visualize: " << filename  << endl;

    const vector<ListElement*>& xList = x_object.getRemappedElements();
    const vector<ListElement*>& yList = y_object.getRemappedElements();

    Grafix png(xList.size() + 4 - (xList.size() % 4), yList.size()); //FIXME is this still necessary?

    const map<int, set<int> > &mat = matrix[MIXED_ORIENT];
    map<int, set<int> >::const_iterator itX = mat.begin();
    map<int, set<int> >::const_iterator endX = mat.end();

    png.setDrawingColor(white);

    // plot all AP not in any cloud in white
    for ( ; itX != endX; itX++) {
        int x = itX->first;
        const set<int> setY = itX->second;

        set<int>::const_iterator itY = setY.begin();
        set<int>::const_iterator endY = setY.end();

        for ( ; itY != endY; itY++) {
            int y = *itY;
            png.putPixel(x, y);
        }
    }

    //draw bounding box all clouds
    png.setDrawingColor(green);

    list<SynthenicCloud*>::const_iterator it=sClouds.begin();
    for (; it!=sClouds.end(); it++) {
        int xmin=(*it)->getBeginX();
        int ymin=(*it)->getBeginY();
        int xmax=(*it)->getEndX();
        int ymax=(*it)->getEndY();
        png.drawBox(xmin,xmax,ymin,ymax);

    }

    //draw bounding box (filteredClouds
    png.setDrawingColor(red);
    for (int i=0; i<filteredSC.size(); i++) {
        int xmin=filteredSC[i]->getBeginX();
        int ymin=filteredSC[i]->getBeginY();
        int xmax=filteredSC[i]->getEndX();
        int ymax=filteredSC[i]->getEndY();
        png.drawBox(xmin,xmax,ymin,ymax);
    }

    //draw AP in clouds
    it=sClouds.begin();
    png.setDrawingColor(yellow);
    for (; it!=sClouds.end(); it++) {
        vector<AnchorPoint>::const_iterator itAP=(*it)->getAPBegin();

        for (; itAP!=(*it)->getAPEnd(); itAP++) {
            int x=itAP->getX();
            int y=itAP->getY();
            png.putPixel(x, y);
        }
    }

    //draw AP in filtered clouds in white
    png.setDrawingColor(white);
    for (int i=0; i<filteredSC.size(); i++) {
        vector<AnchorPoint>::const_iterator it=filteredSC[i]->getAPBegin();

        for (; it!=filteredSC[i]->getAPEnd(); it++) {
            int x=it->getX();
            int y=it->getY();
            png.putPixel(x, y);
        }
    }

    std::string outputName=output_path+filename;
    png.saveCanvasPng(const_cast<char*>(outputName.c_str()));
}
#endif

void GHM::visualizeSynthenicClouds(const std::string& output_path) const
{
    string filename = "SYNTHGHM_"+x_object.getListName()+"_"
                      +y_object.getListName()+".bmp";
    cout << "Visualize: " << filename  << endl;

    const vector<ListElement*>& xList = x_object.getRemappedElements();
    const vector<ListElement*>& yList = y_object.getRemappedElements();

    Grafix bmp(xList.size() + 4 - (xList.size() % 4), yList.size());


    const map<int, set<int> > &mat = matrix[MIXED_ORIENT];
    map<int, set<int> >::const_iterator itX = mat.begin();
    map<int, set<int> >::const_iterator endX = mat.end();

    bmp.setDrawingColor(white);

    // plot all AP not in any cloud in white
    for ( ; itX != endX; itX++) {
        int x = itX->first;
        const set<int> setY = itX->second;

        set<int>::const_iterator itY = setY.begin();
        set<int>::const_iterator endY = setY.end();

        for ( ; itY != endY; itY++) {
            int y = *itY;
            bmp.putPixel(x, y);
        }
    }

    //draw bounding box all good clouds
    list<SynthenicCloud*>::const_iterator it=sClouds.begin();
    bmp.setDrawingColor(green);
    for (; it!=sClouds.end(); it++) {
        const SynthenicCloud& sCloud=*(*it);
        int xmin=sCloud.getBeginX();
        int xmax=sCloud.getEndX();
        int ymin=sCloud.getBeginY();
        int ymax=sCloud.getEndY();
        bmp.drawBox(xmin,xmax,ymin,ymax);
    }

    //draw bounding box filtered clouds
    bmp.setDrawingColor(red);
    for (int i=0; i<filteredSC.size(); i++) {
        int xmin=filteredSC[i]->getBeginX();
        int xmax=filteredSC[i]->getEndX();
        int ymin=filteredSC[i]->getBeginY();
        int ymax=filteredSC[i]->getEndY();
        bmp.drawBox(xmin,xmax,ymin,ymax);
    }

    //draw AP in good clouds
    it=sClouds.begin();
    bmp.setDrawingColor(yellow);
    for (; it!=sClouds.end(); it++) {
        const SynthenicCloud& sCloud=*(*it);
        vector<AnchorPoint>::const_iterator itAP=sCloud.getAPBegin();

        for (; itAP!=sCloud.getAPEnd(); itAP++) {
            int x=itAP->getX();
            int y=itAP->getY();
            bmp.putPixel(x, y);
        }
    }

    //draw AP in filtered clouds
    bmp.setDrawingColor(white);
    for (int i=0; i<filteredSC.size(); i++) {
        vector<AnchorPoint>::const_iterator it=filteredSC[i]->getAPBegin();

        for (; it!=filteredSC[i]->getAPEnd(); it++) {
            int x=it->getX();
            int y=it->getY();
            bmp.putPixel(x, y);
        }
    }

    std::string outputName=output_path+filename;
    bmp.saveCanvasBmp(outputName);
}

#ifdef HAVE_PNG
void GHM::visualizeBaseClustersPNG(const std::string& output_path) const
{
    string filename = "COLGHM__"+x_object.getListName()+"_"
                      +y_object.getListName()+".png";
    cout << "Visualize: " << filename  << endl;

    map<int, set<int> >::iterator itX;
    set<int>::const_iterator itY;

    const vector<ListElement*>& xList = x_object.getRemappedElements();
    const vector<ListElement*>& yList = y_object.getRemappedElements();

    Grafix png(xList.size() + 4 - (xList.size() % 4), yList.size());

    for (int orient = 0; orient < 2; orient++) {
        const map<int, set<int> > &mat = matrix[orient];
        map<int, set<int> >::const_iterator itX = mat.begin();
        map<int, set<int> >::const_iterator endX = mat.end();

        png.setDrawingColor(white);
        // plot the dots which have never been in any cluster
        for ( ; itX != endX; itX++) {
            int x = itX->first;
            const set<int> setY = itX->second;

            set<int>::const_iterator itY = setY.begin();
            set<int>::const_iterator endY = setY.end();

            for ( ; itY != endY; itY++) {
                int y = *itY;
                png.putPixel(x, y);
            }
        }

        vector<BaseCluster*>::const_iterator it;
        it = filteredBC[orient].begin();
        while (it != filteredBC[orient].end()) {
            png.setDrawingColor(red);
            int xmin=(*it)->getLowestX();
            int xmax=(*it)->getHighestX();
            int ymin=(*it)->getLowestY();
            int ymax=(*it)->getHighestY();
            png.drawBox(xmin,xmax,ymin,ymax);

            png.setDrawingColor(white);
            multiset<AnchorPoint>::const_iterator e = (*it)->getAPBegin();
            for ( ; e != (*it)->getAPEnd(); e++) {
                int x = e->getX();
                int y = e->getY();
                png.putPixel(x, y);
            }
            it++;
        }
    }

    //NOTE: code is used to visualize intermediate clusters (during runcollinear)
    /*for (int i=0; i<baseclusters[1].size(); i++)
    {
        png.setDrawingColor(green);
        int xmin=baseclusters[1][i]->getLowestX();
        int xmax=baseclusters[1][i]->getHighestX();
        int ymin=baseclusters[1][i]->getLowestY();
        int ymax=baseclusters[1][i]->getHighestY();
        png.drawBox(xmin,xmax,ymin,ymax);


        png.setDrawingColor(blue);
        multiset<AnchorPoint>::const_iterator f = baseclusters[1][i]->getAPBegin();
        double upL,downL;
        double upR,downR;
        double xL,xR;
        xL=f->getX();
        baseclusters[1][i]->intervalBounds(xL,upL,downL);
        f++;
        for ( ; f != baseclusters[1][i]->getAPEnd(); f++) {
            xR=f->getX();
            baseclusters[1][i]->intervalBounds(xR,upR,downR);
            //png.drawLine(xL,upL,xR,upR);
            //png.drawLine(xL,downL,xR,downR);
            png.drawCircle(xR,upR,1.0);
            png.drawCircle(xR,downR,1.0);

            upL=upR;
            downL=downR;
            xL=xR;
        }

        png.setDrawingColor(yellow);
        multiset<AnchorPoint>::const_iterator e = baseclusters[1][i]->getAPBegin();
        for ( ; e != baseclusters[1][i]->getAPEnd(); e++) {
                int x = e->getX();
                int y = e->getY();
                png.drawCircle(x, y,2.0);
        }



    }

    for (int i=0; i<baseclusters[0].size(); i++)
    {

        png.setDrawingColor(green);
        int xmin=baseclusters[0][i]->getLowestX();
        int xmax=baseclusters[0][i]->getHighestX();
        int ymin=baseclusters[0][i]->getLowestY();
        int ymax=baseclusters[0][i]->getHighestY();
        png.drawBox(xmin,xmax,ymin,ymax);


        png.setDrawingColor(blue);
        multiset<AnchorPoint>::const_iterator f = baseclusters[0][i]->getAPBegin();
        double upL,downL;
        double upR,downR;
        double xL,xR;
        xL=f->getX();
        baseclusters[0][i]->intervalBounds(xL,upL,downL);
        f++;
        for ( ; f != baseclusters[0][i]->getAPEnd(); f++) {
            xR=f->getX();
            baseclusters[0][i]->intervalBounds(xR,upR,downR);
            //png.drawLine(xL,upL,xR,upR);
            //png.drawLine(xL,downL,xR,downR);
            png.drawCircle(xR,upR,1.0);
            png.drawCircle(xR,downR,1.0);
            upL=upR;
            downL=downR;
            xL=xR;
        }



        png.setDrawingColor(yellow);
        multiset<AnchorPoint>::const_iterator e = baseclusters[0][i]->getAPBegin();
        for ( ; e != baseclusters[0][i]->getAPEnd(); e++) {
                int x = e->getX();
                int y = e->getY();
                png.drawPixel(x, y);
        }

    }*/

    //code below is switched of when we want to visualize intermediate clusters
    vector<Multiplicon*>::const_iterator it = multiplicons.begin();
    for (int j=0; j<multiplicons.size(); j++) {
        vector<BaseCluster*> BCs=multiplicons[j]->getBaseClusters();

        for (int i=0; i<BCs.size(); i++) {
            png.setDrawingColor(green);
            int xmin=BCs[i]->getLowestX();
            int xmax=BCs[i]->getHighestX();
            int ymin=BCs[i]->getLowestY();
            int ymax=BCs[i]->getHighestY();
            png.drawBox(xmin,xmax,ymin,ymax);

            png.setDrawingColor(blue);
            multiset<AnchorPoint>::const_iterator f = BCs[i]->getAPBegin();
            double upL,downL;
            double upR,downR;
            int xL,xR;
            xL=f->getX();
            BCs[i]->intervalBounds(xL,upL,downL);
            f++;
            for ( ; f != BCs[i]->getAPEnd(); f++) {
                xR=f->getX();
                BCs[i]->intervalBounds(xR,upR,downR);
                png.drawCircle(xR,upR,1.0);
                png.drawCircle(xR,downR,1.0);
                upL=upR;
                downL=downR;
                xL=xR;
            }

            multiset<AnchorPoint>::const_iterator e = BCs[i]->getAPBegin();
            png.setDrawingColor(yellow);
            for ( ; e != BCs[i]->getAPEnd(); e++) {
                int x = e->getX();
                int y = e->getY();
                png.drawPixel(x, y);
            }
        }
    }

    std::string outputName=output_path+filename;
    png.saveCanvasPng(const_cast<char*>(outputName.c_str()));

}
#endif

void GHM::visualizeBaseClusters(const std::string& output_path) const
{

    string filename = "COLGHM__"+x_object.getListName()+"_"
                      +y_object.getListName()+".bmp";
    cout << "Visualize: " << filename  << endl;

    map<int, set<int> >::iterator itX;
    set<int>::const_iterator itY;

    const vector<ListElement*>& xList = x_object.getRemappedElements();
    const vector<ListElement*>& yList = y_object.getRemappedElements();

    Grafix bmp(xList.size() +  4 - (xList.size() % 4), yList.size());

    for (int orient = 0; orient < 2; orient++) {
        const map<int, set<int> > &mat = matrix[orient];
        map<int, set<int> >::const_iterator itX = mat.begin();
        map<int, set<int> >::const_iterator endX = mat.end();
        bmp.setDrawingColor(white);

        // plot the dots which haven't been in any cluster
        for ( ; itX != endX; itX++) {
            int x = itX->first;
            const set<int> setY = itX->second;

            set<int>::const_iterator itY = setY.begin();
            set<int>::const_iterator endY = setY.end();

            for ( ; itY != endY; itY++) {
                int y = *itY;
                bmp.putPixel(x, y);
            }
        }

        vector<BaseCluster*>::const_iterator it;
        it = filteredBC[orient].begin();
        while (it != filteredBC[orient].end()) {
            bmp.setDrawingColor(red);
            int xmin=(*it)->getLowestX();
            int xmax=(*it)->getHighestX();
            int ymin=(*it)->getLowestY();
            int ymax=(*it)->getHighestY();
            bmp.drawBox(xmin,xmax,ymin,ymax);

            bmp.setDrawingColor(white);
            multiset<AnchorPoint>::const_iterator e = (*it)->getAPBegin();
            for ( ; e != (*it)->getAPEnd(); e++) {
                int x = e->getX();
                int y = e->getY();
                bmp.putPixel(x, y);
            }
            it++;
        }
    }

    vector<Multiplicon*>::const_iterator it = multiplicons.begin();
    for (int j=0; j<multiplicons.size(); j++) {
        vector<BaseCluster*> BCs=multiplicons[j]->getBaseClusters();
        for (int i=0; i<BCs.size(); i++) {
            bmp.setDrawingColor(green);
            int xmin=BCs[i]->getLowestX();
            int xmax=BCs[i]->getHighestX();
            int ymin=BCs[i]->getLowestY();
            int ymax=BCs[i]->getHighestY();
            bmp.drawBox(xmin,xmax,ymin,ymax);

            bmp.setDrawingColor(blue);
            multiset<AnchorPoint>::const_iterator f = BCs[i]->getAPBegin();
            double upL,downL;
            double upR,downR;
            int xL,xR;
            xL=f->getX();
            BCs[i]->intervalBounds(xL,upL,downL);
            f++;
            for ( ; f != BCs[i]->getAPEnd(); f++) {
                xR=f->getX();
                BCs[i]->intervalBounds(xR,upR,downR);
                bmp.drawCircle(xR,upR,1.0);
                bmp.drawCircle(xR,downR,1.0);
                upL=upR;
                downL=downR;
                xL=xR;
            }

            multiset<AnchorPoint>::const_iterator e = BCs[i]->getAPBegin();
            bmp.setDrawingColor(yellow);
            for ( ; e != BCs[i]->getAPEnd(); e++) {
                int x = e->getX();
                int y = e->getY();
                bmp.putPixel(x, y);
            }
        }
    }

    std::string outputName=output_path+filename;
    bmp.saveCanvasBmp(outputName);

}

void GHM::removeAPFromMultiplicons(const vector<Multiplicon*>& mps) {

    int APCounter=0;
    vector<AnchorPoint> APinMps;
    for (int i=0; i<mps.size(); i++) {

        vector<BaseCluster*> BCs=mps[i]->getBaseClusters();

        for (int j=0; j<BCs.size(); j++)
        {
            multiset<AnchorPoint>::const_iterator it=BCs[j]->getAPBegin();

            for (; it!=BCs[j]->getAPEnd(); it++) {
                APinMps.push_back(*it);
                APCounter++;
            }
        }
    }
    removeAddedAnchorPoints(APinMps);
    count_points[MIXED_ORIENT]-=APCounter;
    assert(count_points[MIXED_ORIENT]>=0);
}

void GHM::filterCloudsBinomialDFDR(double probCutoff)
{
    //calculate the random probabilities using the Binomial Distribution
    double APdensity=calculateAPDensity();
    list<SynthenicCloud*>::iterator it = sClouds.begin();
    while (it !=sClouds.end()) {
        (*it)->setRandomProbabilityBinomialD(APdensity);
        it++;
    }

    //sort clouds according to there pvalue
    sClouds.sort(pValueSortClouds);

    int k=0; //number of relevant clouds
    for (it=sClouds.begin(); it!=sClouds.end(); it++ ) {
        if ((*it)->getRandomProbability() > (k+1)*probCutoff/count_points[MIXED_ORIENT])
            break;
        k++;
    }

    // move the filtered clouds to filteredSC and erase them from sClouds
    for (; it!=sClouds.end(); it++) {
        filteredSC.push_back(*it);
        sClouds.erase(it++);
    }

}

void GHM::filterCloudsBinomialDCorrFDR(double probCutoff)
{
    //calculate the random probabilities using the Binomial Distribution
    double APdensity=calculateAPDensity();
    list<SynthenicCloud*>::iterator it = sClouds.begin();
    while (it !=sClouds.end()) {
        (*it)->setRandomProbabilityBinomialDCorr(APdensity);
        it++;
    }

    //sort clouds according to there pvalue
    sClouds.sort(pValueSortClouds);

    int k=0; //number of relevant clouds
    for (it=sClouds.begin(); it!=sClouds.end(); it++) {
        if ((*it)->getRandomProbability() > (k+1)*probCutoff/count_points[MIXED_ORIENT])
            break;
        k++;
    }

    // move the filtered clouds to filteredSC and erase them from sClouds
    for (; it!=sClouds.end(); it++) {
        filteredSC.push_back(*it);
        sClouds.erase(it++);
    }

}
