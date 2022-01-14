#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "DataSet.h"

#include "AnchorPoint.h"
#include "BaseCluster.h"
#include "SynthenicCloud.h"
#include "GeneList.h"
#include "GenePairs.h"
#include "GeneFamily.h"
#include "Multiplicon.h"
#include "Profile.h"
#include "ListElement.h"
#include "Gene.h"
#include "GHMProfile.h"
#include "Settings.h"

#include <cassert>
#include "util.h"
#include "parallel.h"

#include <climits>

#include <pthread.h>

#include "AlignmentDrawer.h"

#include "alignment/GGAligner2.h"

using namespace std;

typedef struct s_triple {
    string genome;
    string genelist;
    string target;

    s_triple(const string A,const string B, const string C) :
        genome(A), genelist(B), target(C) {}

    s_triple() {
        genome = ""; genelist =""; target="";
    }

    bool operator<(const s_triple& A) const {
        string value = genome + "_" + genelist + "_" + target;
        string valueA = A.genome + "_" + A.genelist + "_" + A.target;
        return value<valueA;
    }

} t_triple;
typedef list<string> StrLst;

DataSet::DataSet(const Settings& sett) :
    settings(sett), genepairs(NULL), threads(NULL), nThreads(0), workInProgress(0)
{
    settings.displaySettings();

    Util::startChrono();
    cout << "Creating dataset...";

    const list<ListFile> &listfiles = settings.getListFiles();

    genelists.reserve(listfiles.size());
    list<ListFile>::const_iterator it = listfiles.begin();
    for (int i = 0; it != listfiles.end(); it++, i++) {
        genelists.push_back(new GeneList(it->getListName(),
                                         it->getGenomeName(),
                                         it->getFileName()));
        // give each gene list a unique ID
        genelists.back()->setID(i);
    }

    // give each list element a unique ID
    for (unsigned int i = 0, cnt = 0; i < genelists.size(); i++) {
        geneIDmap[cnt] = genelists[i];
        vector<ListElement*> &elements = genelists[i]->getElements();
        for (unsigned int j = 0; j < elements.size(); j++, cnt++)
            elements[j]->setNumID(cnt);
    }

    cout << "\t\t\tdone. (time: " << Util::stopChrono() << "s)" << endl;
}

const Gene& DataSet::getGene(unsigned int geneID)
{
    // find the first element whose key is greater than geneID
    map<int, GeneList*>::iterator it = geneIDmap.upper_bound(geneID);
    assert(it != geneIDmap.begin());

    // iterator now points to the last element with a key <= geneID
    it--;

    return it->second->getElements()[geneID - it->first]->getGene();
}

void DataSet::deleteEvaluatedMultiplicons()
{
    //delete multiplicons
    vector<Multiplicon*>::const_iterator it = evaluated_multiplicons.begin();
    for ( ; it != evaluated_multiplicons.end(); it++)
        delete (*it);

    evaluated_multiplicons.clear();
}

void DataSet::deleteEvaluatedClouds()
{
    //delete clouds
    vector<SynthenicCloud*>::const_iterator it=clouds.begin();
    for (; it!=clouds.end(); it++)
        delete (*it);

    clouds.clear();
}

DataSet::~DataSet()
{
    for (int i = 0; i < genelists.size(); i++)
        delete genelists[i];

    // delete any remaining evaluated multiplicons
    deleteEvaluatedMultiplicons();

    // delete the genepairs
    delete genepairs;
}

void DataSet::mapGenes()
{
    Util::startChrono();
    if (settings.useFamily() == true) {
        cout << "Mapping gene families..."; cout.flush();
        getGeneFamilies();
    } else {
        cout << "Mapping gene pairs..."; cout.flush();
        getGenePairs();
    }

    cout << "\t\tdone. (time: " << Util::stopChrono() << "s)" << endl;
}

void DataSet::getGenePairs()
{
    // load the gene pairs
    genepairs = new GenePairs(settings.getBlastTable());

    unsigned int count = 0;
    for (unsigned int i = 0; i < genelists.size(); i++) {
        count += genelists[i]->getElementsLength();
    }

    hash_set<string, stringhash> geneID (count);
    for (unsigned int i = 0; i < genelists.size(); i++) {
        vector<ListElement*>& list = genelists[i]->getElements();
        for (unsigned int j = 0; j < list.size(); j++) {
            geneID.insert(list[j]->getGene().getID());
        }
    }

    for (unsigned int i = 0; i < genelists.size(); i++) {
        vector<ListElement*>& list = genelists[i]->getElements();
        for (unsigned int j = 0; j < list.size(); j++) {
            hash_set<string,stringhash>* pairs =
            genepairs->getPairsOf(list[j]->getGene().getID());

            list[j]->getGene().setPairs(*pairs);

            if (pairs == NULL) continue;

            // remove genes from the pairs set, if they do not
            // exist in the genes
            hash_set<string, stringhash>::iterator it;
            vector<string> removeEl;
            for (it = pairs->begin(); it != pairs->end(); it++) {
                if (geneID.find(*it) == geneID.end()) {
                    removeEl.push_back(*it);
                }
            }
            for (unsigned int k = 0; k < removeEl.size(); k ++) {
                pairs->erase(removeEl[k]);
            }
        }
    }
}

void DataSet::getGeneFamilies()
{
    GeneFamily genefamily(settings.getBlastTable());

    for (unsigned int i = 0; i < genelists.size(); i++) {
        vector<ListElement*>& list = genelists[i]->getElements();
        for (unsigned int j = 0; j < list.size(); j++) {

            string ID(list[j]->getGene().getID());

            string fam = genefamily.getFamilyOf(ID);
            list[j]->getGene().setFamily(genefamily.getFamilyOf(ID));
        }
    }
}

void DataSet::remapTandems()
{
    Util::startChrono();
    cout << "Remapping tandem duplicates..."; cout.flush();
    for (unsigned int i = 0; i < genelists.size(); i++)
        genelists[i]->remapTandems(settings.getTandemGap(),
                                   settings.useFamily());

    cout << "\tdone. (time: " << Util::stopChrono() << "s)" << endl;
}

void DataSet::indexToXY(int i, int N, int &x, int &y)
{
    double D = (N-0.5)*(N-0.5)-2.0*(i+1-N);
    assert(D > 0.0);  // if D is negative, input is corrupt

    double r = (N - 0.5) - sqrt(D); // chose solution that is < N

    x = (int)ceil(r);
    y = i - ( N + N*(x-1) -((x-1)*x)/2 ) + x;

    assert ((x >= 0) && (x < N));
    assert ((y >= 0) && (y < N));
}

void DataSet::finishWorkPacket()
{
    int firstItem, nItems;
    while (runGetSomeWork(firstItem, nItems, 0, false))
        (this->*workFunction)(firstItem, nItems, 0, masterTMultiplicons,masterTClouds);

    pthread_mutex_lock(&mpcMutex);
    localMultiplicons.insert(localMultiplicons.end(),
                             masterTMultiplicons.begin(),
                             masterTMultiplicons.end());
    localClouds.insert(localClouds.end(),
                        masterTClouds.begin(),
                        masterTClouds.end());

    masterTMultiplicons.clear();
    masterTClouds.clear();
    pthread_mutex_unlock(&mpcMutex);
}

void DataSet::addMultiplicons(const vector<Multiplicon*>& mplicons)
{
    multiplicons.insert(multiplicons.end(),
                mplicons.begin(), mplicons.end());
}

void DataSet::addClouds(const vector<SynthenicCloud*>& synthClouds)
{
   clouds.insert(clouds.end(),
                synthClouds.begin(), synthClouds.end());
}

bool SortMultiplicons(const Multiplicon *m1, const Multiplicon *m2)
{
    // first sort on number of anchorpoints
    if (m1->getCountAnchorPoints() != m2->getCountAnchorPoints())
        return (m1->getCountAnchorPoints() < m2->getCountAnchorPoints());

    // then sort on dpd size
    if (m1->dpd() != m2->dpd())
        return (m1->dpd() > m2->dpd());

    // then sort on x-GeneList ID
    if (m1->getXObjectID() != m2->getXObjectID())
        return (m1->getXObjectID() < m2->getXObjectID());

    // then sort on y-GeneList ID
    if (m1->getYObjectID() != m2->getYObjectID())
        return (m1->getYObjectID() < m2->getYObjectID());

    // finally sort on lowest x-index
    return (m1->getLowestX() < m2->getLowestX());
}

void DataSet::sortByMultipliconSize(vector<Multiplicon*>& mplicons) const
{
    // a (much) faster sort algoritm than bubblesort
    std::sort(mplicons.begin(), mplicons.end(), SortMultiplicons);
}

bool DataSet::allMasked(const GeneList &list, int begin, int end)
{
    assert(begin < end);
    const vector<ListElement*>& elements = list.getRemappedElements();

    vector<ListElement*>::const_iterator it = elements.begin() + begin;
    for ( ; it != elements.begin() + (end + 1); it++)
        if (!(*it)->isMasked())
            return false;

    return true;
}

void DataSet::prepareOutput()
{
    // create output directory
    const string &outputPath = settings.getOutputPath();
    string command = "mkdir -p ";
    if (system(command.append(outputPath).c_str()) == -1)
        throw runtime_error("Cannot create directory " + outputPath);

    // create the output filenames
    geneFile = settings.getOutputPath();
    geneFile.append("genes.txt");

    multipliconsFile = settings.getOutputPath();
    multipliconsFile.append("multiplicons.txt");

    baseclustersFile = settings.getOutputPath();
    baseclustersFile.append("baseclusters.txt");

    anchorpointsFile = settings.getOutputPath();
    anchorpointsFile.append("anchorpoints.txt");

    mplpairsFile = settings.getOutputPath();
    mplpairsFile.append("multiplicon_pairs.txt");

    segmentsFile = settings.getOutputPath();
    segmentsFile.append("segments.txt");

    listElementsFile = settings.getOutputPath();
    listElementsFile.append("list_elements.txt");

    alignmentFile = settings.getOutputPath();
    alignmentFile.append("alignment.txt");

    synthenicCloudsFile = settings.getOutputPath();
    synthenicCloudsFile.append("clouds.txt");

    cloudAnchorPointsFile = settings.getOutputPath();
    cloudAnchorPointsFile.append("cloudAP.txt");

    // clear the contents of the empty files
    // we can then safely append data to these files
    ofstream ofs;

    switch (settings.getClusterType()){

        case Collinear:

            ofs.open(multipliconsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tgenome_x\tlist_x\tparent\tgenome_y\tlist_y\tlevel\t"
                "number_of_anchorpoints\tprofile_length\tbegin_x\tend_x\t"
                "begin_y\tend_y\tis_redundant" << endl;
            ofs.close();

            ofs.open(baseclustersFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tmultiplicon\tnumber_of_anchorpoints\torientation\t"
                "was_twisted\trandom_probability" << endl;
            ofs.close();

            ofs.open(anchorpointsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tmultiplicon\tbasecluster\tgene_x\tgene_y\tcoord_x\t"
                "coord_y\tis_real_anchorpoint" << endl;
            ofs.close();

            ofs.open(mplpairsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tmultiplicon\t\tgene_x\tgene_y\tcode" << endl;
            ofs.close();

            ofs.open(segmentsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tmultiplicon\tgenome\tlist\tfirst\tlast\torder" << endl;
            ofs.close();

            ofs.open(listElementsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tsegment\tgene\tposition\torientation" << endl;
            ofs.close();

            ofs.open(alignmentFile.c_str(), ios_base::out | ios_base::trunc);
            ofs.close();

            ofs.open(geneFile.c_str(), ios_base::out | ios_base::trunc);
            ofs.close();

        break;

        case Cloud:

            ofs.open(synthenicCloudsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tgenome_x\tlist_x\tgenome_y\tlist_y\tnumber_of_anchorpoints\tcloud_density\tdim_x\tdim_y" << endl;
            ofs.close();

            ofs.open(cloudAnchorPointsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "CloudID\tgene_x\tgene_y\tcoord_x\tcoord_y" << endl;
            ofs.close();

            ofs.open(geneFile.c_str(), ios_base::out | ios_base::trunc);
            ofs.close();

        break;

        case Hybrid:

            ofs.open(multipliconsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tgenome_x\tlist_x\tparent\tgenome_y\tlist_y\tlevel\t"
                "number_of_anchorpoints\tprofile_length\tbegin_x\tend_x\t"
                "begin_y\tend_y\tis_redundant" << endl;
            ofs.close();

            ofs.open(baseclustersFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tmultiplicon\tnumber_of_anchorpoints\torientation\t"
                "was_twisted\trandom_probability" << endl;
            ofs.close();

            ofs.open(anchorpointsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tmultiplicon\tbasecluster\tgene_x\tgene_y\tcoord_x\t"
                "coord_y\tis_real_anchorpoint" << endl;
            ofs.close();

            ofs.open(mplpairsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tmultiplicon\t\tgene_x\tgene_y\tcode" << endl;
            ofs.close();

            ofs.open(segmentsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tmultiplicon\tgenome\tlist\tfirst\tlast\torder" << endl;
            ofs.close();

            ofs.open(listElementsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tsegment\tgene\tposition\torientation" << endl;
            ofs.close();

            ofs.open(alignmentFile.c_str(), ios_base::out | ios_base::trunc);
            ofs.close();

            ofs.open(synthenicCloudsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "id\tgenome_x\tlist_x\tgenome_y\tlist_y\tnumber_of_anchorpoints\tcloud_density\tdim_x\tdim_y" << endl;
            ofs.close();

            ofs.open(cloudAnchorPointsFile.c_str(), ios_base::out | ios_base::trunc);
            ofs << "CloudID\tgene_x\tgene_y\tcoord_x\tcoord_y" << endl;
            ofs.close();

            ofs.open(geneFile.c_str(), ios_base::out | ios_base::trunc);
            ofs.close();
    }

    multipliconID = 1;
    baseclusterID = 1;
    anchorpointID = 1;
    pairID = 1;
    segmentID = 1;
    elementID = 1;
    cloudID = 1;
}

void DataSet::flushCollinear()
{
    ofstream ofs;
    ofstream seg, le;

     /**
    * Creating multiplicons.txt
    */

    ofs.open(multipliconsFile.c_str(), ios::app);
    if (!ofs) {
        cerr << "Error creating file " << multipliconsFile << endl;
    } else {
        for (unsigned int i = 0; i < evaluated_multiplicons.size(); i++) {
            Multiplicon &mpl = *evaluated_multiplicons[i];

            ofs << mpl.getId() << '\t';
            if (mpl.getLevel() == 2) {
                int xID = mpl.getXObjectID();
                GeneList& lX = *genelists[xID];
                ofs << lX.getGenomeName() << '\t';
                ofs << lX.getListName() << '\t';
                ofs << '\t';
            }
            else {
                ofs << '\t';
                ofs << '\t';
                ofs << mpl.getParentID() << '\t';
            }
            int yID = mpl.getYObjectID();
            GeneList& lY = *genelists[yID];
            ofs << lY.getGenomeName() << '\t';
            ofs << lY.getListName() << '\t';
            ofs << mpl.getLevel() << '\t';
            ofs << mpl.getCountAnchorPoints() << '\t';
            ofs << mpl.getProfile()->getSize() << '\t';
            if (mpl.getLevel() == 2) {
                int xID = mpl.getXObjectID();
                GeneList& lX = *genelists[xID];
                ListElement &fstLE = *lX.getRemappedElements()[mpl.getBeginX()];
                ListElement &lstLE = *lX.getRemappedElements()[mpl.getEndX()];
                ofs << fstLE.getGene().getCoordinate() << '\t';
                ofs << lstLE.getGene().getCoordinate() << '\t';
            } else {
                ofs << mpl.getBeginX() << '\t';
                ofs << mpl.getEndX() << '\t';
            }

            ListElement &fstLE = *lY.getRemappedElements()[mpl.getBeginY()];
            ListElement &lstLE = *lY.getRemappedElements()[mpl.getEndY()];
            ofs << fstLE.getGene().getCoordinate() << '\t';
            ofs << lstLE.getGene().getCoordinate() << '\t';
            if ( mpl.getIsRedundant() )
                ofs << "-1";
            else
                ofs << "0";
            ofs << endl;
        }
    }

    ofs.close();

    /**
    * Creating baseclusters.txt
    */

    // give all baseclusters a unique id
    for (unsigned int i = 0; i < evaluated_multiplicons.size(); i++) {
        Multiplicon &mpl = *evaluated_multiplicons[i];
        for (unsigned int j = 0; j < mpl.getBaseClusters().size(); j++)
            mpl.getBaseClusters()[j]->setID(baseclusterID++);
    }

    ofs.open(baseclustersFile.c_str(), ios::app);
    if (!ofs) {
        cerr << "Error creating file " << baseclustersFile << endl;
    } else {
        for (unsigned int i = 0; i < evaluated_multiplicons.size(); i++) {
            Multiplicon &mpl = *evaluated_multiplicons[i];

            for (unsigned int j = 0; j < mpl.getBaseClusters().size(); j++) {
                BaseCluster* basecluster = mpl.getBaseClusters()[j];

                ofs << basecluster->getID() << '\t';
                ofs << evaluated_multiplicons[i]->getId() << '\t';
                ofs << basecluster->getCountAnchorPoints() << '\t';
                if (basecluster->getOrientation())
                    ofs << "+" << '\t';
                else
                    ofs << "-" << '\t';
                if (basecluster->wasTwisted())
                    ofs << "-1" << '\t';
                else
                    ofs << "0" << '\t';
                ofs << basecluster->getRandomProbability();
                ofs << endl;
            }
        }
    }

    ofs.close();


    /**
    * Creating anchorpoints.txt
    */

    ofs.open(anchorpointsFile.c_str(), ios::app);
    if (!ofs) {
        cerr << "Error creating file " << anchorpointsFile << endl;
    } else {
        for (unsigned int i = 0; i < evaluated_multiplicons.size(); i++) {
            Multiplicon &mpl = *evaluated_multiplicons[i];

            for (unsigned int j = 0; j < mpl.getBaseClusters().size(); j++) {
                BaseCluster* basecluster = mpl.getBaseClusters()[j];

                multiset<AnchorPoint>::const_iterator e = basecluster->getAPBegin();
                for ( ; e != basecluster->getAPEnd(); e++) {
                    const AnchorPoint& anchorpoint = *e;

                    ofs << anchorpointID << '\t';
                    ofs << evaluated_multiplicons[i]->getId() << '\t';
                    ofs << basecluster->getID() << '\t';
                    const Gene &geneX = getGene(anchorpoint.getGeneXID());
                    const Gene &geneY = getGene(anchorpoint.getGeneYID());
                    ofs << geneX.getID() << '\t';
                    ofs << geneY.getID() << '\t';
                    ofs << anchorpoint.getX() << '\t';
                    ofs << anchorpoint.getY() << '\t';
                    if (anchorpoint.isRealAnchorPoint())
                        ofs << "-1";
                    else
                        ofs << "0";
                    ofs << endl;

                    anchorpointID++;
                }
            }
        }
    }

    ofs.close();

    /**
    * Creating multiplicon_pairs.txt
    */

    ofs.open(mplpairsFile.c_str(), ios::app);
    if (!ofs) {
        cerr << "Error creating file " << mplpairsFile << endl;
    } else {
        for (unsigned int i = 0; i < evaluated_multiplicons.size(); i++) {
            const Profile *profile = evaluated_multiplicons[i]->getProfile();

            set<Link>::const_iterator it = profile->getHomBegin();
            for ( ; it != profile->getHomEnd(); it++) {
                // code: 0 = unaligned homolog, 1 = aligned homolog,
                //       2 = unaligned AP, 3 = aligned AP
                int code = 0;
                if (it->isAP)
                    code += 2;
                if (it->isAligned)
                    code += 1;
                const Gene &geneX = getGene(it->geneXID);
                const Gene &geneY = getGene(it->geneYID);
                ofs << pairID++ << "\t" << evaluated_multiplicons[i]->getId()
                    << "\t" << geneX.getID() << "\t" << geneY.getID() << "\t"
                    << code << endl;
            }
        }
    }

    ofs.close();

    /**
    * Creating segments.txt and list_elements.txt
    */

    seg.open(segmentsFile.c_str(), ios::app);
    le.open(listElementsFile.c_str(), ios::app);

    if (!seg || !le) {
        if (!seg)
            cerr << "Error creating file " << segmentsFile << endl;
        if (!le)
            cerr << "Error creating file " << listElementsFile << endl;
    } else {
        for (unsigned int i = 0; i < evaluated_multiplicons.size(); i++) {
            Multiplicon &mpl = *evaluated_multiplicons[i];

            try {
                mpl.align(settings.getAlignmentMethod(),
                          settings.getMaxGapsInAlignment());
            } catch(const ProfileException& e) {
                cout << e.what() << endl;
                //continue;
            }

            int order = 0;

            for (unsigned int j = 0; j < mpl.getProfile()->getSegments().size(); j++) {
                const GeneList& segment = *mpl.getProfile()->getSegments()[j];

                Gene* first = NULL;
                Gene* last = NULL;
                unsigned int k = 0;
                while(k < segment.getRemappedElements().size() &&
                    segment.getRemappedElements()[k]->isGap()) {
                    k++;
                }
                if (k < segment.getRemappedElements().size())
                    first = &segment.getRemappedElements()[k]->getGene();
                last = first;

                for (k = 0; k < segment.getRemappedElements().size(); k++) {
                    ListElement* element = segment.getRemappedElements()[k];

                    if (!element->isGap()) {
                        le << elementID << '\t';
                        le << segmentID << '\t';
                        le << element->getGene().getID() << '\t';
                        le << k << '\t';
                        if (element->getOrientation())
                            le << "+";
                        else
                            le << "-";
                        le << endl;

                        if (element->getGene().getCoordinate() < first->getCoordinate()) {
                            first = &element->getGene();
                        }
                        if (element->getGene().getCoordinate() > last->getCoordinate()) {
                            last = &element->getGene();
                        }

                        elementID++;
                    }
                }

                seg << segmentID << '\t';
                seg << evaluated_multiplicons[i]->getId() << '\t';
                seg << segment.getGenomeName() << '\t';
                seg << segment.getListName() << '\t';
                seg << first->getID() << '\t';
                seg << last->getID() << '\t';
                seg << order;
                seg << endl;

                segmentID++;
                order++;
            }
        }
    }

    seg.close();
    le.close();

    if (settings.showAlignedProfiles())
           visualizeAlignedProfiles();
    //NOTE only for debugging purposes
    //  printProfiles();

}

void DataSet::flushClouds()
{
    ofstream ofs;

        /**
        * Creating clouds.txt
        */

        // give all clouds a unique id
        for (unsigned int i = 0; i < clouds.size(); i++)
        {
            clouds[i]->setID(cloudID++);
        }
        ofs.open(synthenicCloudsFile.c_str(), ios::app);
        if (!ofs) {
            cerr << "Error creating file " << synthenicCloudsFile << endl;
        } else {

            for (int i=0; i<clouds.size(); i++) {

                    int xID = clouds[i]->getXObjectID();
                    int yID = clouds[i]->getYObjectID();

                    const GeneList& xlist=*genelists[xID];
                    const GeneList& ylist=*genelists[yID];

                    ofs << clouds[i]->getID()                 <<"\t";
                    ofs << xlist.getGenomeName()              <<"\t";
                    ofs << xlist.getListName()                <<"\t";
                    ofs << ylist.getGenomeName()              <<"\t";
                    ofs << ylist.getListName()                <<"\t";
                    ofs << clouds[i]->getCountAnchorPoints()  <<"\t";
                    ofs << clouds[i]->calculateCloudDensity() <<"\t";
                    ofs << clouds[i]->calculateBoxWidth()     <<"\t";
                    ofs << clouds[i]->calculateBoxHeight()    << endl;
            }
        }
        ofs.close();

        /**
        * Creating cloudAP.txt
        */
        ofs.open(cloudAnchorPointsFile.c_str(), ios::app);
        if (!ofs) {
            cerr << "Error creating file " << cloudAnchorPointsFile << endl;
        } else {

            vector<AnchorPoint>::const_iterator it;

            for (int i=0; i<clouds.size(); i++) {

                int xID = clouds[i]->getXObjectID();
                int yID = clouds[i]->getYObjectID();

                GeneList& xlist=*genelists[xID];
                GeneList& ylist=*genelists[yID];

                it=clouds[i]->getAPBegin();

                for (; it!=clouds[i]->getAPEnd(); it++) {

                    ofs << clouds[i]->getID() << "\t";
                    ofs << xlist.getGeneName(it->getX()) <<"\t";
                    ofs << ylist.getGeneName(it->getY()) <<"\t";
                    ofs << it->getX()         << "\t";
                    ofs << it->getY()         << endl;

                }
            }
        }

        ofs.close();
}

void DataSet::flushOutput()
{
    cout << "Flushing output files...";
    switch (settings.getClusterType()) {
        case Collinear:
            flushCollinear();
            break;
        case Cloud:
            flushClouds();
            break;
        case Hybrid:
            flushCollinear();
            flushClouds();
            break;
    }
    cout << "done." << endl;
}

void DataSet::outputGenes()
{
    assert(!geneFile.empty());

    Util::startChrono();
    cout << "Writing genelists file..."; cout.flush();
    ofstream ofs;
    ofs.open(geneFile.c_str());
    if (!ofs) {
        cerr << "Error creating file " << geneFile << endl;
    } else {
        ofs << "id\tgenome\tlist\tcoordinate\torientation\t"
               "remapped_coordinate\tis_tandem\tis_tandem_representative\t"
               "tandem_representative\tremapped" << endl;

        for (unsigned int i = 0; i < genelists.size(); i++) {
            vector<ListElement*>& list = genelists[i]->getElements();

            for (unsigned int j = 0; j < list.size(); j++) {
                const Gene& gene = list[j]->getGene();
                ofs << gene.getID() << '\t';
                ofs << gene.getGenomeName() << '\t';
                ofs << genelists[i]->getListName() << '\t';
                ofs << gene.getCoordinate() << '\t';
                if (gene.getOrientation())
                    ofs << "+\t";
                else
                    ofs << "-\t";
                ofs << gene.getRemappedCoordinate() << '\t';
                if (gene.isTandem())
                    ofs << "-1" << '\t';
                else
                    ofs << "0" << '\t';
                if (gene.isTandemRepresentative())
                    ofs << "-1" << '\t';
                else
                    ofs << "0" << '\t';
                if (gene.isTandem())
                    ofs << gene.tandemRepresentative().getID() << '\t';
                else
                    ofs << '\t';
                if (gene.isRemapped())
                    ofs << "-1";
                else
                    ofs << "0";
                ofs << endl;
            }
        }
    }

    ofs.close();
    cout << "\t\tdone. (time: " << Util::stopChrono() << "s)" << endl;
}

void DataSet::output()
{
    ofstream ofs;

    /**
     * Creating alignment.txt
     */

    // uncomment to get alignment timings statistics
    /*cout << endl;
    cout << "Min set time : " << minSetTime << endl;
    cout << "Confl time : " << confTime << endl;
    cout << "Init time : " << initConfl << endl;
    cout << "Build compact : " << buildCompact << endl;
    cout << "max flow : " << maxFlow << endl;*/

    if (!settings.getCompareAligners()) return;

    // verbose output on-screen
    cout << endl;
    cout << "\tnumP\tNW\tGG\tRA\tRC\tRAC\tLL\tLLBS\tLS" << endl;
    for (int i = 0; i < NW.getNumLevels(); i++) {
        cout << i << "\t" << NW.getNumProfiles(i) << " &\t" <<
                             NW.getNumAlHom(i) << " &\t" <<
                             GG.getNumAlHom(i)  << " &\t" <<
                             RA.getNumAlHom(i) << " &\t" <<
                             RC.getNumAlHom(i) << " &\t" <<
                             RAC.getNumAlHom(i) << " &\t" <<
                             LL.getNumAlHom(i) << " &\t" <<
                             LLBS.getNumAlHom(i) << " &\t" <<
                             LS.getNumAlHom(i) << " \\\\" << endl;
    }

    cout << endl << " &\t" << NW.getTotalNumProfiles()
                 << " &\t" << NW.getTotalNumAlHom()
                 << " &\t" << GG.getTotalNumAlHom()
                 << " &\t" << RA.getTotalNumAlHom()
                 << " &\t" << RC.getTotalNumAlHom()
                 << " &\t" << RAC.getTotalNumAlHom()
                 << " &\t" << LL.getTotalNumAlHom()
                 << " &\t" << LLBS.getTotalNumAlHom()
                 << " &\t" << LS.getTotalNumAlHom()
                 << endl;

    cout.precision(2);
    cout << endl << " &\t"
                 << " &\t" << NW.getTotalTime()
                 << " &\t" << GG.getTotalTime()
                 << " &\t" << RA.getTotalTime()
                 << " &\t" << RC.getTotalTime()
                 << " &\t" << RAC.getTotalTime()
                 << " &\t" << LL.getTotalTime()
                 << " &\t" << LLBS.getTotalTime()
                 << " &\t" << LS.getTotalTime()
                 << endl;
}

int DataSet::getRemappedGenomeSize(string genome)
{
    int size = 0;

    for (vector<GeneList*>::iterator gl_it = genelists.begin(); gl_it != genelists.end(); gl_it++) {
        string gn = (*gl_it)->getGenomeName();
        if (gn == genome)
        {
            size += (*gl_it)->getRemappedElementsLength();
        }
    }

    return size;
}

int DataSet::getRemappedGeneListSize(string genome, string list)
{
    int size = 0;

    for (vector<GeneList*>::iterator gl_it = genelists.begin(); gl_it != genelists.end(); gl_it++)
    {
        string gn = (*gl_it)->getGenomeName();
        string gl = (*gl_it)->getListName();
        if (gn == genome && gl==list)
        {
            size += (*gl_it)->getRemappedElementsLength();
        }
    }

    return size;
}

void DataSet::getGenomes()
{
    const list<ListFile> &listFiles = settings.getListFiles();

    //get all genomenames
    list<ListFile>::const_iterator it = listFiles.begin();
    for (int i = 0; it != listFiles.end(); it++, i++) {
        string genomename = it->getGenomeName();
        genomes.remove(genomename);
        genomes.push_back(genomename);
    }

    genomes.sort();

}

void DataSet::getDuplicatedPortions(){
    char buffer2 [256];
    if (settings.getOutputPath()[settings.getOutputPath().length() - 1] == '/')
        sprintf(buffer2, "%sduplicated_portions.txt", settings.getOutputPath().c_str());
    else         sprintf(buffer2, "%s/duplicated_portions.txt", settings.getOutputPath().c_str());

    ofstream outLog (buffer2);

    //loop over all genomes and count the % in duplicated regions , ...
    list<string>::iterator git = genomes.begin();
    for (int i=1; git != genomes.end(); git++,i++)
    {
        string gn = *git;
        map<string,StrLst*> listToDuplicates;   //will store all duplicates for a given list
        map<string, int> listToRemappedSize;    //store remapped lengths
        StrLst lists;

        for (vector<GeneList*>::iterator gl_it = genelists.begin(); gl_it != genelists.end(); gl_it++)
        {
            string listGenome = (*gl_it)->getGenomeName();
            if (gn == listGenome)
            {
                string gl = (*gl_it)->getListName();
                lists.push_back((*gl_it)->getListName());
                listToRemappedSize[gl] = (*gl_it)->getRemappedElementsLength();
            }
        }

        for (unsigned int j = 0; j < evaluated_multiplicons.size(); j++) {
            int found_segment = 0; //do we have multiple segments from the required species in this multiplicon ?

            for (unsigned int k = 0; k < evaluated_multiplicons[j]->getProfile()->getSegments().size(); k++) {
                const GeneList& segment = *evaluated_multiplicons[j]->getProfile()->getSegments()[k];

                string segmentg = segment.getGenomeName();

                if (segmentg == gn)
                {
                    found_segment++;
                }

                //found 2 segments of the right genome, this is enough to go through
                if (found_segment > 1)
                {
                    break;
                }
            }

            //if enough segments are present in the genome, continue
            if (found_segment > 1)
            {
                for (unsigned int k = 0; k < evaluated_multiplicons[j]->getProfile()->getSegments().size(); k++) {
                    const GeneList& segment = *evaluated_multiplicons[j]->getProfile()->getSegments()[k];

                    string segmentg = segment.getGenomeName(); //get genome name

                    if (segmentg == gn) //genome is correct & we know the multiplicon is correct so continue
                    {
                        string segmentl = segment.getListName(); //get the listname

                        //now store elements
                        StrLst* pList = listToDuplicates[segmentl];
                        if (pList == NULL)
                        {
                            pList = new StrLst();
                            listToDuplicates[segmentl] = pList;
                        }
                        for (int l = 0; l < segment.getRemappedElements().size(); l++)
                        {
                            ListElement* element = segment.getRemappedElements()[l];

                            pList->remove(element->getGene().getID());
                            if (element->getGene().getID() != "")
                            {
                                pList->push_back(element->getGene().getID());
                            }
                        }
                    }

                }
            }
        }

        outLog << i << "\t" << gn << endl << "\tlist\tsize\tduplicates\tpercentage" << endl;
        int totalSize=0;
        int totalDuplicated=0;
        for(list<string>::iterator genelist_it = lists.begin(); genelist_it != lists.end(); genelist_it++)
        {
            int duplicates = 0;
            string gl = *genelist_it;
            int totalLength = listToRemappedSize[gl];
            float percentage = 0;

            StrLst* duplicatedGenes = listToDuplicates[gl];

            if (duplicatedGenes != NULL)
            {
                duplicates = duplicatedGenes->size();
            }

            percentage = (float)duplicates*100.0f/(float)totalLength;
            totalSize += totalLength;
            totalDuplicated += duplicates;
            outLog << "\t" <<  gl << "\t" << totalLength << "\t" << duplicates << "\t" << percentage << " %"<< endl;
            if (duplicates > totalLength)
            {
                for(list<string>::iterator gene_it = duplicatedGenes->begin(); gene_it != duplicatedGenes->end(); gene_it++)
                {
                    outLog << *gene_it << endl;
                }
            }

        }

        float percentage = (float)totalDuplicated*100.0f/(float)totalSize;
        outLog << "\t" << "total" << "\t" << totalSize << "\t" << totalDuplicated << "\t" << percentage << " %" << endl << endl;
    }
}

void DataSet::sortGeneLists()
{
    // recover from previous run
    indexToList.clear();
    indexToList.reserve(genelists.size());

    // create a mapping weight -> list
    multimap<lluint, uint> weightToList;
    for (uint l = 0; l < genelists.size(); l++) {
        lluint w = genelists[l]->getRemappedElementsLength();
        weightToList.insert(pair<uint, uint>(w, l));
    }

    // create a mapping index -> list
    multimap<lluint, uint>::reverse_iterator it = weightToList.rbegin();
    for ( ; it != weightToList.rend(); it++)
        indexToList.push_back(it->second);
}

void DataSet::getCollinearPortions()
{
    //loop over all genomes and count the % in collinear regions if there is more then one genome in the dataset
    if (genomes.size() > 1)
    {
        map<string,StrLst*> genomeToList;
        map<t_triple, StrLst*> data;

        //store genelists per genomes necessary to generate output
        for (vector<GeneList*>::iterator gl_it = genelists.begin(); gl_it != genelists.end(); gl_it++)
        {
            string gn = (*gl_it)->getGenomeName();
            string gl = (*gl_it)->getListName();

            StrLst* pList = genomeToList[gn];
            if (pList == NULL)
            {
                pList = new StrLst();
                genomeToList[gn] = pList;
            }
            pList->push_back(gl);
        }

        //loop trough all multiplicons
        for (unsigned int i = 0; i < evaluated_multiplicons.size(); i++)
        {
            list<string> genomesInMultiplicon;

            //loop trough all segments and store genomes
            for (unsigned int j = 0; j < evaluated_multiplicons[i]->getProfile()->getSegments().size(); j++) {
                const GeneList& segment = *evaluated_multiplicons[i]->getProfile()->getSegments()[j];

                string genomeName = segment.getGenomeName();

                genomesInMultiplicon.remove(genomeName);
                genomesInMultiplicon.push_back(genomeName);

            }

            //loop trough all segments and store collinearity information
            for (unsigned int j = 0; j < evaluated_multiplicons[i]->getProfile()->getSegments().size(); j++) {
                const GeneList& segment = *evaluated_multiplicons[i]->getProfile()->getSegments()[j];

                string genomeName = segment.getGenomeName();
                string listName = segment.getListName();

                for (list<string>::iterator list_it = genomesInMultiplicon.begin(); list_it != genomesInMultiplicon.end();list_it++)
                {
                    t_triple key;

                    key.genome = genomeName;
                    key.genelist = listName;
                    key.target = *list_it;

                    //if not yet defined define
                    StrLst* pList = data[key];
                    if (pList == NULL)
                    {
                        pList = new StrLst();
                        data[key] = pList;
                    }

                    if (genomeName != *list_it)
                    {
                        for (int k = 0; k < segment.getRemappedElements().size(); k++) {
                            ListElement* element = segment.getRemappedElements()[k];

                            pList->remove(element->getGene().getID());
                            if (element->getGene().getID() != "")
                            {
                                pList->push_back(element->getGene().getID());
                            }
                        }
                    }
                }
            }
        }

        //start writing the output

        char buffer2 [256];
        if (settings.getOutputPath()[settings.getOutputPath().length() - 1] == '/')
            sprintf(buffer2, "%scollinear_portions.txt", settings.getOutputPath().c_str());
        else         sprintf(buffer2, "%s/collinear_portions.txt", settings.getOutputPath().c_str());

        ofstream outLog (buffer2);

        for (list<string>::iterator genome_it = genomes.begin(); genome_it != genomes.end(); genome_it++)
        {
            StrLst* pgeneLists = genomeToList[*genome_it];
            StrLst geneLists = *pgeneLists;

            outLog << *genome_it << endl;
            for (list<string>::iterator second_genome_it = genomes.begin(); second_genome_it != genomes.end(); second_genome_it++)
            {
                if (*genome_it != *second_genome_it)
                {
                    outLog << "list\tsize\t" << *second_genome_it << endl;

                    int totalSize = 0;
                    int totalCollinear = 0;

                    for(list<string>::iterator genelist_it = geneLists.begin(); genelist_it != geneLists.end(); genelist_it++)
                    {
                        t_triple key;
                        key.genome = *genome_it;
                        key.target = *second_genome_it;
                        key.genelist = *genelist_it;
                        int listSize = getRemappedGeneListSize(key.genome, key.genelist);
                        int listCollinearSize = 0;

                        StrLst* collineargenes = data[key];
                        if (collineargenes != NULL)
                        {
                            listCollinearSize = collineargenes->size();
                        }

                        totalSize += listSize;
                        totalCollinear += listCollinearSize;

                        float percentage = (float)listCollinearSize*100.0f/(float)listSize;
                        outLog << key.genelist << "\t" << listSize << "\t" << listCollinearSize << "\t" << percentage << " %" << endl;
                    }
                    float percentage = (float)totalCollinear*100.0f/(float)totalSize;
                    outLog << endl << percentage <<" % of the genome of "<< *genome_it << " is collinear to " << *second_genome_it << "." << endl << endl;

                }
            }
        }
    }
}

void DataSet::statistics()
{
    if (!settings.writeStatistics())
        return;

    cout << "Generating Statistics..."; cout.flush();
    getGenomes(); //this function could be better, genomes should be stored when ini file is read !!!!
    getDuplicatedPortions();
    getCollinearPortions();
    cout << "done." << endl;
}


void DataSet::visualizeAlignedProfiles()
{
    cout << "Visualize AlignedProfiles" << endl;
    char buffer[50];

    string tempfilename=settings.getOutputPath()+"AlignmentMultiplicon";

    for (int i=0; i<evaluated_multiplicons.size(); i++) {
        Multiplicon& multiplicon=*evaluated_multiplicons[i];

        try {
                multiplicon.align(settings.getAlignmentMethod(),
                                  settings.getMaxGapsInAlignment());
        } catch(const ProfileException& e) {
                cout << e.what() << endl;
                continue;
        }

        int id=multiplicon.getId();
        sprintf(buffer,"%i.svg",multiplicon.getId());
        string filename=tempfilename+string(buffer);
        AlignmentDrawer drawer(evaluated_multiplicons[i]->getProfile());
        bool succes=drawer.buildColorMatrix(this,settings.getTandemGap());
        drawer.visualizeConflicts=true;
        if (succes){
            drawer.generateAlignmentSVG(filename);
        }
    }
}

void DataSet::printProfiles()
{
    char buffer[50];
    string tempfilename=settings.getOutputPath()+"PrintMultiplicon";
    for (int i=0; i<evaluated_multiplicons.size(); i++) {
        const Multiplicon& multiplicon=*evaluated_multiplicons[i];

        sprintf(buffer,"%i.data",multiplicon.getId());
        string filename=tempfilename+string(buffer);
        ofstream ofs;
        ofs.open(filename.c_str());

        vector<GeneList*> segments=evaluated_multiplicons[i]->getProfile()->getSegments();
        ofs << "Alignments: " << endl;
        for (int j=0; j<segments.size(); j++) {
            vector<ListElement*>::const_iterator leIt=segments[j]->getLEBegin();
            for (; leIt!=segments[j]->getLEEnd(); leIt++) {

                if ((*leIt)->isGap())
                    ofs << "gap ";
                else
                    ofs << (*leIt)->getGene().getID() << " ";
            }
            ofs << endl;
            ofs << endl;

        }
        ofs << endl;
        ofs << "Links" << endl;

        set<Link> links=evaluated_multiplicons[i]->getHomologs();
        set<Link>::iterator linkIt;
        for (linkIt=links.begin(); linkIt!=links.end(); linkIt++) {
            const Gene& geneX=getGene(linkIt->geneXID);
            const Gene& geneY=getGene(linkIt->geneYID);
            ofs << linkIt->geneXID  << " " << geneX.getID() << "\t" << linkIt->geneYID << " " << geneY.getID() << endl;
        }

        ofs.close();
    }
}

GeneList* DataSet::getGeneList(const string& listName, const string& genomeName) const
{
    for (int i=0; i<genelists.size(); i++){
        if (genelists[i]->getListName().compare(listName)==0){
            if (genelists[i]->getGenomeName().compare(genomeName)==0){
                return genelists[i];
            }
        }
    }
    cerr << "Genelist with this name and genome not found!!" << endl;
    return NULL;
}
