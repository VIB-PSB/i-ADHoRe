#include "Profile.h"
#include "AnchorPoint.h"
#include "BaseCluster.h"
#include "Gene.h"
#include "GHM.h"
#include "ListElement.h"
#include "Multiplicon.h"
#include "alignment/NWAligner.h"
#include "alignment/GGAligner.h"
#include "alignment/GGAligner2.h"
#include "alignment/Aligner.h"
#include "alignment/AlignmentException.h"
#include "alignComp.h"

#include "util.h"

#include <cassert>
#include <iostream>

typedef vector<ListElement* >::const_iterator VecListElementCIt;

using namespace std;

Profile::Profile(const Multiplicon& m, unsigned int _id) :
    multiplicon(m), valid(true), id(_id)
{
    segments = multiplicon.getXSegments();
    unaligned_y_list = multiplicon.getYSegment();

    // create structures to keep track of the y-list permutation of genes
    applyYListPermutation();
    segments.push_back(unaligned_y_list);
}

void Profile::createNodes(const set<Link> &links,
                          const vector<GeneList*> &segments,
                          vector<vector<Node *> >& nodes,
                          bool onlyAP, bool priorityAP) const
{
    // add the nodes and create the node map
    nodes.resize(segments.size());
    vector<map<int, Node*> > nodeMap(segments.size());

    vector<GeneList*>::const_iterator gl = segments.begin();
    for (int i = 0; gl != segments.end(); gl++, i++) {
        vector<ListElement *>::const_iterator e = (*gl)->getLEBegin();
        for (int lePos = 0, nodePos = 0; e != (*gl)->getLEEnd(); e++) {
            ListElement &le = *(*e);
            if (le.isGap()) continue;
            if (!le.hasHomolog) {
                lePos++;
                continue;
            }
            if (onlyAP && (!le.hasAP)) {
                lePos++;
                continue;
            }

            // create the node
            Node *node = new Node(&le, i, nodePos++, lePos++);
            nodes[i].push_back(node);

            // Insert the position of the node in the nodeMap
            int ID = le.getNumID();
            assert(nodeMap[i].find(ID) == nodeMap[i].end());
            nodeMap[i][ID] = nodes[i].back();
        }
    }

    // create the links between the nodes
    set<Link>::const_iterator e = links.begin();
    for ( ; e != links.end(); e++) {
        if (onlyAP && (!e->isAP)) continue;
        // find the nodes based on the geneID
        map<int, Node*>::const_iterator itX, itY;
        itX = nodeMap[e->segmentX].find(e->geneXID);
        itY = nodeMap[e->segmentY].find(e->geneYID);
        assert (itX != nodeMap[e->segmentX].end());
        assert (itY != nodeMap[e->segmentY].end());
        Node &nodeX = *itX->second;
        Node &nodeY = *itY->second;

        // create the links between the nodes
        double weight = 1.0;
        if (priorityAP && e->isAP)
            weight = 1.0;
        if (priorityAP && !e->isAP)
            weight = 0.1;
        nodeX.addLink(nodeY, weight);
        nodeY.addLink(nodeX, weight);

        // safety checks
        const Gene &geneX = nodeX.getListElement()->getGene();
        const Gene &geneY = nodeY.getListElement()->getGene();
        assert(geneX.isPairWith(geneY));
        assert(!nodeX.getListElement()->isGap());
        assert(!nodeY.getListElement()->isGap());
    }

    // safety check, just in case :-)
    for (int i = 0; i < nodes.size(); i++) {
        const vector<Node *> &le = nodes[i];
        for (int j = 0; j < le.size(); j++) {
            assert (le[j]->getNumLinks() > 0);
        }
    }

    if (!priorityAP) return;
}

void Profile::destroyNodes(vector<vector<Node *> >& nodes) const
{
    for (int i = 0; i < nodes.size(); i++)
        for (int j = 0; j < nodes[i].size(); j++)
            delete nodes[i][j];
    nodes.clear();
}

bool Profile::findAnchorPoint(int geneXID, int geneYID) const
{
    int score = 0;

    for (int i = 0; i < segments.size(); i++) {
        for (int j = 0; j < segments[i]->getSize(); j++) {
            if (segments[i]->getRemappedElements()[j]->getNumID() == geneXID)
                score++;
            if (segments[i]->getRemappedElements()[j]->getNumID() == geneYID)
                score++;
            if (score == 2)
                return true;
        }
    }

    return score == 2;
}

void Profile::applyYListPermutation()
{
    const vector<BaseCluster*>& baseclusters = multiplicon.getBaseClusters();
    int size = multiplicon.getEndY() - multiplicon.getBeginY() + 1;
    int offset = multiplicon.getBeginY();

    assert(size > 0);

    // first, determine the orientation of the multiplicon
    map<int, int> BCCenters;

    vector<BaseCluster*>::const_iterator it = baseclusters.begin();
    for ( ; it != baseclusters.end(); it++) {
        BaseCluster &bc = *(*it);
        int centreX = (bc.getHighestX() + bc.getLowestX())/2;
        int centreY = (bc.getHighestY() + bc.getLowestY())/2;
        BCCenters[centreX] = centreY;
    }

    int orientMultiplicons = 0;
    map<int, int>::iterator p1 = BCCenters.begin();
    map<int, int>::iterator p2 = BCCenters.begin();
    for (p2++ ; p2 != BCCenters.end(); p1++, p2++) {
        if (p2->second > p1->second) orientMultiplicons++;
        if (p2->second < p1->second) orientMultiplicons--;
    }

    // flip the complete Y-list if necessary
    if (orientMultiplicons < 0)
        unaligned_y_list->invertSection(0, size-1);

    // now, invert pieces of the Y-list according to the orientation of the BC
    vector <double> invert(size, 0.0);

    for (it = baseclusters.begin(); it != baseclusters.end(); it++) {
        BaseCluster &bc = *(*it);
        double sign = bc.getOrientation() ? 1 : -1;
        // if we already flipped the complete section, change the sign
        if (orientMultiplicons < 0)
            sign *= -1;

        sign *= (double)bc.getCountAnchorPoints() / double(bc.getEndY() - bc.getBeginY());

        for (int i = bc.getBeginY() - offset; i <= bc.getEndY() - offset; i++)
            invert[i] += sign;
    }

    // if we already flipped the complete Y-list, also flip the invert vector
    if (orientMultiplicons < 0)
        reverse(invert.begin(), invert.end());

    int firstI = 0;
    double previous = 1;

    for (int i = 0; i < size; i++) {
        if ((invert[i] < 0) && (previous >= 0))
            firstI = i;
        if ((invert[i] >= 0) && (previous < 0)) { 
            unaligned_y_list->invertSection(firstI, i-1);
        }
        previous = invert[i];
    }

    if (previous < 0)
        unaligned_y_list->invertSection(firstI, size-1);
}

void Profile::checkMasking(bool level2Only)
{
    const vector<ListElement*>& elements = unaligned_y_list->getRemappedElements();

    bool all_masked = false;

    if (!level2Only) {
        all_masked = true;

        vector<ListElement*>::const_iterator it = elements.begin();
        for ( ; it != elements.end(); it++) {
            if (!(*it)->isMasked()) {
                all_masked = false;
                break;
            }
        }
    }

    if (all_masked) {
        throw ProfileException("all elements masked");
    } else {
        // unmask every element
        vector<ListElement*>::const_iterator it = elements.begin();
        for ( ; it != elements.end(); it++)
            (*it)->setMasked(false);
    }
}

void Profile::clearSegments(vector<GeneList*>& dst)
{
    for (int i = 0; i < dst.size(); i++)
        delete dst[i];
    dst.clear();
}

void Profile::deepCopySegments(const vector<GeneList*>& src,
                               vector<GeneList*>& dst)
{
    clearSegments(dst);
    dst.reserve(src.size());
    for (int i = 0; i < src.size(); i++)
        dst.push_back(new GeneList(*src[i], 0, src[i]->getSize()-1));
}

void Profile::compareAligners(AlignScore &NW, AlignScore &GG,
                              AlignScore &RA, AlignScore &RC,
                              AlignScore &RAC, AlignScore &LL,
                              AlignScore &LLBS, AlignScore &LS)
{
    Aligner* aligner = NULL;
    vector<vector<Node*> > nodes;
    int level = segments.size();
    double time = 0.0;

    vector<GeneList*> segmentsCopy;
    int largeInt = (int)1E8;

    try {

    // NW ALIGNER (progressive)
    deepCopySegments(segments, segmentsCopy);
    vector<GeneList*> segmentsNW;
    segmentsNW.clear();
    segmentsNW.push_back(segmentsCopy[0]);
    segmentsNW.back()->removeGaps();
    aligner = new NWAligner(multiplicon.getHomologs(), 0, 0, 1, 1, largeInt);
    for (int i = 1; i < segmentsCopy.size(); i++) {
        segmentsNW.push_back(segmentsCopy[i]);
        segmentsNW.back()->removeGaps();
        Util::startChrono();
        aligner->align(segmentsNW);
        time = Util::stopChrono();
    }
    markAlignedHomologs(segmentsCopy);
    NW.addScore(level, numAP, numAlAP, numHom, numAlHom, segmentsCopy[0]->getSize(), time);
    //cout << "NW: aligned AP: " << numAlAP << "/" << numAP
    //     << ", aligned homologs: " << numAlHom << "/" << numHom << endl;
    //printProfile(segmentsCopy);
    delete aligner;
    clearSegments(segmentsCopy);

    // GG ALIGNER
    deepCopySegments(segments, segmentsCopy);
    createNodes(multiplicon.getHomologs(), segmentsCopy, nodes, false);
    aligner = new GGAligner(largeInt);
    Util::startChrono();
    aligner->align(segmentsCopy);
    time = Util::stopChrono();
    markAlignedHomologs(segmentsCopy);
    GG.addScore(level, numAP, numAlAP, numHom, numAlHom, segmentsCopy[0]->getSize(), time);
    //cout << "GG: aligned AP: " << numAlAP << "/" << numAP
    //     << ", aligned homologs: " << numAlHom << "/" << numHom << endl;
    //printProfile(segmentsCopy);
    delete aligner;
    destroyNodes(nodes);
    clearSegments(segmentsCopy);

    // CH_RA: RANDOM ALIGNER
    deepCopySegments(segments, segmentsCopy);
    createNodes(multiplicon.getHomologs(), segmentsCopy, nodes, false);
    aligner = new GGAligner2(largeInt, nodes, CH_RA);
    Util::startChrono();
    aligner->align(segmentsCopy);
    time = Util::stopChrono();
    markAlignedHomologs(segmentsCopy);
    RA.addScore(level, numAP, numAlAP, numHom, numAlHom, segmentsCopy[0]->getSize(), time);
    //cout << "RA: aligned AP: " << numAlAP << "/" << numAP
    //     << ", aligned homologs: " << numAlHom << "/" << numHom << endl;
    //printProfile(segmentsCopy);
    delete aligner;
    destroyNodes(nodes);
    clearSegments(segmentsCopy);

    // CH_RC: RANDOM ALIGNER WITH CONFLICT SELECTOR
    deepCopySegments(segments, segmentsCopy);
    createNodes(multiplicon.getHomologs(), segmentsCopy, nodes, false);
    aligner = new GGAligner2(largeInt, nodes, CH_RC);
    Util::startChrono();
    aligner->align(segmentsCopy);
    time = Util::stopChrono();
    markAlignedHomologs(segmentsCopy);
    RC.addScore(level, numAP, numAlAP, numHom, numAlHom, segmentsCopy[0]->getSize(), time);
    //cout << "RC: aligned AP: " << numAlAP << "/" << numAP
    //     << ", aligned homologs: " << numAlHom << "/" << numHom << endl;
    //printProfile(segmentsCopy);
    delete aligner;
    destroyNodes(nodes);
    clearSegments(segmentsCopy);

    // CH_RAC: RANDOM ALIGNER WITH ACTIVE CONFLICT SELECTOR
    deepCopySegments(segments, segmentsCopy);
    createNodes(multiplicon.getHomologs(), segmentsCopy, nodes, false);
    aligner = new GGAligner2(largeInt, nodes, CH_RAC);
    Util::startChrono();
    aligner->align(segmentsCopy);
    time = Util::stopChrono();
    markAlignedHomologs(segmentsCopy);
    RAC.addScore(level, numAP, numAlAP, numHom, numAlHom, segmentsCopy[0]->getSize(), time);
    //cout << "RIC: aligned AP: " << numAlAP << "/" << numAP
     //    << ", aligned homologs: " << numAlHom << "/" << numHom << endl;
    //printProfile(segmentsCopy);
    delete aligner;
    destroyNodes(nodes);
    clearSegments(segmentsCopy);

    // LL ALIGNER
    deepCopySegments(segments, segmentsCopy);
    createNodes(multiplicon.getHomologs(), segmentsCopy, nodes, false, true);
    aligner = new GGAligner2(largeInt, nodes, CH_LL);
    Util::startChrono();
    aligner->align(segmentsCopy);
    time = Util::stopChrono();
    markAlignedHomologs(segmentsCopy);
    LL.addScore(level, numAP, numAlAP, numHom, numAlHom, segmentsCopy[0]->getSize(), time);
    //cout << "SKW: aligned AP: " << numAlAP << "/" << numAP
    //     << ", aligned homologs: " << numAlHom << "/" << numHom << endl;
    //printProfile(segmentsCopy);
    delete aligner;
    destroyNodes(nodes);
    clearSegments(segmentsCopy);

    // LLBS ALIGNER
    deepCopySegments(segments, segmentsCopy);
    createNodes(multiplicon.getHomologs(), segmentsCopy, nodes, false, true);
    aligner = new GGAligner2(largeInt, nodes, CH_LLBS);
    Util::startChrono();
    aligner->align(segmentsCopy);
    time = Util::stopChrono();
    markAlignedHomologs(segmentsCopy);
    LLBS.addScore(level, numAP, numAlAP, numHom, numAlHom, segmentsCopy[0]->getSize(), time);
    //cout << "ULF: aligned AP: " << numAlAP << "/" << numAP
    //     << ", aligned homologs: " << numAlHom << "/" << numHom << endl;
    //printProfile(segmentsCopy);
    delete aligner;
    destroyNodes(nodes);
    clearSegments(segmentsCopy);

    // LS ALIGNER
    deepCopySegments(segments, segmentsCopy);
    createNodes(multiplicon.getHomologs(), segmentsCopy, nodes, false, true);
    aligner = new GGAligner2(1E8, nodes, CH_LS);
    Util::startChrono();
    aligner->align(segmentsCopy);
    time = Util::stopChrono();
    markAlignedHomologs(segmentsCopy);
    LS.addScore(level, numAP, numAlAP, numHom, numAlHom, segmentsCopy[0]->getSize(), time);
   // cout << "LS: aligned AP: " << numAlAP << "/" << numAP
   //      << ", aligned homologs: " << numAlHom << "/" << numHom << endl;
    //printProfile(segmentsCopy);
    delete aligner;
    destroyNodes(nodes);
    clearSegments(segmentsCopy);

    } catch (const AlignmentException &e) {
        Util::stopChrono();
        clearSegments(segmentsCopy);
    }
}

void Profile::cutAlignment()
{
    // get the first AP or homolog index of the final segment
    uint first = 0;
    vector<ListElement*>::const_iterator e = segments.back()->getLEBegin();
    for ( ; e != segments.back()->getLEEnd(); e++, first++)
        if (((*e)->hasAlAP) || ((*e)->hasAlHomolog))
            break;

    // get the first AP or homolog index of the final segment
    uint last = segments.back()->getSize() - 1;
    e = segments.back()->getLEEnd();
    for (e--; e != segments.back()->getLEBegin(); e--, last--)
        if ((*e)->hasAlAP || (*e)->hasAlHomolog)
            break;

    if (first >= last)
        throw ProfileException("alignment check failed...cannot cut profile");

    // now actually cut the segments
    for (int i = 0; i < segments.size(); i++) {
        segments[i]->cutAfter(last);
        segments[i]->cutBefore(first);
    }

    //printProfile(segments);
}

void Profile::align(const AlignmentMethod &alignMethod, int maxGaps)
{
    AlignmentMethod method = alignMethod;
    /*if (multiplicon.getLevel() >= 15)
        method = NeedlemanWunsch;*/

    const set<Link> &homologs = multiplicon.getHomologs();

    vector<vector<Node*> > nodes;
    Aligner* aligner = NULL;
    try {
        switch (method) {
        case NeedlemanWunsch:
            aligner = new NWAligner(homologs, 0, 0, 1, 1, maxGaps);
            break;
        case GreedyGraphbased:
            aligner = new GGAligner(maxGaps);
            break;
        case GreedyGraphbased2:
            createNodes(homologs, segments, nodes, false, false);
            aligner = new GGAligner2(maxGaps, nodes, CH_LLBS);
            break;
        case GreedyGraphbased3:
            createNodes(homologs, segments, nodes, true, false);
            aligner = new GGAligner2(maxGaps, nodes, CH_LLBS);
            break;
        case GreedyGraphbased4:
            createNodes(homologs, segments, nodes, false, true);
            aligner = new GGAligner2(maxGaps, nodes, CH_LLBS);
            break;
        }

        //Util::startChrono();
        aligner->align(segments);
        //cout << "Alignment time: " << Util::stopChrono() << endl;
        markAlignedHomologs(segments);
       // cutAlignment();
       // markAlignedHomologs();  // FIXME: quick hack to recount aligned AP, etc.
    }
    catch (const AlignmentException& e) {
        //cout << "Alignment time: " << Util::stopChrono() << endl;
        delete aligner;
        destroyNodes(nodes);
        throw ProfileException(e.what());
    }

    destroyNodes(nodes);
    delete aligner;
    //cout << endl;
    //printProfile(segments);
}

void Profile::markAlignedHomologs(vector<GeneList*> &segments)
{
    // first reset from a possible prevous call to this function
    numAlHom = numAP = numAlAP = 0;
    numAlAPV.clear(); numAlAPV.resize(segments.size(), 0);
    numAlHomV.clear(); numAlHomV.resize(segments.size(), 0);

    for (int i = 0; i < segments.size(); i++) {
        const GeneList &gl = *segments[i];
        vector<ListElement *>::const_iterator e;
        for (e = gl.getLEBegin(); e != gl.getLEEnd(); e++) {
            (*e)->hasAlHomolog = false;
            (*e)->hasAlAP = false;
        }
    }

    const set<Link> &homologs = multiplicon.getHomologs();
    numHom = homologs.size();

    set<Link>::iterator it = homologs.begin();
    for ( ; it != homologs.end(); it++)
        if (it->isAP) numAP++;

    // initialize an iterator for each segment in the segments
    vector<vector<ListElement *>::const_iterator > e(segments.size());
    for (int i = 0; i < segments.size(); i++)
        e[i] = segments[i]->getLEBegin();

    for ( ; e[0] != segments[0]->getLEEnd(); ) {
        // compare all X and Y elements
        for (int i = 0; i < segments.size() - 1; i++) {
            ListElement &elX = *(*e[i]);
            if (!elX.hasHomolog) continue;

            for (int j = i + 1; j < segments.size(); j++) {
                ListElement &elY = *(*e[j]);
                if (!elY.hasHomolog) continue;
                if (!elX.getGene().isPairWith(elY.getGene()))
                    continue;

                // elX and elY are homologs
                elX.hasAlHomolog = true;
                elY.hasAlHomolog = true;
                numAlHom++;

                numAlHomV[i]++;
                numAlHomV[j]++;

                // find the corresponding link to check if they are AP
                set<Link>::iterator it = homologs.find(Link(i, j, elX.getNumID(),
                                                            elY.getNumID()));
                assert (it != homologs.end()); //assert that link still exists!

                if (!it->isAP) continue;
                elX.hasAlAP = true;
                elY.hasAlAP = true;
                numAlAP++;

                numAlAPV[i]++;
                numAlAPV[j]++;
            }
        }

        // advance all iterators by one position
        for (int i = 0; i < segments.size(); i++)
            e[i]++;
    }
}

void Profile::printProfile(const vector<GeneList*>& profile)
{
    for (int i = 0; i < profile.size(); i++) {
        const GeneList &gl = *profile[i];
        vector<ListElement *>::const_iterator e;
        for (e = gl.getLEBegin(); e != gl.getLEEnd(); e++) {
            const ListElement &le = *(*e);
            if (le.isGap())
                cout << "_";
            else if (le.hasAlAP)
                cout << "A";
            else if (le.hasAlHomolog)
                cout << "H";
            else if (le.hasHomolog)
                cout << ".";
            else
                cout << "x";
        }
        cout << endl;
    }
}

double Profile::calculateAlignmentScore(vector<GeneList>& lists)
{
    int gap = 0;
    int match = 1;

    int score = 0;

    int total_genes = 0;

    for (unsigned int j = 0; j < lists[0].getSize(); j++) {
        vector<bool> evaluated (lists.size(), false);
        for (unsigned int i = 0; i < lists.size(); i++) {
            const vector<ListElement*>& elements_i = lists[i].getRemappedElements();
            if (elements_i[j]->isGap()) {
                score += gap;
            }
            else {
                total_genes++;
                bool homolog = false;
                for (unsigned int k = 0; k < lists.size(); k++) {
                    if (k != i && !evaluated[k]) {
                        const vector<ListElement*>& elements_k = lists[k].getRemappedElements();
                        if (!elements_k[j]->isGap() && elements_k[j]->getGene().isPairWith(elements_i[j]->getGene())) {
                            score += match;
                            evaluated[k] = true;
                            homolog = true;
                        }

                        if (homolog) {
                            //extra score
                            score += match;
                        }
                    }
                }
            }
            evaluated[i] = true;
        }
    }

    return (double)score / (double)total_genes / (double)lists.size();
}

void Profile::checkAlignment(int minHomologs)
{
    // create a mapping between the number of aligned homologs and the segment
    multimap<int, int> segHom;
    for (int i = 0; i < segments.size(); i++)
        segHom.insert(pair<int, int>(numAlHomV[i], i));

    vector<GeneList*> newProfile;

    map<int, int>::reverse_iterator it;
    for (it = segHom.rbegin(); it != segHom.rend(); it++) {
        newProfile.push_back(segments[it->second]);
        if (newProfile.size() < 2) continue;

        int countAlHom = 0;
        for (int i = 0; i < newProfile.back()->getSize(); i++) {
            ListElement *elY = newProfile.back()->getRemappedElements()[i];
            if (elY->isGap()) continue;

            bool foundST = false;
            for (int j = 0; j < newProfile.size() - 1; j++) {
                ListElement *elX = newProfile[j]->getRemappedElements()[i];
                if (elX->isGap()) continue;

                if (elX->getGene().isPairWith(elY->getGene())) {
                    countAlHom++;
                    foundST = true;
                    break;
                }

                assert(!foundST);
            }

            // if we found sufficient aligned homologs, get out
            if (countAlHom >= minHomologs)
                break;
        }

        if (countAlHom < minHomologs)
            throw ProfileException("too few aligned homologs");
    }
}
