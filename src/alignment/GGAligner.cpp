#include "GGAligner.h"

#include "../GeneList.h"
#include "../ListElement.h"
#include "Node.h"
#include <climits>

#include "AlignmentException.h"


GGAligner::GGAligner (int gapSize) : gap(gapSize/2) {

}

void GGAligner::align(vector<GeneList*>& unaligned_x_lists)
{
    //search for elements that have equals
    vector<vector<bool> > equals (unaligned_x_lists.size());

    for (unsigned int i = 0; i < unaligned_x_lists.size(); i++) {
        unaligned_x_lists[i]->removeGaps();
        equals[i].resize(unaligned_x_lists[i]->getSize(), false);
    }

    for (unsigned int i = 0; i < unaligned_x_lists.size() - 1; i++) {
        for (unsigned int j = i+1; j < unaligned_x_lists.size(); j++) {
            const vector<ListElement*>& elements_i =
                unaligned_x_lists[i]->getRemappedElements();
            const vector<ListElement*>& elements_j =
                unaligned_x_lists[j]->getRemappedElements();

            for (unsigned int x = 0; x < elements_i.size(); x++) {
                if (elements_i[x]->isGap()) continue;
                for (unsigned int y = 0; y < elements_j.size(); y++) {
                    const Gene& xGene = elements_i[x]->getGene();
                    const Gene& yGene = elements_j[y]->getGene();
                    if (xGene.isPairWith(yGene)) {
                        equals[i][x] = true;
                        equals[j][y] = true;
                    }
                }
            }
        }
    }

    //initialisation of segments vector with nodes and deletion of gaps in original genelists
    segments.resize(unaligned_x_lists.size());
    for (unsigned int i = 0; i < unaligned_x_lists.size(); i++) {
        segments.reserve(unaligned_x_lists[i]->getSize());
        const vector<ListElement*>& elements = unaligned_x_lists[i]->getRemappedElements();
        unsigned int order = 0;
        for (unsigned int j = 0; j < elements.size(); j++) {
            if (!equals[i][j]) continue;
            segments[i].push_back(Node (elements[j], i, order, j));
            order++;
        }
    }

    //creating connections from every node to his equal (homologous) nodes
    for (unsigned int i = 0; i < segments.size() - 1; i++) {
        for (unsigned int j = i+1; j < segments.size(); j++) {
            for (unsigned int x = 0; x < segments[i].size(); x++) {
                for (unsigned int y = 0; y < segments[j].size(); y++) {
                    const Gene& xGene = segments[i][x].getListElement()->getGene();
                    const Gene& yGene = segments[j][y].getListElement()->getGene();
                    if (xGene.isPairWith(yGene)) {
                        segments[i][x].addLink(segments[j][y]);
                        segments[j][y].addLink(segments[i][x]);
                    }
                }
            }
        }
    }

    //bool vector initialisation (false represents a gap, true represents a node)
    alignings.resize(segments.size());

    vector<bool> waiting (segments.size(), false);

    //keep track of positions
    vector<unsigned int> positions (segments.size(), 0);
    vector<unsigned int> new_positions (segments.size(), 0);

    //mark true if position is at the end of its segment
    vector<bool> at_end (segments.size(), false);
    for (uint i = 0; i < segments.size(); i++)
        at_end[i] = (segments[i].size() ==  0);
    //count number of positions at the end of its segment
    unsigned int c_at_end = 0;

    while (c_at_end != segments.size()) {
        for (unsigned int i = 0; i < positions.size(); i++) {
            if (positions[i] == segments[i].size()) continue;
            if (segments[i][positions[i]].hasLinks()) {
                //waiting at anchorpoint
                waiting[i] = true;
            } else {
                //no equals, so just skip node
                new_positions[i] = positions[i] + 1;
                alignings[i].push_back(-1);
            }
        }
        c_at_end = 0;
        for (unsigned int j = 0; j < new_positions.size(); j++) {
            if (new_positions[j] == segments[j].size()) {
                c_at_end++;
                if (at_end[j]) {
                    alignings[j].push_back(false);
                }
                else {
                    at_end[j] = true;
                }
            }
        }

        if (c_at_end == segments.size()) break;

        vector<unsigned int> c_equals (segments.size(), 0);
        vector<vector<unsigned int> > waiting_indexes (segments.size());

        // process the waiting nodes
        for (unsigned int j = 0; j < waiting.size(); j++) {
            if (!waiting[j] || positions[j] >= segments[j].size()) continue;
            c_equals[j] = segments[j][positions[j]].getCountEqualsDistinctSegments();

            for (unsigned int k = 0; k < segments.size(); k++) {
                if (k == j || positions[k] >= segments[k].size()) continue;
                if (segments[k][positions[k]].hasAsEqual(segments[j][positions[j]])) {
                    waiting_indexes[j].push_back(k);
                } else if (c_equals[j] > 1) {
                    const NodeLink* equal = segments[j][positions[j]].getEqual(k);
                    if (equal == NULL) continue;
                    int order1 = equal->getNode()->getOriginalOrder();
                    int order2 = segments[j][positions[j]].getOriginalOrder();
                    if (abs(order1 - order2) > gap)
                        c_equals[j]--;
                }
            }

            set<NodeLink> linksToBeKept;
            if (waiting_indexes[j].size() == c_equals[j]) {
                for (unsigned int k = 0; k < waiting_indexes[j].size(); k++) {
                    linksToBeKept.insert(&segments[waiting_indexes[j][k]][positions[waiting_indexes[j][k]]]);
                }
                linksToBeKept.insert(&segments[j][positions[j]]);

                //the equal nodes
                for (unsigned int k = 0; k < waiting_indexes[j].size(); k++) {
                    segments[waiting_indexes[j][k]][positions[waiting_indexes[j][k]]].removeOtherLinks(linksToBeKept);
                    //check to see if it isnt already changed in previous iteration
                    if (waiting[waiting_indexes[j][k]]) {
                        waiting[waiting_indexes[j][k]] = false;
                        new_positions[waiting_indexes[j][k]] = positions[waiting_indexes[j][k]] + 1;
                        alignings[waiting_indexes[j][k]].push_back(1);
                    }
                }

                //the evaluated node itself
                segments[j][positions[j]].removeOtherLinks(linksToBeKept);
                waiting[j] = false;
                new_positions[j] = positions[j] + 1;
                alignings[j].push_back(1);
            }
        }

        //if all positions are still waiting for advance, make intelligent decision which node to ignore
        unsigned int c_waiting = 0;
        for (unsigned int j = 0; j < waiting.size(); j++) {
            if (waiting[j]) {
                c_waiting++;
            }
        }


        if (c_waiting == waiting.size() - c_at_end) {
            vector<vector<unsigned int> > last (segments.size(), vector<unsigned int> (segments.size()));
            vector<vector<int> > weights (segments.size(), vector<int> (segments.size(), 0));
            for (unsigned int j = 0; j < segments.size(); j++) {
                if (at_end[j]) continue;
                int offset = alignings[j].size() - positions[j];
                for (unsigned int k = positions[j]; k < segments[j].size() && k < positions[j] + 3; k++) {
                    for (unsigned int l = 0; l < segments.size(); l++) {
                        if (j == l) continue;
                        const NodeLink* equal = segments[j][k].getEqual(l);
                        if (equal == NULL) continue;
                        if (k == positions[j]) {
                            weights[j][l] += equal->getNode()->getOrder() + offset;
                            last[j][l] = equal->getNode()->getOrder();
                        }
                        else {
                            if (equal->getNode()->getOrder() < last[j][l]) {
                                weights[j][l] += last[j][l] + offset;
                            }
                            weights[j][l] -= ((int)equal->getNode()->getOrder() + offset);
                            last[j][l] = equal->getNode()->getOrder();
                        }
                    }
                }
            }

            int c_highest = INT_MIN;
            int highestI = 0;

            for (unsigned int j = 0; j < weights.size(); j++) {
                if (at_end[j]) continue;
                int weight = 0;
                for (unsigned int k = 0; k < weights[j].size(); k++) {
                    weight += weights[j][k];
                }
                if (weight > c_highest) {
                    c_highest = weight;
                    highestI = j;
                }
            }

            //advance at equal nodes to reference node, only if corresponding position is at that node
            set<NodeLink> linksToBeKept;
            set<NodeLink>::const_iterator _it = segments[highestI][positions[highestI]].getLinkBegin();
            while (_it != segments[highestI][positions[highestI]].getLinkEnd()) {
                int segment_id = _it->getNode()->getSegmentId();
                if (_it->getNode() == &segments[segment_id][positions[segment_id]]) {
                    linksToBeKept.insert(*_it);
                }
                _it++;
            }
            linksToBeKept.insert(NodeLink(&segments[highestI][positions[highestI]]));

            set<NodeLink>::const_iterator it = segments[highestI][positions[highestI]].getLinkBegin();
            while (it != segments[highestI][positions[highestI]].getLinkEnd()) {
                int segment_id = it->getNode()->getSegmentId();

                if (it->getNode() == &segments[segment_id][positions[segment_id]]) {
                    it->getNode()->removeOtherLinks(linksToBeKept);
                    alignings[segment_id].push_back(1);
                    new_positions[segment_id] = positions[segment_id] + 1;
                    waiting[segment_id] = false;
                }
                else {
                    it->getNode()->removeLink(segments[highestI][positions[highestI]]);
                }

                it++;
            }

            segments[highestI][positions[highestI]].removeOtherLinks(linksToBeKept);
            //advance at reference node
            if (linksToBeKept.size() > 1) {
                alignings[highestI].push_back(1);
            }
            else {
                alignings[highestI].push_back(-1);
            }
            new_positions[highestI] = positions[highestI] + 1;
            waiting[highestI] = false;
        }

        //mark as gaps at which positions are still waiting
        for (unsigned int j = 0; j < waiting.size(); j++) {
            if (waiting[j]) {
                alignings[j].push_back(0);
            }
        }

        //check if all positions are at the end, to stop the iteration
        c_at_end = 0;
        for (unsigned int j = 0; j < new_positions.size(); j++) {
            if (new_positions[j] == segments[j].size()) {
                c_at_end++;
                at_end[j] = true;
            }
        }

        //assign new positions
        for (unsigned int j = 0; j < positions.size(); j++) {
            if (positions[j] != new_positions[j]) {
                positions[j] = new_positions[j];
            }
        }
    }

    vector<unsigned int> position (segments.size(), 0);
    vector<unsigned int> last_position (segments.size(), 0);
    for (unsigned int j = 0; j < alignings[0].size(); j++) {
        for (unsigned int i = 0; i < alignings.size(); i++) {
            if (alignings[i][j] == -1) {
                position[i]++;
                continue;
            }
            if (alignings[i][j] != 1) continue;
            set<NodeLink>::const_iterator ite = segments[i][position[i]].getLinkBegin();
            while (ite != segments[i][position[i]].getLinkEnd()) {
                int segment_id = ite->getNode()->getSegmentId();
                int order = 0;
                int order_ref = 0;

                vector<ListElement*>::const_iterator it_ref = unaligned_x_lists[i]->getRemappedElements().begin();
                while ((*it_ref) != segments[i][position[i]].getListElement()) {
                    it_ref++;
                    order_ref++;
                }
                vector<ListElement*>::const_iterator it = unaligned_x_lists[segment_id]->getRemappedElements().begin();
                while ((*it) != ite->getNode()->getListElement()) {
                    it++;
                    order++;
                }

                if (order > order_ref) {
                    vector<ListElement*>::const_iterator it_prev = it_ref;
                    int order_prev = order_ref;
                    while (position[i] > 0 && it_prev != unaligned_x_lists[i]->getRemappedElements().rend().base()
                            && (*it_prev) != segments[i][last_position[i]].getListElement()) {
                        it_prev--;
                        order_prev--;
                    }
                    if ((order-order_ref) > gap)
                        throw AlignmentException("alignment failed...too many gaps in profile");
                    if (position[i] > 0 && (*it_prev) == segments[i][last_position[i]].getListElement()) {
                        unaligned_x_lists[i]->introduceGaps(order_prev, order_ref, order-order_ref);
                    }
                    else {
                        unaligned_x_lists[i]->introduceGaps(0, order_ref, order-order_ref);
                    }

                    segments[i][position[i]].setOrder(order);
                }
                else if (order < order_ref) {
                    vector<ListElement*>::const_iterator it_prev = it;
                    int order_prev = order;
                    while (position[segment_id] > 0 && it_prev != unaligned_x_lists[segment_id]->getRemappedElements().rend().base()
                            && (*it_prev) != segments[segment_id][last_position[segment_id]].getListElement()) {
                        it_prev--;
                        order_prev--;
                    }
                    if ((order_ref-order) > gap)
                            throw AlignmentException("alignment failed...too many gaps in profile");
                    if (position[segment_id] > 0 && (*it_prev) == segments[segment_id][last_position[segment_id]].getListElement()) {
                        unaligned_x_lists[segment_id]->introduceGaps(order_prev, order, order_ref-order);
                    }
                    else {
                        unaligned_x_lists[segment_id]->introduceGaps(0, order, order_ref-order);
                    }

                    ite->getNode()->setOrder(order_ref);
                }

                ite++;
            }

            last_position[i] = position[i];
            position[i]++;
        }
    }


    //make sure every list is equal in size -> insert gaps at the end
    unsigned int largest = 0;
    for (unsigned int i = 0; i < unaligned_x_lists.size(); i++) {
        if (unaligned_x_lists[i]->getSize() > largest) {
            largest = unaligned_x_lists[i]->getSize();
        }
    }
    for (unsigned int i = 0; i < unaligned_x_lists.size(); i++) {
        if (unaligned_x_lists[i]->getSize() < largest) {
            const vector<ListElement*>& elements = unaligned_x_lists[i]->getRemappedElements();
            unsigned int j = 0;
            if (segments[i].size() > 0)
                while (j < elements.size() && elements[j] != segments[i][segments[i].size()-1].getListElement()) {
                    j++;
                }
            int end = unaligned_x_lists[i]->getSize() - 1;
            if (end < 0)
                end = 0;
            int numGaps = largest - unaligned_x_lists[i]->getSize();
            if (numGaps > gap)
                throw AlignmentException("alignment failed...too many gaps in profile");
            unaligned_x_lists[i]->introduceGaps(j, end, numGaps);
        }
    }
}

int GGAligner::getNumAlignedPoints() const
{
    int nma = 0;

    for (uint i = 0; i < segments.size(); i++)
        for (uint j = 0; j < segments[i].size(); j++)
            nma += segments[i][j].getNumLinks();

    return nma;
}

std::ostream& operator<<(std::ostream& os, const GGAligner& a)
{
    for (uint i = 0; i < a.alignings.size(); i++) {
        for (uint j = 0; j < a.alignings[i].size(); j++) {
            if (a.alignings[i][j] == 1)
                os << "G";
            else if (a.alignings[i][j] == 0)
                os << "_";
            else
                os << ".";
        }
        os << endl;
    }
    os << endl;

    return os;
}
