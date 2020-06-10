#include "Node.h"
#include <climits>
#include <cstdlib>

using namespace std;

// ============================================================================
// NODE CLASS
// ============================================================================

Node::Node (const ListElement* le, unsigned int _segment_id,
            unsigned int _order, unsigned int _original_nr)
        : listelement(le), segment_id(_segment_id),
        order(_order), original_nr(_original_nr) {

}

void Node::addLink(Node& equal, double weight)
{
    // don't add yourself as an equal
    if (&equal == this) return;

    links.insert(NodeLink(&equal, weight));
}

void Node::removeLink(Node& node)
{
    set<NodeLink>::iterator it = links.find(&node);
    if (it == links.end()) return;
    invalidatedLinks.insert(*it);
    links.erase(it);
}

void Node::reinstateLink(Node &node)
{
    set<NodeLink>::iterator it = invalidatedLinks.find(&node);
    if (it == invalidatedLinks.end()) return;
    links.insert(*it);
    invalidatedLinks.erase(it);
}

void Node::removeOtherLinks(set<NodeLink> &linksToBeKept)
{
    set<NodeLink>::iterator it;

    for (it = links.begin(); it != links.end(); ) {
        if (linksToBeKept.find(*it) == linksToBeKept.end()) {
            Node *node = it->getNode();
            node->removeLink(*this);
            it++; // increment the iterator before removing link!
            removeLink(*node);
        } else
            it++;
    }

    for (it = linksToBeKept.begin(); it != linksToBeKept.end(); it++)
        reinstateLink(*it->getNode());
}

const NodeLink* Node::getEqual(unsigned int segment_id) const
{
    int closest = INT_MAX;
    const NodeLink *closestLink = NULL;
    set<NodeLink>::const_iterator it = links.begin();
    for ( ; it != links.end(); it++) {
        Node *node = it->getNode();
        if (node->getSegmentId() != segment_id) continue;
        if (node->order < closest) {
            closest = node->order;
            closestLink = &(*it);
        }
    }

    return closestLink;
}

unsigned int Node::getCountEqualsDistinctSegments() const
{
    set<int> segments;
    set<NodeLink>::iterator e = links.begin();
    for ( ; e != links.end(); e++)
        segments.insert(e->getNode()->segment_id);

    return segments.size();
}
