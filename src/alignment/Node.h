#ifndef NODE_HEADER_H
#define NODE_HEADER_H

#include <set>
#include <map>
#include <vector>

class Node;
class ListElement;

// ============================================================================
// NODELINK CLASS
// ============================================================================

class NodeLink {

public:
    NodeLink(Node *equal_, double weight_ = 1.0) :
            equal(equal_), weight(weight_) { }

    Node *getNode() const {
        return equal;
    }

    double getWeight() const {
        return weight;
    }

private:
    Node *equal;            // pointer to the target node
    double weight;          // weight of the link

    friend bool operator<(const NodeLink &lhs, const NodeLink &rhs);
};

// ============================================================================
// NODE CLASS
// ============================================================================

class Node {

public:
    /**
    * Creates a Node object
    * @param le The ListElement object where the node should be created from
    * @param segment_id An ID for the segment that contains the Listelement
    * @param order An order number of anchorpoints in the initial genelist
    * @param original_number The original order number in the initial
    *                        genelist, including non-anchorpoints
    */
    Node (const ListElement* le, unsigned int segment_id,
          unsigned int order, unsigned int original_number);

    /**
     * Adds a node as an equal node, and creates a link to that node
     * @param equal The node to be added as equal
     * @param weight The weight of the link between the nodes
     */
    void addLink(Node& equal, double weight = 1.0);

    /**
     * Removes a link to the specified node
     * @param node The node where the link to should be removed
     */
    void removeLink(Node& node);

    /**
     * Reinstate a link to the specified node
     * @param node The node where the link to should be reinstated
     */
    void reinstateLink(Node& node);

    /**
     * Remove all links, except for a specified set of links
     * @param nodesToBeKept A set of links that may not be removed
     */
    void removeOtherLinks(std::set<NodeLink> &linksToBeKept);

    /**
     * Returns whether or not the node has links to other equal nodes
     * @return True of false
     */
    bool hasLinks() const {
        return !links.empty();
    }

    /**
     * Returns the total number of links to other equal nodes
     * @return Total number of links to other nodes
     */
    unsigned int getNumLinks() const {
        return links.size();
    }

    /**
     * Returns the total number of invalidated links
     * @return Total number of invalidated links
     */
    unsigned int getNumInvalidatedLinks() const {
        return invalidatedLinks.size();
    }

    /**
     * Returns the total number of links to other nodes on seperate segments
     * @return Total number of links to other nodes
     */
    unsigned int getCountEqualsDistinctSegments() const;

    /**
     * Returns whether or not the node has the specified node as equal
     * @param node The node to test as an equal
     * @return Whether or not the node has the specified node as equal
     */
    bool hasAsEqual(Node& node) const {
        return (links.find(&node) != links.end());
    }

    /**
     * Return the closest node on a specified segment
     * @param segmend_id Specifier for the segment
     */
    const NodeLink* getEqual(unsigned int segment_id) const;

    /**
     * Get an iterator to the first link
     * @return An iterator to the first link
     */
    std::set<NodeLink>::iterator getLinkBegin() const {
        return links.begin();
    }

    /**
     * Get an iterator past the final link
     * @return An iterator past the final link
     */
    std::set<NodeLink>::iterator getLinkEnd() const {
        return links.end();
    }

    /**
     * Sets the order number
     * @param order An order number of this node in a segment
     */
    void setOrder(unsigned int order) {
        this->order = order;
    }

    /**
     * Returns the order number
     * @return The order number of the node
     */
    unsigned int getOrder() const {
        return order;
    }

    /**
     * Returns the original order number of the node in the initial genelist
     * @return The original order number
     */
    unsigned int getOriginalOrder() const {
        return original_nr;
    }

    /**
     * Returns the segment ID
     * @return The segment ID
     */
    unsigned int getSegmentId() const {
        return segment_id;
    }

    /**
     * Returns a pointer to the list element
     * @return A pointer to the list element
     */
    const ListElement* getListElement() const {
        return listelement;
    }

private:
    const ListElement* listelement;
    unsigned int segment_id;
    unsigned int order;
    unsigned int original_nr;
    std::vector<std::set<NodeLink> > links2;
    std::set<NodeLink> links;
    std::set<NodeLink> invalidatedLinks;
};

// ============================================================================
// NODELINK CLASS
// ============================================================================

inline bool operator<(const NodeLink &lhs, const NodeLink &rhs)
{
    if (lhs.equal->getSegmentId() != rhs.equal->getSegmentId())
        return lhs.equal->getSegmentId() < rhs.equal->getSegmentId();
    return (lhs.equal->getOrder() < rhs.equal->getOrder());
}

#endif
