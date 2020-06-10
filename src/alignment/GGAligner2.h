#ifndef __GGAligner2_H
#define __GGAligner2_H

#include "Aligner.h"
#include "Node.h"

#include <set>

extern double minSetTime;
extern double confTime;
extern double initConfl;
extern double maxFlow;
extern double buildCompact;

typedef enum { CH_RA, CH_RC, CH_RAC, CH_LL, CH_LLBS, CH_LS } ConflictHeuristic;

class CLink
{
public:

    CLink(int sIndex_, int eIndex_, int eSeg_, double weight_) :
        sIndex(sIndex_), eIndex(eIndex_), eSeg(eSeg_), flow(0.0),
        weight(weight_), inverseLink(NULL) {};

    double getCf() const {
        return weight - flow;
    }

    int sIndex;
    int eIndex;
    int eSeg;

    double flow;
    double weight;

    CLink * inverseLink;

    friend bool operator<(const CLink &lhs, const CLink &rhs) {
        if (lhs.sIndex != rhs.sIndex)
            return lhs.sIndex < rhs.sIndex;
        if (lhs.eIndex != rhs.eIndex)
            return (lhs.eIndex < rhs.eIndex);
        return (lhs.eSeg < rhs.eSeg);
    }
};

class CNode
{
    public:
        set<CLink> links;
};

class LinkReport
{
public:

    LinkReport(Node* startNode_, Node* endNode_) :
        startNode(startNode_), endNode(endNode_), ownWeight(1.0),
        length(1e8), conflict(false), activeConflict(false), fdirect(0.0),
        fstblocking(0.0), ftsblocking(0.0), linkScore(0.0), estLinkScore(0.0),
        fdmax(0.0), fstmax(0.0), ftsmax(0.0) {}

    Node *startNode;
    Node *endNode;

    double ownWeight;
    double length;

    bool conflict;
    bool activeConflict;

    // flow scores
    double fdirect;
    double fstblocking;
    double ftsblocking;

    // link scores
    double linkScore;
    double estLinkScore;

    // upper limits to the direct, blocking st and blocking ts flow respectively
    double fdmax;
    double fstmax;
    double ftsmax;
};

class GGAligner2 : public Aligner
{

public:
    /**
     * Creates an GGAligner2 object
     * @param gapSize The maximal space between two homologous elements
     * @param nodes The
     */
    GGAligner2 (int gapSize, const vector<vector<Node*> > &nodes,
                ConflictHeuristic conflHeuristic);

    /**
     * Aligns the segments
     * @param segments The vector to aligned (input/output)
     */
    void align(vector<GeneList*>& segments);

private:
    /**
     * Recursively mark nodes that are connected to a node through a blocking
     * path by using a Depth-First search (DFS) in the original graph G
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param node Pointer to the current node
     * @param visited Vector to ensure elementary paths in DFS
     * @param movedLeft True if we've moved left on the current path
     */
    void recMarkConnectivity(const vector<vector<Node*> >& segment,
                             const vector<unsigned int>& position,
                             const Node *node, vector<bool> &visited,
                             bool movedLeft);

    /**
     * Recursively mark nodes that are connected to a node through a direct
     * path by using a Depth-First search (DFS) in the original graph G
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param node Pointer to the current node
     * @param visited Vector to ensure elementary paths in DFS
     */
    void recMarkDirConnectivity(const vector<vector<Node*> >& segment,
                                const vector<unsigned int>& position,
                                const Node *node, vector<bool> &visited);

    /**
     * Recursively mark nodes that are connected to a node through a blocking
     * path by using a Depth-First search (DFS) in the reduced graph G'
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param node Pointer to the current node
     * @param visited Vector to ensure elementary paths in DFS
     * @param movedLeft True if we've moved left on the current path
     */
    bool hasActiveConnectivity(const vector<vector<Node*> >& segment,
                               const vector<unsigned int>& position,
                               const Node *node, vector<bool> &visited,
                               int targetSegment, int targetPos,
                               bool movedLeft);

    /**
     * Make the initial reports of the all links in a conflict situation
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param ACReports Links involved in at least one active conflicts (output)
     * @param RCReports Links involved in at least one active conflicts (output)
     * @param NCReports Links non involved in any conflict (output)
     */
    void makeAllReports(const vector<vector<Node*> >& segment,
                        const vector<unsigned int>& position,
                        vector<LinkReport> &ACReports,
                        vector<LinkReport> &RCReports,
                        vector<LinkReport> &NCReports);

    /**
     * Make the reports of the active links in a conflict situation
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param ACReports Links involved in at least one active conflicts (output)
     */
    void makeACReports(const vector<vector<Node*> >& segment,
                       const vector<unsigned int>& position,
                       vector<LinkReport> &ACReports);

    /**
     * Check whether a link is involved in an active conflict
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param report The link report to fill in (input/output)
     */
    void checkForActiveConflicts(const vector<vector<Node*> >& segment,
                                 const vector<unsigned int>& position,
                                 LinkReport &report);

    /**
     * Check whether a link is involved in a conflict
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param report The link report to fill in (input/output)
     */
    void checkForConflicts(const vector<vector<Node*> >& segment,
                           const vector<unsigned int>& position,
                           LinkReport &report);

    /**
     * Reset connectivity matrix
     */
    void resetConnect() {
        for (uint i = 0; i < connect.size(); i++)// {
          /*  uint j = 0;
            while (connect[i][j] != 0 && (j < searchDepth))
                connect[i][j++] = 0;
        }*/

            std::fill( connect[i].begin(), connect[i].end(), 0);
    }

    /**
     * Find a blocking augmenting path in the compressed graph
     * @param cNodes Nodes of te compressed path
     * @param visited Keep track of path, make sure it is elementary
     * @param nodeSeg Segment of the current node
     * @param nodePos Position of the current node
     * @param flow Current flow in the path
     * @param endSegment Target segment
     * @param endPos End position
     * @param movedRight True if the path is blocking
     */
    double findAugmentingPath(vector<CNode> &cNodes, vector<bool> &visited,
                              int nodeSeg, int nodePos, double flow,
                              int endSegment, int endPos, bool movedRight);

    /**
     * Find a direct augmenting path in the compressed graph
     * @param cNodes Nodes of te compressed path
     * @param visited Keep track of path, make sure it is elementary
     * @param nodeSeg Segment of the current node
     * @param nodePos Position of the current node
     * @param flow Current flow in the path
     * @param endSegment Target segment
     * @param endPos End position
     */
    double findDirAugmentingPath(vector<CNode> &cNodes, vector<bool> &visited,
                                 int nodeSeg, int nodePos, double flow,
                                 int endSegment, int endPos);

    /**
     * Calculate the flow through elementary direct paths
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param report The link report to fill in (output)
     */
    void calcDirectFlow(const vector<vector<Node*> >& segment,
                        const vector<unsigned int>& position,
                        LinkReport &report);

    /**
     * Calculate the st-flow through elementary blocking paths
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param report The link report to fill in (output)
     */
    void calcBlockingSTFlow(const vector<vector<Node*> >& segment,
                            const vector<unsigned int>& position,
                            LinkReport &report);

    /**
     * Calculate the st-flow through elementary blocking paths
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param report The link report to fill in (output)
     */
    void calcBlockingTSFlow(const vector<vector<Node*> >& segment,
                            const vector<unsigned int>& position,
                            LinkReport &report);

    /**
     * Count the link score
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param report The link report to fill in (output)
     */
    void calcLinkScore(const vector<vector<Node*> >& segment,
                       const vector<unsigned int>& position,
                       LinkReport &report);

    /**
     * Try to find a minimal set S among the nodes with index c_i
     * @param segment Matrix of nodes (input)
     * @param position Positions in segments (input)
     * @param isProcessed True of the corresponding node has been processed
     * @param seed Pointer to the node to start searching from
     * @param processed Returns the minimal set (output)
     */
    bool findMinimalSet(const vector<vector<Node*> >& segment,
                        const vector<unsigned int>& position,
                        vector<bool> &isProcessed,
                        Node* seed, set<Node*> &processed);

    /**
     * Remove the link with highest max flow from the segment matrix
     * @param segment Matrix of nodes (input)
     * @param positions Positions in segments (input)
     */
    void removeLowestScoreLink(const vector<vector<Node*> >& segment,
                               const vector<unsigned int>& position);

    /**
     * Remove the link with the lowest estimated score
     * @param segment Matrix of nodes (input)
     * @param positions Positions in segments (input)
     */
    void removeLowestLowerBoundScoreLink(const vector<vector<Node*> >& segment,
                                         const vector<unsigned int>& position);

    /**
     * Remove the 'longest' link
     * @param segment Matrix of nodes (input)
     * @param positions Positions in segments (input)
     */
    void removeLongestLink(const vector<vector<Node*> >& segment,
                           const vector<unsigned int>& position);

    /**
     * Remove a random active conflicting link from the segment matrix
     * @param segment Matrix of nodes (input)
     * @param positions Positions in segments (input)
     */
    void removeRandomACLink(const vector<vector<Node*> >& segment,
                            const vector<unsigned int>& position);

    /**
     * Remove a random conflicting link from the segment matrix
     * @param segment Matrix of nodes (input)
     * @param positions Positions in segments (input)
     */
    void removeRandomCLink(const vector<vector<Node*> >& segment,
                           const vector<unsigned int>& position);

    /**
     * Remove a random link from the segment matrix
     * @param segment Matrix of nodes (input)
     * @param positions Positions in segments (input)
     */
    void removeRandomLink(const vector<vector<Node*> >& segment,
                          const vector<unsigned int>& position);

    /**
     * Return the maximum of two integers
     * @param a First integer
     * @param b Second integer
     * @return The maximum of the two integers
     */
    int max(int a, int b) {
        return (a > b) ? a : b;
    }

    /**
     * Return the maximum of two doubles
     * @param a First double
     * @param b Second double
     * @return The maximum of the two doubles
     */
    double max(double a, double b) {
        return (a > b) ? a : b;
    }

    /**
     * Return the minimum of two doubles
     * @param a First double
     * @param b Second double
     * @return The minimum of the two doubles
     */
    double min(double a, double b) {
        return (a < b) ? a : b;
    }

    /**
     * Return true if two doubles can be considered equal
     * @param a First double
     * @param b Second double
     */
    bool doubleequal(double a, double b) {
        return fabs(a - b) <= 1e-10*max(fabs(a), fabs(b));
    }

    uint maxGaps;        // maximal number of gaps to insert
    uint searchDepth;    // maximum search depth for paths
    int rRobin;         // round robin
    std::vector<std::vector<int> > alignings;
    const vector<vector<Node* > >& nodes;
    void (GGAligner2::*removeLink)(const vector<vector<Node*> >& , const vector<unsigned int>& );

    std::vector<std::vector<int> > connect; // connectivity matrix

    friend std::ostream& operator<<(std::ostream& os, const GGAligner2& a);
};

#endif
