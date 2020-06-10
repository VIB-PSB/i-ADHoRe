#include "GGAligner2.h"

#include "Node.h"
#include <climits>
#include "AlignmentException.h"

#include "../util.h"

#include <cassert>
#include <algorithm>

double minSetTime = 0.0;
double confTime = 0.0;
double initConfl = 0.0;
double maxFlow = 0.0;
double buildCompact = 0.0;

using namespace std;

GGAligner2::GGAligner2 (int gap_, const vector<vector<Node*> > &nodes_,
                        ConflictHeuristic conflHeuristic) :
    maxGaps(gap_), nodes(nodes_)
{
    searchDepth = 30;

    // set the function pointer to the method of link removal
    switch (conflHeuristic) {
        case CH_RA:
            removeLink = &GGAligner2::removeRandomLink;
            break;
        case CH_RC:
            removeLink = &GGAligner2::removeRandomCLink;
            break;
        case CH_RAC:
            removeLink = &GGAligner2::removeRandomACLink;
            break;
        case CH_LL:
            removeLink = &GGAligner2::removeLongestLink;
            break;
        case CH_LLBS:
            removeLink = &GGAligner2::removeLowestLowerBoundScoreLink;
            break;
        case CH_LS:
            removeLink = &GGAligner2::removeLowestScoreLink;
            break;
    };
}

void GGAligner2::recMarkConnectivity(const vector<vector<Node*> >& segment,
                                     const vector<unsigned int>& position,
                                     const Node *node, vector<bool> &visited,
                                     bool movedLeft)
{
    uint nodeSeg = node->getSegmentId();
    uint nodePos = node->getOrder() - position[nodeSeg];

    // snap to the rightmost point if we are beyond searchdepth
    if (nodePos >= searchDepth) {
        nodePos = searchDepth - 1;
        node = segment[nodeSeg][position[nodeSeg] + searchDepth - 1];
        movedLeft = true;
    }

    // code: 1 = direct connect with E, 2 = indirect connect with E
    int code = (movedLeft) ? 2 : 1;

    // if the node has already been marked, get out
    if (connect[nodeSeg][nodePos] >= code) return;

    // if not, mark it now!
    connect[nodeSeg][nodePos] = code;

    // depth-first, so first handle all links originating from this node
    set<NodeLink>::const_iterator it;
    for (it = node->getLinkBegin(); it != node->getLinkEnd(); it++) {
        const Node *link = it->getNode();
        int linkSeg = link->getSegmentId();
        if (visited[linkSeg]) continue;

        visited[linkSeg] = true;
        recMarkConnectivity(segment, position, link, visited, movedLeft);
        visited[linkSeg] = false;
    }

    // then, handle all nodes left of this node
    for (int j = nodePos - 1; j >= 0; j--) {
        const Node *loop = segment[nodeSeg][position[nodeSeg] + j];
        recMarkConnectivity(segment, position, loop, visited, true);
    }
}

void GGAligner2::recMarkDirConnectivity(const vector<vector<Node*> >& segment,
                                        const vector<unsigned int>& position,
                                        const Node *node, vector<bool> &visited)
{
    uint nodeSeg = node->getSegmentId();
    uint nodePos = node->getOrder() - position[nodeSeg];

    // return if we are beyond searchdepth
    if (nodePos >= searchDepth) return;

    // code: 1 = direct connect
    const int code = 1;

    // if the node has already been visited, get out
    if (connect[nodeSeg][nodePos] == code) return;

    // if not, mark it now!
    connect[nodeSeg][nodePos] = code;

    // handle all links originating from this node
    set<NodeLink>::const_iterator it;
    for (it = node->getLinkBegin(); it != node->getLinkEnd(); it++) {
        const Node *link = it->getNode();
        int linkSeg = link->getSegmentId();
        if (visited[linkSeg]) continue;

        visited[linkSeg] = true;
        recMarkDirConnectivity(segment, position, link, visited);
        visited[linkSeg] = false;
    }
}

bool GGAligner2::hasActiveConnectivity(const vector<vector<Node*> >& segment,
                                       const vector<unsigned int>& position,
                                       const Node *node, vector<bool> &visited,
                                       int targetSegment, int targetPos,
                                       bool movedLeft)
{
    uint nodeSeg = node->getSegmentId();
    uint nodePos = node->getOrder() - position[nodeSeg];

    // snap to the rightmost point if we are beyond searchdepth
    if (nodePos >= searchDepth) {
        nodePos = searchDepth - 1;
        node = segment[nodeSeg][position[nodeSeg] + searchDepth - 1];
        movedLeft = true;
    }

    // code: 1 = direct connect with E, 2 = indirect connect with E
    int code = (movedLeft) ? 2 : 1;

    // if the node has already been visited, get out
    if (connect[nodeSeg][nodePos] >= code)
        return false;

    // if not, mark it now!
    connect[nodeSeg][nodePos] = code;

    // if we reached the destination, no need to look further
    if (nodeSeg == targetSegment) {
        assert(visited[nodeSeg]);
        if (nodePos > targetPos)
            return true;
        if (nodePos == targetPos)
            return (code == 2);
        return false;
    }

    // depth-first, so first handle all links originating from this node
    set<NodeLink>::const_iterator it;
    for (it = node->getLinkBegin(); it != node->getLinkEnd(); it++) {
        const Node *link = it->getNode();
        int linkSeg = link->getSegmentId();
        int linkPos = link->getOrder() - position[linkSeg];
        // ensure the link is active!
        if ((linkPos != 0) && (nodePos != 0)) continue;
        if (visited[linkSeg]) continue;

        visited[linkSeg] = true;
        if (hasActiveConnectivity(segment, position, link, visited,
                                  targetSegment, targetPos, movedLeft))
            return true;
        visited[linkSeg] = false;
    }

    // then, handle all nodes left of this node
    for (int j = nodePos - 1; j >= 0; j--) {
        const Node *loop = segment[nodeSeg][position[nodeSeg] + j];
        if (hasActiveConnectivity(segment, position, loop, visited,
                                  targetSegment, targetPos, true))
            return true;
    }

    return false;
}

void GGAligner2::checkForActiveConflicts(const vector<vector<Node*> >& segment,
                                         const vector<unsigned int>& position,
                                         LinkReport &report)
{
    Node &startNode = *report.startNode;
    Node &endNode = *report.endNode;

    uint startSegment = startNode.getSegmentId();
    uint endSegment = endNode.getSegmentId();
    uint endPos = endNode.getOrder() - position[endSegment];

    vector<bool> visited(segment.size(), false);

    bool tooFar = endPos >= searchDepth;

    // if the endNode is beyond the search depth, it's no use to search for
    // active st-conflicts.  Assume they are present and exit.
    if (tooFar) {
        report.activeConflict = true;
        return;
    }

    // check for active st-conflicts
    resetConnect();
    visited[endSegment] = true;
    if (hasActiveConnectivity(segment, position, &endNode, visited, startSegment, 0, false)) {
        report.activeConflict = true;
        return;
    }
    visited[endSegment] = false;

    // check for active ts-conflicts
    resetConnect();
    visited[startSegment] = true;
    if (hasActiveConnectivity(segment, position, &startNode, visited, endSegment, endPos, false)) {
        report.activeConflict = true;
        return;
    }
    visited[startSegment] = false;
}

void GGAligner2::checkForConflicts(const vector<vector<Node*> >& segment,
                                   const vector<unsigned int>& position,
                                   LinkReport &report)
{
    Node &startNode = *report.startNode;
    Node &endNode = *report.endNode;

    uint startSegment = startNode.getSegmentId();
    uint endSegment = endNode.getSegmentId();
    uint endPos = endNode.getOrder() - position[endSegment];

    vector<bool> visited(segment.size(), false);

    bool tooFar = endPos >= searchDepth;

    // if the endNode is beyond the search depth, it's no use to search for
    // st-conflicts.  Assume they are present and exit.
    if (tooFar) {
        report.conflict = true;
        return;
    }

    // check for active st-conflicts
    resetConnect();
    visited[endSegment] = true;
    recMarkConnectivity(segment, position, &endNode, visited, false);
    visited[endSegment] = false;

    if (connect[startSegment][0] == 2) {
        report.conflict = true;
        return;
    }

    // check for ts-conflicts
    resetConnect();
    visited[startSegment] = true;
    recMarkConnectivity(segment, position, &startNode, visited, false);
    visited[startSegment] = false;

    if (connect[endSegment][endPos] == 2)
        report.conflict = true;
}

void GGAligner2::makeACReports(const vector<vector<Node*> >& segment,
                               const vector<unsigned int>& position,
                               vector<LinkReport> &ACReports)
{
    Util::startChrono();

    for (unsigned int i = 0; i < segment.size(); i++) {
        if (position[i] >= segment[i].size()) continue;
        Node &startNode = *segment[i][position[i]];

        // calculate the start node flow (= sum of all start node link weights)
        double startNodeFlow = 0.0;
        set<NodeLink>::const_iterator li;
        for (li = startNode.getLinkBegin(); li != startNode.getLinkEnd(); li++)
            startNodeFlow += li->getWeight();

        // loop over all possible end nodes
        set<NodeLink>::const_iterator it = startNode.getLinkBegin();
        for ( ; it != startNode.getLinkEnd(); it++) {
            Node &endNode = *it->getNode();
            uint endSegment = endNode.getSegmentId();
            uint endIndex = endNode.getOrder();

            // don't count aligned links twice
            if (endIndex == position[endSegment])
                if (endSegment < i)
                    continue;

            LinkReport report(&startNode, &endNode);
            // check if the link is involved in an active conflict
            checkForActiveConflicts(segment, position, report);
            if (!report.activeConflict)
                continue;

            // calculate the end node flow
            double endNodeFlow = 0.0;
            for (li = endNode.getLinkBegin(); li != endNode.getLinkEnd(); li++)
                endNodeFlow += li->getWeight();

            report.ownWeight = it->getWeight();
            report.length = endIndex - position[endSegment];

            // calculate the upper limit to the direct score
            report.fdmax = min(startNodeFlow, endNodeFlow);

            // calculate the upper limit to the blocking st-flow
            for (uint j = position[endSegment]; j < endIndex; j++) {
                Node &node = *segment[endSegment][j];

                for (li = node.getLinkBegin(); li != node.getLinkEnd(); li++)
                    report.fstmax += li->getWeight();
            }
            report.fstmax += endNodeFlow;

            // calculate the upper limit to the blocking ts-flow
            report.ftsmax = startNodeFlow - it->getWeight();

            // calculate the estimated score (lower bound)
            report.estLinkScore = it->getWeight() - max(report.fstmax, report.ftsmax);

            ACReports.push_back(report);
        }
    }

    initConfl += Util::stopChrono();
}

void GGAligner2::makeAllReports(const vector<vector<Node*> >& segment,
                                const vector<unsigned int>& position,
                                vector<LinkReport> &ICReports,
                                vector<LinkReport> &RCReports,
                                vector<LinkReport> &NCReports)
{
    Util::startChrono();

    for (unsigned int i = 0; i < segment.size(); i++) {
        if (position[i] >= segment[i].size()) continue;
        Node &startNode = *segment[i][position[i]];

        // calculate the start node flow (= sum of all start node link weights)
        double startNodeFlow = 0.0;
        set<NodeLink>::const_iterator li;
        for (li = startNode.getLinkBegin(); li != startNode.getLinkEnd(); li++)
            startNodeFlow += li->getWeight();

        // loop over all possible end nodes
        set<NodeLink>::const_iterator it = startNode.getLinkBegin();
        for ( ; it != startNode.getLinkEnd(); it++) {
            Node &endNode = *it->getNode();
            uint endSegment = endNode.getSegmentId();
            uint endIndex = endNode.getOrder();

            // don't count aligned links twice
            if (endIndex == position[endSegment])
                if (endSegment < i)
                    continue;

            // calculate the end node flow
            double endNodeFlow = 0.0;
            for (li = endNode.getLinkBegin(); li != endNode.getLinkEnd(); li++)
                endNodeFlow += li->getWeight();

            // make a report
            LinkReport report(&startNode, &endNode);
            report.ownWeight = it->getWeight();
            report.length = endIndex - position[endSegment];

            // calculate the upper limit to the direct score
            report.fdmax = min(startNodeFlow, endNodeFlow);

            // calculate the upper limit to the blocking st-flow
            for (uint j = position[endSegment]; j < endIndex; j++) {
                Node &node = *segment[endSegment][j];

                for (li = node.getLinkBegin(); li != node.getLinkEnd(); li++)
                    report.fstmax += li->getWeight();
            }
            report.fstmax += endNodeFlow;

            // calculate the upper limit to the blocking ts-flow
            report.ftsmax = startNodeFlow - it->getWeight();

            // calculate the estimated score (lower bound)
            report.estLinkScore = it->getWeight() - max(report.fstmax, report.ftsmax);

            // check if the link is involved in an active conflict
            checkForActiveConflicts(segment, position, report);
            if (report.activeConflict)
                ICReports.push_back(report);
            else {
                // check if the link is involved in a conflict
                checkForConflicts(segment, position, report);

                if (report.conflict)
                    RCReports.push_back(report);
                else
                    NCReports.push_back(report);
            }
        }
    }

    initConfl += Util::stopChrono();
}

double GGAligner2::findAugmentingPath(vector<CNode> &cNodes, vector<bool> &visited,
                                      int nodeSeg, int nodePos, double flow,
                                      int endSegment, int endPos, bool movedRight)
{
    // if we reach the endSegment, stop searching
    if (nodeSeg == endSegment) {
        if ((movedRight) || (nodePos < endPos))
            return flow;
        else
            return -1.0;  // we found a direct path
    }

    // depth-first, so first handle all links originating from this node
    set<CLink>::iterator it = cNodes[nodeSeg].links.lower_bound(CLink(nodePos, 0, 0, 0.0));
    for ( ; it != cNodes[nodeSeg].links.end(); it++) {
        CLink &link = const_cast<CLink&>(*it);

        // don't visit the same segment twice
        if (visited[link.eSeg]) continue;

        // don't follow a saturated link
        double newFlow = min(flow, link.getCf());
        if (newFlow <= 0.0) continue;

        // check if we've moved right
        if (it->sIndex > nodePos)
            movedRight = true;

        visited[link.eSeg] = true;
        double result = findAugmentingPath(cNodes, visited, link.eSeg, link.eIndex,
                                           newFlow, endSegment, endPos, movedRight);
        visited[link.eSeg] = false;

        if (result > 0.0) {
            link.flow += result;
            link.inverseLink->flow -= result;
            return result;
        }
    }

    // no augmenting path found
    return -1.0;
}

double GGAligner2::findDirAugmentingPath(vector<CNode> &cNodes,
                                         vector<bool> &visited, int nodeSeg,
                                         int nodePos, double flow,
                                         int endSegment, int endPos)
{
    // if we reach the endSegment, stop searching
    if (nodeSeg == endSegment) {
        if (nodePos == endPos)
            return flow;
        else
            return -1.0;
    }

    set<CLink>::iterator it = cNodes[nodeSeg].links.lower_bound(CLink(nodePos, 0, 0, 0.0));
    for ( ; it != cNodes[nodeSeg].links.end(); it++) {
        CLink &link = const_cast<CLink&>(*it);

        // don't visit the same segment twice
        if (visited[link.eSeg]) continue;

        // don't follow a saturated link
        double newFlow = min(flow, link.getCf());
        if (newFlow <= 0.0) continue;

        // don't move right (these links shouldn't be present anyway)
        if (link.sIndex > nodePos)
            return -1.0;

        // recursively follow the link
        visited[link.eSeg] = true;
        double result = findDirAugmentingPath(cNodes, visited, link.eSeg,
                                              link.eIndex, newFlow,
                                              endSegment, endPos);
        visited[link.eSeg] = false;

        if (result > 0.0) {
            link.flow += result;
            link.inverseLink->flow -= result;
            return result;
        }
    }

    // no augmenting path found
    return -1.0;
}

void GGAligner2::calcDirectFlow(const vector<vector<Node*> >& segment,
                                const vector<unsigned int>& position,
                                LinkReport &report)
{
    Node &startNode = *report.startNode;
    Node &endNode = *report.endNode;

    uint startSegment = startNode.getSegmentId();
    uint endSegment = endNode.getSegmentId();
    uint endIndex = endNode.getOrder();
    uint endPos = endIndex - position[endSegment];

    if (endPos >= searchDepth) {
        report.fdirect = report.ownWeight;
        return;
    }

    vector<bool> visited(segment.size(), false);
    vector<CNode> cNodes(segment.size());

    // 1) calculate the direct st-flow
    Util::startChrono();

    resetConnect();
    visited[startSegment] = true;
    recMarkDirConnectivity(segment, position, &startNode, visited);
    visited[startSegment] = false;

    // build the compact graph
    for (uint i = 0; i < segment.size(); i++) {
        for (uint j = 0; j < searchDepth; j++) {
            if (connect[i][j] == 0) continue;
            Node &node = *segment[i][position[i]+j];

            set<NodeLink>::const_iterator it;
            for (it = node.getLinkBegin(); it != node.getLinkEnd(); it++) {
                const Node &link = *it->getNode();
                uint linkSeg = link.getSegmentId();
                uint linkPos = link.getOrder() - position[linkSeg];

                if (linkPos >= searchDepth) continue;
                if (linkSeg < i) continue; // don't add links twice
                if (connect[linkSeg][linkPos] == 0) continue;

                pair<set<CLink>::iterator, bool> P1, P2;
                P1 = cNodes[i].links.insert(CLink(j, linkPos, linkSeg, it->getWeight()));
                P2 = cNodes[linkSeg].links.insert(CLink(linkPos, j, i, it->getWeight()));

                if (P1.second) {
                    assert(P2.second);
                    CLink &L1 = const_cast<CLink&>(*P1.first);
                    CLink &L2 = const_cast<CLink&>(*P2.first);
                    L1.inverseLink = &L2;
                    L2.inverseLink = &L1;
                }
            }
        }
    }
    buildCompact += Util::stopChrono();

    Util::startChrono();
    while (true) {
        visited[startSegment] = true;
        double extraFlow = findDirAugmentingPath(cNodes, visited, startSegment,
                                                 0, 1e8, endSegment, endPos);
        if (extraFlow <= 0) {
            maxFlow += Util::stopChrono();
            return;
        }

        report.fdirect += extraFlow;
    }
}

void GGAligner2::calcBlockingSTFlow(const vector<vector<Node*> >& segment,
                                    const vector<unsigned int>& position,
                                    LinkReport &report)
{
    Node &startNode = *report.startNode;
    Node &endNode = *report.endNode;

    int startSegment = startNode.getSegmentId();
    int endSegment = endNode.getSegmentId();
    int endIndex = endNode.getOrder();
    int endPos = endIndex - position[endSegment];

    vector<bool> visited(segment.size(), false);
    vector<CNode> cNodes(segment.size());

    // 1) blocking st-flow
    Util::startChrono();

    resetConnect();
    visited[endSegment] = true;
    recMarkConnectivity(segment, position, &endNode, visited, false);
    visited[endSegment] = false;

    // build the compact graph
    for (uint i = 0; i < segment.size(); i++) {
        for (uint j = 0; j < searchDepth; j++) {
            if (connect[i][j] == 0) break;
            Node &node = *segment[i][position[i]+j];

            set<NodeLink>::const_iterator it;
            for (it = node.getLinkBegin(); it != node.getLinkEnd(); it++) {
                const Node &link = *it->getNode();
                uint linkSeg = link.getSegmentId();
                uint linkPos = link.getOrder() - position[linkSeg];
                if (linkPos >= searchDepth) continue;

                if (linkSeg < i) continue; // don't add links twice
                if (connect[linkSeg][linkPos] == 0) continue;

                pair<set<CLink>::iterator, bool> P1, P2;
                P1 = cNodes[i].links.insert(CLink(j, linkPos, linkSeg, it->getWeight()));
                P2 = cNodes[linkSeg].links.insert(CLink(linkPos, j, i, it->getWeight()));

                if (P1.second) {
                    assert(P2.second);
                    CLink &L1 = const_cast<CLink&>(*P1.first);
                    CLink &L2 = const_cast<CLink&>(*P2.first);
                    L1.inverseLink = &L2;
                    L2.inverseLink = &L1;
                }
            }
        }
    }
    buildCompact += Util::stopChrono();

    Util::startChrono();
    while (true) {
        visited[startSegment] = true;
        double extraFlow = findAugmentingPath(cNodes, visited, startSegment, 0,
                                              1e8, endSegment, endPos, false);
        if (extraFlow <= 0) {
            maxFlow += Util::stopChrono();
            return;
        }

        report.fstblocking += extraFlow;
    }
}

void GGAligner2::calcBlockingTSFlow(const vector<vector<Node*> >& segment,
                                    const vector<unsigned int>& position,
                                    LinkReport &report)
{
    Node &startNode = *report.startNode;
    Node &endNode = *report.endNode;

    int startSegment = startNode.getSegmentId();
    int endSegment = endNode.getSegmentId();
    int endIndex = endNode.getOrder();
    int endPos = endIndex - position[endSegment];

    vector<bool> visited(segment.size(), false);
    vector<CNode> cNodes(segment.size());

    // 1) blocking ts-flow
    Util::startChrono();

    resetConnect();
    visited[startSegment] = true;
    recMarkConnectivity(segment, position, &startNode, visited, false);
    visited[startSegment] = false;

    // build the compact graph
    for (uint i = 0; i < segment.size(); i++) {
        for (uint j = 0; j < searchDepth; j++) {
            if (connect[i][j] == 0) break;
            Node &node = *segment[i][position[i]+j];

            set<NodeLink>::const_iterator it;
            for (it = node.getLinkBegin(); it != node.getLinkEnd(); it++) {
                const Node &link = *it->getNode();
                uint linkSeg = link.getSegmentId();
                uint linkPos = link.getOrder() - position[linkSeg];
                if (linkPos >= searchDepth) continue;

                if (linkSeg < i) continue; // don't add links twice
                if (connect[linkSeg][linkPos] == 0) continue;

                pair<set<CLink>::iterator, bool> P1, P2;
                P1 = cNodes[i].links.insert(CLink(j, linkPos, linkSeg, it->getWeight()));
                P2 = cNodes[linkSeg].links.insert(CLink(linkPos, j, i, it->getWeight()));

                if (P1.second) {
                    assert(P2.second);
                    CLink &L1 = const_cast<CLink&>(*P1.first);
                    CLink &L2 = const_cast<CLink&>(*P2.first);
                    L1.inverseLink = &L2;
                    L2.inverseLink = &L1;
                }
            }
        }
    }
    buildCompact += Util::stopChrono();

    Util::startChrono();
    while (true) {
        visited[endSegment] = true;
        double extraFlow = findAugmentingPath(cNodes, visited, endSegment, endPos,
                                              1e8, startSegment, 0, false);
        if (extraFlow <= 0) {
            maxFlow += Util::stopChrono();
            return;
        }

        report.ftsblocking += extraFlow;
    }
}

void GGAligner2::calcLinkScore(const vector<vector<Node*> >& segment,
                               const vector<unsigned int>& position,
                               LinkReport &report)
{
    calcDirectFlow(segment, position, report);
    calcBlockingSTFlow(segment, position, report);
    calcBlockingTSFlow(segment, position, report);

    report.linkScore = report.fdirect - abs(report.fstblocking - report.ftsblocking);
}

bool SortReports(const LinkReport& lhs, const LinkReport& rhs)
{
    return (lhs.estLinkScore < rhs.estLinkScore);
}

void GGAligner2::removeLowestScoreLink(const vector<vector<Node*> >& segment,
                                       const vector<unsigned int>& position)
{
    vector<LinkReport> ACReports, RCReports, NCReports;
    makeACReports(segment, position, ACReports);

    if (ACReports.empty())
        makeAllReports(segment, position, ACReports, RCReports, NCReports);

    vector<LinkReport> &reports = ACReports.empty() ? ACReports : ACReports;
    assert(!reports.empty());

    std::sort(reports.begin(), reports.end(), SortReports);

    LinkReport *worstReport = &reports.front();
    calcLinkScore(segment, position, *worstReport);
    int worstSeg = worstReport->startNode->getSegmentId();

    vector<LinkReport>::iterator rep = reports.begin();
    for (rep++; rep != reports.end(); rep++) {

        // branch and bound
        if (!doubleequal(rep->estLinkScore, worstReport->linkScore))
            if (rep->estLinkScore > worstReport->linkScore)
                continue;

        calcLinkScore(segment, position, *rep);

        // if the report is (significantly) better than the worst, continue
        if (!doubleequal(rep->linkScore, worstReport->linkScore))
            if (rep->linkScore > worstReport->linkScore)
                continue;

        // if the report is equal, let round robin decide
        int seg = rep->startNode->getSegmentId();
        if (doubleequal(rep->linkScore, worstReport->linkScore)) {
            int dist1 = (seg - rRobin  + nodes.size()) % nodes.size();
            if (dist1 == 0) dist1 += nodes.size();
            int dist2 = (worstSeg - rRobin  + nodes.size()) % nodes.size();
            if (dist2 == 0) dist2 += nodes.size();
            if (dist1 >= dist2) continue;
        }

        worstReport = &(*rep);
        worstSeg = seg;
    }

    rRobin = worstSeg;
    worstReport->startNode->removeLink(*worstReport->endNode);
    worstReport->endNode->removeLink(*worstReport->startNode);
}

void GGAligner2::removeLowestLowerBoundScoreLink(const vector<vector<Node*> >& segment,
                                                 const vector<unsigned int>& position)
{
    vector<LinkReport> ACReports, RCReports, NCReports;
    makeACReports(segment, position, ACReports);

    if (ACReports.empty())
        makeAllReports(segment, position, ACReports, RCReports, NCReports);

    vector<LinkReport> reports = (ACReports.empty()) ? ACReports : ACReports;
    assert(!reports.empty());

    LinkReport *worstReport = &reports.front();
    int worstSeg = worstReport->startNode->getSegmentId();

    vector<LinkReport>::iterator rep;
    for (rep = reports.begin(); rep != reports.end(); rep++) {
        // if the report is better than the worst, continue
        if (rep->estLinkScore > worstReport->estLinkScore) continue;

        // if the report is equal, let round robin decide
        int seg = rep->startNode->getSegmentId();
        if (rep->estLinkScore == worstReport->estLinkScore) {
            int dist1 = (seg - rRobin  + nodes.size()) % nodes.size();
            if (dist1 == 0) dist1 += nodes.size();
            int dist2 = (worstSeg - rRobin  + nodes.size()) % nodes.size();
            if (dist2 == 0) dist2 += nodes.size();
            if (dist1 >= dist2) continue;
        }

        worstReport = &(*rep);
        worstSeg = seg;
    }

    rRobin = worstSeg;
    worstReport->startNode->removeLink(*worstReport->endNode);
    worstReport->endNode->removeLink(*worstReport->startNode);
}

void GGAligner2::removeLongestLink(const vector<vector<Node*> >& segment,
                                   const vector<unsigned int>& position)
{
    vector<LinkReport> ACReports, RCReports, NCReports;
    makeACReports(segment, position, ACReports);

    if (ACReports.empty())
        makeAllReports(segment, position, ACReports, RCReports, NCReports);

    vector<LinkReport> reports = (ACReports.empty()) ? ACReports : ACReports;
    assert(!reports.empty());

    LinkReport *worstReport = &reports.front();
    int worstSeg = worstReport->startNode->getSegmentId();

    vector<LinkReport>::iterator rep;
    for (rep = reports.begin(); rep != reports.end(); rep++) {

        // if the report is better than the worst, continue
        if (rep->length < worstReport->length) continue;

        // if the report is equal, let round robin decide
        int seg = rep->startNode->getSegmentId();
        if (rep->length == worstReport->length) {
            int dist1 = (seg - rRobin  + nodes.size()) % nodes.size();
            if (dist1 == 0) dist1 += nodes.size();
            int dist2 = (worstSeg - rRobin  + nodes.size()) % nodes.size();
            if (dist2 == 0) dist2 += nodes.size();
            if (dist1 >= dist2) continue;
        }

        worstReport = &(*rep);
        worstSeg = seg;
    }

    rRobin = worstSeg;
    worstReport->startNode->removeLink(*worstReport->endNode);
    worstReport->endNode->removeLink(*worstReport->startNode);
}

void GGAligner2::removeRandomACLink(const vector<vector<Node*> >& segment,
                                    const vector<unsigned int>& position)
{
    vector<LinkReport> ACReports, RCReports, NCReports;
    makeACReports(segment, position, ACReports);

    if (ACReports.empty())
        makeAllReports(segment, position, ACReports, RCReports, NCReports);

    vector<LinkReport> reports = (ACReports.empty()) ? ACReports : ACReports;
    assert(!reports.empty());

    // now select a random link
    int random = rand() % reports.size();

    reports[random].startNode->removeLink(*reports[random].endNode);
    reports[random].endNode->removeLink(*reports[random].startNode);
}

void GGAligner2::removeRandomCLink(const vector<vector<Node*> >& segment,
                                   const vector<unsigned int>& position)
{
    vector<LinkReport> ACReports, RCReports, NCReports;
    makeAllReports(segment, position, ACReports, RCReports, NCReports);

    vector<LinkReport> reports = ACReports;
    reports.insert(reports.end(), RCReports.begin(), RCReports.end());
    assert(!reports.empty());

    // now select a random link
    int random = rand() % reports.size();

    reports[random].startNode->removeLink(*reports[random].endNode);
    reports[random].endNode->removeLink(*reports[random].startNode);
}

void GGAligner2::removeRandomLink(const vector<vector<Node*> >& segment,
                                  const vector<unsigned int>& position)
{
    vector<LinkReport> ACReports, RCReports, NCReports;
    makeAllReports(segment, position, ACReports, RCReports, NCReports);

    vector<LinkReport> reports = ACReports;
    reports.insert(reports.end(), RCReports.begin(), RCReports.end());
    reports.insert(reports.end(), NCReports.begin(), NCReports.end());
    assert(!reports.empty());

    // now select a random link
    int random = rand() % reports.size();

    reports[random].startNode->removeLink(*reports[random].endNode);
    reports[random].endNode->removeLink(*reports[random].startNode);
}

bool GGAligner2::findMinimalSet(const vector<vector<Node*> >& segment,
                                const vector<unsigned int>& position,
                                vector<bool> &isProcessed,
                                Node* seed, set<Node*> &processed)
{
    queue<Node*> toProcess;
    toProcess.push(seed);

    while (!toProcess.empty()) {    // while there are still nodes to process
        Node *node = toProcess.front();
        toProcess.pop();
        processed.insert(node);

        set<NodeLink>::const_iterator it;
        for (it = node->getLinkBegin(); it != node->getLinkEnd(); it++) {
            Node *link = it->getNode();
            int linkSeg = link->getSegmentId();
            if (isProcessed[linkSeg])
                return false;
            if (link != segment[linkSeg][position[linkSeg]])
                return false;

            // if the link hasn't been processed: schedule it
            if (processed.find(link) == processed.end())
                toProcess.push(link);
        }
    }

    return true;
}

void GGAligner2::align(vector<GeneList*>& unaligned_x_lists)
{
    // recover from a previous align call
    alignings.clear();
        alignings.resize(nodes.size());
    rRobin = 0;

    connect.clear();
    connect.resize(unaligned_x_lists.size(), vector<int>(searchDepth, 0));

    // remove the gaps in the unaligned list
    for (unsigned int i = 0; i < unaligned_x_lists.size(); i++)
        unaligned_x_lists[i]->removeGaps();

    // initialize auxiliary variables
    vector<bool> atEnd (nodes.size());
    for (uint i = 0; i < nodes.size(); i++)
        atEnd[i] = (nodes[i].size() ==  0);

    vector<bool> isProcessed (nodes.size(), false);
    vector<unsigned int> position (nodes.size(), 0);

    bool done = false;
    while (!done) {
        // make all non-finished segment as "not processed"
        for (unsigned int i = 0; i < position.size(); i++)
            isProcessed[i] = atEnd[i];

        // flag to see if we have processed anything
        bool processedNone = true;

        // try to process the segments
        for (unsigned int i = 0; i < position.size(); i++) {
            if (isProcessed[i]) continue;
            Node &node = *nodes[i][position[i]];

            // the node has no equals
            if (!node.hasLinks()) {
                position[i]++;
                alignings[i].push_back(-1);
                isProcessed[i] = true;
                processedNone = false;
                continue;
            }

            Util::startChrono(); // !
            set<Node*> minimalSet;
            bool foundSet = findMinimalSet(nodes, position, isProcessed,
                                           &node, minimalSet);
            minSetTime += Util::stopChrono(); // !
            if (!foundSet) continue;

            set<Node*>::const_iterator it;
            for (it = minimalSet.begin(); it != minimalSet.end(); it++) {
                int s = (*it)->getSegmentId();
                alignings[s].push_back(1);
                position[s]++;
                isProcessed[s] = true;
            }
            processedNone = false;
        }
        // if we processed nothing, there is a conflict.  We remove a
        // link in the hope to resolve the conflict and try again.
        if (processedNone) {
            Util::startChrono(); // !
            (this->*removeLink)(nodes, position);
            confTime += Util::stopChrono(); // !
            continue;
        }

        // finally, mark the segments which are at their end and
        // push a gap in the nodes that were not processed
        done = true;
        for (unsigned int i = 0; i < position.size(); i++) {
            if (atEnd[i] || (!isProcessed[i]))
                alignings[i].push_back(0);
            if (position[i] == nodes[i].size())
                atEnd[i] = true;
            else
                done = false;
        }
    }

    // uncomment this line for verbose output
    //cout << *this << endl;

    for (unsigned int i = 0; i < position.size(); i++)
        position[i] = 0;

    vector<unsigned int> gaps_inserted(nodes.size(), 0);
    vector<unsigned int> lastOrder(nodes.size(), 0);

    for (unsigned int j = 0; j < alignings[0].size(); j++) {
        for (unsigned int i = 0; i < alignings.size(); i++) {
            if (alignings[i][j] == -1) {
                position[i]++;
                continue;
            }
            if (alignings[i][j] != 1) continue;

            // search for the maximum position
            int refOrder = nodes[i][position[i]]->getOriginalOrder() + gaps_inserted[i];
            int maxOrder = refOrder;
            set<NodeLink>::const_iterator it = nodes[i][position[i]]->getLinkBegin();
            for ( ; it != nodes[i][position[i]]->getLinkEnd(); it++) {
                int segID = it->getNode()->getSegmentId();
                int order = it->getNode()->getOriginalOrder() + gaps_inserted[segID];
                maxOrder = max(maxOrder, order);
            }

            // introduce gaps in reference segment
            unaligned_x_lists[i]->introduceGaps(lastOrder[i], refOrder, maxOrder-refOrder);
            lastOrder[i] = maxOrder;
            gaps_inserted[i] += maxOrder - refOrder;

            it = nodes[i][position[i]]->getLinkBegin();
            for ( ; it != nodes[i][position[i]]->getLinkEnd(); it++) {
                int segID = it->getNode()->getSegmentId();
                uint order = it->getNode()->getOriginalOrder() + gaps_inserted[segID];

                if ((maxOrder - order) > maxGaps) {
                    throw AlignmentException("alignment failed...too many gaps in profile");
                }

                // introduce gaps in other segments
                unaligned_x_lists[segID]->introduceGaps(lastOrder[segID], order, maxOrder-order);
                lastOrder[segID] = maxOrder;
                gaps_inserted[segID] += maxOrder - order;

            }
            position[i]++;
        }
    }

    //make sure every list is equal in size -> insert gaps past the end
    unsigned int largest = 0;
    for (unsigned int i = 0; i < unaligned_x_lists.size(); i++)
        largest = max(largest, (int)unaligned_x_lists[i]->getSize());

    for (unsigned int i = 0; i < unaligned_x_lists.size(); i++) {
        uint size = unaligned_x_lists[i]->getSize();
        if (size == largest) continue;

        if ((largest - size) > maxGaps) {
            throw AlignmentException("alignment failed...too many gaps in profile");
        }

        unaligned_x_lists[i]->introduceGaps(size, size, largest - size);
    }
}

std::ostream& operator<<(std::ostream& os, const GGAligner2& a)
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

    return os;
}
