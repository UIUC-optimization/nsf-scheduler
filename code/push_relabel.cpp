/*****************************************************************************/
/* Author: David R. Morrison                                                 */
/* Date: Oct. 2013                                                           */
/* File: push_relabel.cpp                                                    */
/* Description: Implementation details for the push-relabel network flow     */
/*  algorithm.                                                               */
/*                                                                           */
/* Copyright (c) 2013, 2014 University of Illinois Board of Trustees.        */
/* Published under the University of Illinois/NCSA Open Source License.      */
/* See LICENSE file for more details.                                        */
/*****************************************************************************/
#include "push_relabel.h"
#include "util.h"

#include <vector>
#include <utility>
#include <list>
#include <algorithm>

using namespace std;

namespace PushRelabel
{

// Parameters
int numNodes, source, sink;
int reqStart, reqEnd;
int roomStart, roomEnd;

// Tracking data structures
vector<vector<int>> capacity;
vector<vector<int>> flow;
vector<int> height;
vector<int> excess;
vector<int> nextNode;
list<int> nodeQueue;

// Push flow from u to v
void push(int u, int v)
{
    // Determine the amount of flow we can push from u to v
    int toPush = min(excess[u], capacity[u][v] - flow[u][v]);
    flow[u][v] += toPush;
    flow[v][u] -= toPush;

    // Treat values of MaxInt as positive/negative infinity, so don't update
    if (excess[u] != MaxInt) excess[u] -= toPush;
    if (excess[v] != -MaxInt) excess[v] += toPush;
}

// Increase the height of node u to be one larger than a neighbor with residual
// capacity and the smallest height larger than u's
void relabel(int u)
{
    int minHeight = MaxInt;
    for (int v = 0; v < numNodes; ++v)
    {
        if (capacity[u][v] - flow[u][v] > 0 && height[v] < minHeight)
        {
            minHeight = height[v];
            height[u] = minHeight + 1;
        }
    }
}

// Try to push flow from u to all of its neighbors until all u's excess is
// satisfied.  If u still has excess at the end, increase the height of u and
// start over.
void discharge(int u)
{
    while (excess[u] > 0)
    {
        if (nextNode[u] < numNodes)
        {
            int v = nextNode[u];
            if ((capacity[u][v] - flow[u][v] > 0) && (height[u] > height[v]))
                push(u, v);
            else nextNode[u] += 1;
        }
        else
        {
            relabel(u);
            nextNode[u] = 0;
        }
    }
}

}; // namespace PushRelabel

vector<pair<int, int>> pushRelabel(int numRequests, int numRooms,
        const vector<vector<int>>& validRooms)
{
    using namespace PushRelabel;

    // Set up the graph parameters -- the source is the first vertex, the sink
    // is the last vertex, we list request vertices first and then room
    // vertices next
    numNodes = numRequests + numRooms + 2;
    source = 0;
    sink = numNodes - 1;
    reqStart = 1;
    reqEnd = reqStart + numRequests;
    roomStart = reqEnd;
    roomEnd = roomStart + numRooms;

    // Empty data structures to remove any old information
    for (int i = 0; i < capacity.size(); ++i) {
        capacity[i].clear();
    }
    capacity.clear();
    for (int i = 0; i < flow.size(); ++i) {
        flow[i].clear();
    }
    flow.clear();
    height.clear();
    excess.clear();
    nextNode.clear();
    nodeQueue.clear();

    // Initialize data structures
    capacity.resize(numNodes, vector<int>(numNodes, 0));
    flow.resize(numNodes, vector<int>(numNodes, 0));
    height.resize(numNodes, 0);
    excess.resize(numNodes, 0);
    nextNode.resize(numNodes, 0);

    // Create a bipartite graph connecting requests to rooms, with an
    // artificial source/sink
    // Can only select one room for each request, so capacity of all edges is 1
    for (int i = reqStart; i != reqEnd; ++i)
        capacity[source][i] = 1;
    for (int i = roomStart; i != roomEnd; ++i)
        capacity[i][sink] = 1;
    for (int i = 0; i < validRooms.size(); ++i)
    {
        int reqId = i + reqStart;
        for (int j = 0; j < validRooms[i].size(); ++j)
        {
            // This is where the problem lies. validRooms[i][j] needs to be an
            // index from 0..lRooms.size()-1, but we're putting in the IIDs of
            // the rooms instead.
            int roomId = validRooms[i][j] + roomStart;
            capacity[reqId][roomId] = 1;
        }
    }

    // Create a list of untouched vertices
    for (int i = 0; i < numNodes; ++i)
        if((i != source) && (i != sink))
            nodeQueue.push_back(i);

    // The source is as high as possible, and we push as much flow as we can
    // from it to its neighbors (in this case, just all of the requests)
    height[source] = numNodes;
    excess[source] = MaxInt;
    for (int i = reqStart; i != reqEnd; ++i)
        push(source, i);

    // While there's still some vertex we haven't visited...
    auto iter = nodeQueue.begin();
    while (iter != nodeQueue.end())
    {
        int u = *iter;
        int oldHeight = height[u];

        // Try to satisfy this node's excess by pushing flow to all neighbors
        discharge(u);

        // If we increased this node's height in the last discharge operation,
        // add it back to the front and discharge more flow
        if (height[u] > oldHeight)
        {
            nodeQueue.erase(iter);
            nodeQueue.push_front(u);
            iter = nodeQueue.begin();
        }
        else ++iter;
    }

    // Loop through all room/request pairs to compute the final schedule
    vector<pair<int, int>> assignments;
    for (int i = reqStart; i != reqEnd; ++i)
    {
        bool satisfied = false;
        for (int j = roomStart; j != roomEnd; ++j)
        {
            // A request is scheduled to a room if it has positive flow on it
            if (flow[i][j] > 0)
            {
                if (!satisfied)
                {
                    assignments.push_back({i-reqStart, j-roomStart});
                    satisfied = true;
                }
                else throw ERROR << "Request " << i << " has already been satisfied";
            }
        }
    }

    return assignments;
}

