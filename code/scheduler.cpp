/*****************************************************************************/
/* Author: Jason Sauppe, David Morrison                                      */
/* Date: 2013-08-29                                                          */
/* File: scheduler.cpp                                                       */
/* Description: Runs a scheduling algorithm for NSF panels.                  */
/* Usage: ./scheduler <roomsFile> <datesFile> <panelsFile> [other options]   */
/*                                                                           */
/* Copyright (c) 2013, 2014 University of Illinois Board of Trustees.        */
/* Published under the University of Illinois/NCSA Open Source License.      */
/* See LICENSE file for more details.                                        */
/*****************************************************************************/
#include "scheduler.h"
#include "push_relabel.h"
#include "util.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <utility>
using std::pair;
#include <algorithm>
using std::find;

// Global Variables
vector<Room> rooms;
vector<Date> dates;
vector<Panel> panels;
Config conf;
vector<pair<int,int> > dateIntervals;

// TODO: Debug Greedy Solver
//  - Debug Greedy Solver
//  - Update output to streams instead of file handles

int main(int argc, const char **argv)
{
    parseArgs(argc, argv);

    // Initialize
    initializeRooms();
    initializeDates();
    initializePanels();
    initializeDateIntervals();

    // Set the initial random seed
    //srand48((long) time(NULL));
    srand48(0); // Fixed initial seed for testing

    // DEBUG:
    //printRooms();
    //printDates();
    //printPanels();

    // Prepare output stream and print initial problem information
    initializeOutputStream();
    outputProblemInfo();

    // Run feasibility pre-tests
    feasibilityPreTest();

    // Initialize the SolutionInfo object via the constructor
    SolutionInfo si;

    // Solve
    vector<Assignment> assignments = solve(si);

    // Conduct feasibility post-test on final schedule
    feasibilityPostTest(assignments);

    // Output results
    //printAssignments(assignments);
    if (strlen(conf.scheduleFile) > 0) {
        printSchedule(assignments);
    }
    if (strlen(conf.unscheduledFile) > 0) {
        printUnscheduled(assignments);
    }

    outputSolutionInfo(si);
    finalizeOutputStream();

    return 0;
}

/*****************************************************************************/
/* Options-handling functions                                                */
/*****************************************************************************/
void parseArgs(int argc, const char **argv)
{
    if (argc < 4) {
        printUsageAndExit(argv[0]);
    }
    conf.roomsFile = argv[1];
    conf.datesFile = argv[2];
    conf.panelsFile = argv[3];
    conf.outputFile = "";       // Default, can be overridden by flag
    conf.scheduleFile = "";     // Default, can be overridden by flag
    conf.unscheduledFile = "";  // Default, can be overridden by flag
    // Now handle any additional options after required inputs
    for (int i = 4; i < argc; ++i) {
        if (std::strcmp(*(argv + i), "-o") == 0) {
            ++i; conf.outputFile = *(argv + i);
        } else if (std::strcmp(*(argv + i), "-s") == 0) {
            ++i; conf.scheduleFile = *(argv + i);
        } else if (std::strcmp(*(argv + i), "-u") == 0) {
            ++i; conf.unscheduledFile = *(argv + i);
        } else if (std::strcmp(*(argv + i), "-e") == 0) {
            conf.alg = EXACT;   // Toggles exact algorithm
        } else if (std::strcmp(*(argv + i), "-g") == 0) {
            conf.alg = GREEDY;  // Toggles greedy algorithm
        } else if (std::strcmp(*(argv + i), "-r") == 0) {
            conf.alg = RANDOM;  // Toggles randomized algorithm
        } else if (std::strcmp(*(argv + i), "-n") == 0) {
            ++i; conf.RANDOM_numAttempts = atoi(*(argv + i));
        } else if (std::strcmp(*(argv + i), "-p") == 0) {
            conf.printIntervalInfo = true;
// These options are for testing purposes only...
        } else if (std::strcmp(*(argv + i), "-ic") == 0) {
            conf.ignoreRoomCapacities = true;
        } else if (std::strcmp(*(argv + i), "-i") == 0) {
            conf.ignorePanelReqs = true;
        } else {
            printUsageAndExit(argv[0]);
        }
    }
    return;
}

void printUsageAndExit(const char *argv0)
{
    fprintf(stderr,
        "Usage: %s <rooms file> <dates file> <panels file>\n"
        "       [-o <output file>]      (File to print output info (JSON))\n"
        "       [-s <schedule file>]    (File to print the schedule (CSV))\n"
        "       [-u <unscheduled file>] (File to list unsched. panels (CSV))\n"
        "       [-e]                    (Use exact scheduler (DEFAULT))\n"
        "       [-g]                    (Use greedy scheduler (IN DEVELOPMENT)\n"
        "       [-r]                    (Use random scheduler)\n"
        "       [-n <numAttempts>]      (Number of attempts for random alg.)\n"
        "       [-p]                    (Prints additional solver info)\n",
        argv0);
    fprintf(stderr, "  All input files must be in CSV format\n");
    exit(-1);
}

/*****************************************************************************/
/* Initialization functions                                                  */
/*****************************************************************************/
void initializeRooms()
{
    vector<vector<string> > roomRecords;
    readCSV(conf.roomsFile, roomRecords);
    rooms.resize(roomRecords.size() - 1); // First record is the headers
    int maxNumberOfFeatures = 0;
    for (int i = 1; i < roomRecords.size(); ++i) {
        vector<string> &roomRecord = roomRecords[i];
        if (roomRecord.size() < 3) {
            printf("ERROR: Room record %d is missing information!\n", i-1);
            exit(-1);
        }
        rooms[i-1].IID = i-1;
        rooms[i-1].EID = atoi(roomRecord[0].c_str());
        rooms[i-1].name = roomRecord[1];
        rooms[i-1].capacity = atoi(roomRecord[2].c_str());
        // Remaining elements are binary room features
        for (int j = 3; j < roomRecord.size(); ++j) {
            rooms[i-1].roomFeatures.push_back((roomRecord[j][0] == 'Y'));
        }
        if (maxNumberOfFeatures < rooms[i-1].roomFeatures.size()) {
            maxNumberOfFeatures = rooms[i-1].roomFeatures.size();
        }
    }
    // Ensure all rooms have the same number of features (initializing empty
    // values to false)
    for (int rIID = 0; rIID < rooms.size(); ++rIID) {
        if (rooms[rIID].roomFeatures.size() < maxNumberOfFeatures) {
            rooms[rIID].roomFeatures.resize(maxNumberOfFeatures, false);
        }
    }
    return;
}

void initializeDates()
{
    vector<vector<string> > dateRecords;
    readCSV(conf.datesFile, dateRecords);

    // Process the date labels
    vector<string> &dateLabels = dateRecords[0];
    dates.resize(dateLabels.size() - 1); // First record is list of dates
    for (int j = 1; j < dateLabels.size(); ++j) {
        dates[j-1].IID = j-1;
        dates[j-1].name = dateLabels[j];
        yymmdd temp = parseDate(dates[j-1].name);
        // DEBUG:
        //printf("%s is: %d [%d-%d-%d]\n", dates[j-1].name.c_str(), j-1,
        //        temp.yy, temp.mm, temp.dd);
        dates[j-1].year = temp.yy;
        dates[j-1].month = temp.mm;
        dates[j-1].day = temp.dd;
        std::stringstream ss;
        ss << std::setw(2) << std::setfill('0') << dates[j-1].month << "/";
        ss << std::setw(2) << std::setfill('0') << dates[j-1].day << "/";
        ss << std::setw(2) << std::setfill('0') << dates[j-1].year-2000;
        dates[j-1].fName = ss.str();
    }
    // The first (non-header) row contains the default room settings.
    vector<string> &defaultRoomAvailabilities = dateRecords[1];
    for (int j = 1; j < defaultRoomAvailabilities.size(); ++j) {
        // For each date, if the default availibility is yes, mark all rooms
        // as available
        if (defaultRoomAvailabilities[j][0] == 'Y') {
            for (int rIID = 0; rIID < rooms.size(); ++rIID) {
                dates[j-1].availableRooms.push_back(rIID);
            }
        }
    }
    // Now handle the room-specific availabilities
    for (int i = 2; i < dateRecords.size(); ++i) {
        vector<string> &dateRecord = dateRecords[i];
        int rEID = atoi(dateRecord[0].c_str());
        int rIID = lookupRoomIID(rEID); // Find room's IID
        // For each of the dates, make sure the room's IID appears if it
        // needs to and doesn't appear if it doesn't need to.
        for (int j = 1; j < dateRecord.size(); ++j) {
            bool roomAppearsInDate = false;
            int roomIndex = -1;
            for (int k = 0; k < dates[j-1].availableRooms.size(); ++k) {
                if (dates[j-1].availableRooms[k] == rIID) {
                    roomAppearsInDate = true;
                    roomIndex = k;
                    break;
                }
            }
            if ((roomAppearsInDate) && (dateRecord[j][0] == 'N')) {
                // Need to remove room
                dates[j-1].availableRooms[roomIndex]
                    = dates[j-1].availableRooms.back();
                dates[j-1].availableRooms.pop_back();
            } else if ((!roomAppearsInDate) && (dateRecord[j][0] == 'Y')) {
                // Need to add room
                dates[j-1].availableRooms.push_back(rIID);
            } // else no action needs to be taken
        }
    }
    return;
}

void initializePanels()
{
    vector<vector<string> > panelRecords;
    readCSV(conf.panelsFile, panelRecords);
    panels.resize(panelRecords.size() - 1); // First record is the headers
    for (int i = 1; i < panelRecords.size(); ++i) {
        vector<string> &panelRecord = panelRecords[i];
        if (panelRecord.size() < 6) {
            printf("Panel record %d is missing information.\n", i-1);
            exit(-1);
        }
        panels[i-1].IID = i-1;
        panels[i-1].EID = atoi(panelRecord[0].c_str());
        panels[i-1].name = panelRecord[1];
        panels[i-1].pd = panelRecord[2];
        panels[i-1].numberOfDays = atoi(panelRecord[3].c_str());
        panels[i-1].startingDateIID = lookupDateIID(panelRecord[4].c_str());
        panels[i-1].size = atoi(panelRecord[5].c_str());
        // Remaining elements are binary room features
        for (int j = 6; j < panelRecord.size(); ++j) {
            panels[i-1].roomRequirements.push_back((panelRecord[j][0] == 'Y'));
        }
        // Now ensure that panel has at least the same number of requirements
        // as the features of the rooms
        if (panels[i-1].roomRequirements.size() < rooms[0].roomFeatures.size()){
            panels[i-1].roomRequirements.resize(rooms[0].roomFeatures.size(),
                                                false);
        }
        // Now add panel to dates that it needs
        int sdIID = panels[i-1].startingDateIID;
        if ((sdIID < 0) || (sdIID + panels[i-1].numberOfDays > dates.size())) {
            printf("ERROR: Panel %d does not have a valid starting date\n",
                panels[i-1].EID);
            exit(-1);
        }
        for (int dIID = sdIID; dIID < sdIID+panels[i-1].numberOfDays; ++dIID) {
            dates[dIID].panelRequests.push_back(panels[i-1].IID);
        }
    }
    return;
}

void initializeDateIntervals()
{
    initializeDateIntervalsHelper(0, 0, dates.size() - 1);
    return;
}

void initializeDateIntervalsHelper(int sdIID, int cdIID, int edIID)
{
    if (cdIID == edIID) { // Base Case
        dateIntervals.push_back(std::make_pair(sdIID, edIID));
        return;
    } // else work to do
    bool hasPanelExtendingPastCurDate = false;
    for (int i = 0; i < dates[cdIID].panelRequests.size(); ++i) {
        int pIID = dates[cdIID].panelRequests[i];
        if (panels[pIID].startingDateIID+panels[pIID].numberOfDays-1 > cdIID) {
            hasPanelExtendingPastCurDate = true;
            break;
        }
    }
    if (hasPanelExtendingPastCurDate) {
        initializeDateIntervalsHelper(sdIID, cdIID + 1, edIID);
    } else {
        dateIntervals.push_back(std::make_pair(sdIID, cdIID));
        initializeDateIntervalsHelper(cdIID + 1, cdIID + 1, edIID);
    }
    return;
}

/*****************************************************************************/
/* Output Stream Functions                                                   */
/*****************************************************************************/
void initializeOutputStream()
{
    if (strlen(conf.outputFile) > 0) {
        std::ofstream *ofstream = new std::ofstream();
        ofstream->open(conf.outputFile, std::ofstream::out);
        conf.outstream = ofstream;
    } else {
        conf.outstream = &std::cout;
    }
    return;
}

void outputProblemInfo()
{
    *conf.outstream << "{";
    *conf.outstream << " \n  \"roomsFile\": \"" << conf.roomsFile << "\""
                    << ",\n  \"datesFile\": \"" << conf.datesFile << "\""
                    << ",\n  \"panelsFile\": \"" << conf.panelsFile << "\""
                    << ",\n  \"outputFile\": \"" << conf.outputFile << "\""
                    << ",\n  \"scheduleFile\": \"" << conf.scheduleFile << "\"";
    *conf.outstream << ",\n  \"problemInfo\": {"
                    << " \n    \"numRooms\": " << rooms.size()
                    << ",\n    \"numDates\": " << dates.size()
                    << ",\n    \"numPanels\": " << panels.size()
                    << ",\n    \"numIntervals\": " << dateIntervals.size();
    int numDatesWithPanels = 0;
    int numIntervalsWithPanels = 0;
    int maxPanelsPerDay = 0;
    int numPanelRequests = 0;
    for (int i = 0; i < dateIntervals.size(); ++i) {
        int panelReqsInInterval = 0;
        for (int dIID = dateIntervals[i].first;
                 dIID <= dateIntervals[i].second; ++dIID) {
            panelReqsInInterval += dates[dIID].panelRequests.size();
            if (maxPanelsPerDay < dates[dIID].panelRequests.size()) {
                maxPanelsPerDay = dates[dIID].panelRequests.size();
            }
            if (dates[dIID].panelRequests.size() > 0) { ++numDatesWithPanels; }
        }
        numPanelRequests += panelReqsInInterval;
        if (panelReqsInInterval > 0) { ++numIntervalsWithPanels; }
    }
    *conf.outstream << ",\n    \"numPanelRequests\": " << numPanelRequests
                    << ",\n    \"numDatesWithPanels\": " << numDatesWithPanels
                    << ",\n    \"numIntervalsWithPanels\": "
                    << numIntervalsWithPanels
                    << ",\n    \"maxPanelRequestsPerDay\": " << maxPanelsPerDay
                    << " \n  }";
    return;
}

void outputSolutionInfo(SolutionInfo &si)
{
    // Print stats on best assignment found.
    *conf.outstream << ",\n  \"solutionInfo\": {"
            // Info on Scheduled Panels
                    << " \n    \"numPanels\": " << panels.size()
                    << ",\n    \"scheduledPanelsGlobalUB\": "
                    << si.scheduledPanelsGlobalUB
                    << ",\n    \"numScheduledPanels\": "
                    << si.numScheduledPanels
            // Info on Satisfied Panel Requests
                    << ",\n    \"numPanelReqs\": " << si.numPanelReqs
                    << ",\n    \"satisfiedPanelReqsGlobalUB\": "
                    << si.satisfiedPanelReqsGlobalUB
                    << ",\n    \"numSatisfiedPanelReqs\": "
                    << si.numSatisfiedPanelReqs
            // Info on Interval Completeness and Optimality
                    << ",\n    \"numCompleteIntervals\": "
                    << si.numCompleteIntervals
                    << ",\n    \"numOptimalIntervals\": " << si.numOptIntervals
                    << " \n  }";
    *conf.outstream << "\n}\n\n";
    return;
}

void finalizeOutputStream()
{
    // Close and clean up the output stream
    conf.outstream->flush();
    if (strlen(conf.outputFile) > 0) {
        ((std::ofstream *) conf.outstream)->close();
        delete conf.outstream;
    }
    return;
}

/*****************************************************************************/
/* Feasibility-checking functions                                            */
/*****************************************************************************/
void feasibilityPreTest()
{
    // DEBUG:
    *conf.outstream << ",\n  \"feasibilityPreTest\": ["
                    << " \n    \"Notes\"";
    checkNumberOfPanelsPerDate();
    checkPanelRequirements();
    checkBoundsPerDate();
    *conf.outstream << " \n  ]";
    return;
}

void checkNumberOfPanelsPerDate()
{
    for (int dIID = 0; dIID < dates.size(); ++dIID) {
        if (dates[dIID].availableRooms.size() <
                dates[dIID].panelRequests.size()) {
            *conf.outstream << ",\n    \"INFEAS: Date " << dIID << " ["
                            << dates[dIID].fName << "] has "
                            << dates[dIID].availableRooms.size()
                            << " rooms available for "
                            << dates[dIID].panelRequests.size()
                            << " panel requests\"";
        }
    }
    return;
}

void checkPanelRequirements()
{
    for (int pIID = 0; pIID < panels.size(); ++pIID) {
        int sdIID = panels[pIID].startingDateIID;
        bool atLeastOneRoomSatisfiesReqs = false;
        for (int i = 0; i < dates[sdIID].availableRooms.size(); ++i) {
            int rIID = dates[sdIID].availableRooms[i];
            if ((roomSatisfiesPanelRequirements(rIID, pIID)) &&
                (roomAvailableForRequiredDays(rIID, pIID))) {
                atLeastOneRoomSatisfiesReqs = true; break;
            }
        }
        if (!atLeastOneRoomSatisfiesReqs) {
            *conf.outstream << ",\n    \"INFEAS: No rooms can satisfy the "
                            << "requirements for Panel " << panels[pIID].EID
                            << " [" << panels[pIID].name << "]\"";
        }
    }
    return;
}

void checkBoundsPerDate()
{
    int infeasSum = 0;
    for (int dIID = 0; dIID < dates.size(); ++dIID) {
        vector<Assignment> dateAssignments =
            scheduleSingleDayRequests(dates[dIID].panelRequests,
                                      dates[dIID].availableRooms);
        if (dateAssignments.size() < dates[dIID].panelRequests.size()) {
            infeasSum+=dates[dIID].panelRequests.size()-dateAssignments.size();
            *conf.outstream << ",\n    \"INFEAS: No way to schedule all "
                            << dates[dIID].panelRequests.size()
                            << " panels for Date " << dIID << " ["
                            << dates[dIID].fName << "]; Maximum possible is "
                            << dateAssignments.size() << ".\"";
        }
    }
    if (infeasSum > 0) {
        *conf.outstream << ",\n    \"INFEAS: Sum of Unsatisfiable Panel "
                        << "Requests is " << infeasSum << "\"";
    }
    return;
}

bool roomSatisfiesPanelRequirements(int rIID, int pIID)
{
    if ((!conf.ignoreRoomCapacities) &&
        (rooms[rIID].capacity < panels[pIID].size)) { return false; }
    if (conf.ignorePanelReqs) { return true; } // Stop if ignoring requirements
    for (int i = 0; i < rooms[rIID].roomFeatures.size(); ++i) {
        if (panels[pIID].roomRequirements[i] && !rooms[rIID].roomFeatures[i]) {
            // Room doesn't have feature and panel needs it
            return false;
        }
    }
    return true; // Getting to this point means room satisfies panel reqs
}

bool roomAvailableForRequiredDays(int rIID, int pIID)
{
    int sdIID = panels[pIID].startingDateIID;
    for (int dIID = sdIID; dIID < sdIID + panels[pIID].numberOfDays; ++dIID) {
        // NOTE: This can be optimized by skipping the first day, since we
        // know the room is available because we're processing it.
        bool roomAvailableOnDay = false;
        for (int i = 0; i < dates[dIID].availableRooms.size(); ++i) {
            if (dates[dIID].availableRooms[i] == rIID) {
                roomAvailableOnDay = true;
                break;
            }
        }
        if (!roomAvailableOnDay) { return false; } // Not available
    }
    return true; // Must be true at this point
}

void feasibilityPostTest(vector<Assignment> &assignments)
{
    // DEBUG:
    *conf.outstream << ",\n  \"feasibilityPostTest\": ["
                    << " \n    \"Notes\"";

    vector<vector<Assignment> > schedule = makeSchedule(assignments);

    for (int dIID = 0; dIID < dates.size(); ++dIID) {
        // Initialize vector of all rooms available on this date
        vector<int> unusedRooms;
        for (int i = 0; i < dates[dIID].availableRooms.size(); ++i) {
            unusedRooms.push_back(dates[dIID].availableRooms[i]);
        }
        // Initialize vector of all panels requested for this date
        vector<int> unscheduledPanels;
        for (int i = 0; i < dates[dIID].panelRequests.size(); ++i) {
            unscheduledPanels.push_back(dates[dIID].panelRequests[i]);
        }
        // Now iterate through the assignments for this date, removing rooms
        // from the list of total rooms and panels from the list of panel
        // requests as we see them
        for (int i = 0; i < schedule[dIID].size(); ++i) {
            int rIID = schedule[dIID][i].roomIID;
            int pIID = schedule[dIID][i].panelIID;
            // Check if room satisfies panel requirements
            if (!roomSatisfiesPanelRequirements(rIID, pIID)) {
                *conf.outstream << ",\n    \"ERROR: Room " << rooms[rIID].EID
                                << " does not satisfy the requirements of "
                                << "Panel " << panels[pIID].EID << "\"";
                break;
            }
            // Check if room used in assignment is actually available on
            // this date (and has not been used already)
            int ind = -1;
            for (int j = 0; j < unusedRooms.size(); ++j) {
                if (unusedRooms[j] == rIID) { ind = j; break; }
            }
            if (ind < 0) {
                *conf.outstream << ",\n    \"ERROR: Assignment uses "
                                << "unavailable room (" << panels[pIID].EID
                                << " -> " << rooms[rIID].EID << ")\"";
                break;
            }
            unusedRooms[ind] = unusedRooms.back();
            unusedRooms.pop_back();
            // Check if panel request in assignment was actually for this date
            // (and has not been assigned to another room on the same date)
            ind = -1;
            for (int j = 0; j < unscheduledPanels.size(); ++j) {
                if (unscheduledPanels[j] == pIID) { ind = j; break; }
            }
            if (ind < 0) {
                *conf.outstream << ",\n    \"ERROR: Assignment uses "
                                << "unrequested panel (" << panels[pIID].EID
                                << " -> " << rooms[rIID].EID << ")\"";
                break;
            }
            unscheduledPanels[ind] = unscheduledPanels.back();
            unscheduledPanels.pop_back();
        }
    }
    *conf.outstream << " \n  ]";
    return;
}

/*****************************************************************************/
/* Solver functions                                                          */
/*****************************************************************************/
vector<Assignment> solve(SolutionInfo &si)
{
    // First initialize the totals in the SolutionInfo struct
    si.numPanels = panels.size();
    for (int dIID = 0; dIID < dates.size(); ++dIID) {
        si.numPanelReqs += dates[dIID].panelRequests.size();
    }

    // DEBUG:
    *conf.outstream << ",\n  \"solver\": {";
    clock_t overallStartTime = clock();

    // Solve each of the date intervals separately, combining as we go
    vector<Assignment> assignments;
    if (conf.alg == EXACT) {
        *conf.outstream << " \n    \"algorithm\": \"Exact\""
                        << ",\n    \"details\": ["
                        << " \n      \"Notes\"";
        for (int iIID = 0; iIID < dateIntervals.size(); ++iIID) {
            append(assignments, exactlyScheduleInterval(iIID, si));
        }
    } else if (conf.alg == RANDOM) {
        *conf.outstream << " \n    \"algorithm\": \"Randomized-"
                        << conf.RANDOM_numAttempts << "\""
                        << ",\n    \"details\": ["
                        << " \n      \"Notes\"";
        for (int iIID = 0; iIID < dateIntervals.size(); ++iIID) {
            append(assignments, randomlyScheduleInterval(iIID, si));
        }
    } else if (conf.alg == GREEDY) {
        *conf.outstream << " \n    \"algorithm\": \"Greedy\""
                        << ",\n    \"details\": ["
                        << " \n      \"Notes\"";
        for (int iIID = 0; iIID < dateIntervals.size(); ++iIID) {
            append(assignments, greedilyScheduleInterval(iIID, si));
        }
    } else {
        printf("Unknown algorithm.\n");
        exit(-1);
    }

    double totalTime = ((double) (clock()-overallStartTime)) / CLOCKS_PER_SEC;

    // Update solver info
    si.numScheduledPanels = assignments.size();
    for (int i = 0; i < assignments.size(); ++i) {
        int schedPIID = assignments[i].panelIID;
        si.numSatisfiedPanelReqs += panels[schedPIID].numberOfDays;
    }

    *conf.outstream << " \n    ]"
                    << ",\n    \"totalTime\": " << totalTime
                    << " \n  }";

    return assignments;
}

/*****************************************************************************/
/* Functions for the exact scheduler                                         */
/*****************************************************************************/
vector<Assignment> exactlyScheduleInterval(int iIID, SolutionInfo &si)
{
    int sdIID = dateIntervals[iIID].first;
    int edIID = dateIntervals[iIID].second;

    // Initialize data structure for tracking branching decisions and info
    BranchInfo bi(iIID);
    for (int dIID = sdIID; dIID <= edIID; ++dIID) {
        // Identify the available rooms
        bi.remRooms.push_back(vector<int>());
        for (int i = 0; i < dates[dIID].availableRooms.size(); ++i) {
            bi.remRooms[dIID-sdIID].push_back(dates[dIID].availableRooms[i]);
        }
        // Identify the single-day and multi-day panels
        bi.remSDPanels.push_back(vector<int>());
        for (int i = 0; i < dates[dIID].panelRequests.size(); ++i) {
            int pIID = dates[dIID].panelRequests[i];
            if (panels[pIID].startingDateIID == dIID) {
                bi.numPanelsInInterval += 1;
                bi.numPanelReqsInInterval += panels[pIID].numberOfDays;
                if (panels[pIID].numberOfDays > 1) {
                    bi.remMDPanels.push_back(pIID);
                } else {
                    bi.remSDPanels[dIID-sdIID].push_back(pIID);
                }
            }
        }
    }
    bi.numRemPanels = bi.numPanelsInInterval;
    bi.numRemPanelReqs = bi.numPanelReqsInInterval;
    bi.lookupTables.resize(edIID - sdIID + 1);
    bi.intervalStartTime = clock();

    // Compute an upper bound on the number of panel requests that can be
    // satisfied in this interval, and then use that to compute an upper bound
    // on the number of panels that can be scheduled
    int intervalUB_SatPR = computeSatisfiedPanelReqsBound(bi);
    int intervalUB_SchP = computeScheduledPanelsBound(bi, intervalUB_SatPR);

    si.satisfiedPanelReqsGlobalUB += intervalUB_SatPR;
    si.scheduledPanelsGlobalUB += intervalUB_SchP;

    // Print out information on the bounds
    if (conf.printIntervalInfo) {
        *conf.outstream << ",\n      \"EXACT: Beginning Interval " << bi.iIID
                        << "; Computing bounds...\"";
        *conf.outstream << ",\n      \"       Total Panels: "
                        << bi.numPanelsInInterval << "; UB schedules "
                        << intervalUB_SchP << " panels\"";
        *conf.outstream << ",\n      \"       Total Panel Reqs: "
                        << bi.numPanelReqsInInterval << "; UB satisfies "
                        << intervalUB_SatPR << " panel reqs\"";
    }

    // Construct a schedule
    vector<Assignment> intAssignments = ext_scheduleForInterval(bi);

    // If time limit was hit, try the random solver
    if ((bi.exceededTimeLimit) && (conf.EXACT_useRandomAfterExact) &&
        (intAssignments.size() < bi.numPanelsInInterval)) {
        vector<Assignment> rndAssignments = randomlyScheduleInterval(iIID, si);
        if (rndAssignments.size() > intAssignments.size()) {
            intAssignments.clear();
            append(intAssignments, rndAssignments);
        }
    }

    // Check the termination status of the best solution found
    char termStatus = '?';
    if (intAssignments.size() == bi.numPanelsInInterval) {
        // Interval is complete and therefore optimal
        si.numCompleteIntervals += 1;
        si.numOptIntervals += 1;
        termStatus = 'C';
    } else if (!bi.exceededTimeLimit) {
        // Interval is incomplete, but provably optimal
        si.numOptIntervals += 1;
        termStatus = 'O';
    }

    // DEBUG:
    if (conf.printIntervalInfo) {
        *conf.outstream << ",\n      \"EXACT: Interval " << std::setw(3)
                        << iIID << " [" << std::setw(3) << sdIID << ","
                        << std::setw(3) << edIID << "], #SD=" << std::setw(3)
                        << bi.numPanelsInInterval - bi.remMDPanels.size()
                        << ", #MD=" << std::setw(3) << bi.remMDPanels.size()
                        << ", Scheduled = " << std::setw(3)
                        << intAssignments.size() << "/" << std::setw(3)
                        << bi.numPanelsInInterval
                        << ", Status: " << termStatus << "\"";
    }

    return intAssignments;
}

vector<Assignment> ext_scheduleForInterval(BranchInfo &bi)
{
    int sdIID = dateIntervals[bi.iIID].first;
    int edIID = dateIntervals[bi.iIID].second;

    // Check if we just have single day panels left
    if (bi.remMDPanels.size() == 0) {
        // Just solve each of the individual dates optimally now
        vector<Assignment> sdAssignments;
        for (int dIID = sdIID; dIID <= edIID; ++dIID) {
            int i = dIID - sdIID;
            append(sdAssignments,
                   ext_scheduleSDPanels(dIID, bi.remSDPanels[i],
                                        bi.remRooms[i], bi.lookupTables[i]));
        }
        // Terminal node, so update global incumbent if necessary
        int objValue = sdAssignments.size() + bi.numScheduledPanels;
        if (objValue > bi.bestSolNumScheduledPanels) {
            bi.bestSolNumScheduledPanels = objValue;
            if (conf.printIntervalInfo) {
                *conf.outstream << ",\n      \"EXACT: Updated the incumbent "
                                << "for interval " << bi.iIID << " to "
                                << bi.bestSolNumScheduledPanels
                                << " scheduled panels";
            }
        }
        return sdAssignments;
    } // Else we have a multi-day panel to schedule

    // Track the best set of assignments for this interval
    vector<Assignment> bestIntAssignments;

    // Check for early termination via bounds
    if (conf.EXACT_computeBounds) {
        int partialSolUB_SatPR = computeSatisfiedPanelReqsBound(bi);
        int partialSolUB_SchP =
                computeScheduledPanelsBound(bi, partialSolUB_SatPR);
        // Use the bound on number of scheduled panels for pruning
        if (partialSolUB_SchP <= bi.bestSolNumScheduledPanels) {
            return bestIntAssignments; // Empty vector
        }
    }

    // Identify a random multi-day panel for scheduling
    //int rnd = randomGen.IRandom(0, bi.remMDPanels.size() - 1);
    int rnd = randomInt(0, bi.remMDPanels.size() - 1);
    int pIID = bi.remMDPanels[rnd];
    int psdIID = panels[pIID].startingDateIID;

    // Remove multi-day panel from list
    bi.remMDPanels[rnd] = bi.remMDPanels.back();
    bi.remMDPanels.pop_back();

    // Identify the compatible rooms
    vector<int> potentialRoomIIDs =
        computePotentialRoomsForPanel(pIID, bi.iIID, bi.remRooms);

    // Order the potential rooms to prioritize those most likely to
    // lead to a fully satisfied schedule. Right now we're just using the
    // capacities as sorting criteria, to leave as many high-capacity
    // rooms as we can. Will want to look at incorporating the room
    // features into the sorting process, though.
    if (conf.EXACT_useRoomSort) {
        sort(potentialRoomIIDs.begin(), potentialRoomIIDs.end(), RoomSort);
    }

    // Assume that we can schedule this panel, and update branch tracking vars
    bi.numRemPanels -= 1;
    bi.numRemPanelReqs -= panels[pIID].numberOfDays;
    bi.numScheduledPanels += 1;
    bi.numSatisfiedPanelReqs += panels[pIID].numberOfDays;

    // Now try assigning panel to each possible room
    double time = ((double) (clock() - bi.intervalStartTime)) / CLOCKS_PER_SEC;
    for (int i = 0; i < potentialRoomIIDs.size(); ++i) {
        int rIID = potentialRoomIIDs[i];

        // Assign panel to room
        vector<Assignment> tempAssignments;
        tempAssignments.push_back(Assignment(pIID, rIID));

        // Remove room from vector of remaining rooms
        for (int j = 0; j < panels[pIID].numberOfDays; ++j) {
            int dateInd = psdIID - sdIID + j;
            int rInd = vec_index(bi.remRooms[dateInd], rIID);
            bi.remRooms[dateInd][rInd] = bi.remRooms[dateInd].back();
            bi.remRooms[dateInd].pop_back();
        }

        // Now finish scheduling the remaining panels
        append(tempAssignments, ext_scheduleForInterval(bi));

        // Add room back to vector of remaining rooms
        for (int j = 0; j < panels[pIID].numberOfDays; ++j) {
            int dateInd = psdIID - sdIID + j;
            bi.remRooms[dateInd].push_back(rIID);
        }

        // Update best assignments if this set is better (here better is
        // defined using the number of scheduled panels)
        if (tempAssignments.size() > bestIntAssignments.size()) {
            bestIntAssignments.clear();
            append(bestIntAssignments, tempAssignments);
        }

        // Check for early termination if we successfully scheduled all
        // remaining panels. (No need to try other room assignments for the
        // currently selected panel.)
        if (bestIntAssignments.size() == bi.numRemPanels) { break; }

        // Check for exceeded time limit, in which case we stop here
        time = ((double) (clock() - bi.intervalStartTime)) / CLOCKS_PER_SEC;
        if (time > conf.EXACT_timeLimitOneBranch) {
            bi.exceededTimeLimit = true; break;
        }
    }
    // Now revert the tracking variables after trying to schedule this panel
    bi.numRemPanels += 1;
    bi.numRemPanelReqs += panels[pIID].numberOfDays;
    bi.numScheduledPanels -= 1;
    bi.numSatisfiedPanelReqs -= panels[pIID].numberOfDays;

    // If we were unable to schedule all remaining panels, and in fact
    // failed to schedule two or more panels, then we can try to drop the
    // current multi-day panel and schedule everything else.
    // Only do this if we still have time remaining.
    // (NOTE: bestIntAssignments should have this multi-day panel in it, so we
    // compare to the numRemPanels at the start of this function, not the
    // updated version assuming we scheduled the selected MD panel)
    if ((bestIntAssignments.size() < bi.numRemPanels-1) &&
        (time <= conf.EXACT_timeLimitZeroBranch)) {
        // Ignore this panel, but note that we no longer schedule it
        bi.numRemPanels -= 1;
        bi.numRemPanelReqs -= panels[pIID].numberOfDays;

        vector<Assignment> tempAssignments = ext_scheduleForInterval(bi);

        // Check for better assignment
        if (tempAssignments.size() > bestIntAssignments.size()) {
            bestIntAssignments.clear();
            append(bestIntAssignments, tempAssignments);
        }

        // Now update tracking vars to include this MD panel again
        bi.numRemPanels += 1;
        bi.numRemPanelReqs += panels[pIID].numberOfDays;
    }
    // else we either successfully scheduled everything, or only failed to
    // schedule one panel (in which case dropping this panel wouldn't
    // allow us to improve the total number of panels scheduled)

    // Add multi-day panel back to list
    bi.remMDPanels.push_back(pIID);

    return bestIntAssignments;
}

vector<Assignment> ext_scheduleSDPanels(int dIID, vector<int> &remSDPanels,
                                        vector<int> &remRooms,
                                        rem_rooms_map &lookupTable)
{
    // Optimally schedules the single-day panels for the given date (dIID),
    // using only the rooms given in the remRooms vector.

    // First check if we have this recorded in the memoization table
    if (conf.EXACT_useMemoization) {
        auto pos = lookupTable.find(remRooms);
        if (!(pos == lookupTable.end())) {
            return pos->second; // Return memoized solution
        }
    } // else compute the solution and store it for the first time

    // Construct the optimal schedule
    vector<Assignment> singleDayAssignments =
        scheduleSingleDayRequests(remSDPanels, remRooms);

    // Store assignments in memoization table
    if (conf.EXACT_useMemoization) {
        lookupTable[remRooms] = singleDayAssignments;
    }

    return singleDayAssignments;
}

int computeSatisfiedPanelReqsBound(BranchInfo &bi)
{
    int sdIID = dateIntervals[bi.iIID].first;
    int edIID = dateIntervals[bi.iIID].second;

    int numRemReqsSatisfied = 0;
    for (int dIID = sdIID; dIID <= edIID; ++dIID) {
        // For each date in the interval, create single-day requests for each
        // of the multi-day panels on that date
        vector<int> panelReqsForDate;
        for (int i = 0; i < bi.remMDPanels.size(); ++i) {
            int mdpIID = bi.remMDPanels[i];
            int mdpSDIID = panels[mdpIID].startingDateIID;
            int mdpND = panels[mdpIID].numberOfDays;
            for (int j = 0; j < mdpND; ++j) {
                if (mdpSDIID + j == dIID) {
                    panelReqsForDate.push_back(mdpIID);
                    break;
                }
            }
        }
        // Now add the single day panels for this date
        for (int i = 0; i < bi.remSDPanels[dIID-sdIID].size(); ++i) {
            panelReqsForDate.push_back(bi.remSDPanels[dIID-sdIID][i]);
        }
        // Now schedule all of the remaining ''single day'' panels optimally
        // using the remaining rooms
        vector<Assignment> sdAssignments =
          scheduleSingleDayRequests(panelReqsForDate, bi.remRooms[dIID-sdIID]);
        numRemReqsSatisfied += sdAssignments.size();
    }
    return bi.numSatisfiedPanelReqs + numRemReqsSatisfied;
}

int computeScheduledPanelsBound(BranchInfo &bi, int satisfiedPanelReqsBound)
{
    int sdIID = dateIntervals[bi.iIID].first;
    int edIID = dateIntervals[bi.iIID].second;

    int numRemScheduledPanels = 0;
    int numRemSatisfiedReqs = satisfiedPanelReqsBound;
    for (int dIID = sdIID; dIID <= edIID; ++dIID) {
        // For each of the dates in the interval, assume that all of the
        // single-day panels can be scheduled and remove their contribution
        // to the total satisfied panel requests bound
        numRemScheduledPanels += bi.remSDPanels[dIID-sdIID].size();
        numRemSatisfiedReqs -= bi.remSDPanels[dIID-sdIID].size();
        if (numRemSatisfiedReqs <= 0) {
            // We potentially over-estimated how many panels we could
            // schedule, so make that correction here
            numRemScheduledPanels += numRemSatisfiedReqs;
            break;
        }
    }
    if (numRemSatisfiedReqs > 0) { // Need to handle the multi-day panels
        // First, count the number of multi-day panels with each length
        vector<int> numMultiDayPanelsByLength;
        for (int i = 0; i < bi.remMDPanels.size(); ++i) {
            int lengthMDPi = panels[bi.remMDPanels[i]].numberOfDays;
            while (numMultiDayPanelsByLength.size() <= lengthMDPi - 2) {
                numMultiDayPanelsByLength.push_back(0);
            }
            numMultiDayPanelsByLength[lengthMDPi - 2] += 1;
        }
        // Now start adding panels to the number of scheduled panels while
        // removing their contribution to the total satisfied panel
        // requests bound
        for (int i = 0; i < numMultiDayPanelsByLength.size(); ++i) {
            int curLength = i + 2;
            if (curLength > numRemSatisfiedReqs) { break; }
            for (int j = 0; j < numMultiDayPanelsByLength[i]; ++j) {
                if (curLength <= numRemSatisfiedReqs) {
                    numRemScheduledPanels += 1;
                    numRemSatisfiedReqs -= curLength;
                } else {
                    break;
                }
            }
        }
    }

    return bi.numScheduledPanels + numRemScheduledPanels;
}

/*****************************************************************************/
/* Functions for the randomized scheduler                                    */
/*****************************************************************************/
vector<Assignment> randomlyScheduleInterval(int iIID, SolutionInfo &si)
{
    int sdIID = dateIntervals[iIID].first;
    int edIID = dateIntervals[iIID].second;

    // Compute total panels in interval for early termination option
    int panelsInInterval = 0;
    int numMultiDayPanels = 0;
    for (int dIID = sdIID; dIID <= edIID; ++dIID) {
        for (int i = 0; i < dates[dIID].panelRequests.size(); ++i) {
            if (panels[dates[dIID].panelRequests[i]].startingDateIID == dIID) {
                ++panelsInInterval;
                if (panels[dates[dIID].panelRequests[i]].numberOfDays > 1) {
                    ++numMultiDayPanels;
                }
            }
        }
    }

    // Perform a number of random restarts for this interval
    vector<Assignment> bestIntAssignments;
    for (int i = 0; i < conf.RANDOM_numAttempts; ++i) {
        vector<Assignment> intAssignments = rnd_scheduleForInterval(iIID);
        if (intAssignments.size() > bestIntAssignments.size()) {
            bestIntAssignments.clear();
            append(bestIntAssignments, intAssignments);
        }
        // Check for early termination possibility if we satisfy all panels
        if (bestIntAssignments.size() == panelsInInterval) { break; }
    }

    // DEBUG:
    if (conf.printIntervalInfo) {
        *conf.outstream << ",\n      \" RAND: Interval " << std::setw(3)
                        << iIID << " [" << std::setw(3) << sdIID << ","
                        << std::setw(3) << edIID << "], #SD=" << std::setw(3)
                        << panelsInInterval - numMultiDayPanels
                        << ", #MD=" << std::setw(3) << numMultiDayPanels
                        << ", Scheduled = " << std::setw(3)
                        <<  bestIntAssignments.size() << "/" << std::setw(3)
                        << panelsInInterval << "\"";
    }

    // Done with random scheduling for interval, so return the best assignments
    return bestIntAssignments;
}

vector<Assignment> rnd_scheduleForInterval(int iIID)
{
    int sdIID = dateIntervals[iIID].first;
    int edIID = dateIntervals[iIID].second;

    // Initialize the vectors for this interval
    vector<int> remDateIIDs;
    vector<vector<int> > intRemPanels;
    vector<vector<int> > intRemRooms;
    for (int dIID = sdIID; dIID <= edIID; ++dIID) {
        remDateIIDs.push_back(dIID);
        intRemPanels.push_back(vector<int>());
        intRemRooms.push_back(vector<int>());
        for (int i = 0; i < dates[dIID].panelRequests.size(); ++i) {
            intRemPanels[dIID-sdIID].push_back(dates[dIID].panelRequests[i]);
        }
        for (int i = 0; i < dates[dIID].availableRooms.size(); ++i) {
            intRemRooms[dIID-sdIID].push_back(dates[dIID].availableRooms[i]);
        }
    }

    // DEBUG:
    //printRemPanelsForInterval(iIID, intRemPanels);
    //printRemRoomsForInterval(iIID, intRemRooms);

    // Schedule panels for each date in the interval, one by one
    vector<Assignment> intAssignments;
    while (remDateIIDs.size() > 0) {
        // Select a random date from the ones remaining
        //int rnd = randomGen.IRandom(0, remDateIIDs.size() - 1);
        int rnd = randomInt(0, remDateIIDs.size() - 1);
        int dateIID = remDateIIDs[rnd];

        // Construct a schedule for the date
        append(intAssignments,
               rnd_scheduleForDate(dateIID, iIID, intRemPanels, intRemRooms));

        // Now remove date from further consideration
        remDateIIDs[rnd] = remDateIIDs.back();
        remDateIIDs.pop_back();
    }

    return intAssignments;
}

vector<Assignment> rnd_scheduleForDate(int dIID, int iIID,
                                       vector<vector<int> > &intRemPanels,
                                       vector<vector<int> > &intRemRooms)
{
    int sdIID = dateIntervals[iIID].first;
    int edIID = dateIntervals[iIID].second;

    vector<int> remPanelIIDs;
    for (int i = 0; i < intRemPanels[dIID-sdIID].size(); ++i) {
        remPanelIIDs.push_back(intRemPanels[dIID-sdIID][i]);
    }

    vector<Assignment> dateAssignments;
    while (remPanelIIDs.size() > 0) {
        // Select a random panel from the ones remaining
        //int rnd = randomGen.IRandom(0, remPanelIIDs.size() - 1);
        int rnd = randomInt(0, remPanelIIDs.size() - 1);
        int nextPanelIID = remPanelIIDs[rnd];

        Assignment panelAssignment = rnd_schedulePanel(nextPanelIID, iIID,
                                            intRemPanels, intRemRooms);
        if (panelAssignment.roomIID >= 0) {
            dateAssignments.push_back(panelAssignment);
        } // else panel was not assigned

        // Remove panel IID from further consideration
        remPanelIIDs[rnd] = remPanelIIDs.back();
        remPanelIIDs.pop_back();
    }
    // Any panels still belonging to this date's remPanels vector after this
    // loop finishes are unable to be scheduled. The algorithm may still
    // attempt to schedule a multi-day panel that lies on this day when it
    // selects an adjacent day for scheduling, but it will be unsuccessful
    // then as well.

    return dateAssignments;
}

Assignment rnd_schedulePanel(int pIID, int iIID,
                             vector<vector<int> > &intRemPanels,
                             vector<vector<int> > &intRemRooms)
{
    int sdIID = dateIntervals[iIID].first;
    int edIID = dateIntervals[iIID].second;
    int psdIID = panels[pIID].startingDateIID;

    // Find all rooms that can satisfy this panel request. Note that this
    // panel was selected because one of its requested dates was chosen;
    // however, we do not need to start scheduling it based on that particular
    // date. Instead, we just want to schedule it across the number of days it
    // requires.
    vector<int> potentialRoomIIDs =
        computePotentialRoomsForPanel(pIID, iIID, intRemRooms);

    // If any potential rooms are available, pick one at random and assign it
    // TODO: This can be further optimized by picking the potential room with
    // the smallest set of features (e.g., smallest size), thus ensuring that
    // future panels have a more diverse set of rooms available. That is, we
    // should look at the demand for rooms with various features among the
    // remaining panels and rooms, and use that information to intelligently
    // choose a room for this panel.
    if (potentialRoomIIDs.size() > 0) {
        //int rnd = randomGen.IRandom(0, potentialRoomIIDs.size() - 1);
        int rnd = randomInt(0, potentialRoomIIDs.size() - 1);
        int rIID = potentialRoomIIDs[rnd];

        // Remove panel and rooms from remaining panels and rooms lists for
        // the appropriate days
        for (int i = 0; i < panels[pIID].numberOfDays; ++i) {
            int dateInd = psdIID - sdIID + i;
            int pInd = vec_index(intRemPanels[dateInd], pIID);
            int rInd = vec_index(intRemRooms[dateInd], rIID);
            if ((pInd < 0) || (rInd < 0)) {
                printf("ERROR: Missing pIID or rIID from vector\n"); exit(-1);
            }
            intRemPanels[dateInd][pInd] = intRemPanels[dateInd].back();
            intRemPanels[dateInd].pop_back();
            intRemRooms[dateInd][rInd] = intRemRooms[dateInd].back();
            intRemRooms[dateInd].pop_back();
        }

        return Assignment(pIID, rIID);
    } else { // there are no rooms available that can host this panel
        return Assignment(pIID, -1);
    }
}

/*****************************************************************************/
/* Greedy Solver Functions                                                   */
/*****************************************************************************/
vector<Assignment> greedilyScheduleInterval(int iIID, SolutionInfo &si)
{
    int sdIID = dateIntervals[iIID].first;
    int edIID = dateIntervals[iIID].second;

    // Initialize the available rooms by date
    vector<vector<int> > remRooms;
    for (int dIID = sdIID; dIID <= edIID; ++dIID) {
        remRooms.push_back(vector<int>());
        for (int i = 0; i < dates[dIID].availableRooms.size(); ++i) {
            remRooms[dIID-sdIID].push_back(dates[dIID].availableRooms[i]);
        }
    }

    // Now greedily schedule the panels starting on each date of the interval
    vector<Assignment> assignments;
    for (int dIID = sdIID; dIID <= edIID; ++dIID) {
        // Only look at requests that start on this day
        vector<int> lPanels;
        for (int i = 0; i < dates[dIID].panelRequests.size(); ++i) {
            int pIID = dates[dIID].panelRequests[i];
            if (panels[pIID].startingDateIID == dIID) {
                lPanels.push_back(pIID);
            }
        }
        // Skip this date if we have no panels or no rooms available
        if ((lPanels.size() == 0) || (remRooms[dIID-sdIID].size() == 0)) {
            continue;
        }

        // TODO: DEBUG THIS
        // Need to make sure that rooms are available for the dates in question
        // This is NOT working properly at the moment...
        vector<Assignment> dateAssignments =
            scheduleSingleDayRequests(lPanels, remRooms[dIID-sdIID]);

        // Update the remaining rooms on subsequent days given the assignments
        // for this date
        for (int i = 0; i < dateAssignments.size(); ++i) {
            assignments.push_back(dateAssignments[i]); // Add to list
            int pIID = assignments.back().panelIID;
            int rIID = assignments.back().roomIID;
            for (int l = 0; l < panels[pIID].numberOfDays; ++l) {
                int rrInd = panels[pIID].startingDateIID - sdIID + l;
                // Remove the room from the available rooms on each of the
                // days that the panel requires
                for (int j = 0; j < remRooms[rrInd].size(); ++j) {
                    if (rIID == remRooms[rrInd][j]) {
                        remRooms[rrInd][j] = remRooms[rrInd].back();
                        remRooms[rrInd].pop_back();
                        break;
                    }
                }
            }
        }
    }
    return assignments;
}

/*****************************************************************************/
/* Helper Functions for Solving Strategies                                   */
/*****************************************************************************/
// This function assumes that $lPanels is a list of one-day panels to schedule,
// and $lRooms is the list of rooms that are available for scheduling
vector<Assignment> scheduleSingleDayRequests(const vector<int>& lPanels,
                                             const vector<int>& lRooms)
{
    // NOTE: I get a glibc error when this check is omitted.
    if ((lPanels.size() == 0) || (lRooms.size() == 0)) {
        return vector<Assignment>();
    }

    // Compute the list of available rooms that are valid for each request
    // validRooms[i] is a list of rooms that meet the requirements for panel i
	//
	// TODO - there may be some rooms that are invalid for all requests.  Right
	// now we're just passing these in to push_relabel.  It's probably more
	// efficient to not pass them in, but it makes indexing more complicated,
	// and I don't see this being a big bottleneck...
    vector<vector<int>> validRooms(lPanels.size());
    for (int i = 0; i < lPanels.size(); ++i)
    {
        for (int j = 0; j < lRooms.size(); ++j)
        {
            // Internally, the push-relabel code uses room IDs
            // [0..lRooms.size()) and request IDs [0..lPanels.size()), so just
            // push_back indices
            if (roomSatisfiesPanelRequirements(lRooms[j], lPanels[i]))
                validRooms[i].push_back(j);
        }
    }

    // Run push-relabel to compute an optimal schedule
    vector<pair<int, int>> localAssignments =
            pushRelabel(lPanels.size(), lRooms.size(), validRooms);

    // Convert from push-relabel IDs to internal IDs
    vector<Assignment> assignments;
    for (int i = 0; i < localAssignments.size(); ++i) {
        assignments.push_back(Assignment(lPanels[localAssignments[i].first],
                                         lRooms[localAssignments[i].second]));
    }

    return assignments;
}

vector<int> computePotentialRoomsForPanel(int pIID, int iIID,
                                          vector<vector<int> > &intRemRooms)
{
    int sdIID = dateIntervals[iIID].first;
    int edIID = dateIntervals[iIID].second;
    int psdIID = panels[pIID].startingDateIID;

    // Look at all the rooms that are still available on the first day of
    // this panel. If a room satisfies the size and feature requirements of
    // the panel, and if that room is still available on all of the days the
    // panel needs, then it can be added to the list of potential rooms.
    vector<int> potentialRoomIIDs;
    for (int i = 0; i < intRemRooms[psdIID-sdIID].size(); ++i) {
        int remRoomIID = intRemRooms[psdIID-sdIID][i];
        if (roomSatisfiesPanelRequirements(remRoomIID, pIID)) {
            // Check if the room is available for subsequent days
            bool roomIsAvailable = true;
            for (int j = 1; j < panels[pIID].numberOfDays; ++j) {
                int dateInd = psdIID - sdIID + j;
                if (!contains(intRemRooms[dateInd], remRoomIID)) {
                    roomIsAvailable = false; break;
                }
            }
            if (roomIsAvailable) {
                potentialRoomIIDs.push_back(remRoomIID);
            }
        }
    }
    return potentialRoomIIDs;
}

/*****************************************************************************/
/* Print / output functions                                                  */
/*****************************************************************************/
void printRooms()
{
    printf("*** Rooms ***\n");
    for (int rIID = 0; rIID < rooms.size(); ++rIID) {
        printf("Room %d [%s], Cap=%d, Features: ", rooms[rIID].EID,
            rooms[rIID].name.c_str(), rooms[rIID].capacity);
        for (int i = 0; i < rooms[rIID].roomFeatures.size(); ++i) {
            printf("%d", (rooms[rIID].roomFeatures[i]) ? 1 : 0);
        }
        printf("\n");
    }
    return;
}

void printDates()
{
    printf("*** Dates ***\n");
    for (int dIID = 0; dIID < dates.size(); ++dIID) {
        printf("Date %d [%s], #RoomsAvail=%d, #PanelReqs=%d\n",
            dIID, dates[dIID].fName.c_str(),
            dates[dIID].availableRooms.size(),
            dates[dIID].panelRequests.size());
    }
    return;
}

void printPanels()
{
    printf("*** Panels ***\n");
    for (int pIID = 0; pIID < panels.size(); ++pIID) {
        printf("Panel %d [%s], [%s], Len=%d, StartDateID=%d, "
               "Size=%d, Requirements: ", panels[pIID].EID,
            panels[pIID].name.c_str(), panels[pIID].pd.c_str(),
            panels[pIID].numberOfDays, panels[pIID].startingDateIID,
            panels[pIID].size);
        for (int i = 0; i < panels[pIID].roomRequirements.size(); ++i) {
            printf("%d", (panels[pIID].roomRequirements[i]) ? 1 : 0);
        }
        printf("\n");
    }
    return;
}

void printAssignments(vector<Assignment> &assignments)
{
    printf("*** Assignments ***\n");
    for (int i = 0; i < assignments.size(); ++i) {
        int rIID = assignments[i].roomIID;
        int pIID = assignments[i].panelIID;
        int sdIID = panels[pIID].startingDateIID;
        printf("Panel %d [%s] scheduled in Room %d [%s] for day(s) ",
                panels[pIID].EID, panels[pIID].name.c_str(),
                rooms[rIID].EID, rooms[rIID].name.c_str());
        for (int dIID = sdIID; dIID < sdIID+panels[pIID].numberOfDays; ++dIID){
            printf("%d [%s], ", dIID, dates[dIID].fName.c_str());
        }
        printf("\n");
    }
    return;
}

void printSchedule(vector<Assignment> &assignments)
{
    // TODO: Revise this to print output using streams
    vector<vector<Assignment> > schedule = makeSchedule(assignments);
    FILE *FOS = fopen(conf.scheduleFile, "w");

    fprintf(FOS, "Room Name,");
    for (int dIID = 0; dIID < dates.size(); ++dIID) {
        //fprintf(FOS, "%d [%s],", dIID, dates[dIID].fName.c_str());
        fprintf(FOS, "%s,", dates[dIID].fName.c_str());
    }
    fprintf(FOS, "\n");

    fprintf(FOS, "Total Panels,");
    for (int dIID = 0; dIID < dates.size(); ++dIID) {
        fprintf(FOS, "%d,", schedule[dIID].size());
    }
    fprintf(FOS, "\n");

    for (int rIID = 0; rIID < rooms.size(); ++rIID) {
        //fprintf(FOS, "\"%d [%s]\",", rooms[rIID].EID, rooms[rIID].name.c_str());
        fprintf(FOS, "\"%s\",", rooms[rIID].name.c_str());
        for (int dIID = 0; dIID < dates.size(); ++dIID) {
            for (int i = 0; i < schedule[dIID].size(); ++i) {
                int assnRoomIID = schedule[dIID][i].roomIID;
                if (assnRoomIID == rIID) {
                    int pIID = schedule[dIID][i].panelIID;
                    //fprintf(FOS, "\"%d [%s]\"", panels[pIID].EID,
                    //                           panels[pIID].name.c_str());
                    //fprintf(FOS, "\"%d\"", panels[pIID].EID);
                    fprintf(FOS, "\"%s\"", panels[pIID].name.c_str());
                    break;
                }
            }
            fprintf(FOS, ",");
        }
        fprintf(FOS, "\n");
    }
    fclose(FOS);
    return;
}

void printUnscheduled(vector<Assignment> &assignments)
{
    FILE *FOU = fopen(conf.unscheduledFile, "w");
    fprintf(FOU, "Unscheduled,\n");
    for (int pIID = 0; pIID < panels.size(); ++pIID) {
        bool foundPIID = false;
        // TODO: Make this search more efficient...
        for (int i = 0; i < assignments.size(); ++i) {
            if (assignments[i].panelIID == pIID) { foundPIID = true; break; }
        }
        if (!foundPIID) {
            fprintf(FOU, "%d,%s,%s,%d,%s,%d",
                panels[pIID].EID, panels[pIID].name.c_str(),
                panels[pIID].pd.c_str(), panels[pIID].numberOfDays,
                dates[panels[pIID].startingDateIID].fName.c_str(),
                panels[pIID].size);
            for (int j = 0; j < panels[pIID].roomRequirements.size(); ++j) {
                fprintf(FOU, ",%c", panels[pIID].roomRequirements[j]?'Y':'N');
            }
            fprintf(FOU, "\n");
        }
    }
    fprintf(FOU, "\n");
    fclose(FOU);
    return;
}

vector<vector<Assignment> > makeSchedule(vector<Assignment> &assignments)
{
    vector<vector<Assignment> > schedule;
    schedule.resize(dates.size());
    for (int i = 0; i < assignments.size(); ++i) {
        int rIID = assignments[i].roomIID;
        int pIID = assignments[i].panelIID;
        int sdIID = panels[pIID].startingDateIID;
        // For each of the days that the panel is scheduled
        for (int dIID = sdIID; dIID < sdIID+panels[pIID].numberOfDays; ++dIID){
            // Add the <room,panel> pair to the schedule
            schedule[dIID].push_back(Assignment(pIID, rIID));
        }
    }
    return schedule;
}

void printIntervalRemRooms(int iIID, vector<vector<int> > &intRemRooms)
{
    int sdIID = dateIntervals[iIID].first;
    int edIID = dateIntervals[iIID].second;
    for (int i = 0; i < intRemRooms.size(); ++i) {
        printf(" Rem Rooms for Date IID %4d: ", sdIID + i);
        for (int j = 0; j < intRemRooms[i].size(); ++j) {
            printf("%d, ", intRemRooms[i][j]);
        }
        printf("\n");
    }
    return;
}

void printIntervalRemPanels(int iIID, vector<vector<int> > &intRemPanels)
{
    int sdIID = dateIntervals[iIID].first;
    int edIID = dateIntervals[iIID].second;
    for (int i = 0; i < intRemPanels.size(); ++i) {
        printf("Rem Panels for Date IID %4d: ", sdIID + i);
        for (int j = 0; j < intRemPanels[i].size(); ++j) {
            printf("%d, ", intRemPanels[i][j]);
        }
        printf("\n");
    }
    return;
}

/*****************************************************************************/
/* Miscellaneous helper functions                                            */
/*****************************************************************************/
// Generates uniform random integer i with a <= i <= b
int randomInt(int a, int b)
{
    return ((int) floor(drand48() * (b - a + 1.0))) + a ;
}

/*****************************************************************************/
/* I/O helper functions                                                      */
/*****************************************************************************/
void readCSV(const char *inFile, vector<vector<string> > &records)
{
    std::ifstream ifs;
    ifs.open(inFile);
    if (!ifs) {
        printf("Unable to open file: %s\n", inFile);
        exit(-1);
    }
    while (ifs.good()) {
        string record;
        getline(ifs, record);
        if (record.find_first_not_of(' ') != std::string::npos) {
            records.push_back(vector<string>());
            std::stringstream ss(record);
            while (ss.good()) {
                string entry;
                getline(ss, entry, ',');
                // Check if entry starts with a double-quote, in which case we
                // need to be careful about commas inside the string
                if ((entry.length() > 0) && (entry[0] == '"')) {
                    if (entry[entry.length() - 1] == '"') { // Ends with "
                        entry = entry.substr(1, entry.length()-2);
                    } else { // Quoted string with commas inside
                        string rem;
                        getline(ss, rem, '"');
                        entry = entry.substr(1, entry.length()-1) + "," + rem;
                        getline(ss, rem, ','); // Remove the extra comma
                    }
                }
                records.back().push_back(entry);
            }
        }
    }
    ifs.close();
    return;
}

yymmdd parseDate(string dateString)
{
    // TODO: Expand on the formats accepted here...
    // DEBUG:
    //printf("Parsing date: %s\n", dateString.c_str()); fflush(stdout);
    yymmdd temp;
    if (dateString.find('/') != string::npos) {
        // Assume format is "MM/DD/YY" or "M/D/Y"
        int ind1 = dateString.find('/');
        temp.mm = atoi(dateString.substr(0, ind1).c_str());
        string rem = dateString.substr(ind1+1);
        int ind2 = rem.find('/');
        temp.dd = atoi(rem.substr(0, ind2).c_str());
        temp.yy = atoi(rem.substr(ind2+1).c_str());
        if (temp.yy < 100) { temp.yy += 2000; }
    } else { // Assume format is "DoW, MMM DD, YY" or "MMM DD, YY"
        std::stringstream ss(dateString);
        string entry;
        getline(ss, entry, ',');
        int sInd = 0;
        if (entry.length() <= 3) { // Ignore day of week, if it is present
            getline(ss, entry, ',');
            sInd = 1;
        }
        string month = entry.substr(sInd, sInd + 2);
        string day = entry.substr(5);
        getline(ss, entry);
        string year = entry.substr(1,3);
        temp.yy = atoi(year.c_str()) + 2000;
        temp.mm = monthNameToNumber(month);
        temp.dd = atoi(day.c_str());
    }
    // DEBUG:
    //printf("%s = [%d|%d|%d]\n", dateString.c_str(),
    //          temp.yy, temp.mm, temp.dd);
    return temp;
}

int monthNameToNumber(string monthName)
{
    if      (monthName.compare("Jan") == 0) { return  1; }
    else if (monthName.compare("Feb") == 0) { return  2; }
    else if (monthName.compare("Mar") == 0) { return  3; }
    else if (monthName.compare("Apr") == 0) { return  4; }
    else if (monthName.compare("May") == 0) { return  5; }
    else if (monthName.compare("Jun") == 0) { return  6; }
    else if (monthName.compare("Jul") == 0) { return  7; }
    else if (monthName.compare("Aug") == 0) { return  8; }
    else if (monthName.compare("Sep") == 0) { return  9; }
    else if (monthName.compare("Oct") == 0) { return 10; }
    else if (monthName.compare("Nov") == 0) { return 11; }
    else if (monthName.compare("Dec") == 0) { return 12; }
    else                                    { return -1; }
}

int lookupDateIID(string dateString)
{
    yymmdd temp = parseDate(dateString);
    // TODO: This can be optimized further using a binary search type
    // strategy, if dates are assumed to be ordered.
    for (int dIID = 0; dIID < dates.size(); ++dIID) {
        if ((dates[dIID].year == temp.yy) &&
            (dates[dIID].month == temp.mm) &&
            (dates[dIID].day == temp.dd)) {
            //printf("%s matches [%d|%d|%d] (date %d)\n", dateString.c_str(),
            //        dates[i].year, dates[i].month, dates[i].day, i);
            return dIID;
        }
    }
    // If we get here, date wasn't found
    printf("ERROR: No date found matching %s\n", dateString.c_str());
    exit(-1);
    return -1;
}

int lookupRoomIID(int rEID)
{
    for (int rIID = 0; rIID < rooms.size(); ++rIID) {
        if (rooms[rIID].EID == rEID) { return rIID; }
    }
    // If we get to here, room wasn't found
    printf("ERROR: No room found with ID %d\n", rEID);
    exit(-1);
    return -1;
}
