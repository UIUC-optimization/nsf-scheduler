/*****************************************************************************/
/* Author: Jason Sauppe, David Morrison                                      */
/* Date: 2013-08-29                                                          */
/* File: scheduler.h                                                         */
/* Description: The header file for the scheduling program.                  */
/*                                                                           */
/* Copyright (c) 2013, 2014 University of Illinois Board of Trustees.        */
/* Published under the University of Illinois/NCSA Open Source License.      */
/* See LICENSE file for more details.                                        */
/*****************************************************************************/
#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <string>
using std::string;
#include <vector>
using std::vector;
#include <unordered_map>
using std::unordered_map;
#include <iostream>
#include <utility>

enum Algorithm { EXACT, RANDOM, GREEDY };

/*****************************************************************************/
/* Structs used in scheduler.cpp                                             */
/*****************************************************************************/
struct Config {
    Config() {
        // Default Settings
        printIntervalInfo = false;
        ignoreRoomCapacities = false;
        ignorePanelReqs = false;
        attemptShiftScheduling = false;
        checkStartDateDoubleBooking = false;
        panelsIncludeOrganizer = false;

        alg = EXACT;

        EXACT_computeBounds = true;
        EXACT_useMemoization = true;
        EXACT_useRoomSort = true;
        EXACT_useRandomAfterExact = true;
        EXACT_useSDExactAfterExact = false;
        EXACT_timeLimitOneBranch = 5;
        EXACT_timeLimitZeroBranch = 10;

        RANDOM_numAttempts = 50;
    }

    // TODO: Maybe change these to strings?
    // Input and output files
    const char *roomsFile;
    const char *datesFile;
    const char *organizersFile;
    const char *panelsFile;
    const char *outputFile;
    const char *scheduleFile;
    const char *unscheduledFile;

    // Output stream
    std::ostream *outstream;

    // Other global settings
    bool printIntervalInfo;
    bool ignoreRoomCapacities;
    bool ignorePanelReqs;
    bool attemptShiftScheduling;
    bool checkStartDateDoubleBooking;
    bool panelsIncludeOrganizer;

    // Algorithm Settings
    Algorithm alg;
    bool EXACT_computeBounds;
    bool EXACT_useMemoization;
    bool EXACT_useRoomSort;
    bool EXACT_useRandomAfterExact;
    bool EXACT_useSDExactAfterExact;
    double EXACT_timeLimitOneBranch;
    double EXACT_timeLimitZeroBranch;
    int RANDOM_numAttempts;
};

// TODO: Figure out a better way to handle dates
struct yymmdd {
    int yy;
    int mm;
    int dd;
};

struct Assignment {
    Assignment(int pIID, int rIID) : panelIID(pIID), roomIID(rIID) {}
    int panelIID;
    int roomIID;
};

struct ShiftedAssignment {
    ShiftedAssignment(int pIID, int rIID, int s) :
        panelIID(pIID), roomIID(rIID), shift(s) {}
    int panelIID;
    int roomIID;
    int shift;
};

struct Room {
    int IID;    // Internal ID used by program, equals index in vector
    int EID;    // External ID used in spreadsheets - can be any integer
    string name;

    int capacity;
    vector<bool> roomFeatures; // Binary vector of room features
    //vector<int> roomFeatures; // Binary vector of room features
};

struct Date {
    int IID;    // Internal ID used by program, equals index in vector
    string name;
    int year;
    int month;
    int day;
    string fName; // Formatted name / label for the date as MM/DD/YY

    // These are static after initialization (changes should be made locally)
    vector<int> availableRooms;
    vector<int> panelRequests;
};

struct Organizer {
    int IID;
    string name;
    vector<string> organizerKeys;
};

struct Panel {
    int IID;    // Internal ID used by program, equals index in vector
    int EID;    // External ID used in spreadsheets - can be any integer
    string name;
    string directorate;
    string organizer;   // Optional, not used by program
    string organizerKey;
    int numberOfDays;
    int startingDateIID;
    int size;
    vector<bool> roomRequirements; // Binary vector of room requirements

    // Auxiliary data that is not part of the input file
    int organizerIID;
};

/*****************************************************************************/
/* Details for the unordered map.                                            */
/*****************************************************************************/
namespace std
{
	template<> struct hash<vector<int>>
	{
		size_t operator() (vector<int> const &remRooms) const
		{
			size_t hash_val = 0;
			for (int i = 0; i < remRooms.size(); ++i) {
				hash_val += 71*remRooms[i];
			}
			return hash_val;
		}
	};
}

class RemRoomEquals
{
public:
    bool operator() (const vector<int> &rr1, const vector<int> &rr2) const
    {
        if (rr1.size() != rr2.size()) { return false; }
        for (int i = 0; i < rr1.size(); ++i) {
            int rr1ID = rr1[i];
            bool found = false;
            for (int j = 0; j < rr2.size(); ++j) {
                if (rr1[i] == rr2[j]) { found = true; break; }
            }
            if (!found) { return false; }
        }
        return true;
    }
};

typedef unordered_map<vector<int>, vector<Assignment>, std::hash<vector<int>>,
                      RemRoomEquals> rem_rooms_map;

/*****************************************************************************/
/* Details for sorting rooms for the exact solver                            */
/*****************************************************************************/
extern vector<Room> rooms; // Defined in scheduler.cpp
struct roomSort {
    bool operator() (int i, int j) {
        return (rooms[i].capacity < rooms[j].capacity);
    }
} RoomSort;

/*****************************************************************************/
/* Structs for tracking solver info and branching decisions.                 */
/*****************************************************************************/
struct SolutionInfo {
    SolutionInfo() {
        numPanels = 0;
        numPanelReqs = 0;

        scheduledPanelsGlobalUB = 0;
        satisfiedPanelReqsGlobalUB = 0;

        numScheduledPanels = 0;
        numSatisfiedPanelReqs = 0;

        numOptIntervals = 0;
        numCompleteIntervals = 0;

        numShiftScheduledPanels = 0;
        numShiftSatisfiedPanelReqs = 0;
    }
    // Technically we don't need these here, but it makes output a bit easier
    int numPanels;
    int numPanelReqs;

    int scheduledPanelsGlobalUB;
    int satisfiedPanelReqsGlobalUB;

    int numScheduledPanels;
    int numSatisfiedPanelReqs;

    int numOptIntervals;
    int numCompleteIntervals;

    int numShiftScheduledPanels;
    int numShiftSatisfiedPanelReqs;
    vector<std::pair<int, int> > numScheduledByShift;
};

struct BranchInfo {
    BranchInfo(int branch_iIID) {
        iIID = branch_iIID;
        numPanelsInInterval = 0;
        numPanelReqsInInterval = 0;

        numScheduledPanels = 0;
        numSatisfiedPanelReqs = 0;
        bestSolNumScheduledPanels = 0;
        bestSolNumSatisfiedPanelReqs = 0;

        exceededTimeLimit = false;
    }

    int iIID;
    int numPanelsInInterval;
    int numPanelReqsInInterval;

    int numRemPanels;
    int numRemPanelReqs;
    int numScheduledPanels;
    int numSatisfiedPanelReqs;

    int bestSolNumScheduledPanels;
    int bestSolNumSatisfiedPanelReqs;

    clock_t intervalStartTime;
    bool exceededTimeLimit;

    vector<vector<int> > remSDPanels;
    vector<int> remMDPanels;
    vector<vector<int> > remRooms;
    vector<rem_rooms_map> lookupTables;
};

/*****************************************************************************/
/* Functions for scheduler.cpp                                               */
/*****************************************************************************/
void parseArgs(int argc, const char **argv);
void printUsageAndExit(const char *argv0);

// Initialization functions
void initializeRooms();
void initializeDates();
void initializeOrganizers();
void initializePanels();
void initializeDateIntervals();
void initializeDateIntervalsHelper(int sdIID, int cdIID, int edIID);

// Output Stream Functions
void initializeOutputStream();
void outputProblemInfo();
void outputSolutionInfo(SolutionInfo &si);
void finalizeOutputStream();

// Feasibility-checking functions
void feasibilityPreTest();
void checkNumberOfPanelsPerDate();
void checkPanelRequirements();
void checkBoundsPerDate();
void checkOrganizerConflicts();

bool roomSatisfiesPanelRequirements(int rIID, int pIID);
bool roomAvailableForRequiredDays(int rIID, int pIID);

void feasibilityPostTest(vector<Assignment> &assignments);
bool checkRoomsForPanels(vector<Assignment> &assignments);
void checkForDuplicatePanels(vector<Assignment> &assignments);
void checkRoomUsage(vector<vector<Assignment> > &schedule);
void checkOrganizerDoubleBooking(vector<vector<Assignment> > &schedule);

// Solver functions
vector<Assignment> solve(SolutionInfo &si);

vector<Assignment> exactlyScheduleInterval(int iIID, SolutionInfo &si);
vector<Assignment> ext_scheduleForInterval(BranchInfo &bi);
vector<Assignment> ext_scheduleSDPanels(int dIID, vector<int> &remSDPanels,
                                        vector<int> &remRooms,
                                        rem_rooms_map &lookupTable);

int computeSatisfiedPanelReqsBound(BranchInfo &bi);
int computeScheduledPanelsBound(BranchInfo &bi, int satisfiedPanelReqsBound);

vector<Assignment> randomlyScheduleInterval(int iIID, SolutionInfo &si);
vector<Assignment> rnd_scheduleForInterval(int iIID);
vector<Assignment> rnd_scheduleForDate(int dIID, int iIID,
                                       vector<vector<int> > &intRemPanels,
                                       vector<vector<int> > &intRemRooms);
Assignment rnd_schedulePanel(int pIID, int iIID,
                             vector<vector<int> > &intRemPanels,
                             vector<vector<int> > &intRemRooms);

vector<Assignment> grd_scheduleForDate(int dIID, int iIID,
                                       vector<vector<int> > &intRemPanels,
                                       vector<vector<int> > &intRemRooms);

// Helper functions
vector<Assignment> scheduleSingleDayRequests(const vector<int>& lPanels,
                                             const vector<int>& lRooms);
vector<int> computePotentialRoomsForPanel(int pIID, int iIID,
                                          vector<vector<int> > &intRemRooms);
vector<int> computePotentialRoomsForPanel(int pIID, int shift,
                                    vector<vector<bool> > &roomAvailabilities);

// Additional post-processing functions
vector<ShiftedAssignment> shiftScheduleRemainder(SolutionInfo &si,
                                            vector<Assignment> &assignments);
ShiftedAssignment shiftSchedulePanel(int pIID, int shift,
                                vector<vector<bool> > &roomAvailabilities,
                                vector<vector<bool> > &organizerAvailabilities);


// Print / output functions
void printRooms();
void printDates();
void printPanels();

void printAssignments(vector<Assignment> &assignments);
void printSchedule(vector<Assignment> &assignments);
void printUnscheduled(vector<Assignment> &assignments);
vector<vector<Assignment> > makeSchedule(vector<Assignment> &assignments);

void printIntervalRemPanels(int iIID, vector<vector<int> > &intRemPanels);
void printIntervalRemRooms(int iIID, vector<vector<int> > &intRemRooms);

// Miscellaneous helper functions
int randomInt(int a, int b);

// I/O helper functions
void readCSV(const char *inFile, vector<vector<string> > &records);
yymmdd parseDate(string dateString);
int monthNameToNumber(string monthName);
int lookupDateIID(string dateString);
int lookupRoomIID(int rEID);

#endif // SCHEDULER_H

