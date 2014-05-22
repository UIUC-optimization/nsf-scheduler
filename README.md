nsf-scheduler
===============================================================================

Panel Scheduler for NSF.


Contents
------------------------------------------------------------------------------
  1. General Information
  2. Code Overview
  3. Additional Constraints
  4. Installation and Usage
  5. Licensing Information
  6. Acknowledgements


Code Overview
===============================================================================

Inputs
-------------------------------------------------------------------------------
The NSF Scheduler is designed to schedule a set of panels to available rooms
across a specified time period. The inputs to the scheduling code are:

  * Rooms: A list of rooms that are available for scheduling purposes. Each
    room has an integer ID, a name, a capacity, and then a list of boolean room
    features (e.g., whether or not the room has whiteboards) that can be used
    to differentiate between rooms.

  * Dates: A list of the dates in the scheduling period, along with a listing
    of when rooms are available in that time period. Default availabilities can
    be used to specify when all rooms are available (e.g., weekdays) or
    unavailable (e.g., weekends, holidays), and room-specific availabilities
    can be used to override the default behavior.

  * Organizers: A list of panel organizers. Each organizer has a name (not
    necessarily unique) and a set of one or more unique organizer keys that
    are used to associate panels with specific organizers.

  * Panels: A list of the panels in the scheduling period. Each panel has an
    integer ID, a name, an associated program directorate, an organizer key
    (which must either belong to one of the organizers or be empty if the
    panel has no specific organizer), a number of days for the panel (1, 2, or
    3), a starting date, a size (number of panel members), and a list of
    boolean panel requirements that correspond to the room features (e.g.,
    whether or not the panel needs whiteboards).

Example inputs are located at `data/small/`. Note that all input files must be
specified in .csv format.


Constraints and Objective
-------------------------------------------------------------------------------
The NSF scheduler attempts to assign panels to rooms subject to the following
constraints:

  * Panels must start on the specified start date; they cannot be moved to
    another date (during the first phase of scheduling).

  * Multi-day panels must be assigned to the same room for each subsequent day.

  * No two panels can be assigned to the same room on the same date.

  * Panels can only be assigned to rooms that:
    * Are available for the duration of the panel;
    * Have sufficient capacity to hold all panel members;
    * Possess the features that the panel requires (e.g., if a panel requires
      whiteboards, it can only be assigned to a room that has whiteboards).

The scheduler may leave a panel unscheduled if it cannot find an available room
for that panel. The objective is to maximize the number of scheduled panels
(note that this objective treats all panels equally, regardless of how many
days the panels require).


Algorithm
-------------------------------------------------------------------------------

The algorithm does the following:

  1. Process input data and check for inconsistencies (e.g., a panel starting
    on a date not given in the input data, incomplete specifications of room
    features or panel requirements). The program will exit with an error
    message if inconsistencies are identified.

  2. Associate with each panel a set of panel requests, with each request
    corresponding to one of the dates in which the panel needs a room (e.g., a
    single-day panel has one associated panel request for a room on the
    starting date of the panel; a three-day panel has three associated panel
    requests, one on the starting date, one on the subsequent date, and one on
    the date two days after the starting date). Each date then has an
    associated set of panel requests.

  3. Separate the range of dates into a set of intervals such that the requests
    from each panel are never split into two separate intervals. (As a simple
    example, week-long intervals typically satisfy this requirement because
    panels are generally not scheduled over a weekend.) Separating the dates in
    this way ensures that the panels in one interval can be scheduled
    independently of the panels in the other intervals.

  4. Run a feasibility pre-test on the input data, identifying panel requests
    that cannot be satisfied (e.g., no room is available on the required dates
    that meets the size and feature requirements of the panel). Infeasibilities
    are noted but do not prevent the algorithm from constructing a schedule.

  5. Run one of the solvers (listed in the next subsection) to schedule the
    panels in each interval. Combine the schedules for each interval into a
    master schedule.

  6. Run a feasibility post-test on the resulting schedule to ensure that no
    constraints were violated (this is mainly used for debugging purposes).

  7. Optional: Attempt to schedule any unscheduled panels on alternate dates
    within the same week, while ensuring that rescheduling a panel does not
    cause a panel organizer to be double-booked on any date.

  8. Output the schedule and solution information to the appropriate files.


Solvers
-------------------------------------------------------------------------------

There are currently three different solvers built into the scheduling code.

  * Exact: This solver uses a branch-and-bound approach to enumerate the
    possible schedules for an interval in order to find an optimal schedule. It
    operates by branching on all possible assignments of multi-day panels to
    rooms until only single-day panels are left, and then it uses a
    push-relabel maximum flow algorithm to optimally schedule the single-day
    panels subject to the assignments for the multi-day panels. Bounding is
    accomplished by relaxing the requirement that a multi-day panel be
    scheduled in the same room on every single date.

    For performance reasons, a time limit is imposed when solving each
    interval, with the best identified solution being returned once the time
    limit is reached. Because of this, the solver may not be able to verify
    optimality of the schedule in some cases.

  * Random: This solver selects a date at random within the date range and then
    attempts to schedule all panels that have a panel request for that date.
    This is done by selecting panel requests for that date at random and
    assigning the corresponding panel to one of the remaining available rooms
    chosen at random (the panel is left unscheduled if no available rooms
    remain).

  * Greedy: This solver processes each date interval by randomly selecting a
    date from within the interval and then using the push-relabel algorithm to
    schedule as many of the panel requests from that date as possible. This
    process is repeated until all dates within an interval are examined.


Outputs
-------------------------------------------------------------------------------

There are three primary outputs from the scheduling algorithm:

  * Solution information in JSON format, including the number of panels
    scheduled and information about potential infeasibilities in the input
    data. Also contains information about how many panels were successfully
    rescheduled on a different starting date during the week.

  * A .csv file displaying the assignments of panels to rooms in a matrix
    format with columns indexed by date, rows indexed by room, and entries
    indicating the panel that is assigned to the room on the given date (an
    empty cell indicates that the room is not used on that date).

  * A .csv file listing the panels from the input data that were unable to be
    scheduled.


Additional Constraints
===============================================================================

This section describes how to include additional constraints in the scheduling
process without changing the code. This is done by making appropriate
modifications to the inputs or by post-processing the output in some way.


Reserving Rooms for Future Use
-------------------------------------------------------------------------------

To ensure that a specific room remains available for potential future use, it
is necessary to identify the dates when that room should be free, and then mark
the room as unavailable for those dates in the input file. The constructed
schedule is then guaranteed to leave that room open on the specified dates.


Restricting Access to Rooms
-------------------------------------------------------------------------------

The implementation of room features and panel requirements in the scheduling
algorithm allows for panels to be scheduled in any room as long as that room
has the features that the panel requires. However, it may be desirable to
restrict access to rooms to only certain panels (e.g., only panels from a
specific program directorate can be scheduled in a room). One way to accomplish
this is to create directorate-specific room features. For example, the feature
"G:PD1" can be used to associate rooms with Program Directorate 1. Then a panel
that needs to be scheduled in a room associated with Directorate 1 will have
the requirement "G:PD1".

However, this set of features alone would still allow other panels to be
scheduled in directorate-specific rooms. To remedy this, an additional feature
"G:GEN" can be used to indicate whether a room is available for general use.
Panels that do not need to be in a directorate-specific room would then have a
requirement that they be scheduled in a "G:GEN" room. For this to work, every
panel must state that it either needs to be scheduled in a directorate-specific
room (e.g., "G:PD1" = Y, "G:GEN" = N) or it needs to be scheduled in a general
room ("G:PD1" = N, "G:GEN" = Y). If a panel does not belong to one of these
categories, then it can still be scheduled in any room. This would make it
necessary to decide up front whether a panel should be scheduled in a
directorate-specific room or in a general room.

In many cases it may be preferable to allow the scheduler to decide whether a
panel should be assigned to a directorate-specific room or a general room. This
can be accomplished by introducing one additional room feature, "G:GEN+PD1",
which corresponds to rooms that are available for both general use and for
directorate-specific use. From the panel requirements perspective, a panel with
the "G:GEN+PD1" requirement means that the panel should be scheduled in either
a general room or a room associated with Directorate 1.

This can be summarized as follows:

  * Introduce three new room features: "G:GEN", "G:PD1", and "G:GEN+PD1"
  * Rooms that are available for general use possess the features "G:GEN" and
    "G:GEN+PD1"
  * Rooms that are restricted to Directorate 1 possess the features "G:PD1" and
    "G:GEN+PD1"
  * Panels that must be assigned to a general room have the requirement "G:GEN"
  * Panels that must be assigned to a room restricted to Directorate 1 have the
    requirement "G:PD1"
  * Panels that can be assigned to either a general room or a room restricted
    to Directorate 1 have the requirement "G:GEN+PD1"

The above process can be extended to incorporate directorate-specific rooms for
any number of directorates.


Installation and Usage
===============================================================================

The latest version of the code can be found in the
[nsf-scheduler](https://github.com/UIUC-optimization/nsf-scheduler) repository
on Github.


Requirements
-------------------------------------------------------------------------------

The project requires a C++ compiler and the STL. The program is executed from
the command-line.


Installation
-------------------------------------------------------------------------------

The Makefile in the `code` directory can be used to compile the code on most
Linux systems. For Windows, the Visual Studio 2012 project files can be used to
build the code (this may require some modifications).


Running the Code
-------------------------------------------------------------------------------

To run the code, simply execute the following in the `code/` directory:

    ./scheduler <roomsFile> <datesFile> <organizersFile> <panelsFile>

With the example data, this would be:

    ./scheduler ../data/small/rooms.csv ../data/small/dates.csv \
                ../data/small/organizers.csv ../data/small/panels.csv

The following additional arguments are optional:

  * Specify the file for printing solution information (JSON format) (Defaults
    to stdout):

        -o <outputFile>

  * Specify the schedule file for printing the matrix of panel assignments by
    date and room (CSV format) (Will not print if no file is specified):

        -s <scheduleFile>

  * Specify the file to printing the list of unscheduled panels (CSV format)
    (Will not print if no file is specified):

        -u <unscheduledFile>

  * Use the exact solver (selected by default):

        -e

  * Use the greedy solver:

        -g

  * Use the random solver:

        -r

  * Specify the number of random restarts to use with the random and greedy
    solvers:

        -n <numRandomRestarts>

  * Print additional information during the scheduling process:

        -p

  * Allow the scheduler to shift any unscheduled panels to alternate dates
    within the same week in order to reschedule them:

        -ss


Licensing Information
===============================================================================

This project is released under the University of Illinois/NCSA Open Source
License. See the LICENSE file for more details.


Acknowledgements
===============================================================================

This code was written by David R. Morrison and Jason J. Sauppe while studying
under Dr. Sheldon H. Jacobson in the Department of Computer Science at the
University of Illinois at Urbana-Champaign.

This work has been supported in part by the National Science Foundation
Graduate Research Fellowship Program (DGE-1144245). Any opinion, findings, and
conclusions or recommendations expressed in this material are those of the
author(s) and do not necessarily reflect the views of the National Science
Foundation.

