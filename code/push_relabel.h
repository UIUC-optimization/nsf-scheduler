/*****************************************************************************/
/* Author: David R. Morrison                                                 */
/* Date: Oct. 2013                                                           */
/* File: push_relabel.h                                                      */
/* Description: An implementation of the push-relabel algorithm for solving  */
/*  the network flow problem.  This algorithm is adapted from the sample     */
/*  code posted on Wikipedia.                                                */
/*                                                                           */
/* Copyright (c) 2013, 2014 University of Illinois Board of Trustees.        */
/* Published under the University of Illinois/NCSA Open Source License.      */
/* See LICENSE file for more details.                                        */
/*****************************************************************************/
#ifndef PUSH_RELABEL_H
#define PUSH_RELABEL_H

#include <vector>
#include <utility>

std::vector<std::pair<int, int>> pushRelabel(int numRequests, int numRooms,
	const std::vector<std::vector<int>>& validRooms);

#endif // PUSH_RELABEL_H
