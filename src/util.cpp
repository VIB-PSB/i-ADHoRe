/***************************************************************************
 *   Copyright (C) 2006, 2007, 2008 by Jan Fostier                         *
 *   jan.fostier@intec.ugent.be                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "util.h"
#include <sys/time.h>
#include <cassert>

double Util::startTime[MAX_TIMERS];
int Util::currentTimer = 0;

void Util::startChrono()
{
    // make sure we don't use too many timers
    assert(currentTimer < MAX_TIMERS);

    struct timeval tp;
    gettimeofday(&tp, NULL);
    startTime[currentTimer] = (double)tp.tv_sec+(1.e-6)*tp.tv_usec;

    currentTimer++;
}

double Util::stopChrono()
{
    // make sure stopChrono isn't called too often
    assert(--currentTimer >= 0);

    struct timeval tp;
    gettimeofday(&tp, NULL);
    double stopTime = (double)tp.tv_sec+(1.e-6)*tp.tv_usec;

    return (stopTime - startTime[currentTimer]);
}

double Util::getTime()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return (double)tp.tv_sec+(1.e-6)*tp.tv_usec;
}
