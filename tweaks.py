#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright ©2013. The Regents of the University of California (Regents). All
# Rights Reserved. Permission to use, copy, modify, and distribute this
# software and its documentation for educational, research, and
# not-for-profit purposes, without fee and without a signed licensing
# agreement, is hereby granted, provided that the above copyright notice,
# this paragraph and the following two paragraphs appear in all copies,
# modifications, and distributions. Contact The Office of Technology
# Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA
# 94720-1620, (510) 643-7201, for commercial licensing opportunities.
# 
# Created by William J. Holtz, The California Institute for Quantitative
# Biosciences, University of California, Berkeley.
# 
# IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
# SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
# ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
# REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
# HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE
# MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

import j5lib
import numpy
import tempfile
import random
import math
import copy

random.seed(123456789)  # provide a seed to keep runs repeatable


class J5Tweak(object):

    def __init__(self, parameters = None):
        if parameters is None:
            self.parameters = {}
            self.mutate(1)  # initialize parameters to random values
        else:
            self.parameters = parameters

    def __repr__(self):
        d = dictToOrderedString(self.parameters)
        return "j5evo."+self.__class__.__name__ + "(" + d +")"

    def __cmp__(self, other):
        if self.__repr__ > other.__repr__:
            return 1
        if self.__repr__ < other.__repr__:
            return -1
        return 0

    def execute(self, design):
        raise j5Error('J5Tweak.execute()  was not overriden by child class:'+self.__class__.__name__)

    def mutate(self, degree=1):
        # degree should be a float in range [0,1] where 0 is no mutation and
        # 1 is maximal mutation(aka random parameters)
        raise j5Error('J5Tweak.mutate() was not overriden by child class:'+self.__class__.__name__)

class J5TweakSet(object):

    def __init__(self, numTweaks, baseOutputDir):
        self.baseOutputDir = baseOutputDir
        self.updateTweaks([generateRandomTweak() for i in range(0, numTweaks)])

    def __repr__(self):
        # return self.outputDir +' '+str(self.tweaks)
        return str(self.tweaks)

    def __cmp__(self, other):
        if self.tweaks > other.tweaks:
            return 1
        if  self.tweaks < other.tweaks:
            return -1
        return 0


    def setNewOutputDir(self):
        self.outputDir = tempfile.mkdtemp(dir=self.baseOutputDir)

    def updateTweaks(self, tweakList):
        self.tweaks = tweakList
        self.setNewOutputDir()

    def getRandomHeadOfTweakList(self):
        # return between 1 and all head members of list
        return [self.tweaks[i] for i in range(0, random.choice(range(0, len(self.tweaks))))]

    def getRandomTailOfTweakList(self):
        # return between 1 and all tail members of list
        return [self.tweaks[i] for i in range(random.choice(range(0, len(self.tweaks))), len(self.tweaks))]

    def mutate(self, mutationDegree):
        if random.random()< mutationDegree:
           del self.tweaks[random.choice(range(0,len(self.tweaks)))] 
        for t in self.tweaks:
            t.mutate(mutationDegree)
        if random.random()< mutationDegree:
            pos = random.choice(range(0, len(self.tweaks)+1))
            self.tweaks.insert(pos, generateRandomTweak())
        if random.random()< mutationDegree:
            # swap positions of two tweaks
            r1 = random.choice(range(0,len(self.tweaks))) # random index from self.tweaks
            r2 = random.choice(range(0,len(self.tweaks))) # random index from self.tweaks
            self.tweaks[r2], self.tweaks[r1] = self.tweaks[r1], self.tweaks[r2]
        if len(self.tweaks) == 0: # don't let there be empty tweak lists
            self.tweaks.append(generateRandomTweak())
        self.updateTweaks(self.tweaks)


    def sortBySilence(self):
        # move all non-silent tweaks to the start of the list
        # do not change the relative order within the non-silent 
        # tweaks or within the silent tweaks
        nonSilent = []
        silent = []
        nonSilentTweaks = ['ModifyPartEnd']
        for t in self.tweaks:
            if t.__class__.__name__ in nonSilentTweaks:
                nonSilent.append(t)
            else:
                silent.append(t)
        self.tweaks = nonSilent + silent


def sortedOnRepr(l):
    return sorted(l, key=repr)


def uniqueRepr(l):
    if len(l) <= 1:
        return l
    s = sortedOnRepr(l)
    same = [a==b for (a, b) in zip(s, s[1:]+[s[0]])]
    if sum(same) == len(same): # if all equal
        return [s[0]]
    for i in reversed(range(0,len(s))):
        if same[i]:
            del s[i]
    return s


def dictToOrderedString(d):
    out = '{'
    firstLoop = True
    for k in sorted(d.keys()):
        if not firstLoop:
            out = out + ', '
        else:
            firstLoop = False
        out = out + repr(k) + ': '+repr(d[k])
    return out + '}'


def generateRandomTweak():
    possibleTweaks = [SilentJunctionShift, ReducePartsSilently, SplitPart, ForcedRelativeOverhang, Extra5pCPECOverlap, Extra3pCPECOverlap, ModifyPartEnd]
    return random.choice(possibleTweaks)() # return an instance


def mutateTweakSet(mutationDegree, tweakSet):
    tweakSetCopy = copy.deepcopy(tweakSet)
    tweakSetCopy.mutate(mutationDegree)
    return tweakSetCopy,

def areEqualDicts(a,b):
    if set(a.keys()) != set(b.keys()):
        return False
    for k in a.keys():
        if a[k] != b[k]:
            return False
    return True

def areLessThanDicts(a,b):
    # fairly arbitrary way of comparing dictionaries,
    # but more robust than saying equal key set and equal value sets
    # implies equal dictionaries (wrong, because values could be swapped)
    ka = sorted(set(a.keys))
    kb = sorted(set(b.keys))
    if  ka < kb:
        return True
    if ka > kb:
        return False
    for k in ka:
        if a[k] < b[k]:
            return True
    return False

def areEqualTweaks(a,b):
    if a.__class__.__name__ != b.__class__.__name__:
        return False
    return areEqualDicts(a.parameters, b.parameters)

def areLessThanTweaks(a,b):
    if a.__class__.__name__ < b.__class__.__name__:
        return True
    if a.__class__.__name__ > b.__class__.__name__:
        return False
    return areLessThanDicts(a.parameters, b.parameters)

def areEqualTweakSets(a,b):
    if len(a.tweaks) != len(b.tweaks):
        return False
    for (ta, tb) in  zip(A.tweaks, B.tweaks):
        if not areEqualTweaks(ta,tb):
            return False
    return True

def areLessThanTweakSets(a,b):
    if len(a.tweaks) < len(b.tweaks):
        return True
    if len(a.tweaks) > len(b.tweaks):
        return False
    for (ta, tb) in zip(a.tweaks, b.tweaks):
        if areLessThanTweaks(ta,tb):
            return True
    return False

def mateJ5TweakSets(ts1, ts2):
    #copy J5TweakSets before calling this function. They are modified in place
    ts1.updateTweaks(ts1.getRandomHeadOfTweakList()+ts2.getRandomTailOfTweakList())
    ts2.updateTweaks(ts2.getRandomHeadOfTweakList()+ts1.getRandomTailOfTweakList())
    return ts1, ts2


class SinglePartSingleParameter(J5Tweak):

    def mutate(self, degree=1, maxVariance=6, canBeNegative = False):
        # variance of geometic distribution = (1-p)/p^2
        v = float(degree)*float(maxVariance)
        p=(((4*v+1)**0.5)-1)/(2*v)
        if canBeNegative:
            posNeg = random.choice([-1, 1])
        else:
            posNeg = 1
        randomAmount = numpy.random.geometric(p) * posNeg
        if degree == 1:
            self.parameters['position'] = int(
                math.floor(random.random() * 10000))  # this will get wrapped via %
            self.parameters['amount'] = randomAmount
        else:
            if random.random() < degree:
                self.parameters['position'] += random.choice([-1, 1])
            self.parameters['amount'] += randomAmount


class SinglePartVectorParameter(SinglePartSingleParameter):
    # add an end of the DNA strand to SinglePartSingleParameter
    def mutate(self, degree=1, maxVariance=6, canBeNegative=False):
        super(SinglePartVectorParameter, self).mutate(degree, maxVariance, canBeNegative)
        if random.random() < degree:
            self.parameters['end'] = random.choice(['end3p', 'end5p'])


class EndShift(J5Tweak): # for inheritance only. Not to be directly instantiated

    def mutate(self, degree=1):
        if degree == 1:
            self.parameters['fractionOfMaxShift'] = random.random()
            self.parameters['position'] = int(
                math.floor(random.random() * 10000))  # this will get wrapped via %
        if random.random() < degree:
            self.parameters['end'] = random.choice(['end3p', 'end5p'])
        if random.random() < degree:
            self.parameters['position'] += random.choice([-1, 1])
        self.parameters['fractionOfMaxShift'] = random.random() * degree + self.parameters[
            'fractionOfMaxShift'] * (1 - degree)


class ModifyPartEnd(SinglePartVectorParameter):

    def mutate(self, degree=1, maxVariance=20, canBeNegative=True):
        super(ModifyPartEnd, self).mutate(degree, maxVariance, canBeNegative)

    def execute(self, design):
        # end3p and end5p refer to final assembed orientation of part
        part = design.getPartByPosition(self.parameters['position'])
        if self.parameters['end'] != 'end3p' and self.parameters['end'] != 'end5p':
            raise j5Error( "ModifyPartEnd.parameters['end'] != 'end3p' or 'end5p'")
        modifyStart = part.useInSourceDirection() == (self.parameters['end'] == 'end5p')
        part.updateEdgeAndCheckBounds(design, self.parameters['amount'], modifyStart)


class SilentJunctionShift(EndShift):

    def execute(self, design):
        shifts = design.findSilentJunctionShifts()
        part = design.getPartByPosition(self.parameters['position'])
        if self.parameters['end'] == 'end3p':
            direction = 1
            neighborPart = design.getNextPart(part)
        elif self.parameters['end'] == 'end5p':
            direction = -1
            neighborPart = design.getPreviousPart(part)
        else:
            raise j5Error( "SilentJunctionShift.parameters['end'] != 'end3p' or 'end5p'")
        maxShiftDistance = min(shifts[part.position][1+direction],len(design.getPartSequence(neighborPart))-1)
        if self.parameters['end'] == 'end5p': # need to perform shift from part 5 prime of the junction
            part = neighborPart
        distance = direction*int(
            math.ceil(maxShiftDistance * self.parameters['fractionOfMaxShift']))
        design.shiftPartJunction(part, distance)


class ReducePartsSilently(J5Tweak):

    def execute(self, design):
        design.reducePartsSilently()

    def mutate(self, degree=1):
        pass # there are no parameters to mutate


class SplitPart(J5Tweak):

    def execute(self, design):
        position = design.resolveWrappedAssemblyPosition(self.parameters['position'])
        design.breakPart(position)

    def mutate(self, degree=1):
        if degree == 1:
            self.parameters['position'] = int(
                math.floor(random.random() * 10000))  # this will get wrapped via %
        if random.random() < degree:
            self.parameters['position'] += random.choice([-1, 1])


class ForcedRelativeOverhang(SinglePartSingleParameter):

    def execute(self, design):
        part = design.getPartByPosition(self.parameters['position'])
        part.forcedRelativeOverhang = self.parameters['amount']
    
    def mutate(self, degree=1):
        super(ForcedRelativeOverhang, self).mutate(degree, maxVariance = 8, canBeNegative = True) 

class Extra5pCPECOverlap(SinglePartSingleParameter):

    def execute(self, design):
        part = design.getPartByPosition(self.parameters['position'])
        part.extra5pCPECOverlap= self.parameters['amount']

    def mutate(self, degree=1):
        super(Extra5pCPECOverlap, self).mutate(degree, maxVariance = 5, canBeNegative = False) 


class Extra3pCPECOverlap(SinglePartSingleParameter):

    def execute(self, design):
        part = design.getPartByPosition(self.parameters['position'])
        part.extra3pCPECOverlap= self.parameters['amount']

    def mutate(self, degree=1):
        super(Extra3pCPECOverlap, self).mutate(degree, maxVariance = 5, canBeNegative = False)



if __name__ == "__main__":
    print("Python library for evolving j5 designs")

