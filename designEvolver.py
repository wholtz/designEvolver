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

import tempfile
import shutil
import re
import sys
import j5lib
import tweaks
import argparse
import copy
import random
import xmlrpc.client as xr
from deap import algorithms
from deap import base
from deap import creator
from deap import tools

random.seed(123456789)  # provide a seed to keep runs repeatable

parser = argparse.ArgumentParser(description='submit design to j5')
parser.add_argument('--j5URL', type=str,
                    default= 'http://j5.jbei.org/bin/j5_xml_rpc.pl',
                    help='Location of j5 XMLRPC server interface')
parser.add_argument('--username', type=str, required=True,
                    help='Username for logging into the j5 server')
parser.add_argument('--password', type=str, required=True,
                    help='password for logging into the j5 server')
parser.add_argument('--param', type=str, default='j5parameters.csv',
                    help='filename for j5 parameters')
parser.add_argument('--seq', type=str, default='mastersequences.csv',
                    help='filename for sequence name mapping')
parser.add_argument('--zipseq', type=str, default='masterzippedsequences.zip',
                    help='filename for zipped sequence files')
parser.add_argument('--parts', type=str,
                    default='partslist.csv', help='filename for parts list file')
parser.add_argument('--order', type=str, default='targetpartsorder.csv',
                    help='filename for target part order file')
parser.add_argument('--eugene', type=str,
                    default='eugeneruleslist.eug', help='filename for eugene rules')
parser.add_argument('--plasmid', type=str,
                    default='j5_plasmids.csv', help='filename for plasmid list')
parser.add_argument('--oligo', type=str,
                    default='j5_oligos.csv', help='filename for oligo list')
parser.add_argument('--direct', type=str, default='j5_directsyntheses.csv',
                    help='filename for direct syntheses list')
parser.add_argument('--diff', type=str, default='',
                    help='genbank filename to compare final sequence against')
parser.add_argument('--outdir', type=str, default='./DE/',
                    help='directory to write output and temporary files to')
parser.add_argument('--numgen', type=int, default=10,
                    help='number of generations')
parser.add_argument('--maxpop', type=int, default=50,
                    help='maximum population size')
parser.add_argument('--tweaks', type=int, default=30,
                    help='number of tweaks per design')
parser.add_argument('--multiplier', type=float, default=4,
                    help='childern per parent')
parser.add_argument('--serial', dest='serial', action='store_true',
                    help='do not evalute designs in parallel')
parser.set_defaults(serial=False)
parser.add_argument('--nosubmit', dest='nosubmit', action='store_true',
                    help='checks that design can be loaded, but not sent to j5')
parser.set_defaults(nosubmit=False)
parser.add_argument('--verbose', dest='verbose', action='store_true',
                    help='output more information')
parser.set_defaults(verbose=False)
args = parser.parse_args()

maxPopulation = args.maxpop
mu = args.maxpop
lambda_ = int(args.maxpop*args.multiplier)
maxProcesses = min(40, mu+lambda_)
initialTweaksPerList = args.tweaks
generations = args.numgen
mutationRate = 0.5
mutationDegree = 0.3
crossRate = 0.5
hallOfFameLength = 1

inputFiles = { 'param': args.param, 'seq': args.seq, 'zipseq': args.zipseq, 'parts': args.parts, 'order': args.order, 'eugene': args.eugene}

if not args.serial:
    import multiprocessing
    import logging
    if args.verbose:
        logger = multiprocessing.log_to_stderr()
        logger.setLevel(multiprocessing.SUBDEBUG)

def removeDirectories(removeList, exceptionList = []):
    for d in set(removeList) - set(exceptionList):
        shutil.rmtree(d, ignore_errors=True)

def applyTweaks(server, design, inputFiles, individual, retry=0):
    results = None
    initialLengthPlasmid = len(design.getFinalSequence())
    tweakDesign = copy.deepcopy(design)
    individual.sortBySilence() # move non-silent tweaks to start of list
    tweakDesign.breakAllBigParts() # parts could change size due to tweaks, might create violations of max part size
    for t in individual.tweaks:
        t.execute(tweakDesign)
    tweakDesign.breakAllBigParts() # this is a poor way to ensure no big parts were generated...
    positions = tweakDesign.getMaxPosition()
    lengthPlasmid = len(tweakDesign.getFinalSequence())
    if not args.nosubmit:
        try:
            inputDir = tempfile.mkdtemp(dir=args.outdir)
            returnFromServer = server.submitToj5(tweakDesign, inputFiles,
                    inputDir=inputDir, delInputDir=True)
        except xr.Fault as err:
            if args.verbose:
                print('caught xmlrpc.client.Fault on job submission')
                print('submit dir: '+inputDir)
                print(individual)
                print(tweakDesign.parts)
                sys.stdout.flush()
            retryPattern = re.compile("Can't extract zipped sequence file")
            if retryPattern.match(err.faultString) and retry<3:
                if args.verbose:
                    print('caught ioclt fault on job submission, retrying...')
                    sys.stdout.flush()
                return applyTweaks(server, design, inputFiles, individual, retry=retry+1)
            warnings = 8888888
        else:
            results = j5lib.J5Results(returnFromServer, individual.outputDir)
            warnings = results.warnings[0]
    else:
        warnings = 9999999
    return (positions, lengthPlasmid, warnings, len(individual.tweaks))

design = j5lib.J5Design()
design.loadFromFiles(
    args.param, args.seq, args.zipseq, args.parts, args.order, args.eugene,
    args.plasmid, args.oligo, args.direct)
server = j5lib.J5Connection(args.username, args.password, args.j5URL)

if args.diff != '':
    try:
        s = j5lib.readGenbankSequence(args.diff)
    except TypeError:
        print("Error loading genbank file for diff")
        raise
    finalSeq = list(s.values())[0]
else:
    finalSeq = design.getFinalSequence()
startPlasmidCompiled = re.compile(finalSeq.lower())
startPlasmidLength = len(finalSeq)

def isSamePlasmid(p):
    # use global for compiled re, so that it is only calculated once
    p=p.lower()
    return (re.search(startPlasmidCompiled, p+p) != None) and (len(p) == startPlasmidLength)


# three weights are
# (maxPositionNumber, lengthPlasmid, number of warnings, number of tweaks)
creator.create("FitnessMulti", base.Fitness, weights=(-10000, -1, -1000000, -10))
creator.create("Individual", j5evo.J5TweakSet, fitness=creator.FitnessMulti)
toolbox = base.Toolbox()
toolbox.register("individual", creator.Individual, args.tweaks, args.outdir)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("evaluate", applyTweaks, server, design, inputFiles)
toolbox.register("mate", j5evo.mateJ5TweakSets)
toolbox.register("mutate", j5evo.mutateTweakSet, mutationDegree)
toolbox.register("select", tools.selTournament, tournsize=maxPopulation)

if __name__ == '__main__':
    if (args.verbose):
        print(args)
    if not args.serial:
        pool = multiprocessing.Pool(processes=maxProcesses)
        toolbox.register("map", pool.map)
    tempDirectories = []
    pop = toolbox.population(n=lambda_+mu)
    best = tools.HallOfFame(hallOfFameLength)
    history = tools.History()
    toolbox.decorate("mate", history.decorator)
    toolbox.decorate("mutate", history.decorator)
    history.update(pop)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", tools.mean)
    stats.register("std", tools.std)
    stats.register("min", min)
    stats.register("max", max)

    for g in range(0,generations):
        pop = j5evo.uniqueRepr(pop) # remove duplicates
        while len(pop) < 2:
            pop.append(toolbox.individual())
        algorithms.eaMuPlusLambda(pop, toolbox, mu=mu, lambda_=lambda_,
                        cxpb=crossRate, mutpb=mutationRate,
                        ngen=1, stats=stats, halloffame=best,
                        verbose=args.verbose)
        tempDirectories.extend([i.outputDir for i in history.genealogy_history.values()])
        removeDirectories(tempDirectories, [i.outputDir for i in best])
    print('output from best design: '+str(best[0].outputDir))

