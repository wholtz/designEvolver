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

import shutil
import random
import xmlrpc.client as xr
import zipfile
import re
import os
import math
import copy
import base64
import sys
import tempfile
import csv
import time
from Bio.Seq import Seq
from Bio import SeqIO
from xml.dom import minidom

random.seed(123456789)  # provide a seed to keep runs repeatable

def matchLength(a, b):
    # returns the length of equal elements at start of a & b
    length = 0
    for x,y in zip(a, b):
        if x == y:
            length += 1
        else:
            break
    return length


def randomDigits(n):
    return ''.join([str(random.choice(range(0,9))) for i in range(0,n)])


def readCSVFile(filename):
    data = []
    try:
        with open(filename, 'r') as f:
            for l in csv.reader(f):
                data.append(l)
    except IOError:
        sys.stderr.write(filename + ' either not found or not readable\n')
        raise
    return data


def writeCSVFile(filename, data):
    try:
        with open(filename, 'w') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(data)
    except IOError:
        sys.stderr.write('cannot open '+filename+' for writing')
        raise


def readFile(filename):
    try:
        with open(filename, 'rb') as f:
            data = f.read()
    except IOError:
        sys.stderr.write(
            filename + ' either not found or not readable. Treating it as an empty file.\n')
        data = ""
    return data


def writeFile(filename, s):
    try:
        with open(filename, 'wb') as f:
            f.write(decode(s))
    except IOError:
        sys.stderr.write('Error: Could not open ' + filename + ' for writing.')
        raise


def encode(s):
    if s == "":
        return s
    return base64.b64encode(s).decode("utf-8")


def decode(s):
    return base64.b64decode(s.encode("utf-8"))


def stringReverse(s):
    return s[::-1]


def stringToInt(s):
    if s == '':
        return ''
    return int(s)


def matchAny(s, regexList):
    patterns = []
    for r in regexList:
        patterns.append(re.compile(r))
    for p in patterns:
        if re.search(p, s) != None:
            return True
    return False


def findFilesInZip(regex, zipfilename):
    # returns the list of filenames matching regex within zipfilename
    found = []
    with zipfile.ZipFile(zipfilename, 'r') as myzip:
        for f in myzip.infolist():
            pattern = re.compile(regex)
            if pattern.search(f.filename):
                found.append(f.filename)
    return found


def mergeDictionaries(x, y):
    return dict(list(x.items()) + list(y.items()))


def readJBEISeqSequence(filename):
    try:
        jbeiSeq = minidom.parse(filename)
    except:
        raise j5Error('cannot parse JBEI-seq file: '+filename)
    try:
        names = jbeiSeq.getElementsByTagName('seq:name')
        seq = jbeiSeq.getElementsByTagName('seq:sequence')
    except:
        raise j5Error('cannot find required elements in JBEI-seq file: '+filename)
    sequences = {}
    for i in zip(names, seq):
        sequences[str(i[0].childNodes[0].nodeValue)] = str(
            i[1].childNodes[0].nodeValue)
    return sequences


def readGenbankSequence(filename):
    sequences = {}
    for record in SeqIO.parse(filename, 'genbank'):
        sequences[str(record.id)] = str(record.seq)
    return sequences


def readFASTASequence(filename):
    sequences = {}
    for record in SeqIO.parse(filename, 'fasta'):
        sequences[str(record.id)] = str(record.seq)
    return sequences


def readSBOLXMLSequence(filename):
    raise j5Error('read SBOLXML not yet implemented')
    return {}


def unzipInTemp(zipFilename, baseDir=None):
    tempDir = tempfile.mkdtemp(dir=baseDir)
    zipF = zipfile.ZipFile(zipFilename, 'r')
    zipF.extractall(tempDir)
    zipF.close()
    return tempDir


def isInSegment(start, end, query):
    # return True is query is in the closed interval
    # [start,end]. Segment always goes in the increasing direction
    # from start, even if end wraps to a value smaller than start.
    if start <= end: # no wrapping past 1+ on plasmid
        return (start <= query) and (query <= end)
    # wrapping case
    return (start <= query) or (query <= end)

def shortestDistance(A,B, plasmidSize):
    # return the shortest number of basepairs between A and B on a
    # plasmid. For adjacent basepairs, return 1.
    # Plasmid size is number of base pairs, 1 indexed.
    # Two possible distances for each pair. One that crosses +1 and one
    # that does not cross +1
    return min(abs(A-B), plasmidSize-abs(A-B))

def isFirstPositionCloserThanSecond(A, B, C, plasmidSize):
    # on a circular plasmid, is A closer to C than B is to C?
    # return True if A is closer 
    # Plasmid size is number of base pairs, 1 indexed.
    return shortestDistance(A,C, plasmidSize) < shortestDistance(B,C, plasmidSize)

def resolveWrappedSequencePosition(p, maximum):
    # p in 1<=p<=maximum return p,
    # else return corresponding value in that range
    return (p % maximum) + int(p % maximum == 0) * maximum


class j5Error(Exception):
    pass


class J5Connection(object):

    def __init__(self, userName, password, server_url):
        self.userName = userName
        self.password = password
        self.server_url = server_url
        # don't include the server in the class attributes as I think
        # this has pickling problems when used with multiprocessing
        server = xr.Server(self.server_url, allow_none=True)
        result = server.CreateNewSessionId( {'username': self.userName, 'password': self.password})
        self.sessionId = result['j5_session_id']

    def __repr__(self):
        return self.userName + ":" + str(len(self.password)) + ":" + self.server_url

    def submitToj5(self, design, inputFilenames, inputDir=tempfile.mkdtemp(), delInputDir=True, verbose=False, resolution=0):
        # setup the files for design submission
        for filename in  inputFilenames.values():
            shutil.copyfile(filename, inputDir+'/'+os.path.basename(filename))
        paramFilename = inputDir+'/'+inputFilenames['param']
        partsFilename = inputDir+'/'+inputFilenames['parts']
        orderFilename = inputDir+'/'+inputFilenames['order']
        design.writeToFiles(paramFilename, partsFilename, orderFilename)
        design.xmlrpcParameters['encoded_j5_parameters_file'] = encode( readFile(paramFilename))
        design.xmlrpcParameters['encoded_parts_list_file'] = encode( readFile(partsFilename))
        design.xmlrpcParameters['encoded_target_part_order_list_file'] = encode(readFile(orderFilename))
        if delInputDir:
            shutil.rmtree(inputDir, ignore_errors=True)
        designWithAuth = design.getXmlrpcDesign(self.sessionId)
        server = xr.Server(self.server_url, verbose=verbose, allow_none=True)
        return server.DesignAssembly(designWithAuth)


class J5Results(object):

    def __init__(self, result, outputDir):
        self.outputDir = outputDir
        zipFilename = outputDir+'/'+result['output_filename']
        writeFile(zipFilename, result['encoded_output_file'])
        self.unzipj5Results(zipFilename)
        self.plasmidNames = self.getPlasmidNames()
        self.assemblyFiles = [
            self.outputDir+ '/' + p + '.csv' for p in self.plasmidNames]
        self.warnings = self.getWarningCounts()

    def __repr__(self):
        s = 'zip filename:' + self.zipFilename + '\n'
        s += 'output directory:' + self.outputDir+ '\n'
        s += 'plasmid names:' + str(self.plasmidNames) + '\n'
        s += 'assembly files:' + str(self.assemblyFiles) + '\n'
        s += 'warning counts:' + str(self.warnings) + '\n'
        return s

    def unzipj5Results(self, filename):
        self.zipFilename = filename
        zipF = zipfile.ZipFile(filename, 'r')
        zipF.extractall(self.outputDir)
        zipF.close()
        for filename in os.listdir(self.outputDir):
            zipPattern = re.compile('.zip$')
            if re.search(zipPattern, filename) != None:
                zipF = zipfile.ZipFile(self.outputDir+ '/' + filename)
                zipF.extractall(self.outputDir)
                zipF.close()

    def getWarningCounts(self):
        # these are strings from warnings that we aren't going to count against
        # designs
        excludedWarnings = [
            'file is empty',
            'homologous sequence repeat',
            'is longer than the 36 bp Primer3 Tm calculation limit.'
        ]
        countList = []
        for f in self.assemblyFiles:
            data = readCSVFile(f)
            warnPattern = re.compile("Warning: ")
            count = 0
            for l in data:
                for i in l:
                    if re.search(warnPattern, i) != None:
                        if not matchAny(i, excludedWarnings):
                            count += 1
            countList.append(count)
        return countList

    def getPlasmidNames(self):
        genbank = findFilesInZip(".gb$", self.zipFilename)
        if genbank == []:
            raise j5Error('No Genbank file found in j5 results')
        return [re.sub('\.gb$', '', g) for g in genbank]

    def getPlasmidGenbankFiles(self):
        return [self.outputDir+'/'+p+'.gb' for p in self.getPlasmidNames()]

class Part(object):

    def __init__(self):
        self.source = ''
        self.reverseComplement = False
        self.start = 1  # lowest index is 1
        self.end = 1
        self.forwardDirection = True
        self.forcedAssemblyStrategy = ''
        self.forcedRelativeOverhang = 0
        self.directSynthesisFirewall = False
        self.extra5pCPECOverlap = 0
        self.extra3pCPECOverlap = 0
        self.position = 0  # zero index ordering of parts, from targetpartsorder.csv

    def __repr__(self):
        s  = '\n'
        s += 'pos:' + str(self.position) +  '\t'
        s += 'RC:' + str(self.reverseComplement) + ' '
        s += 's:' + str(self.start) + '\t'
        s += 'e:' + str(self.end) + '\t'
        s += 'dir:' + str(self.forwardDirection) + ' '
        s += 'FAS:' + self.forcedAssemblyStrategy + '\t'
        s += 'FRO:' + str(self.forcedRelativeOverhang) + '\t'
        s += 'DSF:' + str(self.directSynthesisFirewall) + '\t'
        s += 'e5p:' + str(self.extra5pCPECOverlap) + '\t'
        s += 'e3p:' + str(self.extra3pCPECOverlap) + '\t'
        s +=  self.source + '\n'
        return s

    def setForcedRelativeOverhang(self, amount):
        self.forcedRelativeOverhang = amount

    def setExtra5pCPECOverlap(self, amount):
        self.extra5pCPECOverlap = amount

    def setExtra3pCPECOverlap(self, amount):
        self.extra3pCPECOverlap = amount

    def useInSourceDirection(self):
        return self.reverseComplement != self.forwardDirection

    def getEdge(self, forStart):
        if forStart:
            return self.start
        else:
            return self.end

    def setEdge(self, forStart, value):
        if forStart:
            self.start = value
        else:
            self.end = value

    def getBound(self, forLower, forStart):
        # forLower is True when requesting lower bound 
        # and False when requesting upper bound
        # forStart is True for start related bounds
        # and False for end related bounds
        # if self.useInSourceDirection() == forLower:
        if forLower:
            if forStart:
                return self.maxStart
            else:
                return self.minEnd
        else:
            if forStart:
                return self.minStart
            else:
                return self.maxEnd

    def updateEdgeAndCheckBounds(self, design, evalDelta, isStart):
        # add evalDelta to self.start or end and check that the result falls
        # within self.minStart/End and self.maxStart/End. If not in this
        # segment, set start to be equal to the violated bound
        sourceLength = len(design.sequences[self.source])
        lowerBound = self.getBound(forLower=True, forStart=isStart)
        upperBound = self.getBound(forLower=False, forStart=isStart)
        notWrapped = self.getEdge(isStart)+evalDelta # current self edge plus delta
        currentPlusDelta = resolveWrappedSequencePosition(notWrapped, sourceLength)
        if isInSegment(lowerBound, upperBound, currentPlusDelta):
            self.setEdge(isStart, currentPlusDelta)
        elif isFirstPositionCloserThanSecond(lowerBound, upperBound, currentPlusDelta, sourceLength):
            self.setEdge(isStart, lowerBound)
        else:
            self.setEdge(isStart, upperBound)
        if isStart:
            final = self.start
        else:
            final = self.end
        # print("updatedEdgeAndCheckBound - isStart:"+str(isStart)+" lowerBound:"+str(lowerBound)+" upperBound:"+str(upperBound)+" currentPlusDelta:"+str(currentPlusDelta)+" final:"+str(final))


class J5Design(object):

    def __init__(self, maxPartLength=5000):
        self.parts = {}
        self.sequences = {}
        self.xmlrpcParameters = {}
        self.maxPartLength = maxPartLength
        self.sequencesTempDir = '.'


    def __repr__(self):
        return self.parts.__repr__() + self.sequences.__repr__()


    def loadFromFiles(self, paramFilename, seqFilename, zipSeqFilename, partsFilename, orderFilename, eugeneFilename, plasmidFilename, oligoFilename, directFilename, method='SLIC/Gibson/CPEC'):
        self.xmlrpcParameters['reuse_j5_parameters_file'] = False
        self.xmlrpcParameters['encoded_j5_parameters_file'] = encode(
            readFile(paramFilename))
        self.xmlrpcParameters['reuse_sequences_list_file'] = False
        self.xmlrpcParameters['encoded_sequences_list_file'] = encode(
            readFile(seqFilename))
        self.xmlrpcParameters['reuse_zipped_sequences_file'] = False
        self.xmlrpcParameters['encoded_zipped_sequences_file'] = encode(
            readFile(zipSeqFilename))
        self.xmlrpcParameters['reuse_parts_list_file'] = False
        self.xmlrpcParameters['encoded_parts_list_file'] = encode(
            readFile(partsFilename))
        self.xmlrpcParameters['reuse_target_part_order_list_file'] = False
        self.xmlrpcParameters['encoded_target_part_order_list_file'] = encode(
            readFile(orderFilename))
        self.xmlrpcParameters['reuse_eugene_rules_list_file'] = False
        self.xmlrpcParameters['encoded_eugene_rules_list_file'] = encode(
            readFile(eugeneFilename))
        self.xmlrpcParameters['reuse_master_plasmids_file'] = False
        self.xmlrpcParameters['encoded_master_plasmids_file'] = encode(
            readFile(plasmidFilename))
        self.xmlrpcParameters[
            'master_plasmids_list_filename'] = plasmidFilename
        self.xmlrpcParameters['reuse_master_oligos_file'] = False
        self.xmlrpcParameters['encoded_master_oligos_file'] = encode(
            readFile(oligoFilename))
        self.xmlrpcParameters['master_oligos_list_filename'] = oligoFilename
        self.xmlrpcParameters['reuse_master_direct_syntheses_file'] = False
        self.xmlrpcParameters['encoded_master_direct_syntheses_file'] = encode(
            readFile(directFilename))
        self.xmlrpcParameters[
            'master_direct_syntheses_list_filename'] = directFilename
        self.xmlrpcParameters['assembly_method'] = method
        self.sequencesTempDir = unzipInTemp(zipSeqFilename)
        self.getSequences(seqFilename)
        shutil.rmtree(self.sequencesTempDir, ignore_errors=True)
        self.getParts(partsFilename)
        self.getTargetOrder(orderFilename)


    def writeToFiles(self, paramFilename, partsFilename, orderFilename):
        self.writeParameters(paramFilename)
        self.writeParts(partsFilename)
        self.writeOrder(orderFilename)


    def writeParameters(self, paramFilename):
        pass # not yet needed/implemented


    def writeParts(self, partsFilename):
        data = [['Part Name','Part Source (Sequence Display ID)','Reverse Compliment?','Start (bp)','End (bp)']]
        for k, p in self.parts.items():
            data.append([k, p.source, p.reverseComplement, p.start, p.end])
        writeCSVFile(partsFilename, data)


    def writeOrder(self, orderFilename):
        data = [['(>Bin) or Part Name', 'Direction', 'Forced Assembly Strategy?', 'Forced Relative Overhang Position?', 'Direct Synthesis Firewall?', "Extra 5' CPEC overlap bps", "Extra 3' CPEC overlap bps"]]
        for k in self.getKeysByPosition():
            p = self.parts[k]
            if p.forwardDirection:
                forward = 'forward'
            else:
                forward = 'reverse'
            data.append([k, forward, p.forcedAssemblyStrategy, p.forcedRelativeOverhang, p.directSynthesisFirewall, p.extra5pCPECOverlap, p.extra3pCPECOverlap])
        writeCSVFile(orderFilename, data)


    def getXmlrpcDesign(self, sessionId):
        temp = copy.copy(self.xmlrpcParameters)
        temp['j5_session_id'] = sessionId
        return temp


    def getSequences(self, filename):
        data = readCSVFile(filename)
        for l in data[1:]:  # skip header line
            if l[1] == 'jbei-seq':
                readFun = readJBEISeqSequence
            elif l[1] == 'Genbank':
                readFun = readGenbankSequence
            elif l[1] == 'FASTA':
                readFun = readFASTASequence
            elif l[1] == 'SBOLXML':
                readFun = readSBOLXMLSequence
            else:
                raise j5Error(
                    'Unknown sequence file format, ' + l[1] + ' in ' + filename)
            self.sequences = mergeDictionaries(
                self.sequences, readFun(self.sequencesTempDir + '/' + l[0]))


    def getParts(self, partsFilename):
        data = readCSVFile(partsFilename)
        for l in data[1:]:  # skip header line
            self.parts[l[0]] = Part()
            self.parts[l[0]].keyName = l[0]
            self.parts[l[0]].source = str(l[1])
            self.parts[l[0]].reverseComplement = (l[2].upper() == 'TRUE')
            if l[3] == '':
                l[3] = 1
            self.parts[l[0]].start = int(l[3])
            self.parts[l[0]].minStart = self.parts[l[0]].start
            self.parts[l[0]].maxStart = self.parts[l[0]].start
            if l[4] == '':
                l[4] = len(l[1])
            self.parts[l[0]].end = int(l[4])
            self.parts[l[0]].minEnd = self.parts[l[0]].end
            self.parts[l[0]].maxEnd = self.parts[l[0]].end


    def getTargetOrder(self, orderFilename):
        # Currently not doing combinatorial assemblies and ignoring
        # bins. If a partID equals 'MAX_.*' then do not increment the
        # position index and do not create a new part, instead add it to 
        # the last part added as a maximum DNA range.
        data = readCSVFile(orderFilename)
        position = 0
        lastPart = None
        for l in data[1:]:  # skip header line
            if l[0][0] == '>': # bin header row
                continue
            if l[0][:3].upper() == 'MAX':
                try:
                    lastPart.maxStart = self.parts[l[0]].start
                    lastPart.maxEnd= self.parts[l[0]].end
                except NameError:
                    raise j5Error('MAX sequence set before base sequence was defined')
                del self.parts[l[0]] # no longer need this part definition
                continue
            self.storeTargetOrderRecord(l, position)
            lastPart = self.parts[l[0]]
            position = position + 1


    def storeTargetOrderRecord(self, line, position):
        self.parts[line[0]].forwardDirection = (line[1].upper() == 'FORWARD')
        self.parts[line[0]].forcedAssemblyStrategy = line[2]
        self.parts[line[0]].forcedRelativeOverhang = stringToInt(line[3])
        self.parts[line[0]].directSynthesisFirewall = (line[4].upper() == 'TRUE')
        self.parts[line[0]].extra5pCPECOverlap = stringToInt(line[5])
        self.parts[line[0]].extra3pCPECOverlap = stringToInt(line[6])
        self.parts[line[0]].position = position


    def getMaxPosition(self):
        # Positions are zero indexed, so a 2 bin assembly returns 1
        return max([self.parts[k].position for k in self.parts.keys()])


    def getKeysByPosition(self):
        sortedKeyTuples = sorted([(k, self.parts[k].position) for k in self.parts.keys()], key=lambda tup:tup[1])
        return list(zip(*sortedKeyTuples))[0]


    def getPartsByPosition(self):
        return [self.parts[k] for k in self.getKeysByPosition()]


    def getPartByPosition(self, position):
        return self.getPartsByPosition()[self.resolveWrappedAssemblyPosition(position)]

    def getSequenceSlice(self, seq, start, end):
        if start <= end:
            return seq[start - 1:end]
        # else it wraps around plasmid, through position N:1 junction
        return (seq + seq)[start - 1:(len(seq) + end)]


    def getPartSequence(self, part):
        #returns cropped sequence in direction used in final design
        uncroppedSequence = self.sequences[part.source]
        croppedSequence = self.getSequenceSlice(uncroppedSequence, part.start, part.end)
        if not part.useInSourceDirection():
            return str(Seq(croppedSequence).reverse_complement())
        return croppedSequence


    def getRestOfPlasmidSequence(self, part):
        uncroppedSequence = self.sequences[part.source]
        croppedSequence = self.getSequenceSlice(uncroppedSequence, part.end+1, part.start-1)
        if not part.useInSourceDirection():
            return str(Seq(croppedSequence).reverse_complement())
        return croppedSequence


    def getOrderedSequences(self):
        # returns a list of sequences, in the order they are assembled
        # the sequences are cropped to desired part boundaries, not full source
        # sequences, and in the direction they will be assembled
        return [self.getPartSequence(p) for p in self.getPartsByPosition()]


    def resolveWrappedAssemblyPosition(self, position, maximum=None):
        if maximum is None:
            maximum = self.getMaxPosition()
        return position % (maximum + 1)


    def getRelativePosition(self, position, offset):
        return self.resolveWrappedAssemblyPosition(position+offset)


    def getNextPosition(self, position):
        return self.getRelativePosition(position, 1)


    def getPreviousPosition(self, position):
        return self.getRelativePosition(position, -1)


    def getRelativePart(self, part, offset):
        return self.getPartByPosition(part.position+offset)


    def getNextPart(self, part):
        return self.getRelativePart(part, 1)


    def getPreviousPart(self, part):
        return self.getRelativePart(part, -1)



    def getFinalSequence(self):
        return ''.join(self.getOrderedSequences())


    def breakPart(self, position):
        # breaks part at position into two parts of about equal size (first part is 1 bp larger for odd part lengths)
        # does not try to pick an optimal splitting point for assembly
        part = self.getPartByPosition(position)
        length = len(self.getPartSequence(part))
        if length < 2:
            return # refuse silently
        for p in self.parts.values():
            if p.position > position:
                p.position += 1
        oldEnd = part.end
        part.end = resolveWrappedSequencePosition(
            int(math.ceil(length / 2)) + part.start-1, len(self.sequences[part.source]))
        newKey = part.keyName + "_"+str(randomDigits(10))
        newPart = self.parts[newKey] = copy.deepcopy(part)
        newPart.keyName = newKey
        if not part.useInSourceDirection(): #need to swap order as this part is backwards
            part.position += 1
        else:
            newPart.position += 1
        newPart.start = resolveWrappedSequencePosition(
            part.end + 1, len(self.sequences[newPart.source]))
        newPart.end = oldEnd


    def breakAllBigParts(self):
        for p in list(self.parts.values()):
            if len(self.getPartSequence(p)) > self.maxPartLength:
                self.breakPart(p.position)


    def reducePartsSilently(self):
        # tries to reduce the number of parts without changing the sequence of the resulting assembly
        # if adjacent parts already are adjoining in a source sequence, then
        # combine them into a single part
        if self.getMaxPosition() < 1:  # need atleast 2 parts do do a reduction
            return
        shifts = self.findSilentJunctionShifts()
        for currentPart, currentShift in zip(self.getPartsByPosition(), shifts):
            nextPart = self.getNextPart(currentPart)
            lengthNextPart = len(self.getPartSequence(nextPart))
            lengthCurrentPart = len(self.getPartSequence(currentPart))
            if (currentShift[2] >= lengthNextPart) and (lengthNextPart + lengthCurrentPart <= self.maxPartLength):
                # grow current part by length of next part on 3p end
                self.modifyPartSize(currentPart, 0, lengthNextPart)
                self.eliminatePart(nextPart) 
                self.reducePartsSilently()
                break

    def findSilentJunctionShifts(self):
        # returns a list of tuples [(extendReverse, matchReverseSequence, extendForward, matchForwardSequence) ... ]
        # returnValue[N] corresponds to parts[x].position=N
        # where extendReverse is the number of bp 5' of the Nth part that could be
        # included in this part without changing the resulting sequence of the assembly
        # and extendForward is the number of bp 3' of the Nth part that could be
        # included in this part without changing the resulting sequence of the
        # assembly
        shifts = []
        for part in self.getPartsByPosition():
            nextPartSequence = self.getPartSequence(self.getNextPart(part))
            previousPartSequence = self.getPartSequence(self.getPreviousPart(part))
            restOfPlasmid = self.getRestOfPlasmidSequence(part)
            matchForwardLength = matchLength(restOfPlasmid, nextPartSequence)
            matchForwardSequence = restOfPlasmid[:matchForwardLength]
            matchReverseLength = matchLength(stringReverse(restOfPlasmid), stringReverse(previousPartSequence))
            matchReverseSequence = stringReverse(restOfPlasmid)[:matchReverseLength]
            displayLen = 10
            shifts.append(
                (matchReverseLength, matchReverseSequence[:displayLen],
                matchForwardLength, matchForwardSequence[:displayLen]))
        return shifts

    def eliminatePart(self, part):
        self.eliminatePosition(part.position)

    def eliminatePosition(self, delPosition):
        # delPosition is the target part position with zero indexed ordering
        for k in list(self.parts.keys()):
            if self.parts[k].position == delPosition:
                del self.parts[k]
        for k in self.parts.keys():
            if self.parts[k].position > delPosition:
                self.parts[k].position -= 1


    def shiftJunction(self, position, distance):
        part = self.getPartByPosition(position)
        shiftPartJunction(part, distance)


    def shiftPartJunction(self, part, distance):
        # position is the target part order, zero indexed
        # distance is how many base pairs to extend (positive) or contract (negative) the 3' part boundary
        # parts is the dictionary of Part objects
        # sequences is dictionary of sequences

        # if position=N and distance=X where X>0, then part at position N grows by X bp on 3' end and
        # part at position N+1%NumParts shrinks by X bp on 5' end
        # elif position=N and distance=X where X<0, then part at position N shrinks by X bp on 3' end and
        # part at position N+1%NumParts grows by X bp on 5' end

        # shifting a junction will change the resulting assembly sequence, except when the adjoining parts
        # each have the same sequence in the region of the junction shift
        self.modifyPartSize(part, 0, distance)
        self.modifyPartSize(self.getNextPart(part), -1 * distance, 0)

    def modifyPartSize(self, part, delta5p, delta3p):
        # position is the target part order, zero indexed
        # delta5p is how many base pairs to extend (positive) or contract (negative) the 5' part boundary
        # delta3p is how many base pairs to extend (positive) or contract (negative) the 3' part boundary
        #          3' and 5' are relative to the part orientation in the final assembly
        # parts is the list of Part objects -- parts is modified by this
        # function.

        # there are two parameters that determine the orientation of a Part
        #        1) Part.reverseComplement
        #        2) Part.forwardDirection
        # a part is used in the same direction as the source material (genbank,jbei-seq,fasta, etc)
        # when Part.reverseComplement != Part.forwardDirection
        if part.useInSourceDirection():
            part.start -= delta5p
            part.end += delta3p
        else:
            part.start -= delta3p
            part.end += delta5p
        # need to check if our position numbers need to wrap around
        sourceLength = len(self.sequences[part.source])
        part.start = resolveWrappedSequencePosition(part.start, sourceLength)
        part.end = resolveWrappedSequencePosition(part.end, sourceLength)


if __name__ == "__main__":
    print("Python library for working with j5")

