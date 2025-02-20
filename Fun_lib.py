#!/usr/bin/env python
import os
import re
import subprocess
import pandas as pd
import numpy as np
import sys
import itertools
from copy import deepcopy
import copy

globalExonBoundaryOffset=100 # There should be gap of 100 nucleotides from the trimmed boundaries
minimumExonLength=200 # new chopped exons length should be atleast 200nt

class Edge(object):
    """Stores information about an edge, including its location
       and the gene/transcript(s) it belongs to.
       Attributes:
           identifier: Accession ID of the edge
           gene: Accession ID of the gene that the edge belongs to
           transcript_ids: Set of transcript accession IDs that the edge 
           belongs to
           chromosome: Chromosome that the transcript is located on 
           (format "chr1")
           start: The start position of the edge with respect to the
           forward strand 
           end: The end position of the edge with respect to the
           forward strand
           strand: "+" if the edge is on the forward strand, and "-" if
           it is on the reverse strand
 
           length: The length of the edge
    """

    def __init__(self, identifier, chromosome, start, end, strand, gene_id,
                 transcript_id, annotations):
        self.chromosome = str(chromosome)
        self.gene_id = gene_id
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.length = abs(self.end - self.start + 1)
        self.annotations = annotations

        self.identifier = str(identifier)
        self.transcript_ids = set()
        if transcript_id != None:
            self.transcript_ids.add(transcript_id)
        self.v1 = None
        self.v2 = None

    def print_edge(self):
        """ Prints a string representation of the edge"""
        print(self.identifier + ": " + self.chromosome + ":" + \
              str(self.start) + "-" + str(self.end))
        print(self.transcript_ids)
        return
    
    def getExonLineByObject(exonObj):
        exonLine=str(exonObj.chromosome+"\t"+exonObj.annotations.get('source')+"\t"+"exon"+"\t"+str(exonObj.start)+"\t"+str(exonObj.end)+"\t.\t"+exonObj.strand+"\t.\t")
        for x in exonObj.annotations:
            exonLine=exonLine+str(x+" \""+exonObj.annotations[x]+"\"; ")
        return exonLine+"\n"

    def create_edge_from_gtf(edge_info):
        """ Creates an edge object using information from a GTF entry
                Args:
                   edge_info: A list containing fields from a GTF file edge entry.
                   Example:   
                   ['chr1', 'HAVANA', 'exon', '11869', '12227', '.', '+', '.', 
                    'gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; 
                    gene_type "transcribed_unprocessed_pseudogene"; 
                    gene_status "KNOWN"; gene_name "DDX11L1"; 
                    transcript_type "processed_transcript"; 
                    transcript_status "KNOWN"; transcript_name "DDX11L1-002"; 
                    edge_number 1; edge_id "ENSE00002234944.1"; level 2; 
                    tag "basic"; transcript_support_level "1"; 
                    havana_gene "OTTHUMG00000000961.2"; 
                    havana_transcript "OTTHUMT00000362751.1";'] 
        """
        description = edge_info[-1]
        start = int(edge_info[3])
        end = int(edge_info[4])
        chromosome = edge_info[0]
        strand = edge_info[6]

        annotations = extract_edge_annotations_from_GTF(edge_info)
        if "exon_id" not in annotations:
            annotations["exon_id"] = "_".join([chromosome, str(start), str(end), strand])
        gene_id = annotations['gene_id']
        transcript_id = annotations['transcript_id']
        edge_id = "_".join([chromosome, str(start), str(end), strand,gene_id])

        if "gene_id" in description:
            gene_id = (description.split("gene_id ")[1]).split('"')[1]
        if "transcript_id" in description:
            transcript_id = (description.split("transcript_id ")[1]).split('"')[1]

        edge = Edge(edge_id, chromosome, start, end, strand, gene_id, transcript_id,
                    annotations)
        return edge

def extract_edge_annotations_from_GTF(tab_fields):
    """ Extracts key-value annotations from the GTF description field
    """

    attributes = {}

    # remove trailing newline and split by semicolon
    description = tab_fields[-1].strip('\n')
    description = description.split(';')

    # Parse description
    for fields in description:
        if fields == "" or fields == " ": continue
        fields = fields.split()
        if fields[0] == '': fields = fields[1:]

        key = fields[0].replace('"', '')
        val = ' '.join(fields[1:]).replace('"', '')
        
        attributes[key] = val

    # Put in placeholders for important attributes (such as gene_id) if they
    # are absent
    if "gene_id" not in attributes:
        attributes["gene_id"] = "NULL"
    if "transcript_id" not in attributes:
        attributes["transcript_id"] = "NULL"

    attributes["source"] = tab_fields[1]

    return attributes

def get_edge_from_db(vertex_info_1, vertex_info_2):
    """ Uses information from a database edge entry to create an edge object.
    """
    if vertex_info_1["edge_id"] != vertex_info_2["edge_id"]:
        raise ValueError('Tried to create edge from endpoints with different IDs')
    edge_id = vertex_info_1["edge_id"]
    chromosome = vertex_info_1['chromosome']
    start = min(vertex_info_1['position'], vertex_info_2['position'])
    end = max(vertex_info_1['position'], vertex_info_2['position']) 
    strand = vertex_info_1['strand']
    gene_id = vertex_info_1['gene_id']

    edge = Edge(edge_id, chromosome, start, end, strand, gene_id, None, None)
    edge.v1 = str(vertex_info_1["vertex_ID"])
    edge.v2 = str(vertex_info_2["vertex_ID"])
    return edge

def create_novel_edge(chromosome, start, end, strand, gene_id, transcript_id, counter):
    """ Creates a novel edge with a unique identifier (obtained using
        counter). Returns the edge object as well as the updated counter.
    """
    counter["edges"] += 1
    curr_novel = counter["edges"]
    edge = Edge(curr_novel, chromosome, start, end, strand, gene_id, transcript_id,
                None)
    return edge

class Gene(object):
    """ Contains high-level information about a gene, such as its identifiers, 
        genomic location, and transcripts. Does not contain exon information.
        Attributes:
            - identifier: Accession ID of gene, i.e. an Ensembl ID. Required.
            - name: Human-readable name of the gene. This attribute can be left 
              empty if the gene does not have an assigned name.
            - chromosome: Chromosome that the gene is located on (format "chr1")
            - start: The start position of the gene with respect to the forward 
              strand (int). Should always be less than or equal to end.
            - end: The end position of the gene with respect to the forward strand 
              (int). Should always be greater than or equal to start.
            - strand: "+" if the gene is on the forward strand, "-" if it is on 
              the reverse strand
            - annotations: a dictionary of miscellaneous annotation categories
              extracted from a GTF
            
    """

    def __init__(self, identifier, chromosome, start, end, strand, annotations):
        start = int(start)
        end = int(end)

        self.identifier = str(identifier)
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.transcripts = {}
        self.length = end - start + 1
        self.annotations = annotations

        if start > end:
            raise ValueError("""Plus strand gene start must be less than or 
                             equal to end.""")

    def set_name(self, name):
        """ Sets the name attribute of the Gene to the provided value.
        """
        self.annotations['name'] = name
        return

    def add_transcript(self, transcript):
        """ Adds a key-value pair (transcript identifier -> Transcript oject)
            to the gene's transcript dictionary
            Args:
                transcript: object of type Transcript. Must overlap with the 
                location of the gene.
        """
        if transcript.start >= self.end or transcript.end <= self.start:
            raise ValueError('Transcript must overlap the gene it is assigned to')
 
        if transcript.gene_id == self.identifier:
            # In order to belong to a gene, the transcript gene_id must match
            transcript_id = transcript.identifier
            self.transcripts[transcript_id] = transcript
        else:
            raise ValueError('Gene ID of transcript must match gene ' + \
                  'in order for assignment to be made.')
        return             


    def print_gene(self):
        """ Print a string representation of the Gene. Good for debugging. """

        if "name" in self.annotations != "":
            # Include name in output if there is one
            print(self.identifier + " (" + self.annotations['name']  + "):")
        else:
            print(self.identifier + ":")

        print("\tLocation: " + self.chromosome + ":" + str(self.start) + "-" + \
              str(self.end) + "(" + self.strand + ")")
        
        # Print transcripts in shorthand 
        for transcript in self.transcripts:
            print("\t Transcript: " + transcript)

        return
   
    
    def getGeneLineByObject(gene):
        geneLine=str(gene.chromosome+"\t"+gene.annotations.get('source')+"\t"+"gene"+"\t"+str(gene.start)\
         +"\t"+str(gene.end)+"\t.\t"+gene.strand+"\t.\t")
        for x in gene.annotations:
            geneLine=geneLine+str(x+" \""+gene.annotations[x]+"\"; ")
        return geneLine+"\n"

    def get_gene_from_gtf(gene_info):
        """ Creates a Gene object from a GTF file entry
            Args:
                gene_info: A list containing fields from a GTF file gene entry.
                Example:
                ['chr1', 'HAVANA', 'gene', '11869', '14409', '.', '+', '.',
                'gene_id "ENSG00000223972.5";
                gene_type "transcribed_unprocessed_pseudogene";
                gene_status "KNOWN"; gene_name "DDX11L1"; level 2;
                havana_gene "OTTHUMG00000000961.2";']
        """
        chromosome = gene_info[0]
        start = int(gene_info[3])
        end = int(gene_info[4])
        strand = gene_info[6]
        annotations = extract_gene_annotations_from_GTF(gene_info)
        if "gene_id" not in gene_info[-1]:
                raise ValueError('GTF entry lacks a gene_id field')
        gene_id = annotations['gene_id']

        gene = Gene(gene_id, chromosome, start, end, strand, annotations)
        return gene

class Transcript(object):
    """Stores information about a gene transcript, including its location
       and constitutive exons.
       Attributes:
           identifier: Accession ID of transcript, i.e. an Ensembl ID. Must
           be unique.
           name: Human-readable name of the transcript. Does not have to be 
           unique
           chromosome: Chromosome that the transcript is located on 
           (format "chr1")
           start: The start position of the transcript with respect to the
           forward strand 
           end: The end position of the transcript with respect to the
           forward strand
           strand: "+" if the transcript is on the forward strand, and "-" if
           it is on the reverse strand
           gene_id: unique ID of the gene that this transcript belongs to
           exons: List of exon objects belonging to this transcript, in sorted
           order.
    """

    def __init__(self, identifier, chromosome, start, end, strand, gene_id, 
                 annotations):

        self.identifier = str(identifier)
        self.gene_id = str(gene_id)

        self.chromosome = str(chromosome)
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.n_exons = 0
        self.exons = []
        self.introns = []
        self.annotations = annotations

    def get_5prime_vertex(self):
        """ Returns ID of 5' end vertex """

        if self.strand == "+":
            return self.exons[0].v1
        if self.strand == "-":
            return self.exons[-1].v2

    def get_3prime_vertex(self):
        """ Returns ID of 5' end vertex """

        if self.strand == "+":
            return self.exons[-1].v2
        if self.strand == "-":
            return self.exons[0].v1

    def get_edge_path(self):
        edges = self.get_all_edges()
        if len(edges) == 0:
            return None
        path = [ x.identifier for x in edges]

        # Must reverse the path if the transcript is on the '-' strand
        if self.strand == "-":
            path = path[::-1]
        return ",".join(path)

    def get_all_edges(self):
        all_edges = []
        for i in range(0,self.n_exons):
            all_edges.append(self.exons[i])
            try:
                all_edges.append(self.introns[i])
            except:
                pass
            
        return all_edges

    def get_length(self):
        """ Computes the length of the transcript by summing the lengths of
            its exons """

        if len(self.exons) == 0:
            raise ValueError('Cannot compute length: Transcript does not ' + \
                             'have any exons')
        
        transcript_length = 0
        for exon in self.exons:
            transcript_length += exon.length
        return transcript_length

    def get_exon_coords(self):
        """ Returns a list of the exon coordinates in order """
        exon_coords = []
        for exon in self.exons:
            exon_coords.append(int(exon.start))
            exon_coords.append(int(exon.end))
        return exon_coords

    def add_exon(self, exon):
        """Adds an exon object to the transcript."""

        if exon.start > exon.end:
            raise ValueError('Exon start (' + str(exon.start) + ') ' + \
                'is supposed to be before the exon end (' + str(exon.end) + ')')

        # Check where in the list the exon should be added
        for i in range(0,len(self.exons)):
            existing_exon = self.exons[i]
            if exon.end < existing_exon.start:
                self.exons = self.exons[0:i] + [exon] + self.exons[i:]
                self.check_exon_validity()
                self.n_exons += 1
                return
        self.exons.append(exon)
        self.check_exon_validity()
        self.n_exons += 1
        return

    def add_intron(self, intron):
        """Adds an edge object to the transcript."""

        if intron.start > intron.end:
            raise ValueError('Intron start (' + str(intron.start) + ')' + \
                'is supposed to be before the intron end (' + str(intron.end) + ')')

        # Check where in the list the intron should be added
        for i in range(0,len(self.introns)):
            existing_intron = self.introns[i]
            if intron.end < existing_intron.start:
                self.introns = self.introns[0:i] + [intron] + self.introns[i:]
                return
        self.introns.append(intron)
        return
                    
    def check_exon_validity(self):
        """ The transcript's exons are valid if:
            1) Exons are in sorted order (ascending)
            2) Exon bounds do not exceed transcript start and end
            3) Exons are all on the appropriate chromosome
            If these conditions are violated, this function raises an error.
        """
        prev = 0
        for exon in self.exons:
            if exon.chromosome != self.chromosome:
                raise ValueError('Invalid exon in transcript ' + \
                      self.identifier + ': wrong chromosome')
            if exon.start < self.start or exon.end > self.end:
                print("self.start: " + str(self.start))
                print("self.end: " + str(self.end))
                print("exon.start: " + str(exon.start))
                print("exon.end: " + str(exon.end))
                raise ValueError('Invalid exon in transcript ' + \
                      self.identifier + ': (' + str(exon.start) + "-" + \
                      str(exon.end) + \
                      ') is located beyond start or end of transcript')
            if exon.start <= prev:
                # This error would indicate a TALON bug rather than user error,
                # so we shouldn't see it. 
                raise ValueError('Exons of transcript ' + \
                      self.identifier + ' are not stored in ascending order.')
            prev = exon.end
        return

    def get_introns(self):
        """
        Computes introns based on the exon list
        """
        exon_coords = self.get_exon_coords()
        intron_list = []

        i = 1
        while (i < len(exon_coords) - 1):
            j = i + 1

            intron_list.append(exon_coords[i] + 1)
            intron_list.append(exon_coords[j] - 1)
            i += 2

        return intron_list


    def print_transcript(self):
        """ Print a string representation of the Transcript. Good for debugging
        """
        transcript_id = self.identifier
        if transcript_id == None:
            transcript_id = "Transcript"

        print("\tLocation: " + self.chromosome + ":" + str(self.start) + "-" + \
              str(self.end) + "(" + self.strand + ")")

        # Print exons
        print("\tExons: " + "\n".join([str(x.start) + "-" + str(x.end) for x in self.exons]))
        return 
    
    def getTranscriptLineByObject(transcriptObj):
        transcriptLine=str(transcriptObj.chromosome+"\t"+transcriptObj.annotations.get('source')+"\t"+"transcript"+"\t"+str(transcriptObj.start)\
     +"\t"+str(transcriptObj.end)+"\t.\t"+transcriptObj.strand+"\t.\t")
        for x in transcriptObj.annotations:
            transcriptLine=transcriptLine+str(x+" \""+transcriptObj.annotations[x]+"\"; ")
        return transcriptLine+"\n"
    

    def get_transcript_from_gtf(transcript_info):
        """ Uses information from a GTF-formatted transcript entry to create a
        Transcript object.
            Args:
                transcript_info: A list containing fields from a GTF file gene 
                entry. Example:

                chr1	HAVANA	transcript	12010	13670	.	+
                .	gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2"; 
                gene_type "transcribed_unprocessed_pseudogene"; 
                gene_status "KNOWN"; gene_name "DDX11L1"; 
                transcript_type "transcribed_unprocessed_pseudogene"; 
                transcript_status "KNOWN"; transcript_name "DDX11L1-001"; 
                level 2; ont "PGO:0000005"; ont "PGO:0000019"; tag "basic"; 
                transcript_support_level "NA"; havana_gene "OTTHUMG00000000961.2"; 
                havana_transcript "OTTHUMT00000002844.2";
        """
        chromosome = transcript_info[0]
        start = int(transcript_info[3])
        end = int(transcript_info[4])
        strand = transcript_info[6]

        if "transcript_id" not in transcript_info[-1]:
                raise ValueError('GTF entry lacks a transcript_id field')
        annotations = extract_transcript_annotations_from_GTF(transcript_info)


        gene_id = annotations['gene_id']
        transcript_id = annotations['transcript_id']

        transcript = Transcript(transcript_id, chromosome, start, end, strand, 
                                gene_id, annotations)

        return transcript

def extract_transcript_annotations_from_GTF(tab_fields):
    """Extracts key-value annotations from the GTF description field"""

    attributes = {}

    # remove trailing newline and split by semicolon
    description = tab_fields[-1].strip("\n")
    description = description.split(";")

    # Parse description
    for fields in description:
        if fields == "" or fields == " ":
            continue
        fields = fields.split()
        if fields[0] == "":
            fields = fields[1:]

        key = fields[0].replace('"', "")
        val = " ".join(fields[1:]).replace('"', "")

        attributes[key] = val

    # Put in placeholders for important attributes (such as gene_id) if they
    # are absent
    if "gene_id" not in attributes:
        attributes["gene_id"] = "NULL"

    attributes["source"] = tab_fields[1]

    return attributes
 


    def get_transcript_from_exon(exon, gene_id, transcript_id):
        """ In rare cases, GTF exons are listed with gene and transcript IDs that
            do not have corresponding entries. In this case, we create a transcript
            for this exon for bookkeeping purposes."""

        name = transcript_id
        chromosome = exon.chromosome
        start = exon.start
        end = exon.end
        strand = exon.strand
        transcript = Transcript(transcript_id, name, None, chromosome, start, end,
                                strand, gene_id)
        return transcript

    def create_novel_transcript(chromosome, start, end, strand, gene_id, counter,
                                 exons, introns):
        """ Creates a novel transcript with a unique identifier (obtained using
            counter). Returns the transcript object as well as the updated counter.
        """
        counter["transcripts"] += 1
        transcript_id = str(counter["transcripts"])

        transcript = Transcript(transcript_id, chromosome, start, end, strand, 
                                gene_id, None)

        for exon in exons:
            transcript.add_exon(exon)
        for intron in introns:
            transcript.add_intron(intron)

        return transcript

def extract_gene_annotations_from_GTF(tab_fields):
    """Parses the description field of a gene GTF in order to organize the 
       information therein into a dictionary.
    """

    attributes = {}

    # remove trailing newline and split by semicolon
    description = tab_fields[-1].strip('\n')
    description = description.split(';')

    # Parse description
    for fields in description:
        if fields == "" or fields == " ": continue
        fields = fields.split()
        if fields[0] == '': fields = fields[1:]

        key = fields[0].replace('"', '')
        val = ' '.join(fields[1:]).replace('"', '')

        attributes[key] = val

    attributes["source"] = tab_fields[1]

    return attributes  

    def get_gene_from_exon(exon, gene_id):
        """ In rare cases, GTF exons are listed with gene and transcript IDs that
            do not have corresponding entries. In this case, we create a gene
            for this exon for bookkeeping purposes."""

        gene_name = gene_id
        chromosome = exon.chromosome
        start = exon.start
        end = exon.end
        strand = exon.strand
        gene = Gene(gene_id, gene_name, None, chromosome, start, end, strand)
        return gene

    def create_novel_gene(chromosome, start, end, strand, counter):
        """ Creates a novel gene with a unique identifier (obtained using
            counter). Returns the gene object as well as the updated counter.
        """
        gene_id = str(counter["genes"] + 1)
        counter["genes"] += 1
        gene = Gene(gene_id, chromosome, start, end, strand, None)
        return gene


def updateLncRNAExonSingleCase(resultingIntervalsDf,lncRNAStart,lncRNAEnd):
    updatedLncRNAStart=0
    updatedLncRNAEnd=0
    resultingExons=[]
    row=resultingIntervalsDf.iloc[0]
    eStart=row['start']
    eEnd=row['end']
    if (lncRNAStart in range(eStart-1,eEnd) and lncRNAEnd in range(eStart,eEnd+1)):
        updatedLncRNAStart=0
        updatedLncRNAEnd=0
    elif lncRNAStart < eStart: # Use the orignal lncRNA start
        updatedLncRNAStart=lncRNAStart
        updatedLncRNAEnd=eStart-globalExonBoundaryOffset
    elif lncRNAEnd > eEnd: # Use the orignal lncRNA end
        updatedLncRNAStart=eEnd+globalExonBoundaryOffset
        updatedLncRNAEnd=lncRNAEnd
    resultingExons.append([updatedLncRNAStart,updatedLncRNAEnd])
    return resultingExons

def getProteinCodingExonsFromIndex(index,df):
    refStart=df.iloc[index]['refStart']
    refEnd=df.iloc[index]['refEnd']
    proteinCodingExons = []
    for i in range(len(refStart)):
        proteinCodingExons.append([refStart[i],refEnd[i]])
    return proteinCodingExons

# merge the overlapping proteincoding exons into a long non coding exons
def getUnionOfMultipleranges(proteinCodingExons):
    #proteinCodingExons = [[7, 10], [11, 13], [11, 15], [14, 20], [23, 39]] # should output [[7, 20], [23, 39]]
    b = []
    for begin,end in sorted(proteinCodingExons):
        if b and b[-1][1] >= begin - 1:
            b[-1][1] = max(b[-1][1], end)
        else:
            b.append([begin, end])
    return b

# Function to print intersecting intervals
def printIntervals(arr1, arr2):
	resultingIntervals = []
    # i and j pointers for arr1
	# and arr2 respectively
	i = j = 0
	
	n = len(arr1)
	m = len(arr2)

	# Loop through all intervals unless one
	# of the interval gets exhausted
	while i < n and j < m:
		
		# Left bound for intersecting segment
		l = max(arr1[i][0], arr2[j][0])
		
		# Right bound for intersecting segment
		r = min(arr1[i][1], arr2[j][1])
		
		# If segment is valid print it
		if l <= r:
			#print('{', l, ',', r, '}')
			resultingIntervals.append([l,r])
            

		# If i-th interval's right bound is
		# smaller increment i else increment j
		if arr1[i][1] < arr2[j][1]:
			i += 1
		else:
			j += 1
	resultingIntervalsDf = pd.DataFrame(resultingIntervals, columns = ['start', 'end'])
	resultingIntervalsDf=resultingIntervalsDf.drop_duplicates() # sometimes same intervals are returned if multiple protein coding exons overlap in the same region
	return resultingIntervalsDf

def getGeneWithTranscriptsAndExonsAsLines(gene):
    geneLine=Gene.getGeneLineByObject(gene)
    geneWithTranscriptsAndExonsLines=geneLine
    #geneWithTranscriptAndExons=geneLine
    for transcriptID in gene.transcripts:
        transcriptObj=gene.transcripts.get(transcriptID)
        transcriptLine=Transcript.getTranscriptLineByObject(transcriptObj)
        geneWithTranscriptsAndExonsLines=geneWithTranscriptsAndExonsLines+transcriptLine
        for exonObj in transcriptObj.exons:
            exonLine=Edge.getExonLineByObject(exonObj)
            geneWithTranscriptsAndExonsLines=geneWithTranscriptsAndExonsLines+exonLine
    return geneWithTranscriptsAndExonsLines


def updateLncRNAExonMultipleCase(resultingIntervalsDf,lncRNAStart,lncRNAEnd):
    updatedLncRNAStart=0
    updatedLncRNAEnd=0
    
    prevStart=0
    prevEnd=0
    resultingMultipleExons = [] 
    ### check if the lncRNA is contained completely inside the overlapping intervals. Meaning we cannot chop it.
    for index, row in resultingIntervalsDf.iterrows():
        eStart=row['start']
        eEnd=row['end']
        if (lncRNAStart in range(eStart-1,eEnd) and lncRNAEnd in range(eStart,eEnd+1)):
            updatedLncRNAStart=0
            updatedLncRNAEnd=0
            resultingMultipleExons.append([updatedLncRNAStart,updatedLncRNAEnd])
            return resultingMultipleExons
    ### If reached here, then chop it off.
    for i in range(len(resultingIntervalsDf)-1):
        currentRow=resultingIntervalsDf.iloc[i]
        nextRow=resultingIntervalsDf.iloc[i+1]
        resultingMultipleExons.append([currentRow['end']+globalExonBoundaryOffset,nextRow['start']-globalExonBoundaryOffset])
    ## Add Start and end if lncRNA is starting and/or ending after the intervals
    minOfStarts=resultingIntervalsDf['start'].min()
    maxOfEnds=resultingIntervalsDf['end'].max()
    if(lncRNAStart<minOfStarts):
         resultingMultipleExons.append([lncRNAStart,minOfStarts-globalExonBoundaryOffset])
    if(lncRNAEnd>maxOfEnds):
         resultingMultipleExons.append([maxOfEnds+globalExonBoundaryOffset,lncRNAEnd])
    
    resultingMultipleExons = sorted(resultingMultipleExons, key=lambda x: x[1])
    return resultingMultipleExons


# Create New Exons. Chop lncRNA exons or assign [0,0] to delete them
def getNewExons(df,minimumExonLength):
    df['proteinCodingExonsMerged'] = np.empty((len(df), 0)).tolist() # For merged protein coding exons
    df['newExonsWithNoLengthLimit'] = np.empty((len(df), 0)).tolist()# add empty column for new exons
    df['newExons'] = np.empty((len(df), 0)).tolist() # add empty column for new exons
    df['numberOfNewExons'] = 0 # add empty column for numberOfNewExons
    for index in range(len(df)):
        proteinCodingExons = getProteinCodingExonsFromIndex(index,df)
        lncRNAStart=df.iloc[index]['queryStart']
        lncRNAEnd=df.iloc[index]['queryEnd']
        lncRNAExon = [ [lncRNAStart , lncRNAEnd ] ]
        # obtain Union of multiple ranges. This will combine all the overlapping protein coding exons 
        # into one long exon (if they have overlapping boundires or just starting 1 nt after each other.)
        proteinCodingExons=getUnionOfMultipleranges(proteinCodingExons)
        df.at[index,'proteinCodingExonsMerged'] = proteinCodingExons
        resultingIntervals=printIntervals(proteinCodingExons, lncRNAExon)
        newExons=[]
        if(len(resultingIntervals)==1):
            #print("-----------------------Single----------------------------")
            newExons=updateLncRNAExonSingleCase(resultingIntervals,lncRNAStart,lncRNAEnd)
        else:
            #print("-------------------------Multiple--------------------------")
            newExons=updateLncRNAExonMultipleCase(resultingIntervals,lncRNAStart,lncRNAEnd)
        df.at[index,'newExonsWithNoLengthLimit'] = newExons
        
        # remove any exons with length less than minimumExonLength
        filteredNewExon=[]
        for exon in newExons:
            start=exon[0]
            end=exon[1]
            if(end-start>=minimumExonLength):
                filteredNewExon.append([start,end])
        
        df.at[index,'newExons'] = filteredNewExon
        df.at[index,'numberOfNewExons'] = len(filteredNewExon)

    df.sort_values(by=['queryGeneId'])
    return df

def antiSenseExonOverlapCompressedBedFile_antiSenseExonsChoppedDfFilePath(antiSenseExonOverlapBedFile, outputDir):
    antiSenseExonOverlapCompressedBedFile = outputDir+"/antiSenseExonOverlapBedFileCompressed.csv"
    antiSenseExonsChoppedDfFilePath = outputDir+"/antiSenseExonsChoppedDf.csv"
    
    antiSenseExonOverlapBed = pd.read_csv(antiSenseExonOverlapBedFile, sep='\t', lineterminator='\n',header=None)
    antiSenseExonOverlapBed.columns = ["queryChr","queryStart","queryEnd","queryGeneId","queryScore","queryStrand","refChr","refStart","refEnd","refGeneId","refScore","refStrand","overlapLength"]
    df=antiSenseExonOverlapBed.groupby(['queryChr', 'queryStart','queryEnd','queryGeneId','queryStrand'], as_index = False).agg({'refChr': list,'refStart': list,'refEnd': list,'refGeneId': list,'refStrand': list,'overlapLength': list,})
    df["exonId"] = df["queryChr"] +"_" + df["queryStart"].astype(str)+"_"+ df["queryEnd"].astype(str)+"_"+ df["queryStrand"]+"_"+ df["queryGeneId"]
    antiSenseExonOverlapCompressedBedFileObj = open(antiSenseExonOverlapCompressedBedFile, "w")
    df.to_csv(antiSenseExonOverlapCompressedBedFileObj, encoding='utf-8',index=False)
    df=getNewExons(df,minimumExonLength)
    antiSenseExonsChoppedDfObj = open(antiSenseExonsChoppedDfFilePath, "w")
    df.to_csv(antiSenseExonsChoppedDfObj, encoding='utf-8',index=False)
    uniqueGeneIdsOverlappingAntisense=df.queryGeneId.unique().tolist()
    return uniqueGeneIdsOverlappingAntisense, df

def getRemainingAntiSenseToBeDeleted(uniqueGeneIdsOverlappingSenseStrand,uniqueGeneIdsOverlappingAntisense,alreadyDeletedInSenseFile):
    print("unique GeneIds Overlapping Sense Strand: ",len(uniqueGeneIdsOverlappingSenseStrand))
    print("unique GeneIds Overlapping Antisense Strand: ",len(uniqueGeneIdsOverlappingAntisense))
    remainingAntiSenseToBeDeleted=set(uniqueGeneIdsOverlappingAntisense) - set(uniqueGeneIdsOverlappingSenseStrand)
    alreadyDeletedInSense = [gene for gene in uniqueGeneIdsOverlappingAntisense if gene in uniqueGeneIdsOverlappingSenseStrand]
    print("Genes that overlap on both sense and anti-sense strand or alreadyDeletedInSense: ",len(alreadyDeletedInSense))
    print("Remaining Anti Sense Genes To Be Deleted: ",len(remainingAntiSenseToBeDeleted))
    if(len(alreadyDeletedInSense)+len(remainingAntiSenseToBeDeleted)!=len(uniqueGeneIdsOverlappingAntisense)):
         sys.exit("Error!!! Already Deleted In Sense and Remaining Anti Sense Genes To Be Deleted are not equal to total \
        number of genes overlapping in antisense",len(alreadyDeletedInSense)+len(remainingAntiSenseToBeDeleted))
    
    print("Success: Already Deleted In Sense and Remaining Anti Sense Genes To Be Deleted are equal to total \
        number of genes overlapping in antisense",len(alreadyDeletedInSense)+len(remainingAntiSenseToBeDeleted))
    ### Write to a text file for a possible future use
    alreadyDeletedInSenseFileObj = open(alreadyDeletedInSenseFile, "w")
    for element in alreadyDeletedInSense:
        alreadyDeletedInSenseFileObj.write(element + "\n")
    alreadyDeletedInSenseFileObj.close()
    return remainingAntiSenseToBeDeleted


### delete all the LNCRNA genes that have overlapping exon on the same strand with protein coding gene. 
#The file is strand_LnRNA_ExonOverlap_ProteinCodingExons
def deleteLncRNAGenesOverlappingSameExonsSameStrand(senseExonOverlapBedFile, LncBookGenes, outputDir):
    senseExonOverlapBed = pd.read_csv(senseExonOverlapBedFile, sep='\t', lineterminator='\n',header=None)
    senseExonOverlapBed.columns = ["queryChr","queryStart","queryEnd","queryGeneId","queryScore","queryStrand","refChr","refStart","refEnd","refGeneId","refScore","refStrand","overlapLength"]
    uniqueGeneIdsOverlappingSenseStrand=set(senseExonOverlapBed.queryGeneId.unique().tolist())
    
    sameStrandDeletedGenesGTFFile=outputDir+"/same_Strand_Deleted_Genes.gtf"
    cleanedSameStrandExonsGTFFile=outputDir+"/cleaned_Same_Strand_Exons.gtf"
    
    # write new GTF files 
    sameStrandDeletedGenesGTFFileObj = open(sameStrandDeletedGenesGTFFile, "w")
    cleanedSameStrandExonsGTFFileObj = open(cleanedSameStrandExonsGTFFile, "w")

    # keep counts to cross check the consistancy of lines in the files
    geneCountToBeDelted=0
    geneCountToBeKept=0
    for gendId in LncBookGenes.keys():
        gene=LncBookGenes.get(gendId)
        geneWithTranscriptsAndExonsLines=getGeneWithTranscriptsAndExonsAsLines(gene)
        if gendId in uniqueGeneIdsOverlappingSenseStrand:
            geneCountToBeDelted=geneCountToBeDelted+1;
            sameStrandDeletedGenesGTFFileObj.write(geneWithTranscriptsAndExonsLines)
        else:
            geneCountToBeKept=geneCountToBeKept+1
            cleanedSameStrandExonsGTFFileObj.write(geneWithTranscriptsAndExonsLines)
    
    if(len(uniqueGeneIdsOverlappingSenseStrand)==geneCountToBeDelted):
        print("Number of genes to be deleted are equal to the deleted ones! Success",len(uniqueGeneIdsOverlappingSenseStrand),"\tand\t",geneCountToBeDelted)
    else:
        sys.exit("-------Error Number of genes to be deleted are equal to the deleted ones--------"+str(len(uniqueGeneIdsOverlappingSenseStrand))+"\tand\t"+str(geneCountToBeDelted))
    return uniqueGeneIdsOverlappingSenseStrand

def read_gtf_file(gtf_file):
    """ Reads gene, transcript, and edge information from a GTF file.
        Args:
            gtf_file: Path to the GTF file
        Returns:
            genes: A dictionary mapping gene IDs to corresponding gene objects
            transcripts: A dictionary mapping gene IDs to corresponding
                   transcript objects
            exons: A dictionary mapping exon IDs to corresponding edge objects
    """
    genes = {}
    transcripts = {}
    exons = {}

    with open(gtf_file) as gtf:
        for line in gtf:
            line = line.strip()

            # Ignore header
            if line.startswith("#"):
                continue

            # Split into constitutive fields on tab
            tab_fields = line.split("\t")
            chrom = tab_fields[0]
            entry_type = tab_fields[2]

            # Entry is a gene
            if entry_type == "gene":
                gene = Gene.get_gene_from_gtf(tab_fields)
                native_id = gene.identifier
                genes[native_id] = gene

            # Entry is a transcript
            elif entry_type == "transcript":
                transcript = Transcript.get_transcript_from_gtf(tab_fields)
                gene_id = transcript.gene_id
                if gene_id in genes:
                    genes[gene_id].add_transcript(transcript)
                native_id = transcript.identifier
                transcripts[native_id] = transcript

            # Entry is an edge
            elif entry_type == "exon":
                exon = Edge.create_edge_from_gtf(tab_fields)
                # This ID is used because of a rare GENCODE bug
                location_exon_id = exon.identifier
                exons[location_exon_id] = exon

                transcript_id = list(exon.transcript_ids)[0]
                gene_id = exon.annotations["gene_id"]

                if location_exon_id not in exons:
                    # Add the new edge to the data structure
                    exons[location_exon_id] = exon
                else:
                    # Update existing exon entry, including its transcript set
                    exon = exons[location_exon_id]
                    exon.transcript_ids.add(transcript_id)

                if transcript_id in transcripts:
                    currTranscript = transcripts[transcript_id]
                    currTranscript.add_exon(exon)

    return genes, transcripts, exons

def writeGeneToTheFile(singleGene,gtfFilePath):
    gtfFileObj = open(gtfFilePath, "a")
    for i in range(len(singleGene)):
        gtfFileObj.write(singleGene[i])
    gtfFileObj.close()

def writeGeneToTheFile(singleGene,gtfFilePath):
    gtfFileObj = open(gtfFilePath, "a")
    for i in range(len(singleGene)):
        gtfFileObj.write(singleGene[i])
    gtfFileObj.close()
    
def validateGene(gene):
    if len(gene.transcripts)==0:
        raise ValueError("Gene : \t",gene.identifier, "\t has no Transcripts")
    if gene.start > gene.end:
        raise ValueError("Gene : \t",gene.identifier ,"\t start > end")
    for transcriptID in list(gene.transcripts):
        transcriptObj=gene.transcripts.get(transcriptID)
        if len(transcriptObj.exons)==0:
            raise ValueError("Transcript : \t",transcriptObj.identifier ,"\t has no exons")
        if transcriptObj.start > transcriptObj.end:
            raise ValueError("Transcript : \t",transcriptObj.identifier, "\t start > end")
        for exon in transcriptObj.exons:
            if exon.start>exon.end:
                raise ValueError("exon : \t",exon.identifier ,"\t start > end")
            if exon.start==0 or exon.end==0:
                raise ValueError("exon : \t",exon.identifier, "\t start or end is zero")
    return True

def validateAndWriteGTFToFile(CleanedSameStrandGenes,cleanedForSenseAndAntisenseFile):
    if os.path.exists(cleanedForSenseAndAntisenseFile):
        os.remove(cleanedForSenseAndAntisenseFile)
    for geneID in CleanedSameStrandGenes.keys():
        gene=CleanedSameStrandGenes.get(geneID)
        if validateGene(gene)==True:
            geneWithTranscriptsAndExonsLines=getGeneWithTranscriptsAndExonsAsLines(gene)
            writeGeneToTheFile(geneWithTranscriptsAndExonsLines,cleanedForSenseAndAntisenseFile)
        else:
            raise ValueError("Gene : \t",gene.identifier, "\t has some problems") 


def updateOrDeleteExonsWithNewExons(gene,exonsToBeUpdatedForGeneDf,minimumExonLength):
    for transcriptID in gene.transcripts:
        transcriptObj=gene.transcripts.get(transcriptID)
        for exonObj in transcriptObj.exons:
            rowFromCsvForExonId=exonsToBeUpdatedForGeneDf.loc[exonsToBeUpdatedForGeneDf['exonId'] == exonObj.identifier]
            if not rowFromCsvForExonId.empty:
                newExons=rowFromCsvForExonId['newExons']
                numberOfNewExons=rowFromCsvForExonId['numberOfNewExons']            
                if numberOfNewExons.all()==0: # is one value in the file but pandas is loading it as series and 0 0
                    # delete the exon from the Transcript. 
                    #print('deleteng in transcriptID\t'+transcriptID+'\texonObj.identifier\t'+exonObj.identifier)
                    transcriptObj.exons = [item for item in transcriptObj.exons if item.identifier != exonObj.identifier]
                    transcriptObj.n_exons=transcriptObj.n_exons-1
                else:
                    firstExon=True;
                    for start,end in itertools.chain.from_iterable(newExons):
                        if firstExon==True:
                            exonObj.start=start
                            exonObj.end=end
                            exonObj.length=abs(exonObj.end - exonObj.start+1)
                            #print('Updating')
                            firstExon=False;
                        else:
                            #print('Adding more exons')
                            # just updated the start,end and id in the copied exon
                            copiedExonObj = copy.deepcopy(exonObj) 
                            copiedExonObj.start=start
                            copiedExonObj.end=end
                            copiedExonObj.length=abs(copiedExonObj.end - copiedExonObj.start+1)
                            copiedExonObj.identifier=copiedExonObj.chromosome +"_" +str(copiedExonObj.start)\
                            +"_" +str(copiedExonObj.end)+"_" +copiedExonObj.strand
                            transcriptObj.add_exon(copiedExonObj) 
                            transcriptObj.n_exons=transcriptObj.n_exons+1

        gene.transcripts[transcriptID]=transcriptObj
    return gene

def updateOrDeleteTranscripts(gene,numberOfTranscriptsDeletedBecauseOfNoExon):
    for transcriptID in list(gene.transcripts):
        transcriptObj=gene.transcripts.get(transcriptID)
        if(transcriptObj.n_exons==0):
            del gene.transcripts[transcriptID] # remove transcript with no exons
            numberOfTranscriptsDeletedBecauseOfNoExon=numberOfTranscriptsDeletedBecauseOfNoExon+1
        else:
            start=sys.maxsize
            end=0
            for exon in transcriptObj.exons:
                if exon.start<start:
                    start=exon.start
                if exon.end>end:
                    end=exon.end
            transcriptObj.start=start
            transcriptObj.end=end
            gene.transcripts[transcriptID]=transcriptObj
    return gene,numberOfTranscriptsDeletedBecauseOfNoExon

def updateOrDeleteGenes(CleanedSameStrandGenes,gene,numberOfGenesDeletedBecauseOfNoTranscripts):
    start=sys.maxsize
    end=0
    if(len(gene.transcripts)==0):
        del CleanedSameStrandGenes[gene.identifier]# remove gene with no transcripts
        numberOfGenesDeletedBecauseOfNoTranscripts=numberOfGenesDeletedBecauseOfNoTranscripts+1
    else:
        for transcriptID in list(gene.transcripts):
            transcriptObj=gene.transcripts.get(transcriptID)
            if transcriptObj.start<start:
                start=transcriptObj.start
            if transcriptObj.end>end:
                end=transcriptObj.end
        gene.start=start
        gene.end=end
        CleanedSameStrandGenes[gene.identifier]=gene  
    return CleanedSameStrandGenes,numberOfGenesDeletedBecauseOfNoTranscripts

def cleaneAntiSenseOverlappingExons(uniqueGeneIdsOverlappingAntisense,CleanedSameStrandGenes,cleanedForSenseAndAntisenseFile,df):
    nrTotalGenes=len(CleanedSameStrandGenes.keys())
    nrGenesProcessed=0
    nrTranscriptsDelOfNoExon=0
    nrGenesDelOfNoTranscripts=0
    for geneID in uniqueGeneIdsOverlappingAntisense:
        gene=CleanedSameStrandGenes.get(geneID)
        exonsToBeUpdatedForGeneDf=df.loc[df['queryGeneId'] == geneID]
        try:
            nrGenesProcessed=nrGenesProcessed+1
            gene=updateOrDeleteExonsWithNewExons(gene,exonsToBeUpdatedForGeneDf,minimumExonLength)
            gene , nrTranscriptsDelOfNoExon=updateOrDeleteTranscripts(gene,nrTranscriptsDelOfNoExon)
            CleanedSameStrandGenes, nrGenesDelOfNoTranscripts=updateOrDeleteGenes(CleanedSameStrandGenes,gene,nrGenesDelOfNoTranscripts)
        except AttributeError as error:
            sys.exit('Gene causing the erros is\t'+str(geneID))
    print('number Of Total Genes before processing\t'+str(nrTotalGenes))
    print('number Of Genes processed\t'+str(nrGenesProcessed))
    print('nr Transcripts Deleted Because Of No Exon\t'+str(nrTranscriptsDelOfNoExon))  
    print('number Of Genes Deleted Because Of No Transcripts\t'+str(nrGenesDelOfNoTranscripts))
    if(nrGenesProcessed!=len(uniqueGeneIdsOverlappingAntisense)):
        sys.exit("-------Error Number of genes processed and written are not equal to initial-")
    validateAndWriteGTFToFile(CleanedSameStrandGenes,cleanedForSenseAndAntisenseFile)

def createExonBedFile(gtfFile,outputDir_filename):
    lines_seen = set() # holds lines already seen
    with open(gtfFile, "r") as a_file:
            for line in a_file:
                if not line.startswith("#"):
                    gtfLine = line.split("\t")
                    if(gtfLine[2]=="exon"):
                        longString=gtfLine[8].split()
                        geneId=longString[1].split(";")[0].replace('"', '')
                        exonLine=gtfLine[0]+"\t"+gtfLine[3]+"\t"+gtfLine[4]+"\t"+ geneId+"\t"+'0'+"\t"+gtfLine[6]+"\n"
                        if exonLine not in lines_seen: # not a duplicate
                            lines_seen.add(exonLine)
    my_file = open(outputDir_filename, "w")
    my_file.writelines(lines_seen)
    my_file.close()


def listToStringWithoutBrackets(list1):
    s=str(list1).replace('[','').replace(']','')
    return str(s).replace(' ','')
def getUpdatedLine(line):
    exonStart=line[8].split(",")
    exonEnd=line[9].split(",")
    strand=line[2]
    start=int(line[3])
    end=int(line[4])
    nummberOfExons=int(line[7])
    while("" in exonStart) :
        exonStart.remove("")
        exonEnd.remove("")
    exonStart = [int(i) for i in exonStart]
    exonEnd = [int(i) for i in exonEnd]

    block_size = []
    zip_object = zip(exonStart, exonEnd)
    for exonStart_i, exonEnd_i in zip_object:
        if(strand=='+'):
            block_size.append(exonEnd_i-exonStart_i)
        else:
            block_size.append(exonEnd_i-exonStart_i)
    block_size=listToStringWithoutBrackets(block_size)
    block_start=exonStart
    for i in range(len(exonStart)):
        if(strand=='+'):
            block_start[i] = exonStart[i]-start
        else:
            block_start[i] =exonStart[i]-start
    block_start=listToStringWithoutBrackets(block_start)
    lineToReplace=str(line[1]+" "+str(start)+" "+str(end)+" "+line[0]+" "+str(0)+" "+strand+" "+str(start)+" "+str(end)+" "+'0'+" "+str(nummberOfExons)+" "+block_size+" "+block_start+"\n")
    return lineToReplace

def createBed12FormateFromgenePredExt(inputFilePath,outputFilePath):
    with open(inputFilePath, 'r') as file:
        # read a list of lines into data
        data = file.readlines()

    for x in range(len(data)):
        line=getUpdatedLine(data[x].split())
        data[x]=line

    # and write everything back
    with open(outputFilePath, 'w') as file:
        file.writelines( data )
