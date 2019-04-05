import time
import psutil
import os
#Created by Ted

class Needleman(object):
	def __init__(self,gapScore='',seq1='',seq2=''):
		"""This constructor accepts the gap score and two sequences to align, and initialises a new Needleman instance."""
		self._TRACEBACK_TOP = 1
		self._TRACEBACK_LEFT = 2
		self._TRACEBACK_DIAG = 3

		# Declaration of class variables
		self._seq1 = seq1.upper()
		self._seq2 = seq2.upper()
		self._gapScore = gapScore
		self._scoringMatrix = []
		self._tracebackMatrix = []
		self._tracebackSequences = [[],[],[]]
		self._alignmentResults = {"Length":0,"Match":0, "Mismatch":0,"Gap":0,"PathScore":0,"MatchPercent":0,"MismatchPercent":0,"GapPercent":0}
		self._runTime = 0
		self._memSpace = 0

		# Declaration of blosum64 matrix for calculating scores of match/mismatch
		self._blosum64ScoringMatrix = {'A': {'A': '4', 'C': '0', 'E': '-1', 'D': '-2', 'G': '0', 'F': '-2', 'I': '-1', 'H': '-2', 'K': '-1', 'M': '-1', 'L': '-1', 'N': '-2', 'Q': '-1', 'P': '-1', 'S': '1', 'R': '-1', 'T': '0', 'W': '-3', 'V': '0', 'Y': '-2'},
		'C': {'A': '0', 'C': '9', 'E': '-4', 'D': '-3', 'G': '-3', 'F': '-2', 'I': '-1', 'H': '-3', 'K': '-3', 'M': '-1', 'L': '-1', 'N': '-3', 'Q': '-3', 'P': '-3', 'S': '-1', 'R': '-3', 'T': '-1', 'W': '-2', 'V': '-1', 'Y': '-2'},
		'E': {'A': '-1', 'C': '-4', 'E': '5', 'D': '2', 'G': '-2', 'F': '-3', 'I': '-3', 'H': '0', 'K': '1', 'M': '-2', 'L': '-3', 'N': '0', 'Q': '2', 'P': '-1', 'S': '0', 'R': '0', 'T': '-1', 'W': '-3', 'V': '-2', 'Y': '-2'},
		'D': {'A': '-2', 'C': '-3', 'E': '2', 'D': '6', 'G': '-1', 'F': '-3', 'I': '-3', 'H': '-1', 'K': '-1', 'M': '-3', 'L': '-4', 'N': '1', 'Q': '0', 'P': '-1', 'S': '0', 'R': '-2', 'T': '-1', 'W': '-4', 'V': '-3', 'Y': '-3'},
		'G': {'A': '0', 'C': '-3', 'E': '-2', 'D': '-1', 'G': '6', 'F': '-3', 'I': '-4', 'H': '-2', 'K': '-2', 'M': '-3', 'L': '-4', 'N': '0', 'Q': '-2', 'P': '-2', 'S': '0', 'R': '-2', 'T': '-2', 'W': '-2', 'V': '-3', 'Y': '-3'},
		'F': {'A': '-2', 'C': '-2', 'E': '-3', 'D': '-3', 'G': '-3', 'F': '6', 'I': '0', 'H': '-1', 'K': '-3', 'M': '0', 'L': '0', 'N': '-3', 'Q': '-3', 'P': '-4', 'S': '-2', 'R': '-3', 'T': '-2', 'W': '1', 'V': '-1', 'Y': '3'},
		'I': {'A': '-1', 'C': '-1', 'E': '-3', 'D': '-3', 'G': '-4', 'F': '0', 'I': '4', 'H': '-3', 'K': '-3', 'M': '1', 'L': '2', 'N': '-3', 'Q': '-3', 'P': '-3', 'S': '-2', 'R': '-3', 'T': '-1', 'W': '-3', 'V': '3', 'Y': '-1'},
		'H': {'A': '-2', 'C': '-3', 'E': '0', 'D': '-1', 'G': '-2', 'F': '-1', 'I': '-3', 'H': '8', 'K': '-1', 'M': '-2', 'L': '-3', 'N': '1', 'Q': '0', 'P': '-2', 'S': '-1', 'R': '0', 'T': '-2', 'W': '-2', 'V': '-3', 'Y': '2'},
		'K': {'A': '-1', 'C': '-3', 'E': '1', 'D': '-1', 'G': '-2', 'F': '-3', 'I': '-3', 'H': '-1', 'K': '5', 'M': '-1', 'L': '-2', 'N': '0', 'Q': '1', 'P': '-1', 'S': '0', 'R': '2', 'T': '-1', 'W': '-3', 'V': '-2', 'Y': '-2'},
		'M': {'A': '-1', 'C': '-1', 'E': '-2', 'D': '-3', 'G': '-3', 'F': '0', 'I': '1', 'H': '-2', 'K': '-1', 'M': '5', 'L': '2', 'N': '-2', 'Q': '0', 'P': '-2', 'S': '-1', 'R': '-1', 'T': '-1', 'W': '-1', 'V': '1', 'Y': '-1'},
		'L': {'A': '-1', 'C': '-1', 'E': '-3', 'D': '-4', 'G': '-4', 'F': '0', 'I': '2', 'H': '-3', 'K': '-2', 'M': '2', 'L': '4', 'N': '-3', 'Q': '-2', 'P': '-3', 'S': '-2', 'R': '-2', 'T': '-1', 'W': '-2', 'V': '1', 'Y': '-1'},
		'N': {'A': '-2', 'C': '-3', 'E': '0', 'D': '1', 'G': '0', 'F': '-3', 'I': '-3', 'H': '1', 'K': '0', 'M': '-2', 'L': '-3', 'N': '6', 'Q': '0', 'P': '-2', 'S': '1', 'R': '0', 'T': '0', 'W': '-4', 'V': '-3', 'Y': '-2'},
		'Q': {'A': '-1', 'C': '-3', 'E': '2', 'D': '0', 'G': '-2', 'F': '-3', 'I': '-3', 'H': '0', 'K': '1', 'M': '0', 'L': '-2', 'N': '0', 'Q': '5', 'P': '-1', 'S': '0', 'R': '1', 'T': '-1', 'W': '-2', 'V': '-2', 'Y': '-1'},
		'P': {'A': '-1', 'C': '-3', 'E': '-1', 'D': '-1', 'G': '-2', 'F': '-4', 'I': '-3', 'H': '-2', 'K': '-1', 'M': '-2', 'L': '-3', 'N': '-2', 'Q': '-1', 'P': '7', 'S': '-1', 'R': '-2', 'T': '-1', 'W': '-4', 'V': '-2', 'Y': '-3'},
		'S': {'A': '1', 'C': '-1', 'E': '0', 'D': '0', 'G': '0', 'F': '-2', 'I': '-2', 'H': '-1', 'K': '0', 'M': '-1', 'L': '-2', 'N': '1', 'Q': '0', 'P': '-1', 'S': '4', 'R': '-1', 'T': '1', 'W': '-3', 'V': '-2', 'Y': '-2'},
		'R': {'A': '-1', 'C': '-3', 'E': '0', 'D': '-2', 'G': '-2', 'F': '-3', 'I': '-3', 'H': '0', 'K': '2', 'M': '-1', 'L': '-2', 'N': '0', 'Q': '1', 'P': '-2', 'S': '-1', 'R': '5', 'T': '-1', 'W': '-3', 'V': '-3', 'Y': '-2'},
		'T': {'A': '0', 'C': '-1', 'E': '-1', 'D': '-1', 'G': '-2', 'F': '-2', 'I': '-1', 'H': '-2', 'K': '-1', 'M': '-1', 'L': '-1', 'N': '0', 'Q': '-1', 'P': '-1', 'S': '1', 'R': '-1', 'T': '5', 'W': '-2', 'V': '0', 'Y': '-2'},
		'W': {'A': '-3', 'C': '-2', 'E': '-3', 'D': '-4', 'G': '-2', 'F': '1', 'I': '-3', 'H': '-2', 'K': '-3', 'M': '-1', 'L': '-2', 'N': '-4', 'Q': '-2', 'P': '-4', 'S': '-3', 'R': '-3', 'T': '-2', 'W': '11', 'V': '-3', 'Y': '2'},
		'V': {'A': '0', 'C': '-1', 'E': '-2', 'D': '-3', 'G': '-3', 'F': '-1', 'I': '3', 'H': '-3', 'K': '-2', 'M': '1', 'L': '1', 'N': '-3', 'Q': '-2', 'P': '-2', 'S': '-2', 'R': '-3', 'T': '0', 'W': '-3', 'V': '4', 'Y': '-1'},
		'Y': {'A': '-2', 'C': '-2', 'E': '-2', 'D': '-3', 'G': '-3', 'F': '3', 'I': '-1', 'H': '2', 'K': '-2', 'M': '-1', 'L': '-1', 'N': '-2', 'Q': '-1', 'P': '-3', 'S': '-2', 'R': '-2', 'T': '-2', 'W': '2', 'V': '-1', 'Y': '7'},
		'X': {'A': '0', 'C': '-2', 'E': '-1', 'D': '-1', 'G': '-1', 'F': '-1', 'I': '-1', 'H': '-1', 'K': '-1', 'M': '-1', 'L': '-1', 'N': '-1', 'Q': '-1', 'P': '-2', 'S': '0', 'R': '-1', 'T': '0', 'W': '-2', 'V': '-1', 'Y': '-1'},
		'Z': {'A': '-1', 'C': '-3', 'E': '4', 'D': '1', 'G': '-2', 'F': '-3', 'I': '-3', 'H': '0', 'K': '1', 'M': '-1', 'L': '-3', 'N': '0', 'Q': '3', 'P': '-1', 'S': '0', 'R': '0', 'T': '-1', 'W': '-3', 'V': '-2', 'Y': '-2'}}
	

	# Using @property to do method 'overloading' for getter/setters
	# Getter method for seq1
	@property
	def seq1(self):
		return self._seq1
	# Setter method for seq1
	@seq1.setter	
	def seq1(self,value):
		self._seq1 = value

	# Getter method for seq2
	@property
	def seq2(self):
		return self._seq2
	# Setter method for seq2
	@seq2.setter	
	def seq2(self,value):
		self._seq2 = value

	# Getter method for gapScore
	@property
	def gapScore(self):
		return self._gapScore
	# Setter method for gapScore
	@gapScore.setter	
	def gapScore(self,value):
		self._gapScore = value

	# Getter method for matrix with scores
	@property
	def scoringMatrix(self):
		return self._scoringMatrix

	# Getter method for matrix with traceback instructions
	@property
	def tracebackMatrix(self):
		return self._tracebackMatrix
	
	# Getter method for aligned sequences
	@property
	def tracebackSequences(self):
		return self._tracebackSequences
	
	@property
	def runTime(self):
		return self._runTime

	@property
	def memSpace(self):
		return self._memSpace

	@property
	def alignmentResults(self):
		return self._alignmentResults
		

	# Method to empty the matrices 
	def resetData(self):
		"""This method reinitialises/resets the matrices and arrays used for the alignment."""
		self._scoringMatrix = []
		self._tracebackMatrix = []
		self._tracebackSequences = [[],[],[]]
		self._alignmentResults = {key: 0 for key in self._alignmentResults}
		self._runTime = 0
		self._memSpace = 0

	# Initialises alignment of sequences
	def alignSequences(self):
		"""This method accepts no parameters, and builds the scoring and trackback matrices for sequence alignment"""
		self._seq1 = self._seq1.upper()
		self._seq2 = self._seq2.upper()
		start = time.time();
		try:
			self.resetData()
			self.buildMatrix(self._scoringMatrix,self._seq1,self._seq2)
			self.buildMatrix(self._tracebackMatrix,self._seq1,self._seq2)
			self.fillMatrix(self._scoringMatrix, self._tracebackMatrix)
			self.buildTraceBack(self._scoringMatrix, self._tracebackMatrix)
			end = time.time()
			self._runTime = end - start
			self.getMemSpace()
			return 1
		except ValueError as error:
			return 0

	# Measure memory space used by entire program
	def getMemSpace(self):
		process = psutil.Process(os.getpid())
		spaceUsed = process.memory_info()[0]
		if spaceUsed / 1024 >= 1:
			spaceUsed = spaceUsed / 1024
			if spaceUsed / 1024>= 1:
				spaceUsed = spaceUsed / 1024
				self._memSpace = str(spaceUsed)+"mb"
			else:
				self._memSpace = str(spaceUsed)+"kb"
		else:
			self._memSpace = str(spaceUsed)+"b"

	# General method to build matrix of breadth a and length b with headers
	def buildMatrix(self,matrix, a,b):
		start = time.time()
		"""This method accepts an empty array and builds a 2D matrix (array of arrays) with header information."""
		for i in range(len(a)+2):
			matrix.append([])
			for j in range(len(b)+2):
				matrix[i].append(0)
				if i is 0 and j > 1:
					matrix[i][j] = b[j-2]
			if i > 1:
				matrix[i][0] = a[i-2]
		matrix[0][0] = ''
		matrix[0][1] = ''	
		matrix[1][0] = ''
	
	# Match character and get blosum62 score
	def matchChar(self,a,b):
		"""This method accepts two characters and matches them against the blosum62 scoring matrix and returns the score."""
		if a not in self._blosum64ScoringMatrix or b not in self._blosum64ScoringMatrix[a]:
			raise ValueError("Invalid sequence character found.")
		return int(self._blosum64ScoringMatrix[a][b])

	# Counts the score to be filled into the respective cell of the matrix
	def countScore(self,matrix,i,j):
		"""This method accepts the scoring matrix together with two index posiitons for the score to be calculated."""
		if i is 1 and j > 1:
			matrix[i][j] = matrix[i][j-1] + self._gapScore
			return self._TRACEBACK_LEFT
		elif i > 1 and j is 1:
			matrix[i][j] = matrix[i-1][j] + self._gapScore
			return self._TRACEBACK_TOP
		else:
			try:
				scoreDiag = matrix[i-1][j-1] + self.matchChar(matrix[0][j],matrix[i][0])
				scoreTop = matrix[i-1][j] + self._gapScore
				scoreLeft = matrix[i][j-1] + self._gapScore
				matrix[i][j] = max(scoreDiag,scoreTop,scoreLeft)
				if matrix[i][j] is scoreDiag:
					return self._TRACEBACK_DIAG
				elif matrix[i][j] is scoreTop:
					return self._TRACEBACK_TOP
				else:
					return self._TRACEBACK_LEFT
			except ValueError as error:
				raise error
			

	# Traverses through the matrices and builds the data within
	def fillMatrix(self,scoring, traceback):
		"""This method traverses through the scoring matrix and fills in the scoring, together with the traceback instructions in the traceback matrix."""
		for i in range(1,len(scoring)):
			for j in range(1,len(scoring[0])):
				if i is 1 and j is 1:
					scoring[i][j] = 0
				else:
					try:
						traceback[i][j] = self.countScore(scoring,i,j)
					except ValueError as error:
						raise error

					

	# Finds the optimal path for sequence alignment and builds the aligned sequence into tracebackSequences
	def buildTraceBack(self,scoring,traceback):
		"""this method accepts the scoring and trackback matrix, and traces the optimal sequence alignment into tracebackSequences."""
		seq1Len = len(traceback)-1
		seq2Len = len(traceback[0])-1
		startVer,startHor = seq1Len,seq2Len
		endVer,endHor = 1,1
		start = traceback[seq1Len][seq2Len]
		matches = 0
		mismatches = 0
		gaps = 0
		self._alignmentResults['PathScore'] = scoring[seq1Len][seq2Len]
		while True:
			if (startVer is endVer and startHor is endHor):
				break
			if traceback[startVer][startHor] is self._TRACEBACK_DIAG:
				self._tracebackSequences[0].insert(0,traceback[startVer][0]) # Sequence 1
				self._tracebackSequences[1].insert(0,"|") # Match Result
				self._tracebackSequences[2].insert(0,traceback[0][startHor]) # Sequence 2
				if (traceback[startVer][0] == traceback[0][startHor]):
					matches += 1
				else:
					mismatches += 1
				startVer -= 1
				startHor -= 1
			elif traceback[startVer][startHor] is self._TRACEBACK_TOP:
				self._tracebackSequences[0].insert(0,traceback[startVer][0]) # Sequence 1
				self._tracebackSequences[1].insert(0," ") # Match Result
				self._tracebackSequences[2].insert(0,"-") # Sequence 2
				gaps += 1
				startVer -= 1
			elif traceback[startVer][startHor] is self._TRACEBACK_LEFT:
				self._tracebackSequences[0].insert(0,"-") # Sequence 1
				self._tracebackSequences[1].insert(0," ") # Match Result
				self._tracebackSequences[2].insert(0,traceback[0][startHor]) # Sequence 2
				gaps += 1
				startHor -= 1
		self._alignmentResults['Length'] = len(self._tracebackSequences[0])
		self._alignmentResults['Match'] = matches;
		self._alignmentResults['MatchPercent'] = (self._alignmentResults['Match'] * 1.00 / self._alignmentResults['Length'] * 100)
		self._alignmentResults['Mismatch'] = mismatches
		self._alignmentResults['MismatchPercent'] = (self._alignmentResults['Mismatch'] * 1.00 / self._alignmentResults['Length'] * 100)
		self._alignmentResults['Gap'] = gaps
		self._alignmentResults['GapPercent'] = (self._alignmentResults['Gap'] * 1.00 / self._alignmentResults['Length'] * 100)
		self._alignmentResults['Identity'] = self._alignmentResults['Match']
		self._alignmentResults['IdentityPercent'] = (self._alignmentResults['Identity'] * 1.00 / self._alignmentResults['Length'] * 100)
		self._alignmentResults['Similarity'] = self._alignmentResults['Match'] + self._alignmentResults['Mismatch']
		self._alignmentResults['SimilarityPercent'] = (self._alignmentResults['Similarity'] * 1.00 / self._alignmentResults['Length'] * 100)

	# Prints the scoring and traceback matrices to console
	def printMatrices(self):
		"""This method prints the scoring matrix and traceback matrix to console."""
		print "\n\nSCORING MATRIX\n"
		print('\n'.join([''.join(['{:>3}'.format(item) for item in row]) for row in self._scoringMatrix]))
		print "\n\nTRACEBACK MATRIX\n"
		print('\n'.join([''.join(['{:>6}'.format(item) for item in row]) for row in self._tracebackMatrix]))

	# Prints the alignment results to console
	def printResults(self):
		"""This method prints the alignment results to console."""
		print "\n\n\n=========================================================="
		print "Alignment results using BLOSUM62 with a gap score of %d" % (self._gapScore)
		print "=========================================================="

		print "\nOptimal Alignment\n"
		print('\n'.join([''.join(['{:>1}'.format(item) for item in row]) for row in self._tracebackSequences]))	

		print "\nSequence 1 Original Length: %d" % (len(self._seq1))
		print "Sequence 2 Original Length: %d" % (len(self._seq2))
		print "Optimal Aligned Sequence Length: %d" % (len(self.tracebackSequences[0]))


		print "\nAlignment Path Score: %d" % (self._alignmentResults['PathScore'])

		print "\nIdentity: %d/%d (%.2f%%)" % (self._alignmentResults['Identity'], self._alignmentResults['Length'], self._alignmentResults['IdentityPercent'])
		print "Similarity: %d/%d (%.2f%%)" % (self._alignmentResults['Similarity'], self._alignmentResults['Length'], self._alignmentResults['SimilarityPercent'])

		print "\nMatch count: %d (%.2f%%)" % (self._alignmentResults['Match'],self._alignmentResults['MatchPercent'])
		print "Mismatch count: %d (%.2f%%)" % (self._alignmentResults['Mismatch'],self._alignmentResults['MismatchPercent'])
		print "Gap count: %d (%.2f%%)" % (self._alignmentResults['Gap'],self._alignmentResults['GapPercent'])

	def printRunTime(self):
		print "\nAlignment runtime:",self._runTime,"seconds"
