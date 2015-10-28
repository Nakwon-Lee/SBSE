import sys
import math
import random
import copy
from operator import attrgetter

evals = 0
budget = 0

#calculate appropriate box index using global index
#param
#pos: global index
#return
#idx: appropriate box index of global position
def calcIdxOfBox(pos):
	
	quota = pos // 9
	remainder = pos % 9
	idx = 0

	if quota < 3:
		if remainder < 3:
			idx = 0
		elif (remainder >= 3) and (remainder < 6):
			idx = 1
		else:
			idx = 2
	elif (quota >= 3) and (quota < 6):
		if remainder < 3:
			idx = 3
		elif (remainder >= 3) and (remainder < 6):
			idx = 4
		else:
			idx = 5
	else:
		if remainder < 3:
			idx = 6
		elif (remainder >= 3) and (remainder < 6):
			idx = 7
		else:
			idx = 8

	return idx

#class for storing value assignments
#this class have set of assigned values for each column, row, and box
#it is used to check whether one value is already assigned in same column, row, or box 
class Sudoku:
	#because the problem is sudoku, num should be the integer nine
	def __init__(self,num):
		self.col = []
		self.row = []
		self.box = []
	
		for i in range(num):
			self.col.append(set())
			self.row.append(set())
			self.box.append(set())

#read file, which represent sudoku problem, and then save the initial value assignment,
#and global index of vacant positions
#param
#filename: the text file which have sudoku problem
#return
#num: list of storing global indices of vacant positions; num_fix: list of storing initial value assignment and their global indices; sudoku: instance of "Sudoku" that is storing initial value assignment
def read_data(filename):
	global sudoku
	lines = open(filename).readlines()
	#because the problem is sudoku, "Sudoku" should receive integer nine
	sudoku = Sudoku(9)
	num_fix = []
	num = []

	#global index indicator
	pos = 0
	#sudoku problem file format - initial value assignment is represented as integer value
	#when the value is bigger than 0 and smaller than 10. period means vacant position.
	#global index is just the order of each valid digit or period. other symbols are ignored.
	for line in lines:
		for i in range(len(line)):
			if line[i].isdigit() and int(line[i]) >= 1 and int(line[i]) <= 9:
				#if the integer value is valid, it stored to "Sudoku" instance to remember
				#that the value is already assigned in such column, row, and box. 
				sudoku.row[pos//9].add(int(line[i]))
				sudoku.col[pos%9].add(int(line[i]))
				sudoku.box[calcIdxOfBox(pos)].add(int(line[i]))
				#store the initial assinged value and its global index
				num_fix.append((pos, int(line[i])))
				pos += 1
			elif line[i] == '.':
				#if we find the period, it means that the current global index is vacant position.
				num.append(pos)
				pos += 1
			else:
				pass

	return num, num_fix, sudoku

#CrossoverZeroFirst - for changable order
#global indices in first parent that have zero (impossible assignment) are assigned with the values at same index in the second parent. And, those indices are assigned first. remaing global indices are filled with first parent.
#not used
class CrossoverZeroFirst:
	def crossover(self, parent_a, parent_b):
		assert(len(parent_a.permutation) == len(parent_b.permutation))
		size = len(parent_a.permutation)

		child_a_perm = []
		child_b_perm = []
		child_a_set = set()
		child_b_set = set()

		count_a = 0
		for i in range(size):
			if parent_a.permutation[i][1] == 0:
				for j in range(size):
					if parent_b.permutation[j][0] == parent_a.permutation[i][0]:
						child_a_perm.append((parent_b.permutation[j][0], parent_b.permutation[j][1]))
						child_a_set.add(parent_b.permutation[j][0])
						count_a += 1
						break

		count_b = 0
		for i in range(size):
			if parent_b.permutation[i][1] == 0:
				for j in range(size):
					if parent_a.permutation[j][0] == parent_b.permutation[i][0]:
						child_b_perm.append((parent_a.permutation[j][0], parent_a.permutation[j][1]))
						child_b_set.add(parent_a.permutation[j][0])
						count_b += 1
						break

		ch_a_idx = count_a
		p_a_idx = 0
		ch_b_idx = count_b
		p_b_idx = 0

		while ch_a_idx < size:
			if parent_a.permutation[p_a_idx][0] not in child_a_set:
				child_a_perm.append(parent_a.permutation[p_a_idx])
				ch_a_idx += 1
			p_a_idx += 1

		while ch_b_idx < size:
			if parent_b.permutation[p_b_idx][0] not in child_b_set:
				child_b_perm.append(parent_b.permutation[p_b_idx])
				ch_b_idx += 1
			p_b_idx += 1

		return child_a_perm, child_b_perm

#CrossoverOrder - for both changable and non changable order
#swap global indices and their valuations btw two changepoints and fill others with non swap parts with retained order
class CrossoverOrder:
    def crossover(self, parent_a, parent_b):
		assert(len(parent_a.permutation) == len(parent_b.permutation))
		size = len(parent_a.permutation)

		cp1 = random.randrange(size)
		cp2 = random.randrange(size)
		while(cp2 == cp1):
			cp2 = random.randrange(size)

		if cp2 < cp1:
			cp1, cp2 = cp2, cp1

		child_a_perm = []
		child_b_perm = []
		child_a_set = set()
		child_b_set = set()

		#swap global indices and their valuations
		for i in range(cp1, cp2 + 1):
			child_a_perm.append(parent_b.permutation[i])
			child_a_set.add(parent_b.permutation[i][0])
			child_b_perm.append(parent_a.permutation[i])
			child_b_set.add(parent_a.permutation[i][0])

		ch_a_idx = 0
		p_a_idx = 0
		ch_b_idx = 0
		p_b_idx = 0

		#fill other global indices with non swap parts with retained order
		while ch_a_idx < cp1:
			if parent_a.permutation[p_a_idx][0] not in child_a_set:
				child_a_perm.insert(ch_a_idx, parent_a.permutation[p_a_idx])
				ch_a_idx += 1
			p_a_idx += 1
		ch_a_idx = cp2 + 1

		while ch_b_idx < cp1:
			if parent_b.permutation[p_b_idx][0] not in child_b_set:
				child_b_perm.insert(ch_b_idx, parent_b.permutation[p_b_idx])
				ch_b_idx += 1
			p_b_idx += 1
		ch_b_idx = cp2 + 1

		while ch_a_idx < size:
			if parent_a.permutation[p_a_idx][0] not in child_a_set:
				child_a_perm.append(parent_a.permutation[p_a_idx])
				ch_a_idx += 1
			p_a_idx += 1

		while ch_b_idx < size:
			if parent_b.permutation[p_b_idx][0] not in child_b_set:
				child_b_perm.append(parent_b.permutation[p_b_idx])
				ch_b_idx += 1
			p_b_idx += 1

		return child_a_perm, child_b_perm
       
#Mutation - for changable order
#swap randomly selected two global indices with retaining their assignment
class Mutation:
    def mutate(self, solution):
        size = len(solution)

        mp1 = random.randrange(size)
        mp2 = random.randrange(size)
        while(mp2 == mp1):
            mp2 = random.randrange(size)
        solution[mp1], solution[mp2] = solution[mp2], solution[mp1]
        return solution

#MutationSwapValueFixOrder - for changeable and non changeable order
#swap valuations of two randomly selected global indices 
class MutationSwapValueFixOrder:
	def mutate(self, solution):
		size = len(solution)

		mp1 = random.randrange(size)
		mp2 = random.randrange(size)
		while(mp2 == mp1):
		    mp2 = random.randrange(size)

		sol1 = solution[mp1][1]
		sol2 = solution[mp2][1]

		solution[mp1] = (solution[mp1][0], sol2)
		solution[mp2] = (solution[mp2][0], sol1)
	
		return solution

#for changeable order
#swap two randomly selected global indices with zero valuation (forcing new assignment) 
class MutationSwapZero:
	def mutate(self, solution):
		size = len(solution)

		mp1 = random.randrange(size)
		mp2 = random.randrange(size)
		while(mp2 == mp1):
		    mp2 = random.randrange(size)
		solution[mp1], solution[mp2] = solution[mp2], solution[mp1]
		solution[mp1] = (solution[mp1][0],0)
		solution[mp2] = (solution[mp2][0],0)
		return solution

#for changeable order
#move zero assigned (impossible valuation) global indices to the first position and randomly valuate them (maximum two zero assigned indices)
class MutationSwapZeroToFirst:
	def mutate(self, solution):
		size = len(solution)

		val_list = [1,2,3,4,5,6,7,8,9]

		key = False

		for i in range(size):
			mp1 = random.randrange(size)
			if solution[mp1][1] == 0:
				solution[mp1], solution[0] = solution[0], solution[mp1]
				solution[0] = (solution[0][0], random.choice(val_list))
				key = True
				break

		if key:
			for i in range(size):
				mp1 = random.randrange(size)
				if solution[mp1][1] == 0 and mp1 != 0:
					solution[mp1], solution[1] = solution[1], solution[mp1]
					solution[1] = (solution[1][0], random.choice(val_list))
					break
		else:
			for i in range(size):
				mp1 = random.randrange(size)
				if solution[mp1][1] == 0:
					solution[mp1], solution[0] = solution[0], solution[mp1]
					solution[0] = (solution[0][0], random.choice(val_list))
					key = True
					break

		return solution

#for both changeable and non changeable order
#assign zero to randomly selected two global indices
class MutationZero:
    def mutate(self, solution):
		size = len(solution)

		mp1 = random.randrange(size)
		mp2 = random.randrange(size)
		while(mp2 == mp1):
			mp2 = random.randrange(size)

		solution[mp1] = (solution[mp1][0],0)
		solution[mp2] = (solution[mp2][0],0)
		return solution

#for both changeable and non changeable order
#assign zero to randomly selected two global indices (one must be in first three, which is combinatorial part)
class MutationZeroCombMustMut:
    def mutate(self, solution):
		size = len(solution)
		size2 = 3
		

		mp1 = random.randrange(size)			
		mp2 = random.randrange(size2)
		
		while(mp2 == mp1):
			mp2 = random.randrange(size2)

		solution[mp1] = (solution[mp1][0],0)
		solution[mp2] = (solution[mp2][0],0)

		return solution

#randomly select two solutions in population and, return more fitted one
#cannot use current setting
#for one level fitness
class BinaryTournament:
    def select(self, population):
		i = random.randrange(len(population))
		j = random.randrange(len(population))
		while i == j:
		   j = random.randrange(len(population))
		
		a = population[i]
		b = population[j]

		if a.fitnessf < b.fitnessf:
		    return a
		else:
			return b

#for three level fitness
#randomly select two solutions in population and, return more fitted one
class BinaryTournamentComb:
    def select(self, population):
		i = random.randrange(len(population))
		j = random.randrange(len(population))
		while i == j:
		   j = random.randrange(len(population))
		
		a = population[i]
		b = population[j]

		if a.fitness < b.fitness:
		    return a
		elif a.fitness > b.fitness:
			return b
		else:
			if a.fitness3 < b.fitness3:
				return a
			elif a.fitness3 > b.fitness3:
				return b
			else:
				if a.fitness2 < b.fitness2:
					return a
				else:
					return b

#class for solution of sudoku
class Solution:
	def __init__(self):
		#save the vacant global indices and their candidate valuations
		self.permutation = []
		self.fitness = 81 #number of zero assigned global indices in permutation
		self.fitness2 = 99999999.0 #sum of frequency of valuations
		self.fitness3 = 99999999.0 #sum of frequency of combinatorial valuations
		self.age = 0 #aging variable

	#randomly generate the best effort solution
	#param
	#sudoku: initial assignment (given assignment from file); num: list of vacant global indices; freq: non changable freqency of valuations; freq_cl: changeable freqency table; combfreq: non changeable combinatorial freqency table; combfreq_cl: changeable combinatorial freqency table
	def genRanSolComb(self, sudoku, num, freq, freq_cl, combfreq, combfreq_cl):
		global evals
		#copy data to avoid data inconsistancy
		sudoku_cl = copy.deepcopy(sudoku)
		num_cl = copy.deepcopy(num)

		#if activate it, order can be changeable
		#random.shuffle(num_cl)

		self.fitness = 0
		self.fitness2 = 0
		self.fitness3 = 0

		#indices for calculating combinatorial frequency
		#first three vacant indices are combined as combinatorial frequency
		#we assume that first three vacants are always assigned by some value 
		ft = [0 for i in range(3)]

		#first three vacant indices are considered earlier to calculate the combinatorial
		#frequency.
		for i in range(3):
			possible_numbers = set([1,2,3,4,5,6,7,8,9])

			#for every vacant indices, which we want to fill, find the possible set of
			#candidate numbers. it can calculated using set operation. from the whole
			#candidate number set, differentiate already assigned numbers using "Sudoku"
			#class. Therefore, the possible number set is all numbers - already assigned numbers in same row - already assigned numbers in same column - already assigned numbers in same box
			possible_numbers = possible_numbers - sudoku_cl.row[num_cl[i]//9] - sudoku_cl.col[num_cl[i]%9] - sudoku_cl.box[calcIdxOfBox(num_cl[i])]

			#if there is no possible number, assign zero and increase fitness
			#only for non changeable order
			#assume there must be a possible number in first three indices
			if len(possible_numbers) == 0:
				assert i >= 3
				self.permutation.append((num_cl[i], 0))
				self.fitness += 1
			#if there are possible numbers, assign one of possible number.
			else:
				a = random.choice(list(possible_numbers))
				if i < 3:
					ft[i] = a
				self.permutation.append((num_cl[i], a))
				#add assigned number to sets to remember what numbers are assinged
				sudoku_cl.row[num_cl[i]//9].add(a)
				sudoku_cl.col[num_cl[i]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(num_cl[i])].add(a)
				#calculate the fitness3 from non modifiable comb frequency table
				if i == 2:
					self.fitness3 = combfreq[ft[0]-1][ft[1]-1][ft[2]-1]

		#fill vacant positions after index 3
		for i in range(3, len(num_cl)):
			possible_numbers = set([1,2,3,4,5,6,7,8,9])

			#find the possible candidate number using set operation
			possible_numbers = possible_numbers - sudoku_cl.row[num_cl[i]//9] - sudoku_cl.col[num_cl[i]%9] - sudoku_cl.box[calcIdxOfBox(num_cl[i])]

			#if there is no possible number, assign zero and increase fitness
			if len(possible_numbers) == 0:
				self.permutation.append((num_cl[i], 0))
				self.fitness += 1
			else:
				#if there are possible numbers, assign one of them
				a = random.choice(list(possible_numbers))
				self.permutation.append((num_cl[i], a))
				#add assigned number to sets to remember what numbers are assigned
				sudoku_cl.row[num_cl[i]//9].add(a)
				sudoku_cl.col[num_cl[i]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(num_cl[i])].add(a)
				assert freq[num_cl[i]][a-1] >= 0
				#calculated fitness2: sum of frequencies of all assigned numbers except first
				#three
				self.fitness2 = self.fitness2 + freq[num_cl[i]][a-1]

		#we need to find the 0 fitness solution so, a solution has lower fitness is better than higher fitness value. However, lower fitness is not an less optimized solution. It also faild solution. Therefore, we need to avoid the local optimum (around the lower fitness but not a exact solution). So, we try to escape the local optimum if the assigned values are repeted many times and the exact solution cannot be found using frequency of appearing values.
		#actually solutions which have low fitness are not in consideration. they are not converged. but high fitness but not exact solutions are considered as local optimum. Therefore, we add frequency by ratio of how fitted. The assigned values in high fitness vut not exact solutions have more frequency and vice versa.

		#calculate frequencies
		#frequency means how the assigned value is appeared at low fitness.
		#lower fitness value, more frequency
		combfreq_cl[ft[0]-1][ft[1]-1][ft[2]-1] += ((81.0-float(self.fitness))/81.0)

		#calculate the comb frequencies
		for i in range(3, len(num_cl)):
			if self.permutation[i][1] != 0:		
				freq_cl[num_cl[i]][self.permutation[i][1]-1] += ((81.0-float(self.fitness))/81.0)

		#increase evals
		evals += 1

	#randomly generate the best effort solution
	#param
	#sudoku: initial assignment (given assignment from file); num: list of vacant global indices; freq: non changable freqency of valuations; freq_cl: changeable freqency table; combfreq: non changeable combinatorial freqency table; combfreq_cl: changeable combinatorial freqency table
	def genRanSol(self, sudoku, num, freq, freq_cl):
		global evals
		#copy data to avoid data inconsistancy
		sudoku_cl = copy.deepcopy(sudoku)
		num_cl = copy.deepcopy(num)

		#if activate it, order can be changeable
		random.shuffle(num_cl)

		self.fitness = 0
		self.fitness2 = 0

		#fill vacant positions
		for i in range(len(num_cl)):
			possible_numbers = set([1,2,3,4,5,6,7,8,9])

			#find the possible candidate number using set operation
			possible_numbers = possible_numbers - sudoku_cl.row[num_cl[i]//9] - sudoku_cl.col[num_cl[i]%9] - sudoku_cl.box[calcIdxOfBox(num_cl[i])]

			#if there is no possible number, assign zero and increase fitness
			if len(possible_numbers) == 0:
				self.permutation.append((num_cl[i], 0))
				self.fitness += 1
			else:
				#if there are possible numbers, assign one of them
				a = random.choice(list(possible_numbers))
				self.permutation.append((num_cl[i], a))
				#add assigned number to sets to remember what numbers are assigned
				sudoku_cl.row[num_cl[i]//9].add(a)
				sudoku_cl.col[num_cl[i]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(num_cl[i])].add(a)
				assert freq[num_cl[i]][a-1] >= 0
				#calculated fitness2: sum of frequencies of all assigned numbers
				self.fitness2 = self.fitness2 + freq[num_cl[i]][a-1]

		#we need to find the 0 fitness solution so, a solution has lower fitness is better than higher fitness value. However, lower fitness is not an less optimized solution. It also faild solution. Therefore, we need to avoid the local optimum (around the lower fitness but not a exact solution). So, we try to escape the local optimum if the assigned values are repeted many times and the exact solution cannot be found using frequency of appearing values.
		#actually solutions which have low fitness are not in consideration. they are not converged. but high fitness but not exact solutions are considered as local optimum. Therefore, we add frequency by ratio of how fitted. The assigned values in high fitness vut not exact solutions have more frequency and vice versa.

		#calculate frequencies
		#frequency means how the assigned value is appeared at low fitness.
		#lower fitness value, more frequency
		for i in range(len(num_cl)):
			if self.permutation[i][1] != 0:		
				freq_cl[num_cl[i]][self.permutation[i][1]-1] += ((81.0-float(self.fitness))/81.0)

		#increase evals
		evals += 1

	#validate new offspring with given value assignements
	#calculate the frequency and fitness of new offspring
	#param
	#sudoku: instance for storing initially assigned values; perm: global indices and valuations of vacant positions; freq: non modifiable frequency table for calculating fitness2; freq_cl: modifiable frequency table for calculating next generation's frequency table; combfreq: non modifiable comb frequency table; combfreq_cl: modifiable comb frequency table
	#return
	def validate(self, sudoku, perm, freq, freq_cl):
		global evals
		sudoku_cl = copy.deepcopy(sudoku)
		
		self.fitness = 0
		self.fitness2 = 0

		#validate the vacant positions
		for i in range(len(perm)):
			possible_numbers = set([1,2,3,4,5,6,7,8,9])
			possible_numbers = possible_numbers - sudoku_cl.row[perm[i][0]//9] - sudoku_cl.col[perm[i][0]%9] - sudoku_cl.box[calcIdxOfBox(perm[i][0])]

			#if there is no possible number, assign zero
			if len(possible_numbers) == 0:
				self.permutation.append((perm[i][0], 0))
				self.fitness += 1
			else:
				a = 0
				#in validation, we already have list of global indices and their valuations
				#therefore, if the already assigned number is in possible candidate numberset,
				#just select it, if not, select one possible number from candidate set
				if perm[i][1] in possible_numbers:
					a = perm[i][1]
					self.permutation.append((perm[i][0], a)) 
				else:
					a = random.choice(list(possible_numbers))
					self.permutation.append((perm[i][0], a))
				sudoku_cl.row[perm[i][0]//9].add(a)
				sudoku_cl.col[perm[i][0]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(perm[i][0])].add(a)
				assert freq[perm[i][0]][a-1] >= 0
				#calculate the fitness2
				self.fitness2 = self.fitness2 + freq[perm[i][0]][a-1]

		#calculate comb frequencies
		#increase frequency as ratio of fitness to represent that less appeared values
		#with same fitness are more closed to optimal solution.
		for i in range(len(perm)):
			if self.permutation[i][1] != 0:		
				freq_cl[perm[i][0]][self.permutation[i][1]-1] += ((81.0-float(self.fitness))/81.0)

		#increase evals
		evals += 1

	#validate new offspring with given value assignements
	#calculate the frequency and fitness of new offspring
	#param
	#sudoku: instance for storing initially assigned values; perm: global indices and valuations of vacant positions; freq: non modifiable frequency table for calculating fitness2; freq_cl: modifiable frequency table for calculating next generation's frequency table; combfreq: non modifiable comb frequency table; combfreq_cl: modifiable comb frequency table
	#return
	def validateComb(self, sudoku, perm, freq, freq_cl, combfreq, combfreq_cl):
		global evals
		sudoku_cl = copy.deepcopy(sudoku)
		
		self.fitness = 0
		self.fitness2 = 0
		self.fitness3 = 0

		ft = [0 for i in range(3)]

		#first three valuations are considered as combination
		for i in range(3):
			possible_numbers = set([1,2,3,4,5,6,7,8,9])
			#find possible candidate numbers using set operation
			possible_numbers = possible_numbers - sudoku_cl.row[perm[i][0]//9] - sudoku_cl.col[perm[i][0]%9] - sudoku_cl.box[calcIdxOfBox(perm[i][0])]

			#if there is no possible number, assign zero
			#but we assume that possible number must be exist in first three valuations
			if len(possible_numbers) == 0:
				assert i >= 3
				self.permutation.append((perm[i][0], 0))
				self.fitness += 1
			else:
				a = 0
				#in validation, we already have list of global indices and their valuations
				#therefore, if the already assigned number is in possible candidate numberset,
				#just select it, if not, select one possible number from candidate set
				if perm[i][1] in possible_numbers:
					a = perm[i][1]
					self.permutation.append((perm[i][0], a)) 
				else:
					a = random.choice(list(possible_numbers))
					self.permutation.append((perm[i][0], a))
				if i < 3:
					ft[i] = a
				sudoku_cl.row[perm[i][0]//9].add(a)
				sudoku_cl.col[perm[i][0]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(perm[i][0])].add(a)
				#calculate the fitness3 
				if i == 2:
					self.fitness3 = combfreq[ft[0]-1][ft[1]-1][ft[2]-1]

		#validate the remaining parts of vacant positions
		for i in range(3, len(perm)):
			possible_numbers = set([1,2,3,4,5,6,7,8,9])
			possible_numbers = possible_numbers - sudoku_cl.row[perm[i][0]//9] - sudoku_cl.col[perm[i][0]%9] - sudoku_cl.box[calcIdxOfBox(perm[i][0])]

			#if there is no possible number, assign zero
			if len(possible_numbers) == 0:
				self.permutation.append((perm[i][0], 0))
				self.fitness += 1
			else:
				a = 0
				#in validation, we already have list of global indices and their valuations
				#therefore, if the already assigned number is in possible candidate numberset,
				#just select it, if not, select one possible number from candidate set
				if perm[i][1] in possible_numbers:
					a = perm[i][1]
					self.permutation.append((perm[i][0], a)) 
				else:
					a = random.choice(list(possible_numbers))
					self.permutation.append((perm[i][0], a))
				sudoku_cl.row[perm[i][0]//9].add(a)
				sudoku_cl.col[perm[i][0]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(perm[i][0])].add(a)
				assert freq[perm[i][0]][a-1] >= 0
				#calculate the fitness2
				self.fitness2 = self.fitness2 + freq[perm[i][0]][a-1]

		#calculate comb frequencies
		#increase frequency as ratio of fitness to represent that less appeared values
		#with same fitness are more closed to optimal solution.
		combfreq_cl[ft[0]-1][ft[1]-1][ft[2]-1] += ((81.0-float(self.fitness))/81.0)
		
		for i in range(3, len(perm)):
			if self.permutation[i][1] != 0:		
				freq_cl[perm[i][0]][self.permutation[i][1]-1] += ((81.0-float(self.fitness))/81.0)

		#increase evals
		evals += 1

	#re-calculate the fitness and frequency of survived solutions
	#param
	#freq: non modifiable frequency table for calculating solution's fitness2; freq_cl: modifiable frequency table for calculating next generation's frequency; combfreq: non modifiable frequency table for calculating solution's fitness3; combfreq_cl: modifiable frequency table for calculating next generation's frequency
	#return
	def reCalcFitComb(self, freq, freq_cl, combfreq, combfreq_cl):
		global evals
		self.fitness2 = 0
		self.fitness3 = 0

		#increase evals
		evals += 1

		#survived solutions have valid valuations therefore, it is not need to find the candidate numbers
		ft = [0 for i in range(3)]

		for i in range(len(self.permutation)):
			assert freq[self.permutation[i][0]][self.permutation[i][1]-1] >= 0
			#calculate the fitness2
			self.fitness2 = self.fitness2 + freq[self.permutation[i][0]][self.permutation[i][1]-1]
			#calculate the fitness3 and calculate the next generation's frequency
			if self.permutation[i][1] != 0:
				if i < 3:
					ft[i] = self.permutation[i][1]
				if i == 2:
					self.fitness3 = combfreq[ft[0]-1][ft[1]-1][ft[2]-1]
					combfreq_cl[ft[0]-1][ft[1]-1][ft[2]-1] += ((81.0-float(self.fitness))/81.0)
				freq_cl[self.permutation[i][0]][self.permutation[i][1]-1] += ((81.0-float(self.fitness))/81.0)

	def reCalcFit(self, freq, freq_cl):
		global evals
		self.fitness2 = 0

		#increase evals
		evals += 1

		#survived solutions have valid valuations therefore, it is not need to find the candidate numbers
		for i in range(len(self.permutation)):
			assert freq[self.permutation[i][0]][self.permutation[i][1]-1] >= 0
			#calculate the fitness2
			self.fitness2 = self.fitness2 + freq[self.permutation[i][0]][self.permutation[i][1]-1]
			#calculate the next generation's frequency
			if self.permutation[i][1] != 0:
				freq_cl[self.permutation[i][0]][self.permutation[i][1]-1] += ((81.0-float(self.fitness))/81.0)

#main genetic algorithm part
#param
#filename: file that have sudoku problem (initial valuation data); pop: population size
def ga(filename, pop):

	print pop, budget

	#read file and store initial valuation and vacant positions
	num, num_fix, sudoku = read_data(filename)

	#table for combinatorial frequency table
	#combfreq = {}
	#for i in range(9):
	#	combfreq[i] = {}
	#	for j in range(9):
	#		combfreq[i][j] = [0 for col in range(9)]

	#table for frequency table
	freq = [[0 for col in range(9)] for row in range(81)]

	#list for population
	population = []

	#assign selection, crossover, and mutation operators
	#selection_op = BinaryTournament()
	selection_op = BinaryTournamentComb()
	crossover_op = CrossoverZeroFirst()
	#crossover_op = CrossoverOrder()
	mutation_op = Mutation()
	#mutation_op = MutationZero()
	#mutation_op = MutationZeroCombMustMut()
	#mutation_op = MutationSwapValueFixOrder()
	#mutation_op = MutationSwapZeroToFirst()
	#mutation_op = MutationSwapZero()

	#list for elitism
	elitism = []
	aging = []

	#copy frequency table to retain frequency table without changing
	freq_cl = copy.deepcopy(freq)
	#combfreq_cl = copy.deepcopy(combfreq)

	#if we use fixed order, only one shuffle is needed
	#random.shuffle(num)

	pop_size = pop
	#generate initial population
	for i in range(pop_size):
		new_individual = Solution()
		new_individual.genRanSol(sudoku, num, freq, freq_cl)
		population.append(new_individual)

	#change frequency table that changed during initial population generation
	freq = freq_cl
	#combfreq = combfreq_cl

	#calculate probability of survival
	#more fitted and younger solution live more
	for i in range(pop_size):
		elitism.append(0.8 * pow( (float(pop_size-i)/float(pop_size)),5 ))

	for i in range(5):
		aging.append( math.sqrt(math.sqrt(0.8 - (float(i) * 0.2) )) )

	#print elitism
	#print aging

	population = sorted(population, key=attrgetter('fitness', 'fitness3', 'fitness2'))
	current_best = population[0]

	generation = 0

	while (evals < budget) and (current_best.fitness != 0):

		#frequency can exceed the limit of integer value therefore, it is needed to decrease all table before occuring excess
		key = False
		for i in range(81):
			for j in range(9):
				if freq[i][j] >= 99999999:
					key = True
					break
			if key == True:
				break

		if key:
			for i in range(81):
				for j in range(9):
						#two strategy, divide by some value for all or just substitue some value for all
						freq[i][j] = freq[i][j] / 2
						#if freq[i][j] - 50000 >= 0:
						#	freq[i][j] -= 50000
						#else:
						#	freq[i][j] = 0

		#for comb frequency
		#key = False
		#for i in range(9):
		#	for j in range(9):
		#		for k in range(9):
		#			if combfreq[i][j][k] >= 99999:
		#				key = True
		#				break
		#		if key == True:
		#			break
		#	if key == True:
		#		break

		#if key:
		#	for i in range(9):
		#		for j in range(9):
		#			for k in range(9):
		#				#combfreq[i][j][k] = combfreq[i][j][k] / 2
		#				if combfreq[i][j][k]- 50000 >= 0:
		#					combfreq[i][j][k] -= 50000
		#				else:
		#					combfreq[i][j][k] = 0
			
		#copy frequency table to retain tabel unchanged
		freq_cl = copy.deepcopy(freq)
		#combfreq_cl = copy.deepcopy(combfreq)

		nextgeneration = []
		
		#apply elitism and aging
		for i in range(pop_size):
			if random.random() < elitism[i]:
				assert population[i].age < 5
				if(random.random() < aging[population[i].age]):
					#re-calculate suvived solution because frequency can change for each generation
					population[i].reCalcFit(freq, freq_cl)
					nextgeneration.append(population[i])
		
		#offspring generation will proceed until population size met
		while len(nextgeneration) < pop_size:
			#select two solutions
		    parent_a = selection_op.select(population)
		    parent_b = selection_op.select(population)
			#crossover them
		    child_a_p, child_b_p = crossover_op.crossover(parent_a, parent_b)
			#mutate offsprings with given probability
		    if random.random() < 0.2:
		        child_a_p = mutation_op.mutate(child_a_p)
		    if random.random() < 0.2:
		        child_b_p = mutation_op.mutate(child_b_p)

			#actual evaluation of generated offsprings
			child_a = Solution()
			child_a.validate(sudoku, child_a_p, freq, freq_cl)
			child_b = Solution()
			child_b.validate(sudoku, child_b_p, freq, freq_cl)

			nextgeneration.append(child_a)
			nextgeneration.append(child_b)

		population = sorted(nextgeneration, key=attrgetter('fitness', 'fitness3', 'fitness2'))
		#printing parts		
		best = population[0]
		temp = copy.deepcopy(best.permutation)
		temp.sort()
		print temp
		#update current_best if this generation's best is more fitted or equal
		if best.fitness <= current_best.fitness:
		    current_best = best

		#printing parts
		lastbest = 0
		for i in range(len(population)):
			population[i].age += 1
			#check last best fitness			
			if best.fitness != population[i].fitness and lastbest == 0:
				lastbest = i
		generation += 1
		print generation, ' lastbest:', lastbest, ' ', population[0].fitness, ' ', population[0].fitness3, ' ', population[0].fitness2, '...', population[len(population)-1].fitness, current_best.fitness

		#update prequency table
		freq = freq_cl
		#combfreq = combfreq_cl

	#combine initial valuation and calculated vacant positions
	result = current_best.permutation + num_fix

	result.sort()

	return result

if __name__ == '__main__':
	pop_size = 200
	budget = 1000000
	for i in range(len(sys.argv)):
		if sys.argv[i] == '-p': #population size
			pop_size = int(sys.argv[i+1])
		elif sys.argv[i] == '-f': #budget limitation
			budget = int(sys.argv[i+1])

	sol = ga(sys.argv[len(sys.argv)-1], pop_size)
	for i in range(81):
		print sol[i][1], ' ',
		if i % 9 == 8:
			print ' '

