import sys
import math
import random
import copy
from operator import attrgetter

evals = 0
budget = 0

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

class Sudoku:
	def __init__(self,num):
		self.col = []
		self.row = []
		self.box = []
	
		for i in range(num):
			self.col.append(set())
			self.row.append(set())
			self.box.append(set())

def read_data(filename):
	global sudoku
	lines = open(filename).readlines()
	sudoku = Sudoku(9)
	num_fix = []
	num = []

	pos = 0
	for line in lines:
		for i in range(len(line)):
			if line[i].isdigit():
				sudoku.col[pos//9].add(int(line[i]))
				sudoku.row[pos%9].add(int(line[i]))
				sudoku.box[calcIdxOfBox(pos)].add(int(line[i]))
				num_fix.append((pos, int(line[i])))
				pos += 1
			elif line[i] == '.':
				num.append(pos)
				pos += 1
			else:
				pass

	return num, num_fix, sudoku

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

		for i in range(cp1, cp2 + 1):
			child_a_perm.append(parent_b.permutation[i])
			child_a_set.add(parent_b.permutation[i][0])
			child_b_perm.append(parent_a.permutation[i])
			child_b_set.add(parent_a.permutation[i][0])

		ch_a_idx = 0
		p_a_idx = 0
		ch_b_idx = 0
		p_b_idx = 0

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
       
class Mutation:
    def mutate(self, solution):
        size = len(solution)

        mp1 = random.randrange(size)
        mp2 = random.randrange(size)
        while(mp2 == mp1):
            mp2 = random.randrange(size)
        solution[mp1], solution[mp2] = solution[mp2], solution[mp1]
        return solution

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

class Solution:
	def __init__(self):
		self.permutation = []
		self.fitness = 81
		self.fitness2 = 99999999.0
		self.fitness3 = 99999999.0
		self.age = 0

	def genRanSol(self, sudoku, num, freq, freq_cl, combfreq, combfreq_cl):
		global evals
		sudoku_cl = copy.deepcopy(sudoku)
		num_cl = copy.deepcopy(num)

		#random.shuffle(num_cl)

		self.fitness = 0
		self.fitness2 = 0
		self.fitness3 = 0

		ft = [0 for i in range(3)]

		for i in range(3):
			possible_numbers = set([1,2,3,4,5,6,7,8,9])
			possible_numbers = possible_numbers - sudoku_cl.col[num_cl[i]//9] - sudoku_cl.row[num_cl[i]%9] - sudoku_cl.box[calcIdxOfBox(num_cl[i])]

			if len(possible_numbers) == 0:
				assert i >= 3
				self.permutation.append((num_cl[i], 0))
				self.fitness += 1
			else:
				a = random.choice(list(possible_numbers))
				if i < 3:
					ft[i] = a
				self.permutation.append((num_cl[i], a))
				sudoku_cl.col[num_cl[i]//9].add(a)
				sudoku_cl.row[num_cl[i]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(num_cl[i])].add(a)
				if i == 2:
					self.fitness3 = combfreq[ft[0]-1][ft[1]-1][ft[2]-1]

		for i in range(3, len(num_cl)):
			possible_numbers = set([1,2,3,4,5,6,7,8,9])
			possible_numbers = possible_numbers - sudoku_cl.col[num_cl[i]//9] - sudoku_cl.row[num_cl[i]%9] - sudoku_cl.box[calcIdxOfBox(num_cl[i])]

			if len(possible_numbers) == 0:
				self.permutation.append((num_cl[i], 0))
				self.fitness += 1
			else:
				a = random.choice(list(possible_numbers))
				self.permutation.append((num_cl[i], a))
				sudoku_cl.col[num_cl[i]//9].add(a)
				sudoku_cl.row[num_cl[i]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(num_cl[i])].add(a)
				assert freq[num_cl[i]][a-1] >= 0
				self.fitness2 = self.fitness2 + freq[num_cl[i]][a-1]

		#calculate frequencies
		combfreq_cl[ft[0]-1][ft[1]-1][ft[2]-1] += ((81.0-float(self.fitness))/81.0)

		for i in range(3, len(num_cl)):
			if self.permutation[i][1] != 0:		
				freq_cl[num_cl[i]][self.permutation[i][1]-1] += ((81.0-float(self.fitness))/81.0)

		evals += 1

	def validate(self, sudoku, perm, freq, freq_cl, combfreq, combfreq_cl):
		global evals
		sudoku_cl = copy.deepcopy(sudoku)
		
		self.fitness = 0
		self.fitness2 = 0
		self.fitness3 = 0

		ft = [0 for i in range(3)]

		for i in range(3):
			possible_numbers = set([1,2,3,4,5,6,7,8,9])
			possible_numbers = possible_numbers - sudoku_cl.col[perm[i][0]//9] - sudoku_cl.row[perm[i][0]%9] - sudoku_cl.box[calcIdxOfBox(perm[i][0])]

			if len(possible_numbers) == 0:
				assert i >= 3
				self.permutation.append((perm[i][0], 0))
				self.fitness += 1
			else:
				a = 0
				if perm[i][1] in possible_numbers:
					a = perm[i][1]
					self.permutation.append((perm[i][0], a)) 
				else:
					a = random.choice(list(possible_numbers))
					self.permutation.append((perm[i][0], a))
				if i < 3:
					ft[i] = a
				sudoku_cl.col[perm[i][0]//9].add(a)
				sudoku_cl.row[perm[i][0]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(perm[i][0])].add(a)
				if i == 2:
					self.fitness3 = combfreq[ft[0]-1][ft[1]-1][ft[2]-1]
					

		for i in range(3, len(perm)):
			possible_numbers = set([1,2,3,4,5,6,7,8,9])
			possible_numbers = possible_numbers - sudoku_cl.col[perm[i][0]//9] - sudoku_cl.row[perm[i][0]%9] - sudoku_cl.box[calcIdxOfBox(perm[i][0])]

			if len(possible_numbers) == 0:
				self.permutation.append((perm[i][0], 0))
				self.fitness += 1
			else:
				a = 0
				if perm[i][1] in possible_numbers:
					a = perm[i][1]
					self.permutation.append((perm[i][0], a)) 
				else:
					a = random.choice(list(possible_numbers))
					self.permutation.append((perm[i][0], a))
				sudoku_cl.col[perm[i][0]//9].add(a)
				sudoku_cl.row[perm[i][0]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(perm[i][0])].add(a)
				assert freq[perm[i][0]][a-1] >= 0
				self.fitness2 = self.fitness2 + freq[perm[i][0]][a-1]

		#calculate frequencies
		combfreq_cl[ft[0]-1][ft[1]-1][ft[2]-1] += ((81.0-float(self.fitness))/81.0)
		
		for i in range(3, len(perm)):
			if self.permutation[i][1] != 0:		
				freq_cl[perm[i][0]][self.permutation[i][1]-1] += ((81.0-float(self.fitness))/81.0)

		evals += 1

	def reCalcFit(self, freq, freq_cl, combfreq, combfreq_cl):
		global evals
		self.fitness2 = 0
		self.fitness3 = 0
		evals += 1

		ft = [0 for i in range(3)]		

		for i in range(len(self.permutation)):
			assert freq[self.permutation[i][0]][self.permutation[i][1]-1] >= 0
			self.fitness2 = self.fitness2 + freq[self.permutation[i][0]][self.permutation[i][1]-1]
			if self.permutation[i][1] != 0:
				if i < 3:
					ft[i] = self.permutation[i][1]
				if i == 2:
					self.fitness3 = combfreq[ft[0]-1][ft[1]-1][ft[2]-1]
					combfreq_cl[ft[0]-1][ft[1]-1][ft[2]-1] += ((81.0-float(self.fitness))/81.0)
				freq_cl[self.permutation[i][0]][self.permutation[i][1]-1] += ((81.0-float(self.fitness))/81.0)

def ga(filename, pop):

	print pop, budget

	num, num_fix, sudoku = read_data(filename)

	combfreq = {}
	for i in range(9):
		combfreq[i] = {}
		for j in range(9):
			combfreq[i][j] = [0 for col in range(9)]

	freq = [[0 for col in range(9)] for row in range(81)]

	population = []
	#selection_op = BinaryTournament()
	selection_op = BinaryTournamentComb()
	#crossover_op = CrossoverZeroFirst()
	crossover_op = CrossoverOrder()
	#mutation_op = Mutation()
	#mutation_op = MutationZero()
	mutation_op = MutationSwapValueFixOrder()
	#mutation_op = MutationSwapZeroToFirst()
	#mutation_op = MutationSwapZero()

	elitism = []
	aging = []

	freq_cl = copy.deepcopy(freq)
	combfreq_cl = copy.deepcopy(combfreq)

	random.shuffle(num)

	pop_size = pop
	for i in range(pop_size):
		new_individual = Solution()
		new_individual.genRanSol(sudoku, num, freq, freq_cl, combfreq, combfreq_cl)
		population.append(new_individual)

	freq = freq_cl
	combfreq = combfreq_cl

	for i in range(pop_size):
		elitism.append(0.7 * pow( (float(pop_size-i)/float(pop_size)),5 ))

	for i in range(5):
		aging.append( math.sqrt(math.sqrt(0.8 - (float(i) * 0.2) )) )

	#print elitism
	#print aging

	population = sorted(population, key=attrgetter('fitness', 'fitness3', 'fitness2'))
	current_best = population[0]

	generation = 0

	while (evals < budget) and (current_best.fitness != 0):

		key = False
		for i in range(81):
			for j in range(9):
				if freq[i][j] >= 99999:
					key = True
					break
			if key == True:
				break

		if key:
			for i in range(81):
				for j in range(9):
						#freq[i][j] = int(freq[i][j] / 2)
						if freq[i][j] - 50000 >= 0:
							freq[i][j] -= 50000
						else:
							freq[i][j] = 0

		key = False
		for i in range(9):
			for j in range(9):
				for k in range(9):
					if combfreq[i][j][k] >= 99999:
						key = True
						break
				if key == True:
					break
			if key == True:
				break

		if key:
			for i in range(9):
				for j in range(9):
					for k in range(9):
						if combfreq[i][j][k]- 50000 >= 0:
							combfreq[i][j][k] -= 50000
						else:
							combfreq[i][j][k] = 0
			

		freq_cl = copy.deepcopy(freq)
		combfreq_cl = copy.deepcopy(combfreq)

		nextgeneration = []
		
		#apply elitism and aging
		for i in range(pop_size):
			if random.random() < elitism[i]:
				assert population[i].age < 5
				if(random.random() < aging[population[i].age]):
					population[i].reCalcFit(freq, freq_cl, combfreq, combfreq_cl)
					nextgeneration.append(population[i])
		
		while len(nextgeneration) < pop_size:
		    parent_a = selection_op.select(population)
		    parent_b = selection_op.select(population)
		    child_a_p, child_b_p = crossover_op.crossover(parent_a, parent_b)
		    if random.random() < 0.4:
		        child_a_p = mutation_op.mutate(child_a_p)
		    if random.random() < 0.4:
		        child_b_p = mutation_op.mutate(child_b_p)

			child_a = Solution()

			child_a.validate(sudoku, child_a_p, freq, freq_cl, combfreq, combfreq_cl)

			child_b = Solution()

			child_b.validate(sudoku, child_b_p, freq, freq_cl, combfreq, combfreq_cl)

			nextgeneration.append(child_a)
			nextgeneration.append(child_b)

		population = sorted(nextgeneration, key=attrgetter('fitness', 'fitness3', 'fitness2'))
		best = population[0]
		temp = copy.deepcopy(best.permutation)
		#temp.sort()
		print temp
		if best.fitness < current_best.fitness:
		    current_best = best

		lastbest = 0
		for i in range(len(population)):
			population[i].age += 1
			#check last best fitness			
			if best.fitness != population[i].fitness and lastbest == 0:
				lastbest = i
		generation += 1
		print generation, ' lastbest:', lastbest, ' ', population[0].fitness, ' ', population[0].fitness3, ' ', population[0].fitness2, '...', population[len(population)-1].fitness, current_best.fitness

		freq = freq_cl
		combfreq = combfreq_cl

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

