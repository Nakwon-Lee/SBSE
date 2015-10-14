import sys
import math
import random
import copy
from operator import attrgetter

evals = 0
budget = 0

class Solution:
	def __init__(self)
		self.permutation = []
		self.fitness = 81
		self.fitness2 = 10000
		self.age = 1

    def genRanSol(self, sudoku, num, freq):		
		sudoku_cl = copy.deepcopy(sudoku)
		num_cl = copy.deepcopy(num)
		
		random.shuffle(num_cl)
		
		self.fitness = 0
		self.fitness2 = 0

		for i in range(len(num_cl)):
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
				fitness2 = fitness2 + freq[num_cl[i]][a-1]
				freq[num_cl[i]][a-1] += 1

	def validate(self, sudoku, perm, freq):
		sudoku_cl = copy.deepcopy(sudoku)
		
		self.fitness = 0
		self.fitness2 = 0

		for i in range(len(perm)):
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
		
				sudoku_cl.col[num_cl[i]//9].add(a)
				sudoku_cl.row[num_cl[i]%9].add(a)
				sudoku_cl.box[calcIdxOfBox(num_cl[i])].add(a)
				fitness2 = fitness2 + freq[num_cl[i]][a-1]
				freq[num_cl[i]][a-1] += 1

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
	def __init__(self,num)
		self.col = []
		self.row = []
		self.box = []
	
		for i in range(num)
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
			child_a_set.add(parent_b.permutation[i])
			child_b_perm.append(parent_a.permutation[i])
			child_b_set.add(parent_a.permutation[i])

		ch_a_idx = 0
		p_a_idx = 0
		ch_b_idx = 0
		p_b_idx = 0

		while ch_a_idx < cp1:
			if parent_a.permutation[p_a_idx] not in child_a_set:
				child_a_perm.insert(ch_a_idx, parent_a.permutation[p_a_idx])
				ch_a_idx += 1
			p_a_idx += 1
		ch_a_idx = cp2 + 1

		while ch_b_idx < cp1:
			if parent_b.permutation[p_b_idx] not in child_b_set:
				child_b_perm.insert(ch_b_idx, parent_b.permutation[p_b_idx])
				ch_b_idx += 1
			p_b_idx += 1
		ch_b_idx = cp2 + 1

		while ch_a_idx < size:
			if parent_a.permutation[p_a_idx] not in child_a_set:
				child_a_perm.append(parent_a.permutation[p_a_idx])
				ch_a_idx += 1
			p_a_idx += 1

		while ch_b_idx < size:
			if parent_b.permutation[p_b_idx] not in child_b_set:
				child_b_perm.append(parent_b.permutation[p_b_idx])
				ch_b_idx += 1
			p_b_idx += 1

		child_a = Solution(child_a_perm)
		child_b = Solution(child_b_perm)

		return child_a, child_b
       
class Mutation:
    def mutate(self, solution):
        size = len(solution.permutation)

        mp1 = random.randrange(size)
        mp2 = random.randrange(size)
        while(mp2 == mp1):
            mp2 = random.randrange(size)
        solution.permutation[mp1], solution.permutation[mp2] = solution.permutation[mp2], solution.permutation[mp1]
        return solution

class BinaryTournament:
    def select(self, population):
        i = random.randrange(len(population))
        j = random.randrange(len(population))
        while i == j:
           j = random.randrange(len(population))
        
        a = population[i]
        b = population[j]
        if a.fitness < b.fitness:
            return a
        else: 
            return b

def ga(filename):
	num, num_fix, sudoku = read_data(filename)

	population = []
	selection_op = BinaryTournament()
	crossover_op = CrossoverOrder()
	mutation_op = Mutation()

	elitism = []
	aging = []

	pop_size = 200
	for i in range(pop_size):
		perm = []
		for j in range(num):
			perm.append(j)
		random.shuffle(perm)            
		new_individual = Solution(perm)
		evaluate(new_individual)
		population.append(new_individual)

	for i in range(pop_size):
		elitism.append(0.9 * pow( (float(pop_size-i)/float(pop_size)),5 ))

	for i in range(5):
		aging.append( math.sqrt(math.sqrt(0.8- (float(i) * 0.2) )) )

	population = sorted(population, key=attrgetter('fitness'))
	current_best = population[0]

	generation = 0

	while evals < budget:
		nextgeneration = []
		
		#apply elitism and aging
		for i in range(pop_size):
			if random.random() < elitism[i]:
				if(random.random() < aging[population[i].age]):
					nextgeneration.append(population[i])
		
		while len(nextgeneration) < pop_size:
		    parent_a = selection_op.select(population)
		    parent_b = selection_op.select(population)
		    child_a, child_b = crossover_op.crossover(parent_a, parent_b)
		    if random.random() < 0.5:
		        child_a = mutation_op.mutate(child_a)
		    if random.random() < 0.5:
		        child_b = mutation_op.mutate(child_b)
		    
		    evaluate(child_a)
		    evaluate(child_b)

		    nextgeneration.append(child_a)
		    nextgeneration.append(child_b)

		    # print child_a
		    # print child_b
		population = sorted(nextgeneration, key=attrgetter('fitness'))
		best = population[0]
		# print generation, best.fitness, best.translate(coded_word)
		if best.fitness < current_best.fitness:
		    current_best = best
		    # print current_best_str
		#print ",".join([str(generation), str(current_best.fitness)])
		for i in range(len(population)):
			population[i].age += 1
		generation += 1

	return current_best

if __name__ == '__main__':
	budget = 500000
	sol = ga(sys.argv[1])
	for i in range(len(sol.permutation)):
		print sol.permutation[i]+1, ',',
	print ' '
	print sol.fitness

