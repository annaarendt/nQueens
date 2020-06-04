

from deap import benchmarks
import random
from Chromosome import Chromosome


class ChromosomeDEAP:
    #crossover-type 0,1, 2 oder 3 (partitially_mapped, position_based, two_point und order_based)
    def __init__(self, maxLength, crossover_type):
        self.mMaxLength = maxLength
        self.mSelected = False
        self.mSelectionProbability = 0.0
        self.mCrossover = crossover_type
        self.mFitness = 0.0

        self.mData = [0] * maxLength
        for i in range(self.mMaxLength):
            #TODO passende range je nach Funktion auswählen
            self.mData[i] = random.randrange(-500,500,1) #schwefel
            #self.mData[i] = random.randrange(-15,30,1) #ackley
            #self.mData[i] = random.randrange(-6, 6, 1)  # himmelblau
            #self.mData[i] = random.randrange(-600, 600, 1)  # griewank
            #self.mData[i] = round(random.uniform(-5.12, 5.12), 2)  # rastrigin
        return

    #TODO beliebiges benchmarkproblem einsetzen
    def compute_fitness(self):
        self.mFitness = benchmarks.schwefel(Chromosome.toArray(self))[0]
        #self.mFitness = benchmarks.ackley(Chromosome.toArray(self))[0]
        #self.mFitness = benchmarks.himmelblau(Chromosome.toArray(self))[0]
        #self.mFitness = benchmarks.griewank(Chromosome.toArray(self))[0]
        #self.mFitness = benchmarks.rastrigin(Chromosome.toArray(self))[0]
        return



