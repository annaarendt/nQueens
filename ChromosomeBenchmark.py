

from deap import benchmarks
import random
from Chromosome import Chromosome


class ChromosomeDEAP:
    #crossover-type 0,1 oder 2 (partitielle_mapped, position_based und order_based)
    def __init__(self, maxLength, crossover_type):
        self.mMaxLength = maxLength
        self.mSelected = False
        self.mSelectionProbability = 0.0
        self.mCrossover = crossover_type
        self.mFitness = 0.0

        self.mData = [0] * maxLength
        for i in range(self.mMaxLength):
            self.mData[i] = random.randrange(-500,500,1)
        return

    #TODO beliebiges benchmarkproblem einsetzen
    def compute_fitness(self):
        self.mFitness = benchmarks.schwefel(Chromosome.toArray(self))[0]
        return



