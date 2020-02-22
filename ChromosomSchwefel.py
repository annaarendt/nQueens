

from deap import benchmarks
import random


class ChromosomeSchwefel:
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


    def compute_fitness(self):
        self.mFitness=benchmarks.schwefel(self.toArray())[0]
        return

    def set_selection_probability(self, probability):
        self.mSelectionProbability = probability
        return

    def get_selection_probability(self):
        return self.mSelectionProbability

    def set_selected(self, isSelected):
        self.mSelected = isSelected
        return

    def get_selected(self):
        return self.mSelected

    def set_fitness(self, score):
        self.mFitness = score
        return

    def get_fitness(self):
        return self.mFitness


    def set_crossover(self, type):
        self.mCrossover = type
        return

    def get_crossover(self):
        return self.mCrossover

    def set_data(self, index, value):
        self.mData[index] = value
        return

    def get_data(self, index):
        return self.mData[index]

    # Gene des Chromosoms als String darstellen
    def toStr(self):
        array = [0]*self.mMaxLength
        for i in range(self.mMaxLength):
            array[i] = self.get_data(i)

        print(str(array))
        return str(array)

# Gene des Chromosoms als Liste darstellen
    def toArray(self):
        array = [0]*self.mMaxLength
        for i in range(self.mMaxLength):
            array[i] = self.get_data(i)

        return array

