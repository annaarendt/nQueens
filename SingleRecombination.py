import sys
import matplotlib.pyplot as plt
from Chromosome import Chromosome
from EA import NQueen1
from Overall_plots import OverallPlots

START_SIZE = 2  # Population size at start.
MAX_EPOCHS = 0  # Arbitrary number of test cycles. EIGENTLICH AUF 1000 GESETZT!
OFFSPRING_PER_GENERATION = 1  # New offspring created per generation. Range: 0 < OFFSPRING_PER_GENERATION < MAX_SELECT.
MINIMUM_SHUFFLES = 8  # For randomizing starting chromosomes
MAXIMUM_SHUFFLES = 20
PBC_MAX = 4  # Maximum Position-Based Crossover points. Range: 0 < PBC_MAX < 8 (> 8 isn't good).

MAX_LENGTH = 6  # chess board width.


class SingleRekomb:
    def __init__(self, startSize, maxEpochs, generation, minShuffles,
                 maxShuffles, pbcMax, maxLength):
        self.mStartSize = startSize
        self.mEpochs = maxEpochs
        self.mOffspringPerGeneration = generation
        self.mPBCMax = pbcMax
        self.mMaxLength = maxLength
        self.mMinimumShuffles = minShuffles
        self.mMaximumShuffles = maxShuffles

        self.epoch = 0
        self.childCount = 0
        self.population = []
        #1. und 2. stelle für Konflikte der Eltern, letzen beiden der Kinder
        self.conflict_array = [0] * 4
        self.mutations = 0

        return


    def do_mating(self, recomb):

        for i in range(self.mOffspringPerGeneration):
            chromA = self.population[0]
            chromB = self.population[1]

            newChromo1 = Chromosome(self.mMaxLength, 0)
            newChromo2 = Chromosome(self.mMaxLength, 0)
            self.population.append(newChromo1)
            newIndex1 = len(self.population) - 1
            self.population.append(newChromo2)
            newIndex2 = len(self.population) - 1

            sys.stdout.write("Parent1 mit ")
            sys.stdout.write(str(chromA.get_fitness()) + " Konflikten: ")
            chromA.toStr() + "\n"
            sys.stdout.write("Parent2 mit ")
            sys.stdout.write(str(chromB.get_fitness()) + " Konflikten: ")
            chromB.toStr() + "\n"

            #TODO diese Zeile ändern, je nach gewünschter Rekombination (bad funktioniert nur als self.bad...)
            NQueen1.position_based_crossover(self, chromA, chromB, newIndex1, newIndex2)

            newChromo1 = self.population[newIndex1]
            # Konsole:
            sys.stdout.write("Kind1 mit ")
            sys.stdout.write(str(newChromo1.get_fitness()) + " Konflikten: ")
            newChromo1.toStr()+"\n"

            newChromo2 = self.population[newIndex2]
            # Konsole:
            sys.stdout.write("Kind2 mit ")
            sys.stdout.write(str(newChromo2.get_fitness()) + " Konflikten: ")
            newChromo2.toStr()+"\n"

            self.childCount += 2
            sys.stdout.write(str(len(self.population)) + " Länge der Population\n")

            self.set_conflict_array(chromA.get_fitness(), chromB.get_fitness(),
                                    newChromo1.get_fitness(), newChromo2.get_fitness())

            # plot für jeden einzelnen Druchgang
            #self.show_conflicts(chromA, chromB, newChromo1, newChromo2)

        return

    def set_conflict_array(self, chromA, chromB, newChromo1, newChromo2):
        self.conflict_array = [chromA, chromB, newChromo1, newChromo2]
        return

    def get_conflict_array(self):
        return self.conflict_array


    def prep_next_epoch(self):
        # Reset flags for selected individuals.
        popSize = len(self.population)
        for i in range(popSize):
            thisChromo = self.population[i]
            Chromosome.set_selected(thisChromo, False)

        return


    def recomb(self, recomb, a, b, c, d):
        testdata = {
            0: NQueen1.partially_mapped_crossover(self, a, b, c, d),
            1: NQueen1.position_based_crossover(self, a, b, c, d),
            2: NQueen1.order_based_crossover(self, a, b, c, d),
            3: NQueen1.bad_recombination(self, a, b, c, d)
        }
        return testdata.get(recomb, "nothing")


    def genetic_algorithm(self, recomb):

        done = False

        while not done:
            popSize = len(self.population)
            for i in range(popSize):
                if self.epoch == self.mEpochs:
                    done = True

            self.do_mating(recomb)

            self.prep_next_epoch()

            self.epoch += 1

        sys.stdout.write("done.\n")

        sys.stdout.write("Completed " + str(self.epoch) + " epochs.\n")
        return

    #Konflikte der Eltern und Kinder werden in einem Plot gezeigt
    def show_conflicts(self, chromP1, chromP2, chromK1, chromK2 ):
        chromosomes=["P1", "P2", "K1", "K2"]
        conflicts=[chromP1.get_fitness(), chromP2.get_fitness(), chromK1.get_fitness(),
                   chromK2.get_fitness()]

        plt.scatter(chromosomes, conflicts, c='b', marker='o')
        plt.xlabel('Chromosomen', fontsize=16)
        plt.ylabel('Konflikte', fontsize=16)
        plt.title('scatter plot - parents vs offsprings', fontsize=20)
        plt.show()


if __name__ == '__main__':
    COUNTER = 0
    END = 500
    p1 = [0] * END
    p2 = [0] * END
    k1 = [0] * END
    k2 = [0] * END
    #TODO rekombination wegen beschriftung austauschen
    #0= partiallymappd, 1= positionbased, 2=orderbased, 3=bad-positionbased
    recomb = 1
    while (COUNTER != END):
        sr1 = SingleRekomb(START_SIZE, MAX_EPOCHS, OFFSPRING_PER_GENERATION, MINIMUM_SHUFFLES, MAXIMUM_SHUFFLES,
                      PBC_MAX, MAX_LENGTH)

        NQueen1.initialize_chromosomes(sr1)
        sr1.genetic_algorithm(recomb)

        p1[COUNTER] = sr1.get_conflict_array()[0]
        p2[COUNTER] = sr1.get_conflict_array()[1]
        k1[COUNTER] = sr1.get_conflict_array()[2]
        k2[COUNTER] = sr1.get_conflict_array()[3]
        #konsole_konflikte je Durchlauf
        str1 = ','.join(str(e) for e in sr1.get_conflict_array())
        print("Konflikte: "+str1)
        sys.stdout.write("COUNTER: " + str(COUNTER)+"\n")

        COUNTER += 1

    OverallPlots.boxplot(p1, p2, k1, k2, recomb, 0)
    OverallPlots.percentage_table(p1, p2, k1, k2, recomb, 0)









