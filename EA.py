import math
import random
import sys
import matplotlib.pyplot as plt
import numpy as np
from Chromosome import Chromosome
from Overall_plots import OverallPlots
from Operations import Operations

START_SIZE = 75  # Population size at start.
MAX_EPOCHS = 1000  # Arbitrary number of test cycles. Default 1000!
MATING_PROBABILITY = 0.7  # Probability of two chromosomes mating. Range: 0.0 < MATING_PROBABILITY < 1.0
MUTATION_RATE = 0.001  # Mutation Rate. Range: 0.0 < MUTATION_RATE < 1.0
#pro mating aber 2 Kinder -> offspring_per_epoche=1 erstellt 2 Kinder von 2 Eltern. wert auf 20 entspricht 40 kinder von
#40 (nicht unbedignt unterschiedlichn) eltern
OFFSPRING_PER_GENERATION = 20  # New offspring created per generation. Range: 0 < OFFSPRING_PER_GENERATION < MAX_SELECT.
MINIMUM_SHUFFLES = 8  # For randomizing starting chromosomes
MAXIMUM_SHUFFLES = 20
PBC_MAX = 4  # Maximum Position-Based Crossover points. Range: 0 < PBC_MAX < 8 (> 8 isn't good).

MAX_LENGTH = 10 # chess board width.


#welcher crossover für welche phase besser
#mutation von corssover und crossover fixen, auf viele durchgänge testen

#eltern nicht ganzwegschmeiße. kinder erezugen random , zsm schmeisen und dann die besten kinder weiternehmen.
#vlt. nur ein kind? crossover dass das beste weitergegeben werden


class NQueen1:
    def __init__(self, startSize, maxEpochs, matingProb, mutationRate, generation, minShuffles,
                 maxShuffles, pbcMax, maxLength):
        self.mStartSize = startSize
        self.mEpochs = maxEpochs
        self.mMatingProbability = matingProb
        self.mMutationRate = mutationRate
        self.mOffspringPerGeneration = generation
        self.mMinimumShuffles = minShuffles
        self.mMaximumShuffles = maxShuffles
        self.mPBCMax = pbcMax
        self.mMaxLength = maxLength

        self.epoch = 0
        self.childCount = 0
        self.nextMutation = 0  # For scheduling mutations.
        self.mutations = 0
        self.population = []

        #counter für verschiedene permutationen
        self.order_based_co = 0
        self.partielle_mapped_co = 0
        self.position_based_co = 0
        self.array_o_b = [0]
        self.array_p_m = [0]
        self.array_p_b = [0]
        self.current_p_m = 0
        self.current_p_b = 0
        self.current_o_b = 0

        #boolean voariable ob minimum gefunden wurde
        self.found_mimimum = 0
        #beste und schlechteste Fitness
        self.worst=0
        self.best=0

        #variablen für fitness-epochen-plot
        self.current_best_fitness = 0
        self.array_fitness = [0]

        return

    def get_maximum(self):
        # Returns an array index.
        maximum = 0
        done = False

        while not done:
            foundNewMaximum = False
            popSize = len(self.population)
            for i in range(popSize):
                if i != maximum:
                    thisChromo = self.population[i]
                    thatChromo = self.population[maximum]
                    # The maximum has to be in relation to the Target.
                    if math.fabs(thisChromo.get_fitness() > thatChromo.get_fitness()):
                        maximum = i
                        foundNewMaximum = True

            if foundNewMaximum == False:
                done = True

        return maximum

    def get_minimum(self):
        # Returns an array index.
        minimum = 0
        done = False

        while not done:
            foundNewMinimum = False
            popSize = len(self.population)
            for i in range(popSize):
                if i != minimum:
                    thisChromo = self.population[i]
                    thatChromo = self.population[minimum]
                    # The minimum has to be in relation to the Target.
                    if math.fabs(thisChromo.get_fitness() < thatChromo.get_fitness()):
                        minimum = i
                        foundNewMinimum = True

            if foundNewMinimum == False:
                done = True

        return minimum

    def initialize_chromosomes(self):
        for i in range(self.mStartSize):
            #crossover des neuen Chromosoms wird bestimmt
            rand = random.randrange(0, 3)

            newChromo = Chromosome(self.mMaxLength, rand)
            self.population.append(newChromo)
            chromoIndex = len(self.population) - 1

            # Randomly choose the number of shuffles to perform.
            shuffles = random.randrange(self.mMinimumShuffles, self.mMaximumShuffles)
            Operations.exchange_mutation(self, chromoIndex, shuffles)

            newChromo = self.population[chromoIndex]
            newChromo.compute_fitness()

        return

    #wenn die eltern verschiedene crossover haben wird das übernommen von dem elternteil bessere fitness hat
    def do_mating(self):
        for i in range(self.mOffspringPerGeneration):
            #print("for schleife in mating")
            parentA = Operations.choose_first_parent(self)
            # Test probability of mating.
            getRand = random.randrange(0, 100)

            if getRand <= self.mMatingProbability * 100:
                #print ("im ersten if, mating")
                parentB = Operations.choose_second_parent(self,parentA)
                #print("parents gefunden, mating")
                #das crossover des Elternteils mit größerer Fitness wird gewählt
                chromA = self.population[parentA]
                chromB = self.population[parentB]

                if chromA.get_fitness() < chromB.get_fitness():
                    type = chromA.get_crossover()
                else:
                    type = chromB.get_crossover()

                if type == 0:
                    newChromo1 = Chromosome(self.mMaxLength, 0)
                    newChromo2 = Chromosome(self.mMaxLength, 0)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    Operations.partially_mapped_crossover(self, chromA, chromB, newIndex1, newIndex2)
                    self.partielle_mapped_co += 1
                    self.current_p_m +=1
                elif type == 1:
                    newChromo1 = Chromosome(self.mMaxLength, 1)
                    newChromo2 = Chromosome(self.mMaxLength, 1)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    # TODO entweder positionbased crossover oder bad_recombination_qu, beides crossover-typ 1
                    Operations.position_based_crossover(self, chromA, chromB, newIndex1, newIndex2)
                    self.position_based_co += 1
                    self.current_p_b += 1
                else:
                    newChromo1 = Chromosome(self.mMaxLength, 2)
                    newChromo2 = Chromosome(self.mMaxLength, 2)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    Operations.order_based_crossover(self, chromA, chromB, newIndex1, newIndex2)
                    self.order_based_co += 1
                    self.current_o_b += 1


                #mutation erstmal ausgeklammert, da nicht wichtig für rekombinationsbeobachtung

                if self.childCount - 1 == self.nextMutation:
                    Operations.exchange_mutation(self, newIndex1, 1)
                    #self.displacement_mutation(newIndex1)
                elif self.childCount == self.nextMutation:
                    Operations.exchange_mutation(self, newIndex2, 1)
                    #self.displacement_mutation(newIndex2)



                newChromo1 = self.population[newIndex1]
                # Konsole:
                #sys.stdout.write("Kind1 mit ")
                #sys.stdout.write(str(newChromo1.get_fitness()) + " Konflikten: ")
                #newChromo1.toStr() + "\n"
                self.childCount += 1

                newChromo2 = self.population[newIndex2]
                # Konsole:
                #sys.stdout.write("Kind2 mit ")
                #sys.stdout.write(str(newChromo2.get_fitness()) + " Konflikten: ")
                #newChromo2.toStr() + "\n"
                self.childCount += 1

                # Schedule next mutation.
                if math.fmod(self.childCount, Operations.math_round(self,1.0 / self.mMutationRate)) == 0:
                    self.nextMutation = self.childCount + random.randrange(0, Operations.math_round(self,1.0 / self.mMutationRate))

        return

    def print_best_solution(self, bestSolution=Chromosome(MAX_LENGTH,type)):
        board = []
        for i in range(self.mMaxLength):
            board.append([""] * self.mMaxLength)
            board[i][bestSolution.get_data(i)] = "Q"

        # Display the board.
        sys.stdout.write("Board:\n")
        for j in range(self.mMaxLength):
            for i in range(self.mMaxLength):
                if board[i][j] == "Q":
                    sys.stdout.write("Q ")
                else:
                    sys.stdout.write(". ")

            sys.stdout.write("\n")

        return

    def genetic_algorithm(self):

        done = False

        self.mutations = 0
        self.nextMutation = random.randrange(0, Operations.math_round(self, 1.0 / self.mMutationRate))

        while not done:
            popSize = len(self.population)
            for i in range(popSize):
                thisChromo = self.population[i]
                if thisChromo.get_fitness() == 0 or self.epoch == self.mEpochs:
                    done = True

            Operations.get_fitness(self)

            Operations.roulette_selection(self)

            self.do_mating()

            Operations.prep_next_epoch(self)

            self.array_p_m.append(self.current_p_m)
            self.array_p_b.append(self.current_p_b)
            self.array_o_b.append(self.current_o_b)
            self.current_p_m = 0
            self.current_p_b = 0
            self.current_o_b = 0

            self.epoch += 1

            #sys.stdout.write(str(self.array_p_m)+" ARRAY_P_M\n")
            #sys.stdout.write(str(self.array_p_b) + " ARRAY_P_B\n")
            #sys.stdout.write(str(self.array_o_b) + " ARRAY_O_B\n")

            # This is here simply to show the runtime status.
            sys.stdout.write("Epoche: " + str(self.epoch) + "\n")

        sys.stdout.write("done.\n")
        #sys.stdout.write(str(max(self.array_p_m)) + " maximale Value.\n")

        if self.epoch != self.mEpochs:
            popSize = len(self.population)
            for i in range(popSize):
                thisChromo = self.population[i]
                if thisChromo.get_fitness() == 0:
                    sys.stdout.write(str(thisChromo.toStr()) +" hat " + str(thisChromo.get_fitness()) + " Konflikte.\n")
                    self.print_best_solution(thisChromo)
                    self.found_mimimum=1
                    break

        sys.stdout.write("Completed " + str(self.epoch) + " epochs.\n")
        sys.stdout.write(
            "Encountered " + str(self.mutations) + " mutations in " + str(self.childCount) + " offspring.\n")
        sys.stdout.write(
            "Encountered " + str(self.partielle_mapped_co) + " partiell-mapped , " +
            str(self.position_based_co) + " positioon-based and " + str(self.order_based_co)
            + " order-based Crossover.\n")

        return


if __name__ == '__main__':
    foundMinimum = 0
    array = [0, 0, 0]
    counter = 0
    minarray = []
    while(counter != 3):
        print("_________________________________________________________________")
        sys.stdout.write("COUNTER: " + str(counter)+"\n")
        nq1 = NQueen1(START_SIZE, MAX_EPOCHS, MATING_PROBABILITY, MUTATION_RATE, OFFSPRING_PER_GENERATION,
                      MINIMUM_SHUFFLES, MAXIMUM_SHUFFLES, PBC_MAX, MAX_LENGTH)

        nq1.initialize_chromosomes()
        nq1.genetic_algorithm()
        #zeige pro Druchlauf permutation-Anzahl und permutationen pro Epoche
        #nq1.show_permutation_amount()
        OverallPlots.show_crossover_per_epoche(nq1)

        array = np.array(array) + Operations.get_best_crossover(nq1)

        foundMinimum += Operations.get_foundMimimum(nq1)
        minarray.append(Operations.get_best(nq1))

        counter += 1

    print(array)
    OverallPlots.show_overall_permutation_amount(array, "nQueens")
    sys.stdout.write("minArray: "+str(minarray)+"\n")
    OverallPlots.show_overall_minima(minarray, "nQueens")
    sys.stdout.write("Gefundene Minima: "+str(foundMinimum)+"\n")
    sys.stdout.write("Median aller Minima: " + str(np.median(minarray)) + "\n")



