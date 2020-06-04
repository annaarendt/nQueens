import math
import random
import sys
import numpy as np
from Chromosome import Chromosome
from ChromosomeBenchmark import ChromosomeDEAP
from EA import NQueen1, MAXIMUM_SHUFFLES, MINIMUM_SHUFFLES
from Overall_plots import OverallPlots
from Operations import Operations

START_SIZE = 75  # Population size at start.
MAX_EPOCHS = 1000  # Arbitrary number of test cycles. Default 1000!
MATING_PROBABILITY = 0.7  # Probability of two chromosomes mating. Range: 0.0 < MATING_PROBABILITY < 1.0
MUTATION_RATE = 0.001  # Mutation Rate. Range: 0.0 < MUTATION_RATE < 1.0
OFFSPRING_PER_GENERATION = 20  # New offspring created per generation. Range: 0 < OFFSPRING_PER_GENERATION < MAX_SELECT.0
PBC_MAX = 4  # Maximum Position-Based Crossover points. Range: 0 < PBC_MAX < 8 (> 8 isn't good).

MAX_LENGTH = 8 # length of chromosom.


class NDEAP:
    def __init__(self, startSize, maxEpochs, matingProb, mutationRate, generation, pbcMax, maxLength):
        NQueen1.__init__(self, startSize, maxEpochs, matingProb, mutationRate, generation, MINIMUM_SHUFFLES,
                 MAXIMUM_SHUFFLES, pbcMax, maxLength)

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
                    # Der Fitness-Betrag wird hier betrachtet.
                    if abs(Chromosome.get_fitness(thisChromo)) > abs(Chromosome.get_fitness(thatChromo)):
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
                    # Der Fitness-Betrag wird hier betrachtet.
                    if abs(Chromosome.get_fitness(thisChromo)) < abs(Chromosome.get_fitness(thatChromo)):
                        minimum = i
                        foundNewMinimum = True

            if foundNewMinimum == False:
                done = True

        return minimum


    def initialize_chromosomes(self):
        for i in range(self.mStartSize):
            #crossover des neuen Chromosoms wird bestimmt
            rand = random.randrange(0, 4)


            newChromo = ChromosomeDEAP(self.mMaxLength, rand)
            self.population.append(newChromo)
            chromoIndex = len(self.population) - 1

            newChromo = self.population[chromoIndex]
            newChromo.compute_fitness()

        return

    #wenn die eltern verschiedene crossover haben wird das übernommen von dem elternteil bessere fitness hat
    def do_mating(self):
        for i in range(self.mOffspringPerGeneration):
            parentA = Operations.choose_first_parent(self)
            # Test probability of mating.
            getRand = random.randrange(0, 100)

            if getRand <= self.mMatingProbability * 100:
                parentB = Operations.choose_second_parent(self, parentA)
                #das crossover des Elternteils mit größerer Fitness wird gewählt
                chromA = self.population[parentA]
                chromB = self.population[parentB]

                if Chromosome.get_fitness(chromA) < Chromosome.get_fitness(chromB):
                    type = Chromosome.get_crossover(chromA)
                else:
                    type = Chromosome.get_crossover(chromB)

                if type == 0:
                    newChromo1 = ChromosomeDEAP(self.mMaxLength, 0)
                    newChromo2 = ChromosomeDEAP(self.mMaxLength, 0)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    Operations.partially_mapped_crossover(self, chromA, chromB, newIndex1, newIndex2)
                    self.partielle_mapped_co += 1
                    self.current_p_m +=1
                elif type == 1:
                    newChromo1 = ChromosomeDEAP(self.mMaxLength, 1)
                    newChromo2 = ChromosomeDEAP(self.mMaxLength, 1)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    # TODO entweder positionbased crossover oder bad_recombination, beides crossover-typ 1
                    #Operations.bad_recombination_deap(chromA, chromB, newIndex1, newIndex2)
                    Operations.position_based_crossover(self, chromA, chromB, newIndex1, newIndex2)
                    self.position_based_co += 1
                    self.current_p_b += 1
                elif type == 2:
                    newChromo1 = ChromosomeDEAP(self.mMaxLength, 2)
                    newChromo2 = ChromosomeDEAP(self.mMaxLength, 2)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    Operations.two_point_crossover(self, chromA, chromB, newIndex1, newIndex2)
                    self.two_point_co += 1
                    self.current_t_p += 1
                else:
                    newChromo1 = ChromosomeDEAP(self.mMaxLength, 3)
                    newChromo2 = ChromosomeDEAP(self.mMaxLength, 3)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    Operations.order_based_crossover(self, chromA, chromB, newIndex1, newIndex2)
                    self.order_based_co += 1
                    self.current_o_b += 1


                if self.childCount - 1 == self.nextMutation:
                    #Operations.new_gene_mutation(self,newIndex1)
                    Operations.exchange_mutation(self, newIndex1, 1)

                elif self.childCount == self.nextMutation:
                    #self.new_gene_mutation(newIndex2)
                    Operations.exchange_mutation(self, newIndex2, 1)


                newChromo1 = self.population[newIndex1]
                self.childCount += 1

                newChromo2 = self.population[newIndex2]
                self.childCount += 1

                # Schedule next mutation.
                if math.fmod(self.childCount, Operations.math_round(self, 1.0 / self.mMutationRate)) == 0:
                    self.nextMutation = self.childCount + random.randrange(0, Operations.math_round(self, 1.0 / self.mMutationRate))

        return

    def genetic_algorithm(self):

        done = False

        self.mutations = 0
        self.nextMutation = random.randrange(0, Operations.math_round(self, 1.0 / self.mMutationRate))

        while not done:
            popSize = len(self.population)
            for i in range(popSize):
                thisChromo = self.population[i]
                if Chromosome.get_fitness(thisChromo) == 0 or self.epoch == self.mEpochs:
                    done = True

            Operations.get_fitness(self)

            Operations.roulette_selection(self)

            self.do_mating()

            Operations.prep_next_epoch(self)

            self.array_p_m.append(self.current_p_m)
            self.array_p_b.append(self.current_p_b)
            self.array_o_b.append(self.current_o_b)
            self.array_t_p.append(self.current_t_p)
            self.current_p_m = 0
            self.current_p_b = 0
            self.current_o_b = 0
            self.current_t_p = 0

            self.epoch += 1

            # This is here simply to show the runtime status.
            sys.stdout.write("Epoche: " + str(self.epoch) + "\n")

        sys.stdout.write("done.\n")

        if self.epoch != self.mEpochs:
            popSize = len(self.population)
            for i in range(popSize):
                thisChromo = self.population[i]
                if Chromosome.get_fitness(thisChromo) == 0:
                    sys.stdout.write(str(Chromosome.toStr(thisChromo))+" hat " + str(Chromosome.get_fitness(thisChromo)) + " fitness.\n")
                    Operations.set_foundMimimum(self, 1)
                    break

        sys.stdout.write("Completed " + str(self.epoch) + " epochs.\n")
        sys.stdout.write(
            "Encountered " + str(self.mutations) + " mutations in " + str(self.childCount) + " offspring.\n")
        sys.stdout.write(
            "Encountered " + str(self.partielle_mapped_co) + " partiell-mapped , " +
            str(self.position_based_co) + " positioon-based and " + str(self.order_based_co)
            + " order-based Crossover." + str(self.two_point_co)+ " two-point Crossover.\n")

        return


if __name__ == '__main__':
    foundMinimum = 0
    array = [0, 0, 0, 0]
    counter = 0
    minarray = []
    while(counter != 100):
        print("_________________________________________________________________")
        sys.stdout.write("COUNTER: " + str(counter)+"\n")
        nSch = NDEAP(START_SIZE, MAX_EPOCHS, MATING_PROBABILITY, MUTATION_RATE, OFFSPRING_PER_GENERATION,
                     PBC_MAX, MAX_LENGTH)

        nSch.initialize_chromosomes()
        nSch.genetic_algorithm()
        #zeige pro Durchlauf permutation-Anzahl und permutationen pro Epoche
        OverallPlots.show_permutation_amount(nSch)
        OverallPlots.show_crossover_per_epoche(nSch)
        OverallPlots.show_fitness_per_epoche(nSch)

        array = np.array(array) + Operations.get_best_crossover(nSch)

        foundMinimum += Operations.get_foundMimimum(nSch)
        minarray.append(Operations.get_best(nSch))

        counter += 1

    print(array)
    #TODO String je nach Problem ändern
    OverallPlots.show_overall_permutation_amount(array, "Schwefel")
    #OverallPlots.show_overall_permutation_amount(array, "Himmelblau")
    #OverallPlots.show_overall_permutation_amount(array, "Griewank")
    #OverallPlots.show_overall_permutation_amount(array, "Rastrigin")

    sys.stdout.write("minArray: "+str(minarray)+"\n")

    OverallPlots.show_overall_minima(minarray, "Schwefel")
    #OverallPlots.show_overall_minima(minarray, "Himmelblau")
    #OverallPlots.show_overall_minima(minarray, "Griewank")
    #OverallPlots.show_overall_minima(minarray, "Rastrigin")



    sys.stdout.write("Gefundene Minima: "+str(foundMinimum)+"\n")
    sys.stdout.write("Median aller Minima: " + str(np.median(minarray)) + "\n")




