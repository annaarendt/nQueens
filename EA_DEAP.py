import math
import random
import sys
import matplotlib.pyplot as plt
import numpy as np
from Chromosome import Chromosome
from ChromosomeBenchmark import ChromosomeDEAP
from EA import NQueen1
from Overall_plots import OverallPlots

START_SIZE = 75  # Population size at start.
MAX_EPOCHS = 1000  # Arbitrary number of test cycles. Default 1000!
MATING_PROBABILITY = 0.7  # Probability of two chromosomes mating. Range: 0.0 < MATING_PROBABILITY < 1.0
MUTATION_RATE = 0.001  # Mutation Rate. Range: 0.0 < MUTATION_RATE < 1.0
#pro mating aber 2 Kinder -> offspring_per_epoche=1 erstellt 2 Kinder von 2 Eltern. wert auf 20 entspricht 40 kinder von
#40 (nicht unbedignt unterschiedlichn) eltern
OFFSPRING_PER_GENERATION = 20  # New offspring created per generation. Range: 0 < OFFSPRING_PER_GENERATION < MAX_SELECT.0
PBC_MAX = 4  # Maximum Position-Based Crossover points. Range: 0 < PBC_MAX < 8 (> 8 isn't good).

MAX_LENGTH = 8 # length of chromosom.


#welcher crossover für welche phase besser
#mutaion von corssover und crossover fixen, auf viele durchgänge testen
#eine schlechte rekombinationsmethode implementieren zum testen

#eltern nicht ganzwegschmeiße. kinder erezugen random , zsm schmeisen und dann die besten kinder weiternehmen.
#vlt. nur ein kind? crossover dass das beste weitergegeben werden


class NDEAP:
    def __init__(self, startSize, maxEpochs, matingProb, mutationRate, generation, pbcMax, maxLength):
        self.mStartSize = startSize
        self.mEpochs = maxEpochs
        self.mMatingProbability = matingProb
        self.mMutationRate = mutationRate
        self.mOffspringPerGeneration = generation
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
            rand = random.randrange(0, 3)

            newChromo = ChromosomeDEAP(self.mMaxLength, rand)
            self.population.append(newChromo)
            chromoIndex = len(self.population) - 1

            newChromo = self.population[chromoIndex]
            newChromo.compute_fitness()

        return

    def bad_recombination(self, chromA, chromB, child1, child2):
            tempArray1 = [0] * self.mMaxLength
            tempArray2 = [0] * self.mMaxLength
            thisChromo = chromA
            thatChromo = chromB
            newChromo1 = self.population[child1]
            newChromo2 = self.population[child2]

            # Choose and sort the crosspoints.
            numPoints = random.randrange(0, self.mPBCMax)  # if PBC_MAX is set any higher than 6 or 8.
            crossPoints = [0] * numPoints
            for i in range(numPoints):
                crossPoints[i] = NQueen1.get_exclusive_random_integer_by_array(self, 0, self.mMaxLength - 1,
                                                                               crossPoints)
            # Get non-chosens from parent 2: Die Zahlen die bei P2 nicht an den ausgewählten Stellen bei P1 stehen werden in
            # einem Array gesammelt (tempArray1)
            k = 0
            for i in range(self.mMaxLength):
                matchFound = False
                for j in range(numPoints):
                    if Chromosome.get_data(thatChromo, i) == Chromosome.get_data(thisChromo, crossPoints[j]):
                        matchFound = True

                if matchFound == False:
                    tempArray1[k] = Chromosome.get_data(thatChromo, i)
                    k += 1
            # Insert chosens into child 1: in Child 1 werden die Zahlen von P1 an gewählter Position gesetzt, Rest
            # freigelassen
            for i in range(numPoints):
                Chromosome.set_data(newChromo1, crossPoints[i], Chromosome.get_data(thisChromo, crossPoints[i]))

            # Fill in non-chosens to child 1.
            k = 0
            for i in range(self.mMaxLength):
                matchFound = False
                for j in range(numPoints):
                    if i == crossPoints[j]:
                        matchFound = True

                if matchFound == False:
                    Chromosome.set_data(newChromo1, i, tempArray1[k])
                    k += 1
            newChromo1.compute_fitness()

            # Kind1 wird rausgehauen und index child2 um eins heragesetzt weil sich population auch wieder ändert
            if Chromosome.get_fitness(newChromo1) < Chromosome.get_fitness(thisChromo) or \
                    Chromosome.get_fitness(newChromo1) < Chromosome.get_fitness(thatChromo):
                sys.stdout.write("Kind1 zu wenig Konflikte: " + str(Chromosome.get_fitness(newChromo2)) + "\n")
                Chromosome.set_fitness(newChromo1, 2000)
                sys.stdout.write("Kind1 wurde auf Fitness 2000 erhöht\n")
            #
            # Get non-chosens from parent 1
            k = 0
            for i in range(self.mMaxLength):
                matchFound = False
                for j in range(numPoints):
                    if Chromosome.get_data(thisChromo, i) == Chromosome.get_data(thatChromo, crossPoints[j]):
                        matchFound = True

                if matchFound == False:
                    tempArray2[k] = Chromosome.get_data(thisChromo, i)
                    k += 1

            # Insert chosens into child 2.
            for i in range(numPoints):
                Chromosome.set_data(newChromo2, crossPoints[i], Chromosome.get_data(thatChromo, crossPoints[i]))

            # Fill in non-chosens to child 2.
            k = 0
            for i in range(self.mMaxLength):
                matchFound = False
                for j in range(numPoints):
                    if i == crossPoints[j]:
                        matchFound = True

                if matchFound == False:
                    Chromosome.set_data(newChromo2, i, tempArray2[k])
                    k += 1

            newChromo2.compute_fitness()

            if Chromosome.get_fitness(newChromo2) < Chromosome.get_fitness(thisChromo) \
                    or Chromosome.get_fitness(newChromo2) < Chromosome.get_fitness(thatChromo):
                sys.stdout.write("Kind2 zu wenig Konflikte: " + str(Chromosome.get_fitness(newChromo2)) + "\n")
                Chromosome.set_fitness(newChromo2, 2000)
                sys.stdout.write("Kind2 wurde auf Fitness 2000 erhöht\n")

            sys.stdout.write(str(crossPoints) + " CrossPoints\n")

            sys.stdout.write("BAD_recombination verwendet.\n")
            return



    #ersetzt eine gene durch eine zufällig Anderen [-500,500,1]
    def new_gene_mutation(self, index):

        thisChromo = self.population[index]
        point1 = random.randrange(0, self.mMaxLength)

        sys.stdout.write("Punkt: "+str(point1)+", Vorher: " + str(Chromosome.toStr(thisChromo) + "\n"))

        rand= random.randrange(-500,500,1)
        Chromosome.set_data(thisChromo, point1, rand)

        sys.stdout.write("New_Gene_Mutation verwendet.\n")
        sys.stdout.write("Nachher: "+str(Chromosome.toStr(thisChromo)+"\n"))
        self.mutations += 1

        return


    #wenn die eltern verschiedene crossover haben wird das übernommen von dem elternteil bessere fitness hat
    def do_mating(self):
        for i in range(self.mOffspringPerGeneration):
            parentA = NQueen1.choose_first_parent(self)
            # Test probability of mating.
            getRand = random.randrange(0, 100)

            if getRand <= self.mMatingProbability * 100:
                parentB = NQueen1.choose_second_parent(self, parentA)
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
                    NQueen1.partially_mapped_crossover(self, chromA, chromB, newIndex1, newIndex2)
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
                    self.bad_recombination(chromA, chromB, newIndex1, newIndex2)
                    #NQueen1.position_based_crossover(self, chromA, chromB, newIndex1, newIndex2)
                    self.position_based_co += 1
                    self.current_p_b += 1
                else:
                    newChromo1 = ChromosomeDEAP(self.mMaxLength, 2)
                    newChromo2 = ChromosomeDEAP(self.mMaxLength, 2)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    NQueen1.order_based_crossover(self, chromA, chromB, newIndex1, newIndex2)
                    self.order_based_co += 1
                    self.current_o_b += 1


                if self.childCount - 1 == self.nextMutation:

                    self.new_gene_mutation(newIndex1)
                elif self.childCount == self.nextMutation:
                    self.new_gene_mutation(newIndex2)


                newChromo1 = self.population[newIndex1]
                sys.stdout.write("Kind1 mit Fitness " + str(Chromosome.get_fitness(newChromo1))+"\n")
                self.childCount += 1

                newChromo2 = self.population[newIndex2]
                sys.stdout.write("Kind2 mit Fitness " + str(Chromosome.get_fitness(newChromo2))+"\n")
                self.childCount += 1

                # Schedule next mutation.
                if math.fmod(self.childCount, NQueen1.math_round(self, 1.0 / self.mMutationRate)) == 0:
                    self.nextMutation = self.childCount + random.randrange(0, NQueen1.math_round(self, 1.0 / self.mMutationRate))

        return

    def genetic_algorithm(self):

        done = False

        self.mutations = 0
        self.nextMutation = random.randrange(0, NQueen1.math_round(self, 1.0 / self.mMutationRate))

        while not done:
            popSize = len(self.population)
            for i in range(popSize):
                thisChromo = self.population[i]
                if Chromosome.get_fitness(thisChromo) == 0 or self.epoch == self.mEpochs:
                    done = True

            NQueen1.get_fitness(self)

            NQueen1.roulette_selection(self)

            self.do_mating()

            NQueen1.prep_next_epoch(self)

            self.array_p_m.append(self.current_p_m)
            self.array_p_b.append(self.current_p_b)
            self.array_o_b.append(self.current_o_b)
            self.current_p_m = 0
            self.current_p_b = 0
            self.current_o_b = 0

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
                    NQueen1.print_best_solution(self, thisChromo)
                    break

        sys.stdout.write("Completed " + str(self.epoch) + " epochs.\n")
        sys.stdout.write(
            "Encountered " + str(self.mutations) + " mutations in " + str(self.childCount) + " offspring.\n")
        sys.stdout.write(
            "Encountered " + str(self.partielle_mapped_co) + " partiell-mapped , " +
            str(self.position_based_co) + " positioon-based and " + str(self.order_based_co)
            + " order-based Crossover.\n")

        return

    def show_fitness_per_epoche(self):

        popSize = len(self.population)
        fitness_array = [0]

        for i in range(popSize):
            thisChromo = self.population[i]
            fitness_array.append(Chromosome.get_fitness(thisChromo))
        array=np.asarray(self.array_fitness)
        x = np.arange(1, self.epoch, 1)
        y1 = array[x]
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt.xlabel('Epochen')
        plt.ylabel('Best Fitness/Population')
        ax.plot(x, y1, color='orange', label ='partielle_mapped')
        plt.show()
        #plt.savefig('crossovers_per_epoche.png')

        return

if __name__ == '__main__':
    array = [0, 0, 0]
    counter = 0
    while(counter != 10):
        print("_________________________________________________________________")
        sys.stdout.write("COUNTER: " + str(counter)+"\n")
        nSch = NDEAP(START_SIZE, MAX_EPOCHS, MATING_PROBABILITY, MUTATION_RATE, OFFSPRING_PER_GENERATION,
                     PBC_MAX, MAX_LENGTH)

        nSch.initialize_chromosomes()
        nSch.genetic_algorithm()
        #zeige pro Durchlauf permutation-Anzahl und permutationen pro Epoche
        NQueen1.show_permutation_amount(nSch)
        NQueen1.show_crossover_per_epoche(nSch)
        nSch.show_fitness_per_epoche()

        array = np.array(array) + NQueen1.get_best_crossover(nSch)
        counter += 1

    print(array)
    OverallPlots.show_overall_permutation_amount(array)



