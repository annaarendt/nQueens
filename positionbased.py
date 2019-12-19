import math
import random
import sys
import matplotlib.pyplot as plt
import numpy as np
from Chromosome import Chromosome

START_SIZE = 2  # Population size at start.
MAX_EPOCHS = 0  # Arbitrary number of test cycles. EIGENTLICH AUF 1000 GESETZT!
OFFSPRING_PER_GENERATION = 1  # New offspring created per generation. Range: 0 < OFFSPRING_PER_GENERATION < MAX_SELECT.
MINIMUM_SHUFFLES = 8  # For randomizing starting chromosomes
MAXIMUM_SHUFFLES = 20
PBC_MAX = 4  # Maximum Position-Based Crossover points. Range: 0 < PBC_MAX < 8 (> 8 isn't good).

MAX_LENGTH = 6  # chess board width.


class NQueen1:
    def __init__(self, startSize, maxEpochs, generation, minShuffles,
                 maxShuffles, pbcMax, maxLength):
        self.mStartSize = startSize
        self.mEpochs = maxEpochs
        self.mOffspringPerGeneration = generation
        self.mMinimumShuffles = minShuffles
        self.mMaximumShuffles = maxShuffles
        self.mPBCMax = pbcMax
        self.mMaxLength = maxLength

        self.epoch = 0
        self.childCount = 0
        self.population = []
        #1. und 2. stelle für Konflikte der Eltern, letzen beiden der Kinder
        self.conflict_array = [0] * 4

        return

    # gibt zufällige Zahl zurück mit obergrenze(high) -> untere Grenze ist 0, darf aber nicht Wert von numberA haben.
    def get_exclusive_random_integer(self, high, numberA):
        done = False
        numberB = 0

        while not done:
            numberB = random.randrange(0, high)
            if numberB != numberA:
                done = True

        return numberB

    # zufällige zahl, die aber nicht im übergebenen Array (arrayA) vorkommt --> TODO Was ist mit der 0?
    def get_exclusive_random_integer_by_array(self, low, high, arrayA):
        done = False
        getRand = 0

        if high != low:
            while not done:
                done = True
                getRand = random.randrange(low, high)
                for i in range(len(arrayA)):
                    if getRand == arrayA[i]:
                        done = False
        else:
            getRand = high

        return getRand

    def math_round(self, inValue):
        if math.modf(inValue)[0] >= 0.5:
            outValue = math.ceil(inValue)
        else:
            outValue = math.floor(inValue)
        return outValue

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
                    if math.fabs(thisChromo.get_conflicts() > thatChromo.get_conflicts()):
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
                    if math.fabs(thisChromo.get_conflicts() < thatChromo.get_conflicts()):
                        minimum = i
                        foundNewMinimum = True

            if foundNewMinimum == False:
                done = True

        return minimum

    # Mutation, die einfach nur 2 Gene des Chromosoms austauscht
    def exchange_mutation(self, index, exchanges):
        i = 0
        done = False

        thisChromo = self.population[index]

        while not done:
            gene1 = random.randrange(0, self.mMaxLength)
            gene2 = self.get_exclusive_random_integer(self.mMaxLength, gene1)

            # Exchange the chosen genes.
            tempData = thisChromo.get_data(gene1)
            thisChromo.set_data(gene1, thisChromo.get_data(gene2))
            thisChromo.set_data(gene2, tempData)

            if i == exchanges:
                done = True

            i += 1
        return

    def initialize_chromosomes(self):
        for i in range(self.mStartSize):
            newChromo = Chromosome(self.mMaxLength, 0)
            self.population.append(newChromo)
            chromoIndex = len(self.population) - 1

            # Randomly choose the number of shuffles to perform.
            shuffles = random.randrange(self.mMinimumShuffles, self.mMaximumShuffles)

            self.exchange_mutation(chromoIndex, shuffles)

            newChromo = self.population[chromoIndex]
            newChromo.compute_conflicts()

            # Konsole:
            newChromo.toStr()
            sys.stdout.write("erstelltes Chromosom (Parent) mit ")
            sys.stdout.write(str(newChromo.get_conflicts()) + " Konflikten\n")
        return

    def get_fitness(self):
        # Lowest errors = 100%, Highest errors = 0%
        popSize = len(self.population)

        # The worst score would be the one with the highest energy, best would be lowest.
        thisChromo = self.population[self.get_maximum()]
        worstScore = thisChromo.get_conflicts()

        # Convert to a weighted percentage.
        thisChromo = self.population[self.get_minimum()]
        bestScore = worstScore - thisChromo.get_conflicts()
        sys.stdout.write("bestScore = worstScore - thisChromo.get_conflicts(): "+str(bestScore)+" = "+str(worstScore)+" - "+str(thisChromo.get_conflicts())+"\n")

        for i in range(popSize):
            thisChromo = self.population[i]
            if(bestScore!=0):
                thisChromo.set_fitness((worstScore - thisChromo.get_conflicts()) * 100.0 / bestScore)
            else:
                thisChromo.set_fitness((worstScore - thisChromo.get_conflicts()) * 0)


        return

    def position_based_crossover(self, chromA, chromB, child1, child2):
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
            crossPoints[i] = self.get_exclusive_random_integer_by_array(0, self.mMaxLength - 1, crossPoints)
        # Get non-chosens from parent 2: Die Zahlen die bei P2 nicht an den ausgewählten Stellen bei P1 stehen werden in
        # einem Array gesammelt (tempArray1)
        k = 0
        for i in range(self.mMaxLength):
            matchFound = False
            for j in range(numPoints):
                if thatChromo.get_data(i) == thisChromo.get_data(crossPoints[j]):
                    matchFound = True

            if matchFound == False:
                tempArray1[k] = thatChromo.get_data(i)
                k += 1
        # Insert chosens into child 1: in Child 1 werden die Zahlen von P1 an gewählter Position gesetzt, Rest
        # freigelassen
        for i in range(numPoints):
            newChromo1.set_data(crossPoints[i], thisChromo.get_data(crossPoints[i]))

        # Fill in non-chosens to child 1.
        k = 0
        for i in range(self.mMaxLength):
            matchFound = False
            for j in range(numPoints):
                if i == crossPoints[j]:
                    matchFound = True

            if matchFound == False:
                newChromo1.set_data(i, tempArray1[k])
                k += 1
        sys.stdout.write(str(crossPoints) + " CrossPoints\n")
        sys.stdout.write("Parent1")
        thisChromo.toStr()
        sys.stdout.write("Parent2")
        thatChromo.toStr()
        sys.stdout.write("Kind1  ")
        newChromo1.toStr()
        # Get non-chosens from parent 1
        k = 0
        for i in range(self.mMaxLength):
            matchFound = False
            for j in range(numPoints):
                if thisChromo.get_data(i) == thatChromo.get_data(crossPoints[j]):
                    matchFound = True

            if matchFound == False:
                tempArray2[k] = thisChromo.get_data(i)
                k += 1

        # Insert chosens into child 2.
        for i in range(numPoints):
            newChromo2.set_data(crossPoints[i], thatChromo.get_data(crossPoints[i]))

        # Fill in non-chosens to child 2.
        k = 0
        for i in range(self.mMaxLength):
            matchFound = False
            for j in range(numPoints):
                if i == crossPoints[j]:
                    matchFound = True

            if matchFound == False:
                newChromo2.set_data(i, tempArray2[k])
                k += 1

        sys.stdout.write("Position-based Crossover verwendet.\n")
        return

    # Verschiebungsmutatation: 2 Punkte und die Gene dazwischen werden im Chromosom verschoben
    def displacement_mutation(self, index):
        length = 0
        tempArray1 = [0] * self.mMaxLength
        tempArray2 = [0] * self.mMaxLength
        thisChromo = self.population[index]

        # Randomly choose a section to be displaced.
        point1 = random.randrange(0, self.mMaxLength)
        # sys.stdout.write(str(point1))
        # sys.stdout.write(str(self.mMaxLength))

        # Generate re-insertion point.
        candidate = self.mMaxLength - (point1 + 2)
        if candidate <= 0:
            candidate = 1
        point2 = self.get_exclusive_random_integer(candidate, point1)

        j = 0
        for i in range(self.mMaxLength):  # Get non-chosen
            if i < point1 or i > point1 + length:
                tempArray1[j] = thisChromo.get_data(i)
                j += 1

        j = 0
        for i in range(point1, point1 + length + 1):  # Get chosen
            tempArray2[j] = thisChromo.get_data(i)
            j += 1

        j = 0
        for i in range(point2, point2 + length + 1):  # Place chosen
            thisChromo.set_data(i, tempArray2[j])
            j += 1

        j = 0
        for i in range(i, self.mMaxLength):  # Place non-chosen
            if i < point2 or i > point2 + length:
                thisChromo.set_data(i, tempArray1[j])
                j += 1

        self.mutations += 1
        sys.stdout.write("Displacement Mutation verwendet.\n")
        return

    def do_mating(self):

        for i in range(self.mOffspringPerGeneration):
            chromA = self.population[0]
            chromB = self.population[1]

            sys.stdout.write(str(chromA.get_fitness()) + " :Fitness ChromA\n")
            sys.stdout.write(str(chromB.get_fitness()) + " :Fitness ChromB\n")

            newChromo1 = Chromosome(self.mMaxLength, 0)
            newChromo2 = Chromosome(self.mMaxLength, 0)
            self.population.append(newChromo1)
            newIndex1 = len(self.population) - 1
            self.population.append(newChromo2)
            newIndex2 = len(self.population) - 1
            self.position_based_crossover(chromA, chromB, newIndex1, newIndex2)

            newChromo1 = self.population[newIndex1]
            newChromo1.compute_conflicts()
            # Konsole:
            newChromo1.toStr()
            sys.stdout.write("Kind1 mit ")
            sys.stdout.write(str(newChromo1.get_conflicts()) + " Konflikten\n"
                             )
            newChromo2 = self.population[newIndex2]
            newChromo2.compute_conflicts()
            # Konsole:
            newChromo2.toStr()
            sys.stdout.write("Kind2 mit ")
            sys.stdout.write(str(newChromo2.get_conflicts()) + " Konflikten\n")

            self.childCount += 2
            sys.stdout.write(str(len(self.population)) + " Länge der Population\n")

            self.set_conflict_array(chromA.get_conflicts(), chromB.get_conflicts(),
             newChromo1.get_conflicts(), newChromo2.get_conflicts())
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
            thisChromo.set_selected(False)

        return

    def print_best_solution(self, bestSolution=Chromosome(MAX_LENGTH, type)):
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

        while not done:
            popSize = len(self.population)
            for i in range(popSize):
                thisChromo = self.population[i]
                if thisChromo.get_conflicts() == 0 or self.epoch == self.mEpochs:
                    done = True

            self.get_fitness()

            self.do_mating()

            self.prep_next_epoch()

            self.epoch += 1

            # This is here simply to show the runtime status.
            sys.stdout.write("Epoche: " + str(self.epoch) + "\n")

        sys.stdout.write("done.\n")

        if self.epoch != self.mEpochs:
            popSize = len(self.population)
            for i in range(popSize):
                thisChromo = self.population[i]
                if thisChromo.get_conflicts() == 0:
                    sys.stdout.write(str(thisChromo.toStr()) + " hat " + str(thisChromo.get_conflicts()) + " Konflikte.\n")
                    self.print_best_solution(thisChromo)

        sys.stdout.write("Completed " + str(self.epoch) + " epochs.\n")
        sys.stdout.write("Done: "+str(done)+"\n")
        return

    #Konflikte der Eltern und Kinder werden in einem Plot gezeigt
    def show_conflicts(self, chromP1, chromP2, chromK1, chromK2 ):
        chromosomes=["P1", "P2", "K1", "K2"]
        conflicts=[chromP1.get_conflicts(), chromP2.get_conflicts(),chromK1.get_conflicts(),
                            chromK2.get_conflicts()]

        plt.scatter(chromosomes, conflicts, c='b', marker='o')
        plt.xlabel('Chromosomen', fontsize=16)
        plt.ylabel('Konflikte', fontsize=16)
        plt.title('scatter plot - parents vs offsprings', fontsize=20)
        plt.show()

def show_overall_conflicts(p1, p2, k1, k2):
    x = np.arange(0, len(p1), 1)
    y1 = p1
    y2 = p2
    y3 = k1
    y4 = k2
    lines = plt.plot(x, y1, x, y2, x, y3, x, y4)
    l1, l2, l3, l4 = lines
    plt.setp(l1, color='lightcoral', alpha=0.5)  # line1 is red -> parent1
    plt.setp(l2, color='red', alpha=0.5)  # line2 orange -> parent2
    plt.setp(l3, color='forestgreen', alpha=0.5)  # line3 is green -> offspring1
    plt.setp(l4, color='limegreen', alpha=0.5)  # line4 is darkgreen -> offspring2

    plt.ylabel("Konflikte")
    plt.xlabel("Durchläufe")
    plt.title('parents and offsprings')
    plt.legend((l1, l2, l3, l4), ('parent1', 'parent2', 'child1', 'child2'))
    plt.show()
    return

def boxplot(p1, p2, k1, k2):
    data_to_plot = [p1, p2, k1, k2]

    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    plt.ylabel("Konflikte")
    plt.title('parents and offsprings')
    ax.boxplot(data_to_plot)
    ax.set_xticklabels(['Parent1', 'Parent2', 'Child1', 'Child2'])
    plt.show()
    return

def percentage_table(p1, p2, k1, k2):

    len_array = len(p1)

    #Anzahl, wie oft Kinder weniger Konflikte ggü. der Eltern haben
    vgl_k1_p1 = 0
    vgl_k1_p2 = 0
    vgl_k2_p1 = 0
    vgl_k2_p2 = 0
    vgl_k1_beide = 0
    vgl_k2_beide = 0

    for i in range(len(p1)):
        #Kind1 weniger Konflikte als Elternteil1 und/oder Elternteil2
        if k1[i] < p1[i]:
            vgl_k1_p1 += 1
            if k1[i] < p2[i]:
                vgl_k1_beide += 1
        elif k1[i] < p2[i]:
            vgl_k1_p2 += 1

        #Kind2 weniger Konflikte als Elternteil1 und/oder Elternteil2
        if k2[i] < p1[i]:
            vgl_k2_p1 += 1
            if k2[i] < p2[i]:
                vgl_k2_beide += 1
        elif k2[i] < p2[i]:
            vgl_k2_p2 += 1

    #Berechnung der Prozente
    k1_vs_p1 = vgl_k1_p1/len_array
    k1_vs_p2 = vgl_k1_p2 / len_array
    k2_vs_p1 = vgl_k2_p1 / len_array
    k2_vs_p2 = vgl_k2_p2 / len_array
    k1_vs_beide = vgl_k1_beide / len_array
    k2_vs_beide = vgl_k2_beide / len_array

    fig = plt.figure(dpi=80)
    title="Durchgänge = "+str(len_array)
    ax = fig.add_subplot(1, 1, 1)
    table_data = [
        ["Kind1 besser als Parent1", str(k1_vs_p1)+"%", vgl_k1_p1 ],
        ["Kind1 besser als Parent2", str(k1_vs_p2)+"%", vgl_k1_p2],
        ["Kind2 besser als Parent1", str(k2_vs_p1)+"%", vgl_k2_p1],
        ["Kind2 besser als Parent2", str(k2_vs_p2)+"%", vgl_k2_p2],
        ["Kind1 besser als beide", str(k1_vs_beide)+"%", vgl_k1_beide],
        ["Kind2 besser als beide", str(k2_vs_beide)+"%", vgl_k2_beide]]
    table = ax.table(cellText=table_data, loc='center', colLabels=[title,"Prozentsatz", "Anzahl"])
    table.set_fontsize(14)
    table.scale(1, 4)
    ax.axis('off')
    plt.show()

    return

if __name__ == '__main__':
    COUNTER = 0
    END = 5000
    p1 = [0] * END
    p2 = [0] * END
    k1 = [0] * END
    k2 = [0] * END
    while (COUNTER != END):
        nq1 = NQueen1(START_SIZE, MAX_EPOCHS, OFFSPRING_PER_GENERATION, MINIMUM_SHUFFLES, MAXIMUM_SHUFFLES,
                      PBC_MAX, MAX_LENGTH)

        nq1.initialize_chromosomes()
        nq1.genetic_algorithm()

        p1[COUNTER]= nq1.get_conflict_array()[0]
        p2[COUNTER] = nq1.get_conflict_array()[1]
        k1[COUNTER] = nq1.get_conflict_array()[2]
        k2[COUNTER] = nq1.get_conflict_array()[3]
        #konsole_konflikte je Durchlauf
        sys.stdout.write("COUNTER: " + str(COUNTER)+"\n")

        str1 = ','.join(str(e) for e in nq1.get_conflict_array())
        print(str1)

        COUNTER += 1


    show_overall_conflicts(p1, p2, k1, k2)
    boxplot(p1, p2, k1, k2)
    percentage_table(p1, p2, k1, k2)








