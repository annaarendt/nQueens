import math
import random
import sys
import matplotlib.pyplot as plt
import numpy as np
from Chromosome import Chromosome

START_SIZE = 75  # Population size at start.
MAX_EPOCHS = 40  # Arbitrary number of test cycles. Default 1000!
MATING_PROBABILITY = 0.7  # Probability of two chromosomes mating. Range: 0.0 < MATING_PROBABILITY < 1.0
MUTATION_RATE = 0.001  # Mutation Rate. Range: 0.0 < MUTATION_RATE < 1.0
MIN_SELECT = 10  # Minimum parents allowed for selection.
MAX_SELECT = 50  # Maximum parents allowed for selection. Range: MIN_SELECT < MAX_SELECT < START_SIZE
#pro mating aber 2 Kinder -> offspring_per_epoche=1 erstellt 2 Kinder von 2 Eltern. wert auf 20 entspricht 40 kinder von
#40 (nicht unbedignt unterschiedlichn) eltern
OFFSPRING_PER_GENERATION = 20  # New offspring created per generation. Range: 0 < OFFSPRING_PER_GENERATION < MAX_SELECT.
MINIMUM_SHUFFLES = 8  # For randomizing starting chromosomes
MAXIMUM_SHUFFLES = 20
PBC_MAX = 4  # Maximum Position-Based Crossover points. Range: 0 < PBC_MAX < 8 (> 8 isn't good).

MAX_LENGTH = 8 # chess board width.


#welcher crossover für welche phase besser
#mutaion von corssover und crossover fixen, auf viele durchgänge testen
#eine schlechte rekombinationsmethode implementieren zum testen

#eltern nicht ganzwegschmeiße. kinder erezugen random , zsm schmeisen und dann die besten kinder weiternehmen.
#vlt. nur ein kind? crossover dass das beste weitergegeben werden


class NQueen1:
    def __init__(self, startSize, maxEpochs, matingProb, mutationRate, minSelect, maxSelect, generation, minShuffles,
                 maxShuffles, pbcMax, maxLength):
        self.mStartSize = startSize
        self.mEpochs = maxEpochs
        self.mMatingProbability = matingProb
        self.mMutationRate = mutationRate
        self.mMinSelect = minSelect
        self.mMaxSelect = maxSelect
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
        return

    #gibt zufällige Zahl zurück mit obergrenze(high) -> untere Grenze ist 0, darf aber nicht Wert von numberA haben.
    def get_exclusive_random_integer(self, high, numberA):
        done = False
        numberB = 0

        while not done:
            numberB = random.randrange(0, high)
            if numberB != numberA:
                done = True

        return numberB

    #zufällige zahl, die aber nicht im übergebenen Array (arrayA) vorkommt --> TODO Was ist mit der 0?
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

    #Mutation, die einfach nur 2 Gene des Chromosoms austauscht
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

        self.mutations += 1
        return

    def initialize_chromosomes(self):
        for i in range(self.mStartSize):
            #crossover des neuen Chromosoms wird bestimmt
            rand = random.randrange(0, 3)

            newChromo = Chromosome(self.mMaxLength, rand)
            self.population.append(newChromo)
            chromoIndex = len(self.population) - 1

            # Randomly choose the number of shuffles to perform.
            shuffles = random.randrange(self.mMinimumShuffles, self.mMaximumShuffles)

            self.exchange_mutation(chromoIndex, shuffles)

            newChromo = self.population[chromoIndex]
            newChromo.compute_conflicts()

        return

    def get_fitness(self):
        """""
        # Lowest errors = 100%, Highest errors = 0%
        popSize = len(self.population)

        # The worst score would be the one with the highest energy, best would be lowest.
        thisChromo = self.population[self.get_maximum()]
        sys.stdout.write(str(thisChromo.get_conflicts()) + " Maximale Konflikte\n")
        worstScore = thisChromo.get_conflicts()

        # Convert to a weighted percentage.
        thisChromo = self.population[self.get_minimum()]
        sys.stdout.write(str(thisChromo.get_conflicts()) + " Minimale Konflikte\n")
        #bestScore ist die bandbreite der Konflikte. Maximale minus minimale Anzahl an Konflikten
        bestScore = worstScore - thisChromo.get_conflicts()
        sys.stdout.write(str(bestScore) + " bestscore\n")

        #die wenigsten Konflikte: Fitness ist 100, die meisten: Fitness ist 0
        for i in range(popSize):
            thisChromo = self.population[i]
            thisChromo.set_fitness((worstScore - thisChromo.get_conflicts()) * 100.0 / bestScore)
        """""

        #fitness entspricht den Konflikten. Niedrige Fitness ist gut, hohe ist schlecht
        popSize = len(self.population)

        for i in range(popSize):
            thisChromo = self.population[i]
            worst = self.population[self.get_maximum()]
            thisChromo.set_fitness(worst.get_conflicts()-thisChromo.get_conflicts())

        sys.stdout.write(str(worst.get_conflicts()) + " Maximale Konflikte\n")


        return

    #Selektion der Eltern mit Roulette-Methode. Je besser (kleiner) Fitness, desto mehr Anteil auf dem Roulette
    def roulette_selection(self):
        genTotal = 0.0
        sumProp = 0.0
        popSize = len(self.population)

        for i in range(popSize):
            thisChromo = self.population[i]
            genTotal += thisChromo.get_fitness()

        sys.stdout.write(str(genTotal) + " GenTotal(Roulette)\n")
        sys.stdout.write(str(popSize) + " popSize\n")

        for i in range(popSize):
            thisChromo = self.population[i]
            sys.stdout.write(str(thisChromo.get_fitness()) + " fitness\n")
            propability = (thisChromo.get_fitness() / genTotal)
            sumProp += propability
            thisChromo.set_selection_probability(propability)
            sys.stdout.write(str(propability) + " propability von chromosom\n")

        sys.stdout.write(str(sumProp) + " alle Props\n")

        for i in range(self.mOffspringPerGeneration):
            rouletteSpin = random.uniform(0, 1)
            j = 0
            selTotal = 0
            done = False
            while not done:
                thisChromo = self.population[j]
                selTotal += thisChromo.get_selection_probability()
                if selTotal >= rouletteSpin:
                    if j == 0:
                        thatChromo = self.population[j]
                    elif j >= popSize - 1:
                        thatChromo = self.population[popSize - 1]
                    else:
                        thatChromo = self.population[j - 1]

                    thatChromo.set_selected(True)
                    done = True
                else:
                    j += 1
        return

    def choose_first_parent(self):
        parent = 0
        done = False

        while not done:
            # Randomly choose an eligible parent.
            parent = random.randrange(0, len(self.population) - 1)
            thisChromo = self.population[parent]
            if thisChromo.get_selected() == True:
                done = True

        return parent

    def choose_second_parent(self, parentA):
        parentB = 0
        done = False

        while not done:
            # Randomly choose an eligible parent.
            parentB = random.randrange(0, len(self.population) - 1)
            if parentB != parentA:
                thisChromo = self.population[parentB]
                if thisChromo.get_selected() == True:
                    done = True

        return parentB

    def partially_mapped_crossover(self, chromA, chromB, child1, child2):
        thisChromo = chromA
        thatChromo = chromB
        newChromo1 = self.population[child1]
        newChromo2 = self.population[child2]

        crossPoint1 = random.randrange(0, self.mMaxLength)
        crossPoint2 = self.get_exclusive_random_integer(self.mMaxLength, crossPoint1)
        if crossPoint2 < crossPoint1:
            j = crossPoint1
            crossPoint1 = crossPoint2
            crossPoint2 = j

        # Copy Parent genes to offspring.
        for i in range(self.mMaxLength):
            newChromo1.set_data(i, thisChromo.get_data(i))
            newChromo2.set_data(i, thatChromo.get_data(i))

        for i in range(crossPoint1, crossPoint2 + 1):
            # // Get the two items to swap.
            item1 = thisChromo.get_data(i)
            item2 = thatChromo.get_data(i)
            pos1 = 0
            pos2 = 0

            # Get the items' positions in the offspring.
            for j in range(self.mMaxLength):
                if newChromo1.get_data(j) == item1:
                    pos1 = j
                elif newChromo1.get_data(j) == item2:
                    pos2 = j

            # Swap them.
            if item1 != item2:
                newChromo1.set_data(pos1, item2)
                newChromo1.set_data(pos2, item1)

            # Get the items'  positions in the offspring.
            for j in range(self.mMaxLength):
                if newChromo2.get_data(j) == item2:
                    pos1 = j
                elif newChromo2.get_data(j) == item1:
                    pos2 = j

            # Swap them.
            if item1 != item2:
                newChromo2.set_data(pos1, item1)
                newChromo2.set_data(pos2, item2)

        sys.stdout.write(str(crossPoint1) + " CrossPoints\n")
        sys.stdout.write(str(crossPoint2) + " CrossPoints\n")
        sys.stdout.write("Parent1")
        thisChromo.toStr()
        sys.stdout.write("Parent2")
        sys.stdout.write(" hat "+str(thatChromo.get_conflicts())+" Konflikte\n")
        thatChromo.toStr()

        sys.stdout.write("Partially-maped Crossover verwendet.\nMit crossovertyp: "+str(thisChromo.get_crossover())+
                         ", "+str(thatChromo.get_crossover())+"\nUnd fitness "+str(thisChromo.get_fitness())+
                         ", "+str(thatChromo.get_fitness())+"\nUnd crossover der neuen: "
                         +str(newChromo1.get_crossover())+", " + str(newChromo2.get_crossover()) + "\n")

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

    # TODO methode klappt nicht ganz, es werden mehrere Boards angezeigt die nichtmal Bedingungen der Lsg entsprechen
    def order_based_crossover(self, chromA, chromB, child1, child2):

        thisChromo = chromA
        thatChromo = chromB
        newChromo1 = self.population[child1]
        newChromo2 = self.population[child2]



        # Choose and sort the crossNumbers
        # points sind die positionen in parent 1 -> dann werden die Zahlen herausgefunden und die Positionen im
        # 2. Parent gesucht
        numPoints = random.randrange(0, self.mPBCMax)  # if PBC_MAX is set any higher than 6 or 8.
        points = [0] * numPoints
        for i in range(numPoints):
            points[i] = self.get_exclusive_random_integer_by_array(0, self.mMaxLength - 1, points)
        points.sort()

        crossNumbers = self.findCrossNumbers(numPoints, points, thisChromo)
        crossNumbers2 = self.findCrossNumbers(numPoints, points, thatChromo)


        # von Parent2 werden Positionen der Zahlen von Parent1 (in crossNumbers) gesucht und in Array cpoints
        # gespeichert.
        # Von den Zahlen die bei P2  den ausgewählten Zahlen in crossNumbers entsprechen werden die Positionen in
        # einem Array gesammelt (points)
        cpoints = self.findCPoints(numPoints, thatChromo, crossNumbers)
        cpoints2 = self.findCPoints(numPoints, thisChromo, crossNumbers2)



        #temporäres Array erzeugen mit Lücken in Positionen von cpoints/cpoint2, danach Lücken mit crossNumbers/
        # crossNumbers2 füllen

        array1 = self.fullfill(numPoints, cpoints, thatChromo, crossNumbers)
        array2 = self.fullfill(numPoints, cpoints2, thisChromo, crossNumbers2)

        for i in range(self.mMaxLength):
            newChromo1.set_data(i, array1[i])

        for i in range(self.mMaxLength):
            newChromo2.set_data(i, array2[i])

        if not points:
            for i in range(self.mMaxLength):
                newChromo1.set_data(i, thatChromo.get_data(i))
                newChromo2.set_data(i, thisChromo.get_data(i))


        sys.stdout.write(str(points) + " Positionen der ausgewähten Zahlen\n")
        sys.stdout.write(str(crossNumbers) + " CrossNumbers: Zahlen im 1. Parent\n")
        sys.stdout.write(str(crossNumbers2) + " CrossNumbers2: Zahlen im 2. Parent\n")
        sys.stdout.write(str(cpoints) + " cpoints im 2.Parent\n")
        sys.stdout.write(str(cpoints2) + " cpoints im 1.Parent\n")
        sys.stdout.write("Parent1")
        thisChromo.toStr()
        sys.stdout.write("Parent2")
        thatChromo.toStr()
        sys.stdout.write("Kind1  ")
        newChromo1.toStr()
        sys.stdout.write("Kind2  ")
        newChromo2.toStr()
        sys.stdout.write("Order-based Crossover verwendet.\n")

        return

    def findCrossNumbers(self, numPoints, points, chromosome):

        crossNumbers = [0] * numPoints
        for i in range(numPoints):
            for j in range(self.mMaxLength):
                if points[i] == j:
                    crossNumbers[i] = chromosome.get_data(points[i])

        return crossNumbers

    def findCPoints(self, numPoints, chromosome, numbers):

        cpoints = [0] * numPoints

        k = 0
        for i in range(self.mMaxLength):
            matchFound = False
            for j in range(numPoints):
                if chromosome.get_data(i) == numbers[j]:
                    matchFound = True

            if matchFound == True:
                #testen ob cpoints leer ist
                if cpoints :
                    # TODO diese if-Abfrage war nur quickfix
                    if k < numPoints:
                        cpoints[k] = i
                        k += 1

        return cpoints

    def fullfill(self, numPoints, cpoints, chromosome, numbers):

        tempArray = [0] * self.mMaxLength

        for i in range(self.mMaxLength):
            matchFound=False
            for j in range(numPoints):
                if i == cpoints[j]:
                    matchFound=True
                else:
                    tempArray[i]=chromosome.get_data(i)

            if matchFound==True:
                tempArray[i] = 0

        # Lücken in tempArray werden mit Zahlen von numbers aufgefüllt
        k = 0
        # hier wird gecheckt ob crossNumbers leer ist, wenn ja wird nichts gemacht
        #TODO exception bearbeiten
        if numbers:
            for i in range(self.mMaxLength):
                if tempArray[i] == 0:
                    if chromosome.get_data(i) != 0:
                        tempArray[i] = numbers[k]
                        k += 1
                    elif 0 in numbers:
                        try:
                            tempArray[i] = numbers[k]
                            k += 1
                        except IndexError:
                            sys.stdout.write(str(k)+" :k, schiefgelaufen.. i: "+str(i)+ "\n")


        return tempArray


    #Verschiebungsmutatation: 2 Punkte und die Gene dazwischen werden im Chromosom verschoben
    def displacement_mutation(self, index):
        length = 0
        tempArray1 = [0] * self.mMaxLength
        tempArray2 = [0] * self.mMaxLength
        thisChromo = self.population[index]

        # Randomly choose a section to be displaced.
        point1 = random.randrange(0, self.mMaxLength)

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

    #wenn die eltern verschiedene crossover haben wird das übernommen von dem elternteil bessere fitness hat
    def do_mating(self):

        for i in range(self.mOffspringPerGeneration):
            parentA = self.choose_first_parent()
            # Test probability of mating.
            getRand = random.randrange(0, 100)

            if getRand <= self.mMatingProbability * 100:
                parentB = self.choose_second_parent(parentA)

                #das crossover des Elternteils mit größerer Fitness wird gewählt
                chromA = self.population[parentA]
                chromB = self.population[parentB]

                if chromA.get_fitness() > chromB.get_fitness():
                    type = chromA.get_crossover()
                else:
                    type = chromB.get_crossover()

                print(chromA.get_fitness())
                print(chromA.get_crossover())
                print(chromB.get_fitness())
                chromB.toStr()

                if type == 0:
                    newChromo1 = Chromosome(self.mMaxLength, 0)
                    newChromo2 = Chromosome(self.mMaxLength, 0)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    self.partially_mapped_crossover(chromA, chromB, newIndex1, newIndex2)
                    self.partielle_mapped_co += 1
                    self.current_p_m +=1
                elif type == 1:
                    newChromo1 = Chromosome(self.mMaxLength, 1)
                    newChromo2 = Chromosome(self.mMaxLength, 1)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    self.position_based_crossover(chromA, chromB, newIndex1, newIndex2)
                    self.position_based_co += 1
                    self.current_p_b += 1
                else:
                    newChromo1 = Chromosome(self.mMaxLength, 2)
                    newChromo2 = Chromosome(self.mMaxLength, 2)
                    self.population.append(newChromo1)
                    newIndex1 = len(self.population) - 1
                    self.population.append(newChromo2)
                    newIndex2 = len(self.population) - 1
                    self.order_based_crossover(chromA, chromB, newIndex1, newIndex2)
                    self.order_based_co += 1
                    self.current_o_b += 1

                if self.childCount - 1 == self.nextMutation:
                    #self.exchange_mutation(newIndex1, 1)
                    self.displacement_mutation(newIndex1)
                elif self.childCount == self.nextMutation:
                    #self.exchange_mutation(newIndex2, 1)
                    self.displacement_mutation(newIndex2)

                newChromo1 = self.population[newIndex1]
                newChromo1.compute_conflicts()
                newChromo2 = self.population[newIndex2]
                newChromo2.compute_conflicts()

                self.childCount += 2

                # Schedule next mutation.
                if math.fmod(self.childCount, self.math_round(1.0 / self.mMutationRate)) == 0:
                    self.nextMutation = self.childCount + random.randrange(0, self.math_round(1.0 / self.mMutationRate))

        return


    def prep_next_epoch(self):
        # Reset flags for selected individuals.
        popSize = len(self.population)
        for i in range(popSize):
            thisChromo = self.population[i]
            thisChromo.set_selected(False)

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
        self.nextMutation = random.randrange(0, self.math_round(1.0 / self.mMutationRate))

        while not done:
            popSize = len(self.population)
            for i in range(popSize):
                thisChromo = self.population[i]
                if thisChromo.get_conflicts() == 0 or self.epoch == self.mEpochs:
                    done = True

            self.get_fitness()

            self.roulette_selection()

            self.do_mating()

            self.prep_next_epoch()

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
                if thisChromo.get_conflicts() == 0:
                    sys.stdout.write(str(thisChromo.toStr())+" hat " + str(thisChromo.get_conflicts()) + " Konflikte.\n")
                    self.print_best_solution(thisChromo)
                    break

        sys.stdout.write("Completed " + str(self.epoch) + " epochs.\n")
        sys.stdout.write(
            "Encountered " + str(self.mutations) + " mutations in " + str(self.childCount) + " offspring.\n")
        sys.stdout.write(
            "Encountered " + str(self.partielle_mapped_co) + " partiell-mapped , " +
            str(self.position_based_co) + " positioon-based and " + str(self.order_based_co)
            + " order-based Crossover.\n")

        return

    #TODO plt.show() nicht mehr auskommentieren
    def show_permutation_amount(self):
        data = [self.partielle_mapped_co, self.position_based_co, self.order_based_co]
        permutations = ["partielle_mapped", "position_based", "order_based"]
        colors = ['orangered', 'blue', 'yellow']
        plt.bar(permutations, data, color=colors)
        plt.ylabel("Anzahl")
        plt.title("Anzahl der genutzten Crossover")
        plt.show();

        return


    def show_crossover_per_epoche(self):
        pmarray=np.asarray(self.array_p_m)
        pbarray=np.asarray(self.array_p_b)
        obarray=np.asarray(self.array_o_b)
        x = np.arange(1, self.epoch, 1)
        #sys.stdout.write(str(x) +" this is it\n")
        #sys.stdout.write(str(pmarray) + " pmarray\n")
        y1 = pmarray[x]
        y2 = pbarray[x]
        y3 = obarray[x]
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt.xlabel('Epochen')
        plt.ylabel('Anzahl')
        p_m, = ax.plot(x, y1, color='orange', label ='partielle_mapped')
        p_b, = ax.plot(x, y2, color='blue', label ='position_based')
        o_b, = ax.plot(x, y3, color='yellow', label ='order_based')
        ax.legend([p_m, p_b, o_b],["partielle_mapped", "position_based", "order_based"])
        plt.show()

        return

    #methode gibt ein numpy-Array zuruek, abhängig, welches crossover das beste in dem jeweiligen
    #durchgang war. erste stelle im array ist partielle-mapped, dann kommt position-based und dann
    #an dritter stelle order-based
    def get_best_crossover(self):
        pb = self.position_based_co
        ob = self.order_based_co
        pm = self.partielle_mapped_co
        if(pm >= pb and pm >= ob):
            array = [1, 0, 0]
            print("PM IS BEST")
        elif(pb >= pm and pb >= ob):
            array = [0, 1, 0]
            print("PB IS BEST")
        else:
            array = [0, 0, 1]
            print("OB IS BEST")


        return np.array(array)




def show_overall_permutation_amount(array):
    permutations = ["partielle_mapped", "position_based", "order_based"]
    colors = ['lightsalmon', 'darkblue', 'gold']
    plt.bar(permutations, array, color=colors)
    plt.ylabel("Anzahl")
    plt.title("Anzahl der besten Crossover über alle Durchläufe")
    plt.show();

    return



if __name__ == '__main__':
    array = [0, 0, 0]
    counter = 0
    while(counter != 3):
        sys.stdout.write("COUNTER: " + str(counter)+"\n")
        nq1 = NQueen1(START_SIZE, MAX_EPOCHS, MATING_PROBABILITY, MUTATION_RATE, MIN_SELECT, MAX_SELECT,
                      OFFSPRING_PER_GENERATION, MINIMUM_SHUFFLES, MAXIMUM_SHUFFLES, PBC_MAX, MAX_LENGTH)

        nq1.initialize_chromosomes()
        nq1.genetic_algorithm()
        #zeige pro Druchlauf permutation-Anzahl und permutationen pro Epoche
        nq1.show_permutation_amount()
        nq1.show_crossover_per_epoche()

        array = np.array(array) + nq1.get_best_crossover()
        counter += 1

    print(array)
    show_overall_permutation_amount(array)



