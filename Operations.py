import math
import random
import sys
import numpy as np
from Chromosome import Chromosome
from Overall_plots import OverallPlots


class Operations:

    def bad_recombination_qu(self, chromA, chromB, child1, child2):
        tempArray1 = [0] * self.mMaxLength
        tempArray2 = [0] * self.mMaxLength
        thisChromo = chromA
        thatChromo = chromB
        newChromo1 = self.population[child1]
        newChromo2 = self.population[child2]

        # Choose and sort the crosspoints.
        numPoints = random.randrange(0, self.mPBCMax)  # wenn PBC_MAX höher als 6 or 8.
        crossPoints = [0] * numPoints
        for i in range(numPoints):
            crossPoints[i] = Operations.get_exclusive_random_integer_by_array(self, 0, self.mMaxLength - 1, crossPoints)
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
        newChromo1.compute_fitness()

        #Kind1 wird rausgehauen und index child2 um eins heragesetzt weil sich population auch wieder ändert
        if newChromo1.get_fitness() < thisChromo.get_fitness() or newChromo1.get_fitness() < thatChromo.get_fitness():
            newChromo1.set_fitness(18)
#
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

        newChromo2.compute_fitness()

        if newChromo2.get_fitness() < thisChromo.get_fitness() or newChromo2.get_fitness() < thatChromo.get_fitness():
            sys.stdout.write("Kind2 zu wenig Konflikte: " + str(newChromo2.get_fitness()) + "\n")
            newChromo2.set_fitness(18)
            sys.stdout.write("Kind2 wurde auf 18 Konlikte erhöht\n")


        sys.stdout.write("BAD_recombination verwendet.\n")
        return

    def two_point_crossover(self, chromA, chromB, child1, child2):
        thisChromo = chromA
        thatChromo = chromB
        newChromo1 = self.population[child1]
        newChromo2 = self.population[child2]

        # Sample the 2 crossover points randomly
        cp1 = np.random.randint(self.mMaxLength)
        cp2 = np.random.randint(self.mMaxLength)

        # Avoid that cp1 and cp2 are equal
        while cp1 == cp2:
            cp2 = np.random.randint(self.mMaxLength)

        # Swap if cp1 is bigger than cp2 (otherwise array slicing won't work)
        if cp1 > cp2:
            cp1, cp2 = cp2, cp1

        # Copy Parent genes to offspring.
        for i in range(self.mMaxLength):
            Chromosome.set_data(newChromo1, i, Chromosome.get_data(thisChromo, i))
            Chromosome.set_data(newChromo2, i, Chromosome.get_data(thatChromo, i))

        newChromo1.mData[cp1:cp2], newChromo2.mData[cp1:cp2] = newChromo2.mData[cp1:cp2], newChromo1.mData[cp1:cp2]


        newChromo1.compute_fitness()
        newChromo2.compute_fitness()

        return

    def partially_mapped_crossover(self, chromA, chromB, child1, child2):
        thisChromo = chromA
        thatChromo = chromB
        newChromo1 = self.population[child1]
        newChromo2 = self.population[child2]

        crossPoint1 = random.randrange(0, self.mMaxLength)
        crossPoint2 = Operations.get_exclusive_random_integer(self, self.mMaxLength, crossPoint1)
        if crossPoint2 < crossPoint1:
            j = crossPoint1
            crossPoint1 = crossPoint2
            crossPoint2 = j

        # Copy Parent genes to offspring.
        for i in range(self.mMaxLength):
            Chromosome.set_data(newChromo1, i, Chromosome.get_data(thisChromo, i))
            Chromosome.set_data(newChromo2, i, Chromosome.get_data(thatChromo, i))

        for i in range(crossPoint1, crossPoint2 + 1):
            # // Get the two items to swap.
            item1 = Chromosome.get_data(thisChromo, i)
            item2 = Chromosome.get_data(thatChromo, i)
            pos1 = 0
            pos2 = 0

            # Get the items' positions in the offspring.
            for j in range(self.mMaxLength):
                if Chromosome.get_data(newChromo1, j) == item1:
                    pos1 = j
                elif Chromosome.get_data(newChromo1, j) == item2:
                    pos2 = j

            # Swap them.
            if item1 != item2:
                Chromosome.set_data(newChromo1, pos1, item2)
                Chromosome.set_data(newChromo1, pos2, item1)

            # Get the items'  positions in the offspring.
            for j in range(self.mMaxLength):
                if Chromosome.get_data(newChromo2, j) == item2:
                    pos1 = j
                elif Chromosome.get_data(newChromo2, j) == item1:
                    pos2 = j

            # Swap them.
            if item1 != item2:
                Chromosome.set_data(newChromo2, pos1, item1)
                Chromosome.set_data(newChromo2, pos2, item2)

            newChromo1.compute_fitness()
            newChromo2.compute_fitness()


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
            crossPoints[i] = Operations.get_exclusive_random_integer_by_array(self, 0, self.mMaxLength - 1, crossPoints)
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

        newChromo1.compute_fitness()
        newChromo2.compute_fitness()

        # sys.stdout.write("Position-based Crossover verwendet.\n")
        return

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
            points[i] = Operations.get_exclusive_random_integer_by_array(self, 0, self.mMaxLength - 1, points)
        points.sort()

        crossNumbers = Operations.findCrossNumbers(self, numPoints, points, thisChromo)
        crossNumbers2 = Operations.findCrossNumbers(self, numPoints, points, thatChromo)

        # von Parent2 werden Positionen der Zahlen von Parent1 (in crossNumbers) gesucht und in Array cpoints
        # gespeichert.
        # Von den Zahlen die bei P2  den ausgewählten Zahlen in crossNumbers entsprechen werden die Positionen in
        # einem Array gesammelt (points)
        cpoints = Operations.findCPoints(self, numPoints, thatChromo, crossNumbers)
        cpoints2 = Operations.findCPoints(self, numPoints, thisChromo, crossNumbers2)

        # temporäres Array erzeugen mit Lücken in Positionen von cpoints/cpoint2, danach Lücken mit crossNumbers/
        # crossNumbers2 füllen

        array1 = Operations.fullfill(self, numPoints, cpoints, thatChromo, crossNumbers)
        array2 = Operations.fullfill(self, numPoints, cpoints2, thisChromo, crossNumbers2)

        for i in range(self.mMaxLength):
            Chromosome.set_data(newChromo1, i, array1[i])

        for i in range(self.mMaxLength):
            Chromosome.set_data(newChromo2, i, array2[i])

        if not points:
            for i in range(self.mMaxLength):
                Chromosome.set_data(newChromo1, i, Chromosome.get_data(thatChromo, i))
                Chromosome.set_data(newChromo2, i, Chromosome.get_data(thisChromo, i))

        newChromo1.compute_fitness()
        newChromo2.compute_fitness()

        # sys.stdout.write("Order-based Crossover verwendet.\n")

        return

    def findCrossNumbers(self, numPoints, points, chromosome):

        crossNumbers = [0] * numPoints
        for i in range(numPoints):
            for j in range(self.mMaxLength):
                if points[i] == j:
                    crossNumbers[i] = Chromosome.get_data(chromosome,points[i])

        return crossNumbers

    def findCPoints(self, numPoints, chromosome, numbers):

        cpoints = [0] * numPoints

        k = 0
        for i in range(self.mMaxLength):
            matchFound = False
            for j in range(numPoints):
                if Chromosome.get_data(chromosome,i) == numbers[j]:
                    matchFound = True

            if matchFound == True:
                if cpoints:
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
                    tempArray[i]=Chromosome.get_data(chromosome,i)

            if matchFound==True:
                tempArray[i] = 0

        # Lücken in tempArray werden mit Zahlen von numbers aufgefüllt
        k = 0
        # hier wird gecheckt ob crossNumbers leer ist, wenn ja wird nichts gemacht
        if numbers:
            for i in range(self.mMaxLength):
                if tempArray[i] == 0:
                    if Chromosome.get_data(chromosome,i) != 0:
                        tempArray[i] = numbers[k]
                        k += 1
                    elif 0 in numbers:
                        try:
                            tempArray[i] = numbers[k]
                            k += 1
                        except IndexError:
                            sys.stdout.write(str(k)+" :k, Endphase.. i: "+str(i)+ "\n")


        return tempArray

    def bad_recombination_deap(self, chromA, chromB, child1, child2):
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
                crossPoints[i] = Operations.get_exclusive_random_integer_by_array(self, 0, self.mMaxLength - 1,
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
                Chromosome.set_fitness(newChromo1, 5000)
                sys.stdout.write("Kind1 wurde auf Fitness 5000 erhöht\n")
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
                Chromosome.set_fitness(newChromo2, 5000)
                sys.stdout.write("Kind2 wurde auf Fitness 5000 erhöht\n")

            sys.stdout.write(str(crossPoints) + " CrossPoints\n")

            sys.stdout.write("BAD_recombination verwendet.\n")
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

    def exchange_mutation(self, index, exchanges):
        i = 0
        done = False

        thisChromo = self.population[index]
        sys.stdout.write("Zahl exchange: " + str(exchanges) + "\n")
        sys.stdout.write("Vorher: " + str(Chromosome.toStr(thisChromo) + "\n"))

        while not done:
            gene1 = random.randrange(0, self.mMaxLength)
            gene2 = Operations.get_exclusive_random_integer(self, self.mMaxLength, gene1)

            # Exchange the chosen genes.
            tempData = Chromosome.get_data(thisChromo, gene1)
            Chromosome.set_data(thisChromo, gene1, Chromosome.get_data(thisChromo, gene2))
            Chromosome.set_data(thisChromo, gene2, tempData)

            if i == exchanges:
                done = True

            i += 1
        sys.stdout.write("Nachher: "+str(Chromosome.toStr(thisChromo)+"\n"))
        thisChromo.compute_fitness()

        self.mutations += 1

        return

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

    #ersetzt eine gene durch eine zufällig Anderen [-500,500,1]
    def new_gene_mutation(self, index):

        thisChromo = self.population[index]
        point1 = random.randrange(0, self.mMaxLength)

        rand= random.randrange(-500,500,1)
        Chromosome.set_data(thisChromo, point1, rand)

        sys.stdout.write("New_Gene_Mutation verwendet.\n")
        thisChromo.compute_fitness()
        self.mutations += 1

        return

    #Selektion der Eltern mit Roulette-Methode. Je besser (kleiner) Fitness, desto mehr Anteil auf dem Roulette
    def roulette_selection(self):
        genTotal = 0.0
        sumProp = 0.0
        popSize = len(self.population)

        for i in range(popSize):
            worst = self.population[self.get_maximum()]
            thisChromo = self.population[i]
            #zählt quasi die Fitness aller Chromosomen zusammen (abgezogen von Indivuduum mit höchstem Konfliktwert)
            genTotal += (Chromosome.get_fitness(worst) - Chromosome.get_fitness(thisChromo))

        sys.stdout.write(str(genTotal) + " GenTotal(Roulette)\n")
        sys.stdout.write(str(popSize) + " popSize\n")

        for i in range(popSize):
            worst = self.population[self.get_maximum()]

            thisChromo = self.population[i]
            if genTotal != 0:
                probability = ((Chromosome.get_fitness(worst) - Chromosome.get_fitness(thisChromo)) / genTotal)
            else:
                probability = 0
            sumProp += probability
            Chromosome.set_selection_probability(thisChromo, probability)

        sys.stdout.write(str(round(sumProp, 3)) + " alle Props\n")

        for i in range(self.mOffspringPerGeneration):
            rouletteSpin = random.uniform(0, 1)
            j = 0
            selTotal = 0
            done = False
            while not done:
                thisChromo = self.population[j]
                selTotal += Chromosome.get_selection_probability(thisChromo)
                if selTotal >= rouletteSpin:
                    if j == 0:
                        thatChromo = self.population[j]
                    elif j >= popSize - 1:
                        thatChromo = self.population[popSize - 1]
                    else:
                        thatChromo = self.population[j - 1]

                    Chromosome.set_selected(thatChromo, True)
                    done = True
                else:
                    j += 1

                #wenn nur noch Individuen mit 0 Konflikten vorhanden sind, immer die Chromosomen mit jeweiliger zahl
                #der anzahl an Nachwüchsen selected
                if sumProp == 0:
                    Chromosome.set_selected(self.population[i],True)
                    done = True

        return

    def choose_first_parent(self):
        parent = 0
        done = False
        counter = 0

        while not done and counter < 100:
            # Randomly choose an eligible parent.
            counter += 1
            #print(counter)
            parent = random.randrange(0, len(self.population) - 1)
            thisChromo = self.population[parent]
            if Chromosome.get_selected(thisChromo) == True:
                done = True

        return parent

    def choose_second_parent(self, parentA):
        parentB = 0
        done = False

        while not done:
            #print("in choose second parent angekommen")
            #test um zu sehen wieviele Chromosomen in eine rPopulation seleced wurden (wenn nur einer, dann Grund warum
            #sich schleife aufhängt
            probe_counter = 0
            for i in range (len(self.population)):
                if Chromosome.get_selected(self.population[i]) == True:
                    probe_counter += 1

            # Randomly choose an eligible parent.
            parentB = random.randrange(0, len(self.population) - 1)
            if parentB != parentA:
                thisChromo = self.population[parentB]
                if Chromosome.get_selected(thisChromo) == True:
                    done = True
            #im Falle, dass nur ein Chromosom selected wurde (im Falle von alle Konflikte sind 0), wird ein zufälliges
            #Chromosom als zweiten Parent gewählt
                elif probe_counter <= 2:
                    sys.stdout.write("2. parent gefunden durch Schleifenaufhängung\n")
                    done = True

        return parentB

    #funktion wird nicht mehr gebraucht, da konflikte als Fitness dienen. Zeigt nun lediglich den aktuellen
    #Stand der Konflikte am Anfang jeder epoche an.
    def get_fitness(self):
        #fitness entspricht den Konflikten. Niedrige Fitness ist gut, hohe ist schlecht
        popSize = len(self.population)

        for i in range(popSize):
            thisChromo = self.population[i]
            worst = self.population[self.get_maximum()]
            self.worst = Chromosome.get_fitness(worst)
            best = self.population[self.get_minimum()]
            self.best = Chromosome.get_fitness(best)

        #sys.stdout.write("\n"+str(Chromosome.get_fitness(worst)) + ": Maximale Fitness\n")
        #sys.stdout.write(str(Chromosome.get_fitness(best)) + ": Minimale Fitness\n")
        self.current_best_fitness = Chromosome.get_fitness(best)
        self.array_fitness.append(self.current_best_fitness)


        return

    def math_round(self, inValue):
        if math.modf(inValue)[0] >= 0.5:
            outValue = math.ceil(inValue)
        else:
            outValue = math.floor(inValue)
        return outValue

    def prep_next_epoch(self):
        popSize = len(self.population)
        #schlechten Individuen rauslöschen und Zahl der Population wieder auf Startzahl anpassen

        for i in range(popSize):
            while len(self.population) > self.mStartSize:
                #worst=self.population[self.get_maximum()]
                self.population.pop(self.get_maximum())


        popSize = len(self.population)
        # Reset flags for selected individuals.
        for i in range(popSize):
            thisChromo = self.population[i]
            Chromosome.set_selected(thisChromo,False)

        return

    def get_foundMimimum(self):
        return self.found_mimimum

    def set_foundMimimum(self, min):
        self.found_mimimum = min
        return

    def get_best(self):
        return self.best

    #methode gibt ein numpy-Array zuruek, abhängig, welches crossover das beste in dem jeweiligen
    #durchgang war. erste stelle im array ist partielle-mapped, dann kommt position-based und dann
    #an dritter stelle order-based
    def get_best_crossover(self):
        pb = self.position_based_co
        ob = self.order_based_co
        pm = self.partielle_mapped_co
        tp = self.two_point_co
        if(pm >= pb and pm >= ob and pm >= tp):
            array = [1, 0, 0,0]
            print("PM IS BEST")
        elif(pb >= pm and pb >= ob and pb >= tp):
            array = [0, 1, 0,0]
            print("PB IS BEST")
        elif(tp >= pm and tp >= ob and tp >= pb):
            array = [0, 0, 1,0]
            print("TP IS BEST")
        else:
            array = [0, 0, 0,1]
            print("OB IS BEST")


        return np.array(array)




