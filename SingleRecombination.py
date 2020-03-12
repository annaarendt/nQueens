import random
import sys
import matplotlib.pyplot as plt
from Chromosome import Chromosome
from EA import NQueen1

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
        self.mMinimumShuffles = minShuffles
        self.mMaximumShuffles = maxShuffles
        self.mPBCMax = pbcMax
        self.mMaxLength = maxLength

        self.epoch = 0
        self.childCount = 0
        self.population = []
        #1. und 2. stelle für Konflikte der Eltern, letzen beiden der Kinder
        self.conflict_array = [0] * 4
        self.mutations = 0

        return

    def do_mating(self):

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
            sys.stdout.write(str(chromA.get_conflicts()) + " Konflikten: ")
            chromA.toStr() + "\n"
            sys.stdout.write("Parent2 mit ")
            sys.stdout.write(str(chromB.get_conflicts()) + " Konflikten: ")
            chromB.toStr() + "\n"

            #TODO diese Zeile ändern, je nach gewünschter Rekombination (bad funktioniert nur als self.bad...)
            NQueen1.bad_recombination(self, chromA, chromB, newIndex1, newIndex2)

            newChromo1 = self.population[newIndex1]
            # Konsole:
            sys.stdout.write("Kind1 mit ")
            sys.stdout.write(str(newChromo1.get_conflicts()) + " Konflikten: ")
            newChromo1.toStr()+"\n"

            newChromo2 = self.population[newIndex2]
            # Konsole:
            sys.stdout.write("Kind2 mit ")
            sys.stdout.write(str(newChromo2.get_conflicts()) + " Konflikten: ")
            newChromo2.toStr()+"\n"

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


    def genetic_algorithm(self):

        done = False

        while not done:
            popSize = len(self.population)
            for i in range(popSize):
                if self.epoch == self.mEpochs:
                    done = True

            self.do_mating()

            self.prep_next_epoch()

            self.epoch += 1

        sys.stdout.write("done.\n")

        sys.stdout.write("Completed " + str(self.epoch) + " epochs.\n")
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

    #Anzahl, wie oft Kinder weniger Konflikte als Eltern haben
    better_k1_p1 = 0
    better_k1_p2 = 0
    better_k2_p1 = 0
    better_k2_p2 = 0
    better_k1_beide = 0
    better_k2_beide = 0

    worse_k1_p1 = 0
    worse_k1_p2 = 0
    worse_k2_p1 = 0
    worse_k2_p2 = 0
    worse_k1_beide = 0
    worse_k2_beide = 0


    for i in range(len(p1)):
        #Kind1 weniger Konflikte als Elternteil1 und/oder Elternteil2
        if k1[i] < p1[i]:
            better_k1_p1 += 1
            if k1[i] < p2[i]:
                better_k1_beide += 1
        if k1[i] < p2[i]:
            better_k1_p2 += 1

        #Kind2 weniger Konflikte als Elternteil1 und/oder Elternteil2
        if k2[i] < p1[i]:
            better_k2_p1 += 1
            if k2[i] < p2[i]:
                better_k2_beide += 1
        if k2[i] < p2[i]:
            better_k2_p2 += 1

        #Kind1 mehr Konflikte als Elternteil1 und/oder Elternteil2
        if k1[i] > p1[i]:
            worse_k1_p1 += 1
            if k1[i] > p2[i]:
                worse_k1_beide += 1
        if k1[i] > p2[i]:
            worse_k1_p2 += 1

        #Kind2 mehr Konflikte als Elternteil1 und/oder Elternteil2
        if k2[i] > p1[i]:
            worse_k2_p1 += 1
            if k2[i] > p2[i]:
                worse_k2_beide += 1
        if k2[i] > p2[i]:
            worse_k2_p2 += 1

    #Berechnung der Prozente, auf zwei Nachkommastellen ge
    b_k1_vs_p1 = round((better_k1_p1 / len_array * 100),2)
    b_k1_vs_p2 = round((better_k1_p2 / len_array * 100),2)
    b_k2_vs_p1 = round((better_k2_p1 / len_array * 100),2)
    b_k2_vs_p2 = round((better_k2_p2 / len_array * 100),2)
    b_k1_vs_beide = round((better_k1_beide / len_array * 100),2)
    b_k2_vs_beide = round((better_k2_beide / len_array * 100),2)
    w_k1_vs_p1 = round((worse_k1_p1 / len_array * 100),2)
    w_k1_vs_p2 = round((worse_k1_p2 / len_array * 100),2)
    w_k2_vs_p1 = round((worse_k2_p1 / len_array * 100),2)
    w_k2_vs_p2 = round((worse_k2_p2 / len_array * 100),2)
    w_k1_vs_beide = round((worse_k1_beide / len_array * 100),2)
    w_k2_vs_beide = round((worse_k2_beide / len_array * 100),2)

    fig = plt.figure()
    title="Durchgänge = "+str(len_array)
    ax = fig.add_subplot(1, 1, 1)

    table_data = [
        ["Kind1/Parent1", str(b_k1_vs_p1)+"% / "+str(w_k1_vs_p1)+"%", str(better_k1_p1)+" / "+str(worse_k1_p1) ],
        ["Kind1/Parent2", str(b_k1_vs_p2)+"% / "+str(w_k1_vs_p2) + "%", str(better_k1_p2)+" / "+str(worse_k1_p2)],
        ["Kind2/Parent1", str(b_k2_vs_p1)+"% / "+str(w_k2_vs_p1) + "%", str(better_k2_p1)+" / "+str(worse_k2_p1)],
        ["Kind2/Parent2", str(b_k2_vs_p2)+"% / "+str(w_k2_vs_p2) + "%", str(better_k2_p2)+" / "+str(worse_k2_p2)],
        ["Kind1/beide", str(b_k1_vs_beide)+"% / "+str(w_k1_vs_beide) + "%", str(better_k1_beide)+" / "+str(worse_k1_beide)],
        ["Kind2/beide", str(b_k2_vs_beide)+"% / "+str(w_k2_vs_beide) + "%", str(better_k2_beide)+" / "+str(worse_k2_beide)]]

    table = ax.table(cellText=table_data, loc='center', colLabels=[title,"Prozentsatz besser/schlechter", "Anzahl besser/schlechter"])
    table.scale(1.2, 1.5)
    table.set_fontsize(20)
    ax.axis('off')
    plt.show()

    return

if __name__ == '__main__':
    COUNTER = 0
    END = 500
    p1 = [0] * END
    p2 = [0] * END
    k1 = [0] * END
    k2 = [0] * END
    while (COUNTER != END):
        sr1 = SingleRekomb(START_SIZE, MAX_EPOCHS, OFFSPRING_PER_GENERATION, MINIMUM_SHUFFLES, MAXIMUM_SHUFFLES,
                      PBC_MAX, MAX_LENGTH)

        NQueen1.initialize_chromosomes(sr1)
        sr1.genetic_algorithm()

        p1[COUNTER] = sr1.get_conflict_array()[0]
        p2[COUNTER] = sr1.get_conflict_array()[1]
        k1[COUNTER] = sr1.get_conflict_array()[2]
        k2[COUNTER] = sr1.get_conflict_array()[3]
        #konsole_konflikte je Durchlauf
        str1 = ','.join(str(e) for e in sr1.get_conflict_array())
        print("Konflikte: "+str1)
        sys.stdout.write("COUNTER: " + str(COUNTER)+"\n")

        COUNTER += 1

    boxplot(p1, p2, k1, k2)
    percentage_table(p1, p2, k1, k2)








