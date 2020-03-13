import sys
from Overall_plots import OverallPlots
from SingleRecombination import SingleRekomb
from ChromosomeBenchmark import ChromosomeDEAP
from Chromosome import Chromosome
from EA import NQueen1
from EA_DEAP import NDEAP

START_SIZE = 2  # Population size at start.
MAX_EPOCHS = 0  # Arbitrary number of test cycles. EIGENTLICH AUF 1000 GESETZT!
OFFSPRING_PER_GENERATION = 1  # New offspring created per generation. Range: 0 < OFFSPRING_PER_GENERATION < MAX_SELECT.
MINIMUM_SHUFFLES = 8  # For randomizing starting chromosomes
MAXIMUM_SHUFFLES = 20
PBC_MAX = 4  # Maximum Position-Based Crossover points. Range: 0 < PBC_MAX < 8 (> 8 isn't good).
MAX_LENGTH = 6  # chess board width.


class SingleRekomb_bench:
    def __init__(self, startSize, maxEpochs, generation, minShuffles,
                 maxShuffles, pbcMax, maxLength):

        SingleRekomb.__init__(self, startSize, maxEpochs, generation, minShuffles,
                 maxShuffles, pbcMax, maxLength)

        return

    def do_mating(self):

        for i in range(self.mOffspringPerGeneration):
            chromA = self.population[0]
            chromB = self.population[1]

            newChromo1 = ChromosomeDEAP(self.mMaxLength, 0)
            newChromo2 = ChromosomeDEAP(self.mMaxLength, 0)
            self.population.append(newChromo1)
            newIndex1 = len(self.population) - 1
            self.population.append(newChromo2)
            newIndex2 = len(self.population) - 1

            sys.stdout.write("Parent1 mit Fitness ")
            sys.stdout.write(str(Chromosome.get_fitness(chromA))+" : ")
            Chromosome.toStr(chromA) + "\n"
            sys.stdout.write("Parent2 mit Fitness ")
            sys.stdout.write(str(Chromosome.get_fitness(chromB)) + " : ")
            Chromosome.toStr(chromB) + "\n"

            #TODO NQueen.(gewünschte rekombination), bei bad-recombination NDEAP.bad_recombination
            NQueen1.position_based_crossover(self, chromA, chromB, newIndex1, newIndex2)

            newChromo1 = self.population[newIndex1]
            # Konsole:
            sys.stdout.write("Kind1 mit Fitness ")
            sys.stdout.write(str(Chromosome.get_fitness(newChromo1)) + " : ")
            Chromosome.toStr(newChromo1)+"\n"

            newChromo2 = self.population[newIndex2]
            # Konsole:
            sys.stdout.write("Kind2 mit Fitness ")
            sys.stdout.write(str(Chromosome.get_fitness(newChromo2)) + " : ")
            Chromosome.toStr(newChromo2)+"\n"

            self.childCount += 2
            sys.stdout.write(str(len(self.population)) + " Länge der Population\n")

            SingleRekomb.set_conflict_array(self, Chromosome.get_fitness(chromA), Chromosome.get_fitness(chromB),
                                    Chromosome.get_fitness(newChromo1), Chromosome.get_fitness(newChromo2))

            # plot für jeden einzelnen Druchgang
            #self.show_conflicts(chromA, chromB, newChromo1, newChromo2)

        return


    def genetic_algorithm(self):

        done = False

        while not done:
            popSize = len(self.population)
            for i in range(popSize):
                if self.epoch == self.mEpochs:
                    done = True

            self.do_mating()

            SingleRekomb.prep_next_epoch(self)

            self.epoch += 1

        sys.stdout.write("done.\n")

        sys.stdout.write("Completed " + str(self.epoch) + " epochs.\n")
        return

if __name__ == '__main__':
    COUNTER = 0
    END = 500
    p1 = [0] * END
    p2 = [0] * END
    k1 = [0] * END
    k2 = [0] * END
    while (COUNTER != END):
        sr1 = SingleRekomb_bench(START_SIZE, MAX_EPOCHS, OFFSPRING_PER_GENERATION, MINIMUM_SHUFFLES, MAXIMUM_SHUFFLES,
                      PBC_MAX, MAX_LENGTH)

        NDEAP.initialize_chromosomes(sr1)
        sr1.genetic_algorithm()

        p1[COUNTER] = SingleRekomb.get_conflict_array(sr1)[0]
        p2[COUNTER] = SingleRekomb.get_conflict_array(sr1)[1]
        k1[COUNTER] = SingleRekomb.get_conflict_array(sr1)[2]
        k2[COUNTER] = SingleRekomb.get_conflict_array(sr1)[3]
        #konsole_konflikte je Durchlauf
        str1 = ','.join(str(e) for e in SingleRekomb.get_conflict_array(sr1))
        print("Konflikte: "+str1)
        sys.stdout.write("COUNTER: " + str(COUNTER)+"\n")

        COUNTER += 1

    OverallPlots.boxplot(p1, p2, k1, k2)
    OverallPlots.percentage_table(p1, p2, k1, k2)









