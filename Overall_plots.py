import matplotlib.pyplot as plt
import numpy as np
from Chromosome import Chromosome


class OverallPlots:

    def recomb_str(recomb):
        recombis = {
            0: "PMX",
            1: "POS",
            2: "TPX",
            3: "OX",
            4: "mPOS"
            }
        return recombis.get(recomb, "nothing")


    def problem_str(problem):
        problems = {
            0: "nQueens",
            1: "Schwefel",
            2: "Himmelblau",
            3: "Griewank"
        }
        return problems.get(problem, "nothing")


    def show_overall_permutation_amount(array, problem):
        permutations = ["PMX", "POS", "TPX", "OX"]
        colors = ['midnightblue', 'royalblue', 'lightsteelblue', 'lavender']
        plt.bar(permutations, array, color=colors)
        plt.ylabel("Anzahl")
        plt.title("Anzahl der besten Crossover über alle Evolutionen, "+problem)
        plt.show()
        #plt.savefig('overall_permutation_amount.png')

        return

    def show_overall_minima(minarray, problem):
        x = np.arange(0, len(minarray), 1)
        minarray = np.asarray(minarray)
        y = minarray[x]
        plt.xlabel('Durchgänge')
        plt.ylabel('Minimale Fitness')
        plt.plot(x, y, color = 'darkblue')
        plt.title("Minima pro Durchlauf " + problem)
        plt.axis([0, 100, 0, 2500])
        plt.show()
        OverallPlots.show_overall_minbox(minarray, problem)

        return

    def show_overall_minbox(minarray, problem):

        plt.ylabel("Fitness")
        plt.title("Boxplot: Minima pro Durchlauf " + problem)
        plt.boxplot(minarray)
        plt.axis([0, 2, 0, 2500])
        plt.show()

        return

    def boxplot(p1, p2, k1, k2, recomb, problem):
        data_to_plot = [p1, p2, k1, k2]

        plt.figure(1, figsize=(9, 6))
        plt.ylabel("Fitness")
        plt.title('Eltern und Kinder: '+OverallPlots.recomb_str(recomb)+", "+OverallPlots.problem_str(problem))
        plt.boxplot(data_to_plot)
        plt.xticks([1,2,3,4],['Elternteil1', 'Elternteil2', 'Kind1', 'Kind2'])
        plt.show()
        return

    def percentage_table(p1, p2, k1, k2, recomb, problem):

        len_array = len(p1)

        # Anzahl, wie oft Kinder weniger Konflikte als Eltern haben
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
            # Kind1 weniger Konflikte als Elternteil1 und/oder Elternteil2
            if k1[i] < p1[i]:
                better_k1_p1 += 1
                if k1[i] < p2[i]:
                    better_k1_beide += 1
            if k1[i] < p2[i]:
                better_k1_p2 += 1

            # Kind2 weniger Konflikte als Elternteil1 und/oder Elternteil2
            if k2[i] < p1[i]:
                better_k2_p1 += 1
                if k2[i] < p2[i]:
                    better_k2_beide += 1
            if k2[i] < p2[i]:
                better_k2_p2 += 1

            # Kind1 mehr Konflikte als Elternteil1 und/oder Elternteil2
            if k1[i] > p1[i]:
                worse_k1_p1 += 1
                if k1[i] > p2[i]:
                    worse_k1_beide += 1
            if k1[i] > p2[i]:
                worse_k1_p2 += 1

            # Kind2 mehr Konflikte als Elternteil1 und/oder Elternteil2
            if k2[i] > p1[i]:
                worse_k2_p1 += 1
                if k2[i] > p2[i]:
                    worse_k2_beide += 1
            if k2[i] > p2[i]:
                worse_k2_p2 += 1

        # Berechnung der Prozente, auf zwei Nachkommastellen ge
        b_k1_vs_p1 = round((better_k1_p1 / len_array * 100), 2)
        b_k1_vs_p2 = round((better_k1_p2 / len_array * 100), 2)
        b_k2_vs_p1 = round((better_k2_p1 / len_array * 100), 2)
        b_k2_vs_p2 = round((better_k2_p2 / len_array * 100), 2)
        b_k1_vs_beide = round((better_k1_beide / len_array * 100), 2)
        b_k2_vs_beide = round((better_k2_beide / len_array * 100), 2)
        w_k1_vs_p1 = round((worse_k1_p1 / len_array * 100), 2)
        w_k1_vs_p2 = round((worse_k1_p2 / len_array * 100), 2)
        w_k2_vs_p1 = round((worse_k2_p1 / len_array * 100), 2)
        w_k2_vs_p2 = round((worse_k2_p2 / len_array * 100), 2)
        w_k1_vs_beide = round((worse_k1_beide / len_array * 100), 2)
        w_k2_vs_beide = round((worse_k2_beide / len_array * 100), 2)

        plt.figure()
        plt.title(OverallPlots.recomb_str(recomb)+", "+OverallPlots.problem_str(problem))
        title = " Durchgänge = " + str(len_array)

        table_data = [
            ["Kind1/Parent1", str(b_k1_vs_p1) + "% / " + str(w_k1_vs_p1) + "%",
             str(better_k1_p1) + " / " + str(worse_k1_p1)],
            ["Kind1/Parent2", str(b_k1_vs_p2) + "% / " + str(w_k1_vs_p2) + "%",
             str(better_k1_p2) + " / " + str(worse_k1_p2)],
            ["Kind2/Parent1", str(b_k2_vs_p1) + "% / " + str(w_k2_vs_p1) + "%",
             str(better_k2_p1) + " / " + str(worse_k2_p1)],
            ["Kind2/Parent2", str(b_k2_vs_p2) + "% / " + str(w_k2_vs_p2) + "%",
             str(better_k2_p2) + " / " + str(worse_k2_p2)],
            ["Kind1/beide", str(b_k1_vs_beide) + "% / " + str(w_k1_vs_beide) + "%",
             str(better_k1_beide) + " / " + str(worse_k1_beide)],
            ["Kind2/beide", str(b_k2_vs_beide) + "% / " + str(w_k2_vs_beide) + "%",
             str(better_k2_beide) + " / " + str(worse_k2_beide)]]

        table = plt.table(cellText=table_data, loc='center',
                         colLabels=[title, "Prozentsatz besser/schlechter", "Anzahl besser/schlechter"])
        table.scale(1.2, 1.5)
        table.set_fontsize(20)
        plt.axis('off')
        plt.show()

        return

    def show_permutation_amount(self):
        data = [self.partielle_mapped_co, self.position_based_co, self.two_point_co, self.order_based_co]
        permutations = ["PMX", "POS", "TPX", "OX"]
        colors = ['midnightblue', 'royalblue', 'lightsteelblue','lavender']
        plt.bar(permutations, data, color=colors)
        plt.ylabel("Anzahl")
        plt.title("Anzahl der genutzten Crossover")
        plt.show()
        #plt.savefig('permutation_amount.png')

        return

    def show_crossover_per_epoche(self):
        pmarray=np.asarray(self.array_p_m)
        pbarray=np.asarray(self.array_p_b)
        tparray=np.asarray(self.array_t_p)
        obarray=np.asarray(self.array_o_b)
        x = np.arange(1, self.epoch, 1)
        #sys.stdout.write(str(x) +" this is it\n")
        #sys.stdout.write(str(pmarray) + " pmarray\n")
        y1 = pmarray[x]
        y2 = pbarray[x]
        y3 = tparray[x]
        y4 = obarray[x]
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt.xlabel('Epochen')
        plt.ylabel('Anzahl')
        p_m, = ax.plot(x, y1, color='midnightblue', label ='PMX')
        p_b, = ax.plot(x, y2, color='royalblue', label ='POS')
        t_p, = ax.plot(x, y3, color='lightsteelblue', label ='TPX')
        o_b, = ax.plot(x, y4, color='lavender', label ='OX')
        ax.legend([p_m, p_b, t_p, o_b],["partially_mapped", "position_based","two_point", "order_based"])
        plt.show()
        #plt.savefig('crossovers_per_epoche.png')

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
        ax.plot(x, y1, color='cornflowerblue', label ='partially_mapped')
        plt.show()
        #plt.savefig('crossovers_per_epoche.png')

        return
