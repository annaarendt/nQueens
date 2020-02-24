import matplotlib.pyplot as plt

class OverallPlots:

    def show_overall_permutation_amount(array):
        permutations = ["partielle_mapped", "position_based", "order_based"]
        colors = ['lightsalmon', 'darkblue', 'gold']
        plt.bar(permutations, array, color=colors)
        plt.ylabel("Anzahl")
        plt.title("Anzahl der besten Crossover über alle Durchläufe")
        plt.show()
        #plt.savefig('overall_permutation_amount.png')

        return