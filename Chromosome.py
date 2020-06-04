class Chromosome:
    #crossover-type 0,1,2 oder 3 (partitielle_mapped, position_based, two-point und order_based)
    def __init__(self, maxLength, crossover_type):
        self.mMaxLength = maxLength
        self.mFitness = 0.0
        self.mSelected = False
        self.mSelectionProbability = 0.0
        self.mCrossover = crossover_type

        self.mData = [0] * maxLength
        for i in range(self.mMaxLength):
            self.mData[i] = i
        return

    def compute_fitness(self):

        #wenn es mehrere Vorkommen einer Zahl im Chromosom gibt, wird die Fitness sehr schlecht.
        if len(set(self.mData)) == len(self.mData):

            board = []
            conflicts = 0
            dx = [-1, 1, -1, 1]
            dy = [-1, 1, 1, -1]

            for i in range(self.mMaxLength):
                board.append([""] * self.mMaxLength)
                board[i][self.mData[i]] = "Q"

            # Walk through each of the Queens and compute the number of conflicts.
            for i in range(self.mMaxLength):
                x = i
                y = self.mData[i]

                # Check diagonals.
                for j in range(4):
                    tempx = x
                    tempy = y
                    done = False
                    while not done:
                        tempx += dx[j]
                        tempy += dy[j]
                        if (tempx < 0 or tempx >= self.mMaxLength) or (tempy < 0 or tempy >= self.mMaxLength):
                            done = True
                        else:
                            if board[tempx][tempy] == "Q":
                                conflicts += 1

            self.mFitness = conflicts

        #wenn doppete Allele vorkommen ist Fitness schlecht(hoch)
        else: self.mFitness = 30


        return

    def get_fitness(self):
        return self.mFitness

    def set_fitness(self, fitness):
        self.mFitness = fitness

    def set_selection_probability(self, probability):
        self.mSelectionProbability = probability
        return

    def get_selection_probability(self):
        return self.mSelectionProbability

    def set_selected(self, isSelected):
        self.mSelected = isSelected
        return

    def get_selected(self):
        return self.mSelected

    def set_crossover(self, type):
        self.mCrossover = type
        return

    def get_crossover(self):
        return self.mCrossover

    def set_data(self, index, value):
        self.mData[index] = value
        return

    def get_data(self, index):
        return self.mData[index]

    # Gene des Chromosoms als String darstellen
    def toStr(self):
        array = [0]*self.mMaxLength
        for i in range(self.mMaxLength):
            array[i] = Chromosome.get_data(self,i)

        print(str(array))
        return str(array)

        # Gene des Chromosoms als Liste darstellen

    def toArray(self):
        array = [0] * self.mMaxLength
        for i in range(self.mMaxLength):
            array[i] = Chromosome.get_data(self,i)

        return array

