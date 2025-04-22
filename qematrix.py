import numpy as np
import math as m

class QE:
    def __init__(self):
        self.j = complex(0, 1)
        self.TE = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.T1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
        self.T2 = np.array([[0, -self.j, 0], [self.j, 0, 0], [0, 0, 0]])
        self.T3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]])
        self.T4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
        self.T5 = np.array([[0, 0, -self.j], [0, 0, 0], [self.j, 0, 0]])
        self.T6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
        self.T7 = np.array([[0, 0, 0], [0, 0, -self.j], [0, self.j, 0]])
        self.T8 = np.array([[1 / m.sqrt(3), 0, 0], [0, 1 / m.sqrt(3), 0], [0, 0, -2 / m.sqrt(3)]])
        self.Tilist = [self.T1, self.T2, self.T3, self.T4, self.T5, self.T6, self.T7, self.T8]

    def sumMatrix(self, matrixlist):
        """
        Sums all matrices in the given list and returns the resulting matrix.
        """
        for i in range(len(matrixlist)):
            matrixlist[0] += matrixlist[i]
        return matrixlist[0]

    def densityMatrix(self, Ai, Bi, Cij):
        """

        The density matrix is computed using the formula:
        $$
        \\rho = \\frac{1}{9} \sum_{i,j} T_i \otimes T_j + \\frac{1}{3} \sum_i a_i T_i \otimes T_E + \\frac{1}{3} \sum_i b_i T_E \otimes T_i + \sum_{i,j} c_{ij} T_i \otimes T_j
        $$
        """
        term1 = 1 / 9 * np.kron(self.TE, self.TE).astype(complex)
        for i, ai in enumerate(Ai):
            term1 += ai * np.kron(self.TE, self.Tilist[i])
        for i, bi in enumerate(Bi):
            term1 += bi * np.kron(self.Tilist[i], self.TE)
        for i, ci in enumerate(Cij):
            for j, cj in enumerate(ci):
                term1 += cj * np.kron(self.Tilist[i], self.Tilist[j])
        return term1

    def densityMatrixA(self, Ai):
        term1 = 1 / 3 * self.TE.astype(complex)
        for i, ai in enumerate(Ai):
            term1 += 1 / 3 * ai * self.Tilist[i]
        return term1

    def densityMatrixB(self, Bi):
        term1 = 1 / 3 * self.TE.astype(complex)
        for i, bi in enumerate(Bi):
            term1 += 1 / 3 * bi * self.Tilist[i]
        return term1

    def C2(self, Ai, Bi, Cij):
        dm = self.densityMatrix(Ai, Bi, Cij)
        dma = self.densityMatrixA(Ai)
        dmb = self.densityMatrixB(Bi)
        return 2 * max(
            0,
            np.real(np.trace(np.dot(dm, dm)) - np.trace(np.dot(dma, dma))),
            np.real(np.trace(np.dot(dm, dm)) - np.trace(np.dot(dmb, dmb)))
        )
    
    def mean_veclist(self,veclist):
        """
        ## Calculate the Mean of a List of Vectors

        This function computes the mean of a list of vectors by summing all the vectors 
        and dividing by the number of vectors.

        ### Args:
        - **veclist** (*list*): A list of vectors. Each vector should support the `copy`, `+`, 
          and `/` operations, and all vectors in the list should have the same dimensions.

        ### Returns:
        - The mean vector, which has the same dimensions as the input vectors.

        ### Raises:
        - **ZeroDivisionError**: If the input list `veclist` is empty.
        - **AttributeError**: If the elements of `veclist` do not support the required operations.
        """
        ans=veclist[0].copy()
        for vec in veclist[1:]:
            ans+=vec
        return ans/len(veclist)
    def mean_matlist(self, matlist):
        ans=matlist[0].copy()
        for mat in matlist[1:]:
            ans+=mat
        return ans/len(matlist)
    def mean_event(self, Ai_event, Bi_event, Cij_event, nexp, eventperjob):
        Ai_list=[]
        Bi_list=[]
        Cij_list=[]
        if 10<nexp<100:
            Ai_event_temp = np.concatenate([Ai_event] * 10)
            Bi_event_temp = np.concatenate([Bi_event] * 10)
            Cij_event_temp = np.concatenate([Cij_event] * 10)
            
            # Shuffle the concatenated arrays
            shuffle_indices = np.random.permutation(len(Ai_event_temp))
            Ai_event_temp = Ai_event_temp[shuffle_indices]
            Bi_event_temp = Bi_event_temp[shuffle_indices]
            Cij_event_temp = Cij_event_temp[shuffle_indices]
            nexp = nexp * 10
            for i in range(nexp):
                start_idex = i * eventperjob
                end_index = (i+1) * eventperjob
                chuck_Ai = Ai_event_temp[start_idex:end_index].copy()
                chuck_Bi = Bi_event_temp[start_idex:end_index].copy()
                chuck_Cij = Cij_event_temp[start_idex:end_index].copy()
                Ai_list.append(self.mean_veclist(chuck_Ai))
                Bi_list.append(self.mean_veclist(chuck_Bi))
                Cij_list.append(self.mean_matlist(chuck_Cij))
        elif nexp<=10:
            Ai_event_temp = np.concatenate([Ai_event] * 100)
            Bi_event_temp = np.concatenate([Bi_event] * 100)
            Cij_event_temp = np.concatenate([Cij_event] * 100)
            # Shuffle the concatenated arrays
            shuffle_indices = np.random.permutation(len(Ai_event_temp))
            Ai_event_temp = Ai_event_temp[shuffle_indices]
            Bi_event_temp = Bi_event_temp[shuffle_indices]
            Cij_event_temp = Cij_event_temp[shuffle_indices]
            nexp = nexp * 100
            for i in range(nexp):
                start_idex = i * eventperjob
                end_index = (i+1) * eventperjob
                chuck_Ai = Ai_event_temp[start_idex:end_index].copy()
                chuck_Bi = Bi_event_temp[start_idex:end_index].copy()
                chuck_Cij = Cij_event_temp[start_idex:end_index].copy()
                Ai_list.append(self.mean_veclist(chuck_Ai))
                Bi_list.append(self.mean_veclist(chuck_Bi))
                Cij_list.append(self.mean_matlist(chuck_Cij))
        else:
            for i in range(nexp):
                start_idex = i * eventperjob
                end_index = (i+1) * eventperjob
                chuck_Ai = Ai_event[start_idex:end_index].copy()
                chuck_Bi = Bi_event[start_idex:end_index].copy()
                chuck_Cij = Cij_event[start_idex:end_index].copy()
                Ai_list.append(self.mean_veclist(chuck_Ai))
                Bi_list.append(self.mean_veclist(chuck_Bi))
                Cij_list.append(self.mean_matlist(chuck_Cij))
        return Ai_list, Bi_list, Cij_list
