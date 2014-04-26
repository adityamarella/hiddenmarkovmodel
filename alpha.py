#!/usr/bin/env python
"""
Implementation of Foward, Backward algorithm for Hidden Markov Model evaluation and learning. 
 - alpha, beta matrices are incrementally built using dynamic programming
 - one quick way to verify the correctness of the implementation is to see\
         if the values in alpha and beta matrics are the same.
"""

import sys, os


def create_matrix(matrix, trans):
    fp = open(trans)
    for line in fp:
        arr = line.strip().split(' ')
        for item in arr[1:]:
            k,v = item.split(':')
            matrix.add(arr[0], k, float(v))
    fp.close()

class Matrix(object):
"""
    In addition to maintaining the state transitions as an adjacency list,\
    this class also assigns indices to the states. The main purpose of this\
    class though is to provide a uniform implementation agnostic interface for\
    Transition, Emission matrices. For ex: the state transitions could be maintained\
    in an adjacency matrix depending on whether the state transitions are sparse or dense
"""

    def __init__(self):
        self.matrix = {}
        #maintained in the order of insertion
        self.labels = []

    def add(self, frm, to, prob):
        self.matrix.setdefault(frm, []).append((to, prob))
        self.labels.append(frm)

    def index(self, label):
        return self.labels.index(label)

    def get(self, i, j):
        a = self.labels[i]
        b = self.labels[j]

        arr = self.matrix[a]
        for item in enumerate(arr):
            if item[0] == b:
                return item[1]
        return 0

class Forward(object):
"""
@input: transition, emission matrices and prior probabilities of the initial state
@output: alpha matrix

Note: The transition and emission matrices are instances of Matrix class defined above. 
"""

    def __init__(self, t, e, prior):
        self.t = t
        self.e = e
        self.prior = prior

    def create_alpha_matrix(self):
        

def main():

    dev = sys.argv[1]
    trans = sys.argv[2]
    emit = sys.argv[3]
    prior_file = sys.argv[4]

    #read prior probs from file
    prior = []
    fp = open(prior_file)
    for line in fp:
        a,b = line.strip().split(" ")
        prior.append((a,float(b)))
    fp.close()

    #read emission, transition matrices from file
    tmatrix = Matrix()
    ematrix = Matrix()
    create_matrix(tmatrix, trans)
    create_matrix(ematrix, emit)

    #Forward algorithm
    f = Forward(tmatrix, ematrix, prior)
    f.create_alpha_matrix()

    #Backward algorithm

if __name__=='__main__':
    main()
