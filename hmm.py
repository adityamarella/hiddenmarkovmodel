#!/usr/bin/env python
"""
Implementation of Foward, Backward algorithm for Hidden Markov Model evaluation and learning. 
 - alpha, beta matrices are incrementally built using dynamic programming
 - one quick way to verify the correctness of the implementation is to see\
         if the values in alpha and beta matrics are the same.
"""

import sys, os
import math
from math import log1p,exp

def log_sum(left,right):
	if right < left:
		return left + log1p(exp(right - left))
	elif left < right:
		return right + log1p(exp(left - right));
	else:
		return left + log1p(1)

class Emission(object):

    def __init__(self, state_labels):
        self.states =   {label:i for i,label in enumerate(state_labels)}
        self.emission = [{} for i in state_labels]

    def add(self, state_label, outcome_label, outcome_prob):
        try:
            sidx = self.states[state_label]
        except KeyError:
            return

        self.emission[sidx][outcome_label] = outcome_prob

    def getk(self, state_idx, outcome_label):
        return self.emission[state_idx][outcome_label]
    
class Transition(object):
    """
    In addition to maintaining the state transitions as an adjacency matrix,\
    this class also assigns indices to the states. The main purpose of this\
    class though is to provide a uniform implementation agnostic interface for\
    Transition matrix. For ex: the state transitions could be maintained\
    in an adjacency matrix depending on whether the state transitions are sparse or dense

    """

    def __init__(self, state_labels):
        self.states = state_labels
        self.labels = {label:i for i,label in enumerate(state_labels)}
        self.matrix = [[0]*len(state_labels) for _ in state_labels]

    def add(self, frm, to, prob):
        self.matrix[self.labels[frm]][self.labels[to]] = prob

    def getk(self, i, j):
        return self.matrix[i][j]


class Forward(object):
    """
    @input: transition, emission matrices and prior probabilities of the initial state
    @output: alpha matrix

    Note: The transition and emission matrices are instances of Matrix class defined above. 
    """

    def __init__(self, t, e, prior):
        self.transition = t
        self.emission = e
        self.prior = prior

    def create_alpha_matrix(self, O):
        """
        Computes the probability that the seq of outcomes O_{0}, O_{2}..., O_{t-1} is generated at time t\
            and that the system is in state s_{i} at time t.
        @params t(time or outcome), i(state i)
        """
       
        states = self.transition.states
        self.alpha = [[0]*len(O) for _ in states]
        #initialization
        for i in xrange(len(states)):
            self.alpha[i][0] = math.log(self.prior[i][1]) + math.log(self.emission.getk(i, O[0]))

        for t in xrange(1, len(O)):
            for i in xrange(len(states)):
                for j in xrange(0, len(states)):
                    if j == 0:
                        sigma = self.alpha[j][t-1] + math.log(self.transition.getk(j, i))
                    else:
                        sigma = log_sum(sigma, self.alpha[j][t-1] + math.log(self.transition.getk(j, i)) )
                self.alpha[i][t] = math.log(self.emission.getk(i, O[t])) + sigma

class Backward(object):
    """
    @input: transition, emission matrices and prior probabilities of the initial state
    @output: beta matrix

    Note:  
    """
    def __init__(self, t, e, prior):
        self.transition = t
        self.emission = e
        self.prior = prior

    def create_beta_matrix(self, O):
        """
        Computes the probability of generating outcomes O_{t+1}, O_{t+2}..., O_{T} after t\
            provided the system is in state s_{i} at time t.
        @params t(time or outcome), i(state i)
        """

        states = self.transition.states
        self.beta = [[0]*len(O) for _ in states]
        T = len(O) - 1
        #initialization
        for i in xrange(len(states)):
            self.beta[i][T] = 0

        for t in xrange(T-1, -1, -1):
            for i in xrange(len(states)):
                for j in xrange(0, len(states)):
                    if j == 0:
                        sigma = self.beta[j][t+1] + math.log(self.transition.getk(i, j)) + math.log(self.emission.getk(j, O[t+1]))
                    else:
                        sigma = log_sum(sigma, self.beta[j][t+1] + math.log(self.transition.getk(i, j)) + math.log(self.emission.getk(j, O[t+1])) )
                self.beta[i][t] = sigma

