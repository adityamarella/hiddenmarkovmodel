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

    def setk(self, state_idx, outcome_label, prob):
        self.emission[state_idx][outcome_label] = prob
    
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

    def setk(self, i, j, prob):
        self.matrix[i][j] = prob


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
                    value = self.alpha[j][t-1] + math.log(self.transition.getk(j, i))
                    if j == 0:
                        sigma = value
                    else:
                        sigma = log_sum(sigma, value)

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
                    value = self.beta[j][t+1] +\
                            math.log(self.transition.getk(i, j)) +\
                            math.log(self.emission.getk(j, O[t+1]))
                    if j == 0:
                        sigma = value 
                    else:
                        sigma = log_sum(sigma, value)
                self.beta[i][t] = sigma

class Viterbi(object):
    """
    """

    def __init__(self, t, e, prior):
        self.transition = t
        self.emission = e
        self.prior = prior

    def most_likely_transition_path(self, O):
        """
        """
       
        states = self.transition.states
        self.VP = [[0]*len(O) for _ in states]
        self.q = [['']*len(O) for _ in states]

        #initialization
        for i in xrange(len(states)):
            self.VP[i][0] = math.log(self.prior[i][1]) + math.log(self.emission.getk(i, O[0]))
            self.q[i][0] = states[i]

        for t in xrange(1, len(O)):
            for i in xrange(len(states)):
                k = 0
                maxvp = float("-inf")
                for j in xrange(0, len(states)):
                    value = self.VP[j][t-1] + math.log(self.transition.getk(j, i)) + math.log(self.emission.getk(i, O[t]))
                    if maxvp < value:
                        maxvp = value
                        k = j
                self.VP[i][t] = maxvp
                self.q[i][t] = self.q[k][t-1] + '-' + states[i]
        
        k = 0
        maxvp = float("-inf")
        T = len(O) - 1
        for i in xrange(len(states)):
            if maxvp < self.VP[i][T]:
                maxvp = self.VP[i][T]
                k = i
        return self.q[k][T].split("-")

class BaumWelch(object):
    """
    Estimation-Maximization method
    """

    def __init__(self, t, e, prior, training_sentences):
        """
        initialized to some random values sampled from a uniform distribution
        """
        self.transition = t
        self.emission = e
        self.prior = prior

        self.training = training_sentences

    def run(self):
        prevll = float("-inf")
        for i in xrange(20):
            ll = self.avgll()
            if ll - prevll < 0.1:
                break
            self.maximize()
            self.estimate()
            print ll
            prevll = ll

    def avgll(self):
        states = self.transition.states
        M = len(self.training)
        s = 0
        for m in xrange(M):
            O = self.training[m]
            f = Forward(self.transition, self.emission, self.prior)
            f.create_alpha_matrix(O)
            sigma = 0
            for i in xrange(len(states)):
                if i == 0:
                    sigma = f.alpha[i][len(O)-1]
                else:
                    sigma = log_sum(sigma, f.alpha[i][len(O)-1])
            s = s + sigma
        return s/M

    def maximize(self):

        states = self.transition.states

        self.gamma = [[] for _ in self.training]
        self.xi = [[] for _ in self.training]

        for m,O in enumerate(self.training):
            
            self.gamma[m] = [[0]*len(O) for _ in xrange(len(states))]
            self.xi[m] = [ [[0]*len(O) for _ in xrange(len(states))] for _ in xrange(len(states))]

            self.f = Forward(self.transition, self.emission, self.prior)
            self.f.create_alpha_matrix(O)
            
            self.b = Backward(self.transition, self.emission, self.prior)
            self.b.create_beta_matrix(O)
            
            for t in xrange(len(O)):
                for i in xrange(len(states)):
                    for j in xrange(len(states)):
                        value = self.f.alpha[j][t] + self.b.beta[j][t]
                        if j == 0:
                            sigma = value 
                        else:
                            sigma = log_sum(sigma, value)
                        
                    for j in xrange(len(states)):
                        if t!=len(O)-1:
                            self.xi[m][i][j][t] = self.f.alpha[j][t]\
                                    + self.b.beta[j][t]\
                                    + math.log(self.transition.getk(i,j))\
                                    + math.log(self.emission.getk(j,O[t+1]))

                            self.xi[m][i][j][t] = self.xi[m][i][j][t] - sigma
                                 
                    self.gamma[m][i][t] = self.f.alpha[i][t] + self.b.beta[i][t]

                    self.gamma[m][i][t] = math.exp(self.gamma[m][i][t] - sigma)

    def estimate(self):
        states = self.transition.states
        N = len(states)
        M = len(self.training)

        for i in xrange(N):
            for m in xrange(M):
                if m == 0:
                    sigma = self.gamma[m][i][0]
                else:
                    sigma = sigma + self.gamma[m][i][0]
            self.prior[i] = (states[i], math.exp(sigma - math.log(M)))
        
        #new transition probs
        for i in xrange(N):
            denominator = 0
            for j in xrange(N):
                sigma = 0
                for m in xrange(M):
                    O = self.training[m]
                    for t in xrange(len(O)):
                        if sigma == 0:
                            sigma = self.xi[m][i][j][t]
                        else:
                            sigma = log_sum(sigma, self.xi[m][i][j][t])

                self.transition.setk(i, j, sigma)
                if denominator == 0:
                    denominator = sigma
                else:
                    denominator = log_sum(denominator, sigma)

            for j in xrange(N):
                self.transition.setk(i,j, math.exp(self.transition.getk(i,j) - denominator))

        #new emission probs
        outcomes = self.emission.emission[0].keys()
        for i in xrange(N):
            for k in xrange(len(outcomes)):
                sigma = 0
                denominator = 0
                for m in xrange(M):
                    O = self.training[m]
                    for t in xrange(len(O)):
                        if outcomes[k] == O[t]:
                            if sigma == 0:
                                sigma = self.gamma[m][i][t]
                            else:
                                sigma = log_sum(sigma, self.gamma[m][i][t])

                        if denominator == 0:
                            denominator = self.gamma[m][i][t]
                        else:
                            denominator = log_sum(denominator, self.gamma[m][i][t])

                self.emission.setk(i, k, math.exp(sigma - denominator))

