#!/usr/bin/env python

import sys, os
import math
from hmm import *

def find_state_labels(tfile, efile):
    labels = [[],[]]
    files = [tfile, efile]
    for i, file_path in enumerate(files):
        fp = open(file_path)
        for line in fp:
            arr = line.strip().split(' ')
            labels[i].append(arr[0])
        fp.close()

    #check if labels are consistent
    #first check the uniqueness and length of the arr
    for item in labels:
        if len(item) != len(set(item)):
            break
    else:
        
        first = labels[0]
        ll = len(first)
        for item in labels[1:]:
            if ll != len([i for i, j in zip(item, first) if i == j]):
                break
        else:
            return first

    #None is returned when the tranisition and emission files are not consistent
    #This should not happen; terminate if this happens
    return None

def find_outcomes(efile):
    outcomes = set([])
    fp = open(efile)
    for line in fp:
        arr = line.strip().split(' ')
        for item in arr[1:]:
            outcomes.add(item.split(':')[0])
    fp.close()
    return list(outcomes)

def create_matrix(matrix, file_path):
    fp = open(file_path)
    for line in fp:
        arr = line.strip().split(' ')
        for item in arr[1:]:
            k,v = item.split(':')
            matrix.add(arr[0], k, float(v))
    fp.close()


def main():

    dev_file_path = sys.argv[1]
    transition_file_path = sys.argv[2]
    emission_file_path = sys.argv[3]
    prior_file_path = sys.argv[4]

    #read prior probs from file
    prior = []
    fp = open(prior_file_path)
    for line in fp:
        a,b = line.strip().split(" ")
        prior.append((a,float(b)))
    fp.close()

    #find state labels
    state_labels = find_state_labels(transition_file_path, emission_file_path)
    assert state_labels
    assert [j for i,j in zip(prior, state_labels) if i[0] == j]
    
    #find outcomes
    outcome_labels = find_outcomes(emission_file_path)

    #read emission, transition matrices from file
    tmatrix = Transition(state_labels)
    ematrix = Emission(state_labels)
    create_matrix(tmatrix, transition_file_path)
    create_matrix(ematrix, emission_file_path)

    
    fp = open(dev_file_path)
    for line in fp: 
        outcomes = line.strip().split(' ')
        #Forward algorithm
        f = Forward(tmatrix, ematrix, prior)
        f.create_alpha_matrix(outcomes)
        for i in xrange(len(state_labels)):
            if i == 0:
                sigma = f.alpha[i][len(outcomes)-1]
            else:
                sigma = log_sum(sigma, f.alpha[i][len(outcomes)-1])
        
        print sigma
        
    fp.close()

if __name__=='__main__':
    main()
