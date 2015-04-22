from pymix import mixture
import argparse
from random import random
import numpy as np
import matplotlib.pyplot as plt


def plot_histogram(labels, allele_freq):
    allele_freq = np.array(allele_freq)
    color = ['b','r','g','y']    
    
    label_alphabet = np.unique(labels)

    if len(label_alphabet) > 4:
        print "XXX Not enough colors specified to plot clustering. XXX"
        raise RuntimeError


    cind= []
    for u in label_alphabet:
        cind.append(np.where(labels==u)[0])

    plt.figure()
    for i,c in enumerate(cind):
        plt.hist(allele_freq[c])
    
    plt.title('Cluster Plot')
    plt.xlabel('Allele Frequency')
    plt.ylabel('Frequency')
    plt.savefig('Cluster_Hist.jpeg')   
    


def mixture_mode(allele_freq, max_components):
    data = mixture.DataSet()
    data.fromList(allele_freq)
    
    distributions = []
    for i in xrange(max_components):
        mean = random()
        std = random()
        distributions.append(mixture.NormalDistribution(mean, std))

    total_models = []
    for i in xrange(max_components):
        weights = list(np.repeat(1.0/(i+1), i+1))
        components = distributions[0:i+1] 
        model = mixture.MixtureModel(i+1, weights, components)
    
        model.EM(data, 40, 0.1)
        print
        print model
        print '------------------------------'
        total_models.append(model)
    
    model_selections = mixture.modelSelection(data, total_models)
    best_model = total_models[model_selections[1].index(min(model_selections[1]))]
    labels = best_model.classify(data, silent=1)
    
    return best_model, labels, allele_freq


def read_vcf(vcf_file):
    allele_freq = []
    for line in vcf_file:
        if line[0] != '#':
            line = line.strip().split('\t')
            AF = line[8].split(';')[0].split('=')[1]
            allele_freq.append(float(AF))
    return allele_freq

def main():
    '''
    main function
    takes arguments from the cammandline
    '''
    parser=argparse.ArgumentParser()

    parser.add_argument('--vcf', help='vcf file')
    parser.add_argument('--max', help='max number of sub-clones', type=int, default=10)
    args=parser.parse_args()
    
    vcf_file = open(args.vcf, mode='r')
    max_components = args.max
    allele_freq = read_vcf(vcf_file)
    best_model, labels, allele_freqs = mixture_mode(allele_freq, max_components)
    print '-----------------------'
    print 'BEST MODEL:'
    print best_model
    plot_histogram(labels, allele_freq)
    
if __name__ == '__main__':
    main()