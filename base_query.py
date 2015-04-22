import pysam
import argparse
from collections import Counter
from math import exp
from scipy.stats import binom

def chernoff(base_record, baseq_record):
    
    m = float(sum([pow(10.0, (-q/10.0)) for q in baseq_record.values()]))
    sorted_bases = base_record.most_common(2)
    alt_count = sorted_bases[-1][1]
    alt_allele = sorted_bases[-1][0]
    total_count = sum(base_record.values())
    alt_freq = alt_count / float(total_count)
    d = (alt_count/m) - 1
    
    try:
        prob_se = pow((exp(d) / pow((1+d), (1+d))), m)
    except OverflowError:
        prob_se = 1e-15
    
    prob_germ = binom.cdf(alt_count, total_count, 0.5)

    return prob_se, alt_allele, alt_freq, prob_germ


def query_seq(samfile, chrom, start, end, filt, summary):
    
    fmt = "{base_call}\t{base_count}\t"
    fmt += "{mean_baseq}\t{mean_mapq}\t"
    fmt += "{fwd_strand}\t{rvs_strand}"
    
    fmt2 = "{alt_allele}\t{alt_freq}\t"
    fmt2 += "{se_pval}\t{germ_pval}"
    
    fmt3 = "{pos}\t{A}\t{T}\t{G}\t{C}"
    if summary:
        print fmt3.replace("{","").replace("}","")
    
    count = 0
    shift = 1
    if start == end:
        end += 1
        shift = 2

    for pileupcolumn in samfile.pileup(chrom, start-1, end):
        pos = pileupcolumn.pos
        if pos >= start-1 and pos <= end-shift:
            count +=1
            base_record = Counter()
            baseq_record = Counter()
            mapq_record = Counter()
            strand_record = Counter()
            if not summary:
                print ("\ncoverage at base %s = %s" %
                            (pos+1, pileupcolumn.n))
            for pileupread in pileupcolumn.pileups:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                baseq = pileupread.alignment.query_qualities[pileupread.query_position]
                mapq = pileupread.alignment.mapping_quality
                strand = pileupread.alignment.is_reverse
		second = "False"
		if pileupread.alignment.is_secondary:
			second = "True"
                
                if filt:
                    if mapq >= 10 and baseq >= 10 and pileupcolumn.n >= 5:
                        base_record[base]+=1
                        baseq_record[base]+=baseq
                        mapq_record[base]+=mapq
                        strand_record[base]+=strand
                        if not summary:
                            print ('\tbase = %s\tbase quality = %s\tmapping quality = %s' %
                                    (base, baseq, mapq))
                else:
                    base_record[base]+=1
                    baseq_record[base]+=baseq
                    mapq_record[base]+=mapq
                    strand_record[base]+=strand
                    if not summary:
                        print ('\tbase = %s\tbase quality = %s\tmapping quality = %s\tsecondary = %s' %
                                (base, baseq, mapq, second))
            if not summary:
                if filt:
                    if len(base_record) > 0:
                        print
                        print fmt.replace("{","").replace("}","")
                        for base_call in base_record:
                            base_count = base_record[base_call]
                            mean_baseq = baseq_record[base_call]/base_count
                            mean_mapq = mapq_record[base_call]/base_count
                            fwd_strand = base_record[base_call] - strand_record[base_call]
                            rvs_strand = strand_record[base_call]
                            print fmt.format(**locals())
                        se_pval, alt_allele, alt_freq, germ_pval = chernoff(base_record, baseq_record)
                        print
                        print fmt2.replace("{","").replace("}","")
                        print fmt2.format(**locals())
                else:
                    print
                    print fmt.replace("{","").replace("}","")
                    for base_call in base_record:
                        base_count = base_record[base_call]
                        mean_baseq = baseq_record[base_call]/base_count
                        mean_mapq = mapq_record[base_call]/base_count
                        fwd_strand = base_record[base_call] - strand_record[base_call]
                        rvs_strand = strand_record[base_call]
                        print fmt.format(**locals())
                    se_pval, alt_allele, alt_freq, germ_pval = chernoff(base_record, baseq_record)
                    print
                    print fmt2.replace("{","").replace("}","")
                    print fmt2.format(**locals())
            
            if summary:
                A=T=G=C = 0
                for base_call in base_record:
                    if base_call == "A":
                        A+=base_record[base_call]
                    elif base_call == "T":
                        T+=base_record[base_call]
                    elif base_call == "G":
                        G+=base_record[base_call]
                    elif base_call == "C":
                        C+=base_record[base_call]
                print fmt3.format(**locals())
                         
    if not summary:
        print
        print count, "positions found"
                         

def main():
    parser=argparse.ArgumentParser()

    parser.add_argument('--bam', help='input bam file', required=True)
    parser.add_argument('--chrom', help='query chromosome', required=True)
    parser.add_argument('--start', help='query start', required=True, type=int)
    parser.add_argument('--end', help='query end', required=True, type=int)
    parser.add_argument('--filter', help='flag to filter all regions', action="store_true",
                        default=False)
    parser.add_argument('--summary', help='only print summary stuff', action="store_true",
                        default=False)

    args=parser.parse_args()
    
    samfile = pysam.Samfile(args.bam)
    chrom = args.chrom
    start = args.start
    end = args.end
    filt = args.filter
    summary = args.summary
    
    query_seq(samfile, chrom, start, end, filt, summary)
    
if __name__ == "__main__":
    main()
