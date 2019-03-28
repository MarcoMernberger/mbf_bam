import mbf_bam
import pysam
from _common import calculate_duplicate_distribution

import time

def do_time(what, callable):
    start = time.time()
    c = callable()
    stop = time.time()
    print("%s\t%.2fs" % (what ,stop - start))
    return c

def pure_python():
    s = pysam.Samfile('../../results/lanes/HEK293t_scrambled_siRNA_ERR2223563/aligned/STAR/Homo_sapiens_94/HEK293t_scrambled_siRNA_ERR2223563.bam')
    c = 0
    for ref in s.references:
        ll = s.get_reference_length(ref)
        for r in s.fetch(ref, 0, ll):
            c += 1
    return c



#print(mbf_bam.calculate_duplicate_distribution('../mbf_align/tests/sample_data/rnaseq_spliced.bam'))

#20.15s
#single core
#do_time('rustsc', lambda:  mbf_bam.calculate_duplicate_distribution(
    #'../../results/lanes/HEK293t_scrambled_siRNA_ERR2223563/aligned/STAR/Homo_sapiens_94/HEK293t_scrambled_siRNA_ERR2223563.bam'))


#21.94s
#do_time('python', pure_python)

fn = '../../results/lanes/HEK293t_scrambled_siRNA_ERR2223563/aligned/STAR/Homo_sapiens_94/HEK293t_scrambled_siRNA_ERR2223563.bam'
#fn = '../mbf_align/tests/sample_data/rnaseq_spliced.bam'
#3.97


def in_cython():
    s = pysam.Samfile(fn)
    return calculate_duplicate_distribution(s)

python_value = {
    1: 4370276, 
    2: 1116266, 
    3: 421485, 
    4: 207372, 
    5: 118289, 
    6: 74848, 
    7: 51486, 
    8: 37223, 
    9: 27946, 
}
#python_value = do_time('cython', in_cython)

rust_value = do_time('rustmc', lambda:  mbf_bam.calculate_duplicate_distribution(
    fn))

if python_value != rust_value:
    import pprint
    print("Discrepancy")
    print('python')
    for x in list(range(20)) + [9999999]:
        p =  python_value.get(x,-1)
        r = rust_value.get(x,-1)
        print(x,p, r, r==p , sep="\t")

def sum(x):
    s = 0
    if 9999999 in x:
        del x[9999999]
    for k,v in x.items():
        s += k * v
    return s
print('python sum', sum(python_value))
print('rust sum', sum(rust_value))

