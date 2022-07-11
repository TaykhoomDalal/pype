import os
import argparse
import pandas as pd
import subprocess

def main():
    parser = argparse.ArgumentParser(description='Extract the data from the input file')
    parser.add_argument('-v', '--variants', help='List of variants', action='append', required=True)
    parser.add_argument('-s', '--samples', help='List of samples', action='append', required=True)
    parser.add_argument('-b', '--bfiles', help = 'Directory of bfiles', required=True)
    parser.add_argument('-g', '--geno', help = 'Directory of extracted genotypes', required=True)
    parser.add_argument('-k', '--keep', help ='Whether the samples are in / out of sample', action='store_true')
    parser.add_argument('-t', '--threads', help = 'Number of threads', required=True, type=int)
    parser.add_argument('-m', '--memory', help = 'Memory to use (MB)', required=True, type=int)
    args = parser.parse_args()

    # Parse the command line arguments
    variant_list = args.variants
    sample_IDs = args.samples
    bfiles_dir = args.bfiles
    geno_dir = args.geno
    keep = args.keep
    threads = args.threads
    memory = args.memory
    
    operation_on_samples = ''
    if keep:
        operation_on_samples = '--keep'
    else:
        operation_on_samples = '--remove'

    os.system('module load plink2')

    for input_file, sample_file in zip(variant_list, sample_IDs):

        # Determine the age predictor being used
        age_Predictor = input_file.split('/')[-1].split('.')[0].split('_')[0]

        # read the file of variants
        variants = pd.read_csv(input_file, sep='\t', header=0)
        chr = set(variants['CHR'].values)
        rsIDs = variants['rsID'].values

        temp_file = geno_dir + '/' + age_Predictor + '_temp_rsIDs_file.txt'

        with open(temp_file, 'w+') as f:
            f.write('\n'.join(rsIDs))

        for i in chr:

            # need to make sure that the directory used only contains bfiles and no other file, and that they all have the same prefix (before extension)
            file = ""
            for filename in os.listdir(bfiles_dir):
                if filename.startswith('chr' + str(i) + '_'): # every file in the directory starts with 'chr#_', where # is the chromosome
                    file = filename.split('.')[0] # get the file name without the extension
            
            command = " ".join(["plink2", "--bfile", bfiles_dir + '/' + file, 
                            "--extract", temp_file, 
                            operation_on_samples, sample_file,
                            "--threads", str(threads),
                            "--memory", str(memory),
                            "--make-bed", 
                            "--out", geno_dir + '/' + age_Predictor + '_' + file])

            os.system('module load plink2 && ' + command)
            print('chr' + str(i) + ' done')

        os.remove(temp_file)
       

if __name__ == '__main__':
    main()