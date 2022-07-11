import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Create raw geno files')
    parser.add_argument('-g', '--geno_dir', help='Geno directory', required=True)
    parser.add_argument('-r', '--raw_dir', help='Raw directory', required=True)
    
    args = parser.parse_args()
    geno_dir = args.geno_dir
    raw_dir = args.raw_dir

    # get a list of the files in the geno_dir
    geno_files = [os.path.join(geno_dir, f.split('.bed')[0]) for f in os.listdir(geno_dir) if f.endswith('.bed')]

    for geno_file in geno_files:
        
        geno_file_prefix = geno_file.split('/')[-1].split('.')[0]
        
        # create a raw file for the genotype data
        command = " ".join(["plink2", "--export", 'A', 
                            "--bfile", geno_file, 
                            "--out", raw_dir + '/' + geno_file_prefix])

        os.system('module load plink2 && ' + command)

if __name__ == '__main__':
    main()