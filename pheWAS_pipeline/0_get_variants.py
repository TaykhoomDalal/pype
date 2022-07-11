import argparse
import pandas as pd

def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Pass the input variant file')
    parser.add_argument('-i', '--input', help='Url of file', required=True)
    parser.add_argument('-p', '--p_val', help='P-value field name', required=True)
    parser.add_argument('-c', '--chrom', help='Chromosome field name', required=True)
    parser.add_argument('-r', '--rsid', help='RSID field name', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)

    args = parser.parse_args()
    input = args.input # downloaded input files from FUMA https://fuma.ctglab.nl/browse/488 and 487
    output = args.output
    p_val = args.p_val
    chrom = args.chrom
    rsid = args.rsid

    # load the input file
    variants = pd.read_csv(input, sep="\t", usecols = [chrom, rsid, p_val]).sort_values(by=p_val)

    # output fields
    cols = {chrom: 'CHR', rsid: 'rsID'}

    # only retain the variants with a P_BOLT_LMM_INF value smaller than or equal to 5e-8 also drop the pvalue column, 
    # rename the columns, and reorder them to match our other files (not necessary, but makes it easier to compare)
    variants = variants[variants[p_val] <= 5e-8].rename(columns=cols).drop(p_val, axis=1)

    # IDK if this makes sense - Update: it doesn't, but im saving this for now
    #variants.rsID = variants.rsID.str.replace(r'.*:', '', regex = True).str.replace(r"(\d+)_.*",r"rs\1", regex = True)

    # write the variants to the output file
    variants.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()