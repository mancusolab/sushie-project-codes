import argparse as ap
import gzip
import logging
import os
import sys

from intervaltree import Interval, IntervalTree

__log__ = "make.annot.encode"

def get_logger(name):
    logger = logging.getLogger(name)
    if not logger.handlers:
        # Prevent logging from propagating to the root logger
        logger.propagate = 0
        console = logging.StreamHandler()
        logger.addHandler(console)

        log_format = "[%(asctime)s - %(levelname)s] %(message)s"
        date_format = "%Y-%m-%d %H:%M:%S"
        formatter = logging.Formatter(fmt=log_format, datefmt=date_format)
        console.setFormatter(formatter)

    return logger

def data_reducer(lhs, rhs):
    if rhs != "" and rhs != None:
        return f"{lhs},rhs"
    else:
        return lhs


def main(args):
    parser = ap.ArgumentParser()
    parser.add_argument("annot")
    parser.add_argument("snps")
    parser.add_argument("chrom")
    parser.add_argument("-s", "--split", default="")
    parser.add_argument("-l", "--long", action="store_true", default=False)
    parser.add_argument("-f", "--filter", action="store_true", default=False)
    parser.add_argument("-d", "--drop", nargs="+", type=str, default=None)
    parser.add_argument("-o", "--output", type=ap.FileType("w"), default=sys.stdout)

    args = parser.parse_args(args)

    log = get_logger(__log__)
    log.setLevel(logging.INFO)
    anno_chr = f"chr{args.chrom}"

    # create interval tree to manage annotations
    log.info(f"Constructing interval tree for annotations on chromosome {args.chrom}.")
    tree = IntervalTree()

    # chr1	10033	10250	EH38D4327497	EH38E2776516	pELS
    with gzip.open(args.annot, "rt") as annot:
        for line in annot:
            row = line.split()
            chrom = row[0]

            if "chr" not in chrom:
                chrom = "chr" + chrom

            if chrom != anno_chr:
                continue

            start = int(row[1])
            end = int(row[2])

            if args.split != "":
                tmp_annots = row[-1].split(args.split)
            else:
                tmp_annots = [row[-1]]

            if args.drop is not None:
                tmp_annots = [x for x in tmp_annots if x not in args.drop]

            if len(tmp_annots) == 0:
                continue

            l_annots = ",".join(tmp_annots)
            tree[start:end] = l_annots

    # tree.merge_overlaps(data_reducer=data_reducer)
    log.info(f"Completed constructing interval tree for annotations on chromosome {args.chrom}.")

    log.info("Annotating variant data")
    # CHR SNP CM BP ALT REF 
    with open(args.snps, "rt") as snps:
        args.output.write("\t".join(["CHR", "SNP", "CM", "BP",  "ALT", "REF", "ANNOT"]) + os.linesep)
        for line in snps:
            row = line.split()
            chrom = row[0]
            if chrom != args.chrom:
                continue

            bp = int(row[3])
            l_annots = ",".join(entry.data for entry in tree[bp])
            row = row + [l_annots]

            if args.filter and len(l_annots) == 0:
                continue

            if args.long:
                tmp_split = row[-1].split(",")
                for idx in range(len(tmp_split)):
                    row[-1] = tmp_split[idx]
                    args.output.write("\t".join(row) + os.linesep)
            else:
                args.output.write("\t".join(row) + os.linesep)

    log.info(f"Completed annotating variant data on chromosome {args.chrom}.")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
