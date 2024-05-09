#!/usr/bin/env python

# import required libraries
import argparse
import sys
import subprocess
from pathlib import Path
import pandas as pd

# parse args
parser = argparse.ArgumentParser()

parser.add_argument("-q", "--query_file", help="sequencing read file, can be gzipped")
parser.add_argument("-f", "--fasta_db_file", help="gene of interest diamond database file")
parser.add_argument("-d", "--db_file", help="outgroup diamond database file")
parser.add_argument("-V", "--version", action="store_true", help="show script version and exit")
parser.add_argument("-k", "--max_target_seqs", default=500, help="number of DIAMOND hits to consider before sorting on score, default 500")
parser.add_argument("-t", "--threads", default=1, help="number of threads to use for DIAMOND, default 1")


# parser.add_argument("-o", "--out_file", help="Define output file")

args = parser.parse_args()

# check existence of input files
if args.version:
    print("python script to run diamond blastx, avoiding the max-target-seqs problem. Version 0.1")
    sys.exit(0)

if not args.query_file:
    print("input query (read) file is required, specify with '-q' or '--query_file'")
    sys.exit(0)
query = Path(args.query_file)

if not args.fasta_db_file:
    print("outgroup DIAMOND DB file is required, specify with '-d' or '--db_file'")
    sys.exit(0)
gene_db = Path(args.fasta_db_file)

if not args.db_file:
    print("outgroup DIAMOND DB file is required, specify with '-d' or '--db_file'")
    sys.exit(0)
og_db = Path(args.db_file)

if not query.is_file():
    print("provided input query file does not exist")
    sys.exit(0)

if not gene_db.is_file():
    print("provided gene-db fasta file does not exist")
    sys.exit(0)

if not og_db.is_file():
    print("provided outgroup DIAMOND DB file does not exist")
    sys.exit(0)

# check that OG_DB is diamond db (extension)
# check that gene DB is fasta (first char)

# defining the output directory and files based on query + db
extensions = "".join(query.suffixes)
query_filename = query.name
query_base = str(query_filename).removesuffix(extensions)

outdir = query_base+"__"+gene_db.stem+"__read_search"
out = query_base+"__"+gene_db.stem+"__result"
out_og = query_base+"__"+og_db.stem+"__result"

gene_dmnd_db = outdir+"/"+gene_db.stem+"_dmnd"
unfiltered_1st = outdir+"/"+out
filtered_1st  = outdir+"/"+out+"_clean"
unfiltered_2nd = outdir+"/"+out_og
filtered_2nd  = outdir+"/"+out_og+"_clean"
bsr_file = outdir+"/"+out+"_bsr"

subprocess.run(["mkdir", outdir], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# defining the options for DIAMOND not specified by script args
max_target = str(args.max_target_seqs)
blast_col=["qseqid", "sseqid", "pident", "qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "score"]
scoring = ["--matrix", "blosum45", "--gapopen", "14", "--gapextend", "2", "--masking", "0", "--min-score", "10", "--comp-based-stats", "0"]
performance = ["-b", "12", "-c", "1"]

# run 1st diamond search
subprocess.run(["diamond", "makedb", "--in", gene_db, "-d", gene_dmnd_db, "-p", str(args.threads)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

run_diamond_1 = ["diamond", "blastx", "-q", query.resolve(), "-d", gene_dmnd_db, "-o", unfiltered_1st, "-p", str(args.threads), "-k", max_target, "--outfmt", "6", *blast_col, *scoring, *performance]
subprocess.run(run_diamond_1, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# dereplicate hits, keeping highest scoring one
hits_df = pd.read_csv(unfiltered_1st, sep="\t", names=blast_col)
derep_df = hits_df.sort_values('score', ascending=False).drop_duplicates(['qseqid'])
derep_df.to_csv(filtered_1st, sep="\t", index=False)
subprocess.run(["rm", unfiltered_1st])

# retrieve reads
subprocess.run(["blast_based_read_lookup_new.pl", filtered_1st, query.resolve(), filtered_1st+".faa"], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# run 2nd diamond search
run_diamond_og = ["diamond", "blastx", "-q", filtered_1st+".faa", "-d", og_db.resolve(), "-o", unfiltered_2nd, "-p", str(args.threads), "-k", max_target, "--outfmt", "6", *blast_col, *scoring, *performance]
subprocess.run(run_diamond_og, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

# dereplicate hits, keeping highest scoring one
hits_df_og = pd.read_csv(unfiltered_2nd, sep="\t", names=blast_col)
derep_df_og = hits_df_og.sort_values('score', ascending=False).drop_duplicates(['qseqid'])
derep_df_og.to_csv(filtered_2nd, sep="\t", index=False)
subprocess.run(["rm", unfiltered_2nd])

# calculate bsr
db_score = derep_df[["qseqid", "sseqid", "score", "pident"]].rename(columns={'sseqid': 'sseqid_db', 'score': 'score_db', 'pident': 'pident_db'})
og_score = derep_df_og[["qseqid", "sseqid", "score", "pident"]].rename(columns={'sseqid': 'sseqid_og', 'score': 'score_og', 'pident': 'pident_og'})
bsr_df = db_score.merge(og_score, on="qseqid")
bsr_df["bsr"] = bsr_df["score_db"] / bsr_df["score_og"]
bsr_df.to_csv(bsr_file, sep="\t", index=False)
