#!/usr/bin/env python3
import sys
import gffutils
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

if len(sys.argv) != 4:
    sys.exit("Usage: gff3_cds_to_gbk.py genome.fna cds.gff3 output.gbk")

genome_fa, gff3, out_gbk = sys.argv[1:]

# Load genome FASTA
genome = SeqIO.to_dict(SeqIO.parse(genome_fa, "fasta"))

# Build in-memory GFF database
db = gffutils.create_db(
    gff3,
    dbfn=":memory:",
    force=True,
    keep_order=True,
    merge_strategy="merge",
    sort_attribute_values=True,
)

# Reset features and add required annotation
for rec in genome.values():
    rec.features = []
    rec.annotations["molecule_type"] = "DNA"

seen = set()  # (seqid, start0, end, strand)
kept = skipped_dups = skipped_missing_seq = 0

for cds in db.features_of_type("CDS", order_by=("seqid", "start", "end")):
    if cds.seqid not in genome:
        skipped_missing_seq += 1
        continue

    start0 = cds.start - 1   # GFF is 1-based inclusive
    end = cds.end            # Biopython end is exclusive
    strand = 1 if cds.strand == "+" else -1 if cds.strand == "-" else 0

    key = (cds.seqid, start0, end, strand)
    if key in seen:
        skipped_dups += 1
        continue
    seen.add(key)

    qualifiers = {}
    if "ID" in cds.attributes and cds.attributes["ID"]:
        qualifiers["locus_tag"] = [cds.attributes["ID"][0]]
    if "Parent" in cds.attributes and cds.attributes["Parent"]:
        qualifiers["parent"] = [cds.attributes["Parent"][0]]

    genome[cds.seqid].features.append(
        SeqFeature(
            FeatureLocation(start0, end, strand=strand),
            type="CDS",
            qualifiers=qualifiers,
        )
    )
    kept += 1

SeqIO.write(genome.values(), out_gbk, "genbank")
print(
    f"Wrote CDS={kept}; skipped_dups={skipped_dups}; "
    f"skipped_missing_seq={skipped_missing_seq} -> {out_gbk}"
)
