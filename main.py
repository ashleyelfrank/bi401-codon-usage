from pathlib import Path
from collections import Counter
import csv
data_file=(
    Path(__file__).resolve().parent.parent/"data"
    /"ncbi_dataset"
    /"ncbi_dataset"
    /"data"
    /"GCA_000001215.4"
    /"cds_from_genomic.fna"
)
print(data_file.exists())
with data_file.open("r") as fasta_file: 
    first_line = fasta_file.readline()
    print(first_line)
with data_file.open("r") as fasta_file:
    first_line = fasta_file.readline()
    second_line = fasta_file.readline()
    print(second_line)
part1="ATGTCG"
part2="TGCAAG"
combined=part1 + part2
print(combined)
line=">gene_A"
print(line.startswith(">"))
print(first_line.startswith(">"))
print(second_line.startswith(">"))
with data_file.open("r") as fasta_file:
    header=fasta_file.readline().strip()
    sequence= ""
    for line in fasta_file:
        line=line.strip()
        if line.startswith(">"):
            break
        sequence += line
print(header)
print(sequence)
print(len(sequence))
codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
print(codons[0:10])
from collections import Counter
codon_counts = Counter(codons)
print(codon_counts.most_common(10))
print("Sequence length:", len(sequence))
print("Number of codons:", len(codons))
print("Sum of codon counts:", sum(codon_counts.values()))
start=first_line.find("[gene=")
print(start)
gene_start= start + len("[gene=")
gene_end= first_line.find("]", gene_start)
gene_name= first_line[gene_start:gene_end]
print(gene_name)
gene_results ={}
gene_results[gene_name] = codon_counts
print(gene_results[gene_name].most_common(10))
import csv
record_count=0
with data_file.open("r") as fasta_file:
    for line in fasta_file:
        if line.startswith(">"):
            record_count += 1
print("Total FASTA records:", record_count)
codon_usage={
    "ABW0931.1": {
        "gene":"CG17636",
        "codons":{"ATG": 16, "GAA": 24}
    }
}
from pathlib import Path
from collections import Counter
def process_record(header, sequence):
    start = header.find("[gene=")
    if start == -1:
        gene_name= None
    else:
        gene_start= start + len("[gene=")
        gene_end= header.find("]", gene_start)
        gene_name= header[gene_start:gene_end]
    codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
    return {
    "gene": gene_name,
    "codons": Counter(codons)
    }
codon_usage = {}
record_count = 0
header = None
sequence_parts = []
with data_file.open("r") as fasta_file:
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                record_id = header[1:].split()[0] 
                codon_usage[record_id] = process_record(header, "".join(sequence_parts))
                record_count += 1
            header = line
            sequence_parts = []
        else:
            sequence_parts.append(line)
if header is not None:
        record_id = header[1:].split()[0]
        codon_usage[record_id] = process_record(header, "".join(sequence_parts))
        record_count += 1
print("Reords processed:", record_count)
print("Records stored:", len(codon_usage))
print ("First 3 record IDs:", list(codon_usage.keys())[:3])
print ("Last record ID:", list(codon_usage.keys())[-1])
print ("Unique record IDs:", len(set(codon_usage.keys())))
missing_genes= sum(1 for record in codon_usage.values() if record["gene"] is None)
print("Records missing gene annoations:", missing_genes)
bases = "ATCG"
all_codons=[
    first + second + third
    for first in bases
    for second in bases
    for third in bases
]
print("Number of standard codons:", len(all_codons))
print("Data file exists:", data_file.exists())
output_file=Path(__file__).resolve().parent.parent/"data/codon_usage.csv" 
with output_file.open("w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["record_id", "gene_name"] + all_codons)
    for record_id, record_data in codon_usage.items():
        row = [record_id, record_data["gene"]]
        for codon in all_codons:
            row.append(record_data["codons"].get(codon, 0))
        writer.writerow(row)
print("CSV export complete:", output_file)