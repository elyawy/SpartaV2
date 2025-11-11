from itaxotools.mafftpy import ginsi
from pathlib import Path

# Create test FASTA
with open("test.fasta", "w") as f:
    f.write(">seq1\nACGTACGT\n>seq2\nACGTACGT\n")

# Provide an output file
result = ginsi(input=Path("test.fasta"), output=Path("output.fasta"))

print(f"Type of result: {type(result)}")
print(f"Result: {result}")

# Read the output file
if Path("output.fasta").exists():
    print("\nOutput file contents:")
    print(Path("output.fasta").read_text())