# DNA to Protein Translator + Alignment core logic
# (Same logic as appp2.py, with the CLI loop guarded so this module
# can be safely imported by app.py without launching the terminal menu.)

genetic_code = {
    "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
    "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
    "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
    "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
    "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
    "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
    "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
    "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
    "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
    "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
    "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
    "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
    "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
    "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
    "TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_",
    "TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W",
}

def dna_to_protein(dna_seq):
    protein = ""
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        amino_acid = genetic_code.get(codon.upper(), "X")
        protein += amino_acid
    return protein


def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    score = [[0]*(m+1) for _ in range(n+1)]

    for i in range(1, n+1):
        score[i][0] = i * gap
    for j in range(1, m+1):
        score[0][j] = j * gap

    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up = score[i-1][j] + gap
            left = score[i][j-1] + gap
            score[i][j] = max(diag, up, left)

    align1, align2 = "", ""
    i, j = n, m
    while i > 0 and j > 0:
        current = score[i][j]
        if current == score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1; j -= 1
        elif current == score[i-1][j] + gap:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1

    while i > 0:
        align1 = seq1[i-1] + align1
        align2 = "-" + align2
        i -= 1
    while j > 0:
        align1 = "-" + align1
        align2 = seq2[j-1] + align2
        j -= 1

    return align1, align2, score[n][m]


def smith_waterman(seq1, seq2, match=3, mismatch=-3, gap=-2):
    n, m = len(seq1), len(seq2)
    score = [[0]*(m+1) for _ in range(n+1)]
    max_score, max_pos = 0, (0, 0)

    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up = score[i-1][j] + gap
            left = score[i][j-1] + gap
            score[i][j] = max(0, diag, up, left)
            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    align1, align2 = "", ""
    i, j = max_pos
    while i > 0 and j > 0 and score[i][j] > 0:
        current = score[i][j]
        if current == score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1; j -= 1
        elif current == score[i-1][j] + gap:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2, max_score


def run_cli():
    """The original terminal menu, preserved but no longer run on import."""
    while True:
        print("DNA information tool")
        print("choose your operation a or b")
        print("a. DNA to Protein")
        print("b. DNA/Protein Alignment")
        option = input("Enter your choice: ")

        if option == "a":
            dna_sequence = input("Enter a DNA sequence: ")
            protein_sequence = dna_to_protein(dna_sequence)
            print("DNA:", dna_sequence)
            print("Protein:", protein_sequence)

        elif option == "b":
            seqA = input("Enter sequence A: ")
            seqB = input("Enter sequence B: ")
            print("which alignment would you like to use? pick 1 or 2")
            print("1. Needleman-Wunsch (Global)")
            print("2. Smith-Waterman (Local)")
            choice = int(input("Enter your choice: "))
            if choice == 1:
                print("Needleman-Wunsch (Global):")
                nw_align = needleman_wunsch(seqA, seqB)
                print(nw_align[0], "\n", nw_align[1], "\nScore:", nw_align[2])
            elif choice == 2:
                print("Smith-Waterman (Local):")
                sw_align = smith_waterman(seqA, seqB)
                print(sw_align[0], "\n", sw_align[1], "\nScore:", sw_align[2])
            else:
                print("Invalid option")
        else:
            print("Invalid option")

        print("do you want to continue? y/n")
        pick = input("enter your choice: ")
        if pick == "n":
            break
        elif pick != "y":
            print("Invalid option")


if __name__ == "__main__":
    run_cli()
