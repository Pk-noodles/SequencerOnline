from flask import Flask, render_template, request
from dna_core import dna_to_protein, needleman_wunsch, smith_waterman

app = Flask(__name__, template_folder=".")

# Amino acid property groups -> CSS class used for color-coding in the UI
AA_CLASS = {}
for aa in "AVLIMFWPG":
    AA_CLASS[aa] = "aa-nonpolar"
for aa in "STCYNQ":
    AA_CLASS[aa] = "aa-polar"
for aa in "KRH":
    AA_CLASS[aa] = "aa-positive"
for aa in "DE":
    AA_CLASS[aa] = "aa-negative"
AA_CLASS["_"] = "aa-stop"
AA_CLASS["X"] = "aa-unknown"


def build_codon_trace(dna_sequence, protein_sequence):
    """Pairs each codon with its amino acid + a color-coding class."""
    trace = []
    for idx, aa in enumerate(protein_sequence):
        codon = dna_sequence[idx*3:idx*3+3].upper()
        trace.append({
            "codon": codon,
            "aa": aa,
            "cls": AA_CLASS.get(aa, "aa-unknown"),
        })
    return trace


def build_alignment_pairs(align1, align2):
    """Per-character classification for the alignment viewer: match / mismatch / gap."""
    pairs = []
    for a, b in zip(align1, align2):
        if a == "-" or b == "-":
            cls = "align-gap"
        elif a == b:
            cls = "align-match"
        else:
            cls = "align-mismatch"
        pairs.append({"a": a, "b": b, "cls": cls})
    return pairs


@app.route("/", methods=["GET"])
def index():
    return render_template("index.html", active_tab="translate")


@app.route("/translate", methods=["POST"])
def translate():
    dna_sequence = request.form.get("dna_sequence", "").strip()
    error = None
    protein_sequence = None
    codon_trace = None

    if not dna_sequence:
        error = "Enter a DNA sequence first."
    else:
        protein_sequence = dna_to_protein(dna_sequence)
        codon_trace = build_codon_trace(dna_sequence, protein_sequence)

    return render_template(
        "index.html",
        active_tab="translate",
        dna_sequence=dna_sequence,
        protein_sequence=protein_sequence,
        codon_trace=codon_trace,
        translate_error=error,
    )


@app.route("/align", methods=["POST"])
def align():
    seq_a = request.form.get("seq_a", "").strip()
    seq_b = request.form.get("seq_b", "").strip()
    algorithm = request.form.get("algorithm", "nw")
    error = None
    align1 = align2 = score = None
    alignment_pairs = None

    if not seq_a or not seq_b:
        error = "Enter both sequence A and sequence B."
    else:
        if algorithm == "sw":
            align1, align2, score = smith_waterman(seq_a, seq_b)
        else:
            align1, align2, score = needleman_wunsch(seq_a, seq_b)
        alignment_pairs = build_alignment_pairs(align1, align2)

    return render_template(
        "index.html",
        active_tab="align",
        seq_a=seq_a,
        seq_b=seq_b,
        algorithm=algorithm,
        align1=align1,
        align2=align2,
        score=score,
        alignment_pairs=alignment_pairs,
        align_error=error,
    )


if __name__ == "__main__":
    app.run(debug=True)
