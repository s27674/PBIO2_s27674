from Bio import Entrez, SeqIO
import pandas as pd, matplotlib.pyplot as plt, time

def fetch(email, key, taxid, a, b, n):
    Entrez.email, Entrez.api_key = email, key
    r = Entrez.read(Entrez.efetch(db="taxonomy", id=taxid, retmode="xml"))[0]
    print(f"{r['ScientificName']} (TaxID: {taxid})")
    s = Entrez.read(Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", usehistory="y", retmax=0))
    w, q, c, out, i = s["WebEnv"], s["QueryKey"], int(s["Count"]), [], 0
    while i < c and len(out) < n:
        h = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", retstart=i, retmax=100, webenv=w, query_key=q)
        for r in SeqIO.parse(h, "genbank"):
            l = len(r.seq)
            if a <= l <= b: out.append({"accession": r.id, "length": l, "description": r.description})
            if len(out) >= n: break
        i += 100; time.sleep(0.4)
    return out

def report(d, t):
    df = pd.DataFrame(d).sort_values("length", ascending=False)
    df.to_csv(f"{t}_report.csv", index=False)
    plt.plot(df["accession"], df["length"], marker='o')
    plt.xticks(rotation=90, fontsize=7); plt.tight_layout()
    plt.savefig(f"{t}_plot.png")

if __name__ == "__main__":
    e = input("Email: "); k = input("API key: "); t = input("TaxID: ")
    a, b = map(int, input("Min/Max length: ").split()); n = int(input("Max records: "))
    d = fetch(e, k, t, a, b, n)
    report(d, t) if d else print("No matches.")
