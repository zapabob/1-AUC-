import tkinter as tk
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem import MACCSkeys

def calculate_similarity(smiles1, smiles2, method="morgan"):
    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)

        if mol1 is None or mol2 is None:
            raise ValueError("Invalid SMILES string")

        if method == "morgan":
            fp1 = AllChem.GetMorganFingerprint(mol1, 2)
            fp2 = AllChem.GetMorganFingerprint(mol2, 2)
        elif method == "maccs":
            fp1 = MACCSkeys.GenMACCSKeys(mol1)
            fp2 = MACCSkeys.GenMACCSKeys(mol2)
        else:
            raise ValueError("Invalid method: choose 'morgan' or 'maccs'")

        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        return similarity
    except Exception as e:
        print(f"Error: {e}")
        return None

def calculate_and_display_similarity():
    smiles1 = entry_amphetamine.get()
    smiles2 = entry_candidate.get()
    similarity = calculate_similarity(smiles1, smiles2, method="morgan")
    if similarity is not None:
        label_result.config(text=f"類似性: {similarity:.4f}")
    else:
        label_result.config(text="類似性: エラー")

def update_dataframe():
    amphetamine_smiles = entry_amphetamine.get()
    candidates_smiles = text_candidates.get("1.0", tk.END).splitlines()
    similarities = [calculate_similarity(amphetamine_smiles, candidate_smiles) for candidate_smiles in candidates_smiles]

    data = {"Candidate SMILES": candidates_smiles, "Similarity": similarities}
    df = pd.DataFrame(data)
    print(df)

root = tk.Tk()
root.title("類似性計算")

label_amphetamine = tk.Label(root, text="アンフェタミンのSMILES")
entry_amphetamine = tk.Entry(root)
label_candidate = tk.Label(root, text="候補化合物のSMILES")
entry_candidate = tk.Entry(root)
button_calculate = tk.Button(root, text="類似性を計算する", command=calculate_and_display_similarity)
label_result = tk.Label(root, text="類似性: ")

label_candidates = tk.Label(root, text="複数の候補化合物のSMILES（1行に1つ）")
text_candidates = tk.Text(root, height=10, width=40)
button_update_dataframe = tk.Button(root, text="データフレームを更新する", command=update_dataframe)

label_amphetamine.grid(row=0, column=0, sticky="e")
entry_amphetamine.grid(row=0, column=1)
label_candidate.grid(row=1, column=0, sticky="e")
entry_candidate.grid(row=1, column=1)

