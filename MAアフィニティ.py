import requests
import tkinter as tk
from rdkit import Chem
from selenium import webdriver
import pandas as pd


def get_iupac_from_psychonaut(substance):
    url = f'https://www.psychonautwiki.org/api.php?action=browsebysubject&subject={substance}&format=json'
    response = requests.get(url)
    data = response.json()
    iupac = data['query']['pages'][0]['fulltext'][0]['content'][0]['iupac']
    return iupac


def iupac_to_smiles(iupac):
    mol = Chem.MolFromIUPACName(iupac)
    smiles = Chem.MolToSmiles(mol)
    return smiles


def get_adme_data(smiles):
    driver = webdriver.Chrome()
    driver.get('http://www.swissadme.ch/')
    driver.find_element_by_name('SMILES').send_keys(smiles)
    driver.find_element_by_name('qsar').click()
    driver.find_element_by_link_text('Export QSAR results as CSV').click()
    driver.quit()
    adme_data = pd.read_csv('QSAR_results.csv', delimiter=';')
    return adme_data


def get_ki_data(adme_data):
    ki_data = {
        'DAT': adme_data.loc[adme_data['Target'] == 'DAT', 'Ki (nM)'].iloc[0],
        'NAT': adme_data.loc[adme_data['Target'] == 'NAT', 'Ki (nM)'].iloc[0],
        'SERT': adme_data.loc[adme_data['Target'] == 'SERT', 'Ki (nM)'].iloc[0]
    }
    return ki_data


def button_click():
    iupac = iupac_entry.get()
    smiles = iupac_to_smiles(iupac)
    adme_data = get_adme_data(smiles)
    ki_data = get_ki_data(adme_data)
    result_text.configure(state='normal')
    result_text.delete(1.0, tk.END)
    result_text.insert(tk.END, f'DAT Ki: {ki_data["DAT"]} nM\n'
                      f'NAT Ki: {ki_data["NAT"]} nM\n'
                      f'SERT Ki: {ki_data["SERT"]} nM')
    result_text.configure(state='disabled')


substance = 'Methamphetamine'
iupac = get_iupac_from_psychonaut(substance)
print(iupac)

root = tk.Tk()
root.title('Monoamine Receptor Affinity')

iupac_entry = tk.Entry(root, width=40)
iupac_entry.pack(pady=10)

button = tk.Button(root, text='Get Affinities', command=button_click)
button.pack(pady=10)

result_text = tk.Text(root, width=50, height=10, state='disabled')
result_text.pack(pady=10)

root.mainloop()

