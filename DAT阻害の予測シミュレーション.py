import pandas as pd
import tkinter as tk
from tkinter import ttk
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import threading
from deap import base, creator, tools, algorithms

# 適応度クラスを定義
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

# 遺伝アルゴリズムの初期化
def init_ga():
    toolbox = base.Toolbox()
    toolbox.register("attr_bool", random.randint, 0, 1)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=100)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    return toolbox

# 評価関数を定義
def eval_one_max(individual):
    return sum(individual),

# 遺伝アルゴリズムのメインループ
def ga_main_loop(toolbox):
    toolbox.register("evaluate", eval_one_max)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=3)

    population = toolbox.population(n=100)
    ngen = 100
    for gen in range(ngen):
        offspring = algorithms.varAnd(population, toolbox, cxpb=0.5, mutpb=0.1)
        fits = toolbox.map(toolbox.evaluate, offspring)
        for fit, ind in zip(fits, offspring):
            ind.fitness.values = fit
        population = toolbox.select(offspring, k=len(population))

    top_ind = tools.selBest(population, k=1)[0]
    return top_ind

# GUIアプリケーション
class AffinityPredictorApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("DAT Affinity Predictor")
        self.geometry("400x200")
        self.create_widgets()

    def create_widgets(self):
        self.lbl_iupac = ttk.Label(self, text="Enter IUPAC name:")
        self.lbl_iupac.grid(column=0, row=0)

        self.entry_iupac = ttk.Entry(self)
        self.entry_iupac.grid(column=1, row=0)

        self.btn_predict = ttk.Button(self, text="Predict Affinity", command=self.predict_affinity)
        self.btn_predict.grid(column=1, row=1)

        self.lbl_result = ttk.Label(self, text="Predicted Affinity:")
        self.lbl_result.grid(column=0, row=2)

        self.txt_result = tk.Text(self, width=30, height=300)
