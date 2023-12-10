#!/bin/python3
#-*- coding:utf-8 -*-

################################# IMPORTATION DES MODULES #################################

import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict, List, TextIO, Tuple
from sklearn.manifold import TSNE
from sklearn.metrics import classification_report, confusion_matrix, ConfusionMatrixDisplay
from config import GlobalValue, transformer

################################# IMPORTATION DES MODULES #################################


################################### VARIABLES GLOBALES ####################################

GVALUES: GlobalValue = GlobalValue()
CHERIF: Dict[str, List[str]] = GVALUES.get_values("primers")["CHERIF"]
CHERIF_RANGES: Dict[str, Tuple[int, int]] = GVALUES.get_values("ranges")["CHERIF"]
ENCODER: Dict[str, int] = GVALUES.get_values("encoder")
CLF_FILE: str = GVALUES.get_values("clf_file")
SPECIES_AMPLICONS_FILE: str = "./datasets/fakes_species_amplicons.csv"
RANDOMS_AMPLICONS_FILE: str = "./datasets/fakes_unidentified_amplicons.csv"
TEST_SET_FILE: str = "./datasets/Sd_testSet_BC06.csv"
VALIDATION_SET_FILE: str = "./datasets/Sd_testSet_BC12.csv"
RANDOM_STATE: int = 42
DATA_SIZE: int =  2000

################################### VARIABLES GLOBALES ####################################


####################################### FONCTIONS #########################################

def exec_cmd(cmd: str) -> None:
	"""
	Exécute une commande système.

	Args:
		cmd (str): La commande à exécuter.

	Returns:
		None
	"""
	sp.run(args=[cmd], shell=False)

def readFasta(filename: str) -> List[str]:
	"""
	Récupère les séquences à partir d'un fichier au format FASTA.

	Args:
		filename (str): Le chemin du fichier FASTA à lire.

	Returns:
		List[str]: Une liste contenant les séquences du fichier FASTA.
	"""
	# Initialiser les variables 
	sequences: List[str] = []
	fasta: TextIO[str] = ''
	current_line: str = ''
	current_seq: str = ''

	# Ouvrir le fichier
	with open(file=filename, mode='r') as fasta:
		# Récupérer la ligne courrante
		current_line = fasta.readline().rstrip()

		# Initialiser la séquence courrante
		current_seq = ''

		# Tant que l'on n'a pas parcouru tout le fichier
		while len(current_line) != 0:
			# Si la première lettre de la ligne n'est pas un nucléotide
			if current_line[0] not in ('>', '\n', '\t', '\r', '', ' '):
				# Ajouter la ligne à la séquence courrante
				current_seq += current_line
			
			# Sinon
			else:
				# Si on n'est pas à la fin du fichier
				if len(current_seq) != 0:
					# Ajouter la séquence courrante à la liste des séquences
					sequences.append(current_seq)

					# Réinitialiser la séquence courrante
					current_seq = ''
			
			# Récupérer la nouvelle ligne courrante
			current_line = fasta.readline().rstrip()

			# S'il n'y a plus de lignes dans le fichier, mais que la séquence courrante contient une séquence
			if len(current_line) == 0 and len(current_seq) != 0:
				# Ajouter la séquence à la liste des séquences
				sequences.append(current_seq)

	# Retourner les séquences sous forme d'une liste triée
	return sorted(sequences)

def init_dataFrame(sequences: List[str], n_features: int) -> pd.DataFrame:
	"""
	Initialise un DataFrame à partir d'une liste de séquences.

	Args:
		sequences (List[str]): La liste des séquences à inclure dans le DataFrame.
		n_features (int): Le nombre de caractéristiques (colonnes) à créer pour chaque séquence.

	Returns:
		pd.DataFrame: Un DataFrame pandas contenant les séquences en tant que lignes et les caractéristiques en tant que colonnes.
	"""
	# Retourner le dataframe
	return pd.DataFrame(data=[transformer(sequence=sequence, n_features=n_features) for sequence in sequences])

def encode_sp(x: str) -> int:
	"""
	Transforme les labels selon l'encodage défini par le dictionnaire 'ENCODER'.
	
	Args:
		x (str): La chaîne de caractères à encoder.
		
	Returns:
		int: La valeur encodée correspondante. Si la chaîne n'est pas trouvée dans
			 le dictionnaire, la fonction retourne 0 par défaut.
	"""
	# Retourner l'encodage de la valeur
	return ENCODER[x]

def plot_data(
		X: pd.DataFrame,
		y: pd.Series,
		n_sample: int = 5000,
		figsize: Tuple[int, int] = (20, 20),
		palette: List[str] = ["turquoise", "crimson"],
		save: str = None) -> None:
	"""
	Plot les données en utilisant t-SNE.

	Args:
		X (pd.DataFrame): DataFrame des caractéristiques.
		y (pd.Series): Series des étiquettes.
		n_sample (int): Nombre d'échantillons à utiliser pour t-SNE. (Par défaut = 5000).
		figsize (Tuple[int, int]): Taille de la figure. (Par défaut = (20, 20)).
		palette (List[str]): Palette de couleurs pour les classes. (Par défaut = ["turquoise", "crimson"]).
		save (str): Si fourni, enregistre la figure dans le fichier spécifié. (Par défaut = None).

	Returns:
		None: La fonction affiche le plot mais ne renvoie rien.
	"""
	# Effectuer le TSNE
	tsne = TSNE(n_components=2, init="pca", learning_rate="auto", random_state=RANDOM_STATE).fit_transform(X=X.iloc[:n_sample], y=y.iloc[:n_sample])
	tsne = pd.concat(objs=[pd.DataFrame(data=tsne), y.iloc[:n_sample].reset_index(drop=True)], axis=1)
	tsne.columns = ['x', 'y', "specie"]

	# Créer la figure le TSNE
	plt.figure(figsize=figsize)

	# Pour chaque espèce
	for i, specie in enumerate(np.unique(tsne["specie"]).tolist(), start=1):
		# Plotter l'espèce
		plt.subplot(4, 4, i)
		tmp = tsne["specie"].apply(lambda x: specie if x == specie else -1)
		sns.scatterplot(data=tsne, x='x', y='y', hue=tmp, palette=palette)

	# Sauvegarder la figure
	if save != None:
		plt.savefig(save)

	# Afficher la figure
	plt.show()

def plot_cm(
		y_true: pd.DataFrame,
		y_pred: pd.Series,
		figsize: Tuple[int, int] = (14, 8),
		cmap: str = "magma_r",
		save: str = None) -> None:
	"""
	Plot une matrice de confusion.

	Args:
		y_true (pd.DataFrame): DataFrame des vraies étiquettes.
		y_pred (pd.Series): Series des étiquettes prédites.
		figsize (Tuple[int, int]): Taille de la figure. (Par défaut = (14, 8)).
		cmap (str): Colormap à utiliser pour le plot. (Par défaut = "magma_r").
		save (str): Si fourni, enregistre la figure dans le fichier spécifié. (Par défaut = None).

	Returns:
		None: La fonction affiche le plot mais ne renvoie rien.
	"""
	# Calculer la matrice
	cm = confusion_matrix(y_true=y_true, y_pred=y_pred)
	cm = ConfusionMatrixDisplay(confusion_matrix=cm)
	
	# Créer la figure
	_, ax = plt.subplots(figsize=figsize)

	# Plotter la matrice
	cm.plot(ax=ax, cmap=cmap)

	# Sauvegarder la figure
	if save != None:
		plt.savefig(save)

	# Afficher la figure
	plt.show()

def plot_classReport(classifier: object,
					 X: pd.DataFrame,
					 y: pd.Series,
					 figsize: Tuple[int, int] = (14, 8),
					 cmap: str = "magma_r",
					 save: str = None) -> None:
	"""
	Affiche les résultats de classification d'un modèle.

	Args:
		classifier (object): Le modèle de classification entraîné.
		X (pd.DataFrame): DataFrame des caractéristiques.
		y (pd.Series): Series des étiquettes.
		figsize (Tuple[int, int]): Taille de la figure. (Par défaut = (14, 8)).
		cmap (str): Colormap à utiliser pour le plot. (Par défaut = "magma_r").
		save (str): Si fourni, enregistre la figure dans le fichier spécifié. (Par défaut = None).

	Returns:
		None: La fonction affiche le plot mais ne renvoie rien.
	"""
	# Effectuer les prédictions sur les données
	y_pred = classifier.predict(X=X)

	# Afficher le rapport de classification
	print(classification_report(y_true=y, y_pred=y_pred, zero_division=0))

	# Afficher la matrice de confusion
	plot_cm(y_true=y, y_pred=y_pred, figsize=figsize, cmap=cmap, save=save)

####################################### FONCTIONS #########################################