#!/bin/python3
#-*- coding:utf-8 -*-

################################# IMPORTATION DES MODULES #################################

import os
import pandas as pd
from typing import List, Dict, Tuple
from tqdm import tqdm
from reads import Reads
from config import GlobalValue, transformer

################################# IMPORTATION DES MODULES #################################


################################### VARIABLES GLOBALES ####################################

GVALUES: GlobalValue = GlobalValue()
CHERIF: Dict[str, List[str]] = GVALUES.get_values("primers")["CHERIF"]
CHERIF_RANGES: Dict[str, Tuple[int, int]] = GVALUES.get_values("ranges")["CHERIF"]
FORWARDS: str = CHERIF["forwards"]
REVERSES: str = CHERIF["reverses"]
MIN_LENGTH: int = CHERIF_RANGES[0]
MAX_LENGTH: int = CHERIF_RANGES[0]
PATH: str = "../01_clf_ra/fastqs"
FASTQS: List[str] = os.listdir(path=PATH)
N_BARCODES: int = len(FASTQS)
BARCODES: List[int] = [i + 1 for i in range(N_BARCODES)]

################################### VARIABLES GLOBALES ####################################


########################################## MAIN ###########################################

def main() -> None:
	"""
	Fonction principale pour la vectorisation des barcodes contenus dans 'PATH'.

	Args:
		Aucun.

	Returns:
		None.
	"""
	# Boucle sur les barcodes
	for barcode, name in zip(BARCODES, FASTQS):
		# Afficher le barcode actuel
		print(f"\n\nProcessing barcode : {barcode}")

		# Charger le barcode dans un objet reads
		reads_object = Reads(filename=f"{PATH}/{name}")

		# Traitement des primers
		reads_object.cutPrimers(forwards=FORWARDS, reverses=REVERSES, seuil_length=MIN_LENGTH)
		
		# Filtrage des tailles
		reads_object.filter_size(inf=MIN_LENGTH, sup=MAX_LENGTH)

		# Récupérer les séquences
		sequences = reads_object.getColumn("sequence")

		# Supprimer 'reads_object'
		del reads_object

		# Initialisation de la liste qui contiendra les sequences vectorisées
		df = [None for _ in range(len(sequences))]

		# Boucle sur les séquences
		for i, sequence in tqdm(iterable=enumerate(sequences), desc="Vectorisation", unit=" Reads"):
			# Vectorisation de la séquence
			sequence_vectorized = transformer(sequence=sequence, n_features=MAX_LENGTH)
			
			# Ajout de la séquence vectorisée à la liste
			df[i] = sequence_vectorized

		# Création du DataFrame final labelisé
		df = pd.DataFrame(data=df)

		# Sauvegarde du DataFrame dans un fichier CSV
		filename = f"./datasets/Sd_testSet_BC0{barcode}.csv" if barcode <= 9 else f"./datasets/Sd_testSet_BC{barcode}.csv"
		df.to_csv(path_or_buf=filename, index=False)
		
		# Supprimer 'df'
		del df

########################################## MAIN ###########################################

if __name__ == "__main__":
	main()