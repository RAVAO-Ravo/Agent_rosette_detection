#!/bin/python3
#-*- coding:utf-8 -*-

################################# IMPORTATION DES MODULES #################################

import os
import pandas as pd
import subprocess as sp
from typing import List
from tqdm import tqdm
from reads import Reads
from utilitaire import CHERIF, CHERIF_RANGES, transformer

################################# IMPORTATION DES MODULES #################################


################################### VARIABLES GLOBALES ####################################

FORWARDS: str = CHERIF["forwards"]
REVERSES: str = CHERIF["reverses"]
MIN_LENGTH: int = CHERIF_RANGES[0]
MAX_LENGTH: int = CHERIF_RANGES[1]
PATH: str = "fastqs"
FASTQS: List[str] = os.listdir(path=PATH)
N_BARCODES: int = len(FASTQS)
BARCODES: List[int] = [i + 1 for i in range(N_BARCODES)]
TMP: str = "tmp.fastq"
CMD: List[str] = ["minimap2 -a -x map-ont -k 28 species/Sphaerothecum_destruens.fasta tmp.fastq | samtools view -S -F 4 -q 60"]

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
		if barcode in [6, 12]:
			# Afficher le barcode actuel
			print(f"\n\nProcessing barcode : {barcode}")

			# Charger le barcode dans un objet reads
			reads_object = Reads(filename=f"{PATH}/{name}")

			# Traitement des primers
			reads_object.cutPrimers(forwards=FORWARDS, reverses=REVERSES, seuil_length=MIN_LENGTH)
			
			# Filtrage des tailles
			reads_object.filter_size(inf=MIN_LENGTH, sup=MAX_LENGTH)

			# Enregistrement dans un fichier temporaire
			reads_object.write(TMP)

			# Supprimer 'reads_object'
			del reads_object

			# Mapping sur la séquence de référence de l'Agent Rosette
			sequences = sp.run(args=CMD, shell=True, capture_output=True, text=True).stdout.strip()

			# Si l'output est une chaîne non-vide
			if sequences != '':
				# Récupérer les séquences appartenant à l'Agent Rosette
				sequences = [line.split('\t')[9] for line in sequences.split('\n')]

			# Initialisation de la liste qui contiendra les sequences vectorisées
			df = [None for _ in range(len(sequences))]

			# Boucle sur les séquences
			for i, sequence in tqdm(iterable=enumerate(sequences), desc="Vectorisation", unit=" Reads"):
				# Vectorisation de la séquence
				sequence_vectorized = transformer(sequence=sequence, n_features=MAX_LENGTH)
				
				# Ajout de la séquence vectorisée aux tableaux
				df[i] = sequence_vectorized

			# Création du DataFrame final labelisé
			df = pd.DataFrame(data=df)
			df[MAX_LENGTH] = "Sphaerothecum_destruens"

			# Sauvegarde du DataFrame dans un fichier CSV
			filename = f"./datasets/Sd_testSet_BC0{barcode}.csv" if barcode <= 9 else f"./datasets/Sd_testSet_BC{barcode}.csv"
			df.to_csv(path_or_buf=filename, index=False)
			
			# Supprimer 'df'
			del df

	# Supprimer le fichier temporaire
	os.remove(path=TMP)

########################################## MAIN ###########################################

if __name__ == "__main__":
	main()