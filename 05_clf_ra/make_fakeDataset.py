#!/bin/python3
#-*- coding:utf-8 -*-

################################# IMPORTATION DES MODULES #################################

import os
import pandas as pd
import random as rd
from typing import List, Tuple
from argparse import ArgumentParser, HelpFormatter
from utilitaire import CHERIF, CHERIF_RANGES, RANDOM_STATE, init_dataFrame, readFasta, exec_cmd 

################################# IMPORTATION DES MODULES #################################


################################# ARGUMENTS DU PROGRAMME ##################################

parser = ArgumentParser(
				description="Ce programme permet de générer un dataset de données (de séquençages) simulées, pour les espèces du projet eDNA.",
				formatter_class=lambda prog: HelpFormatter(prog=prog, max_help_position=100, width=200)
			)
parser.add_argument("-n", "--n_sequences", type=int, default=10000, help="number of sequences to generate for each specie")
parser.add_argument("--min_error", type=float, default=0.0, help="minimum error rate when generating a sequence")
parser.add_argument("--max_error", type=float, default=0.15, help="maximum error rate when generating a sequence")
args = parser.parse_args()

################################# ARGUMENTS DU PROGRAMME ##################################


################################### VARIABLES GLOBALES ####################################

FORWARD: str = CHERIF["forwards"][0]
REVERSE: str = CHERIF["reverses"][0]
MIN_LENGTH: int = CHERIF_RANGES[0]
MAX_LENGTH: int = CHERIF_RANGES[1]
N_SEQUENCES: int = args.n_sequences
MIN_ERROR: float = args.min_error
MAX_ERROR: float = args.max_error
PATH: str = "species"
SPECIES: List[str] = sorted(os.listdir(path=PATH))
RES_DIRECTORY: str = "datasets"
AMPLICON_POS: List[Tuple[int, int]] = [(595, 709), (583, 702),
									   (542, 654), (566, 680),
									   (585, 699), (554, 675),
									   (564, 686), (490, 604),
									   (566, 680), (577, 693)]

################################### VARIABLES GLOBALES ####################################


################################## UNIDENTIFIED SPECIES ###################################

def random_sequence() -> str:
	"""
	Génère un amplicon aléatoire d'ADN.

	Args:
		None

	Returns:
		str: Séquence d'ADN générée.
	"""
	# Générer la séquence entre les amorces (primers)
	size_seq = rd.randrange(start=MIN_LENGTH, stop=MAX_LENGTH, step=1) - len(FORWARD) - len(REVERSE)

	# Retourner l'amplicon
	return FORWARD + ''.join([rd.choice("ATCG") for _ in range(size_seq)]) + REVERSE

def generateUnidentifiedAmplicons(outputFile: str) -> None:
	""" 
	Génère le dataset des séquences 'non-identifiées'.

	Args:
		outputFile (str): Le chemin du fichier de sortie où le dataset sera enregistré.

	Returns:
		None
	"""
	# Créer le répertoire résultat s'il n'existe pas
	if RES_DIRECTORY not in os.listdir(path='.'):
		exec_cmd(cmd=f"mkdir {RES_DIRECTORY}")

	# Initialiser un set pour s'assurer que chaque séquence est unique
	sequences = set()

	# Tant que l'on n'a pas le nombre de séquences souhaités
	while len(sequences) != N_SEQUENCES:
		# Créer une séquence aléatoire et ajouter la au set de séquences
		sequences.add(random_sequence())
	
	# Transformer le set en dataframe
	dataset = init_dataFrame(sequences=sorted(sequences), n_features=MAX_LENGTH).assign(species="unidentified")
	dataset.columns = [i for i in range(MAX_LENGTH + 1)]
	dataset.to_csv(path_or_buf=outputFile, index=False)

################################## UNIDENTIFIED SPECIES ###################################


################################### SIMULATED SPECIES #####################################

def addError(sequence: str) -> str:
	""" 
	Génère une séquence nucléotidique avec des erreurs, à partir d'une séquence modèle.

	Args:
		sequence (str): La séquence nucléotidique modèle à partir de laquelle générer la séquence avec erreurs.

	Returns:
		str: La séquence nucléotidique avec des erreurs.
	"""
	# Définir le taux d'erreur
	error_rate = rd.randrange(start=round(MIN_ERROR * 100), stop=round(MAX_ERROR * 100), step=1) / 100
	
	# Récupérer sa longueur
	len_seq = len(sequence)
	
	# Définir le nombre d'erreurs que contiendra la séquence
	n_errors = round(len_seq * error_rate)

	# Définir la position des erreurs
	errors = [0 if i < n_errors else 1 for i in range(len_seq)]
	rd.shuffle(errors)

	# Initialiser la futur séquence finale
	new_sequence = list(sequence)

	# Initialiser le seuil à zéro
	seuil = 0.0

	# Initialiser le tableau des indels
	pos_indel = []

	# Pour chaque nucléotide de la séquence
	for i in range(len_seq):
		# S'il est sur une position d'erreur
		if errors[i] == 0:
			# Générer un seuil
			seuil = rd.random()

			# Si le seuil est entre 0.5 et 1.0
			if 0.5 <= seuil < 1.0:
				# Si le seuil est entre 0.75 et 1.0 
				if 0.75 <= seuil < 1.0:
					# Ajouter une insertion à la liste des indels
					pos_indel.append((rd.choice(seq="ATCG"), i))

				# Sinon
				else:
					# Ajouter une délétion à la liste des indels
					pos_indel.append(('D', i))
			
			# Sinon
			else:
				# Remplacer le par un autre nucléotide
				if new_sequence[i] == 'A':
					new_sequence[i] = rd.choice(seq="TCG")
				elif new_sequence[i] == 'T':
					new_sequence[i] = rd.choice(seq="ACG")
				elif new_sequence[i] == 'C':
					new_sequence[i] = rd.choice(seq="ATG")
				else:
					new_sequence[i] = rd.choice(seq="ATC")
	
	# Transforme la liste en string
	new_sequence = ''.join(new_sequence)
	len_forwards = len(FORWARD)
	len_reverses = len(REVERSE)

	# Pour chaque indel dans la liste des indels
	for indel, pos in pos_indel:
		# S'il s'agit d'une insertion
		if indel in ('A', 'T', 'C', 'G'):
			# Si la taille de la nouvelle séquence n'excède pas la limite
			if len(new_sequence) + 1 < MAX_LENGTH - len_forwards - len_reverses:
				# Ajouter l'insertion à la séquence
				new_sequence = new_sequence[:pos] + indel + new_sequence[pos:]

		# Sinon
		else:
			# Ajouter la délétion
			new_sequence = new_sequence[:pos] + new_sequence[pos + 1:]

	# Retourner la nouvelle séquence
	return new_sequence

def gnrErronedAmplicons(sequence: str) -> List[str]:
	""" 
	Génère une liste de séquences avec des erreurs, à partir d'une séquence modèle.

	Args:
		sequence (str): La séquence nucléotidique modèle à partir de laquelle générer les séquences avec erreurs.

	Returns:
		List[str]: Une liste de séquences nucléotidiques résultantes avec des erreurs.
	"""
	# Retourner la liste des séquences
	return [FORWARD + addError(sequence=sequence) + REVERSE for _ in range(N_SEQUENCES)]

def generateSpeciesAmplicons(outputFile: str) -> None:
	""" 
	Crée les données de séquençage factices.

	Args:
		outputFile (str): Le chemin du fichier de sortie où les données de séquençage seront enregistrées.

	Returns:
		None
	"""
	# Initialiser les variables
	dataset = pd.DataFrame(columns=[i for i in range(MAX_LENGTH)])
	specie = ''
	start = 0
	end = 0
	sequence = ''
	sequences = []
	
	# Initialiser une graine aléatoire pour la reproductibilité des résultats
	rd.seed(RANDOM_STATE)

	# Pour chaque espèce
	for i, specie in enumerate(SPECIES):
		# Récupérer les position pour générer l'amplicon de PCR
		start, end = AMPLICON_POS[i][0], AMPLICON_POS[i][1]

		# Récupérer l'amplicon
		sequence = readFasta(filename=f"{PATH}/{specie}")[0][start:end]
		
		# Générer les séquences factices
		sequences = gnrErronedAmplicons(sequence=sequence)

		# Ajouter les séquences au dataframe
		sequences = init_dataFrame(sequences=sequences, n_features=MAX_LENGTH).assign(species=specie.split(sep='.')[0])
		sequences.columns = [i for i in range(MAX_LENGTH + 1)]
		dataset = pd.concat(objs=[dataset, sequences], axis=0).reset_index(drop=True)

	# Sauver le dataset de séquences factices
	dataset.to_csv(path_or_buf=outputFile, index=False)

################################### SIMULATED SPECIES #####################################

if __name__ == "__main__":
	# Définir la graine aléatoire
	rd.seed(RANDOM_STATE)

################################## UNIDENTIFIED SPECIES ###################################

	generateUnidentifiedAmplicons(outputFile=f"{RES_DIRECTORY}/fakes_unidentified_amplicons.csv")

################################## UNIDENTIFIED SPECIES ###################################


################################### SIMULATED SPECIES #####################################

	generateSpeciesAmplicons(outputFile=f"{RES_DIRECTORY}/fakes_species_amplicons.csv")

################################### SIMULATED SPECIES #####################################