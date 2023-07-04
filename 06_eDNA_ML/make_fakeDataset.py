#!/bin/python3
#-*- coding:utf-8 -*-


################################# IMPORTATION DES MODULES #################################

import os
import pandas as pd
import random as rd
import typing as tp
from argparse import ArgumentParser, HelpFormatter
from utilitaire import cherif, init_dataFrame, readFasta, exec_cmd

################################# IMPORTATION DES MODULES #################################


################################# ARGUMENTS DU PROGRAMME ##################################

formatter: tp.Callable[[str], HelpFormatter] = lambda prog : HelpFormatter(prog=prog, max_help_position=100, width=200)
parser = ArgumentParser(description="Description: Ce programme permet de générer un dataset de données (de séquençages) simulées, pour les espèces du projet eDNA.",
						formatter_class=formatter)
parser.add_argument("-r", "--random_state", type=int, default=42, help="value of the seed for the reproductibility")
parser.add_argument("-n", "--n_sequences", type=int, default=10000, help="number of sequences to generate for each specie")
parser.add_argument("-m", "--min_length", type=int, default=150, help="minimum length of sequences")
parser.add_argument("--max_length", type=int, default=160, help="minimum length of sequences")
parser.add_argument("--min_error", type=float, default=0.0, help="minimum error rate when generating a sequence")
parser.add_argument("--max_error", type=float, default=0.15, help="maximum error rate when generating a sequence")
args = parser.parse_args()

################################# ARGUMENTS DU PROGRAMME ##################################


################################### VARIABLES GLOBALES ####################################

random_state: int = args.random_state
n_sequences: int = args.n_sequences
min_length: int = args.min_length
max_length: int = args.max_length
min_error: float = args.min_error
max_error: float = args.max_error
reverses: str =  cherif["reverses"][0]
forwards: str = cherif["forwards"][0]
path: str = "./species/"
species: tp.List[str] = sorted(os.listdir(path=path))
amplicon_pos: tp.List[tp.Tuple[int, int]] = [(595, 709),
											 (583, 702),
											 (542, 654),
											 (566, 680),
											 (585, 699),
											 (554, 675),
											 (564, 686),
											 (490, 604),
											 (566, 680),
											 (577, 693)]
res_directory: str = "dataset"

################################### VARIABLES GLOBALES ####################################


################################## UNIDENTIFIED SPECIES ###################################

def random_sequence() -> str : 

	""" Génère un amplicon aléatoire d'ADN. """

	# Initialiser les variables globales
	global forwards
	global reverses
	global min_length
	global max_length

	# Générer la séquence entre les primers
	size_seq = rd.randrange(start=min_length, stop=max_length, step=1)-len(forwards)-len(reverses)

	# Retourner l'amplicon
	return forwards+''.join([rd.choice("ATCG") for _ in range(0, size_seq, 1)])+reverses


def generateUnidentifiedAmplicons(outFile: str) -> None :

	""" Génére le dataset des séquences 'non-identifiées'. """

	# Initialiser les variables globales
	global res_directory
	global random_state
	global max_length
	
	# Créer le répertoire résultat s'il n'existe pas
	if res_directory not in os.listdir(path='.') :
		exec_cmd(cmd=f"mkdir {res_directory}")

	# Initialiser un set pour s'assurer que chaque séquence est unique
	sequences = set()

	# Tant que l'on n'a pas le nombre de séquences souhaités
	while len(sequences) != n_sequences :

		# Créer une séquence aléatoire et ajouter la au set de séquences
		sequences.add(random_sequence())
	
	# Transformer le set en dataframe
	dataset = init_dataFrame(sequences=sorted(sequences), n_features=max_length).assign(species="unidentified")
	dataset.columns = [i for i in range(0, max_length+1, 1)]
	dataset.to_csv(path_or_buf=outFile, index=False)

################################## UNIDENTIFIED SPECIES ###################################


################################### SIMULATED SPECIES #####################################

def addError(sequence: str) -> str :

	""" Génère une séquence nucléotidique avec erreurs, à partir d'une séquence modèle. """

	# Initialiser variables globales
	global max_error
	global min_error
	global max_length
	global forwards
	global reverses

	# Définir le taux d'erreur
	error_rate = rd.randrange(start=round(min_error*100), stop=round(max_error*100), step=1)/100
	
	# Récupérer sa longueur
	len_seq = len(sequence)
	
	# Définir le nombre d'erreurs que contiendra la séquence
	n_errors = round(len_seq*error_rate)

	# Définir la position des erreurs
	errors = [0 if i < n_errors else 1 for i in range(0, len_seq, 1)]
	rd.shuffle(errors)

	# Initialiser la futur séquence finale
	new_sequence = list(sequence)

	# Initialiser le seuil à zéro
	seuil = 0.0

	# Initialiser le tableau des indels
	pos_indel = []

	# Pour chaque nucléotide de la séquence
	for i in range(0, len_seq, 1) :

		# S'il est sur une position d'erreur
		if errors[i] == 0 :
			
			# Générer un seuil
			seuil = rd.random()

			# S'il le seuil est entre 0.5 et 1.0
			if 0.5 <= seuil < 1.0 :
				
				# Si le seuil est entre 0.75 et 1.0 
				if 0.75 <= seuil < 1.0 :
					
					# Ajouter une insertion à la liste des indels
					pos_indel.append((rd.choice(seq="ATCG"), i))

				# Sinon
				else :

					# Ajouter une délétion à la liste des indels
					pos_indel.append(('D', i))
			
			# Sinon
			else :

				# Remplacer le par un autre nucléotide
				if new_sequence[i] == 'A' :
					new_sequence[i] = rd.choice(seq="TCG")
				elif new_sequence[i] == 'T' :
					new_sequence[i] = rd.choice(seq="ACG")
				elif new_sequence[i] == 'C' :
					new_sequence[i] = rd.choice(seq="ATG")
				else :
					new_sequence[i] = rd.choice(seq="ATC")
	
	# Transforme la liste en string
	new_sequence = ''.join(new_sequence)
	len_forwards = len(forwards)
	len_reverses = len(reverses)

	# Pour chaque indel dans la liste des indels
	for indel, pos in pos_indel :
		
		# S'il s'agit d'une insertion
		if indel in ('A', 'T', 'C', 'G') :
			
			# Si la taille de la nouvelle séquence n'excède pas la limite
			if len(new_sequence)+1 < max_length-len_forwards-len_reverses :

				# Ajouter l'insertion à la séquence
				new_sequence = new_sequence[:pos]+indel+new_sequence[pos:]
		
		# Sinon
		else :
			
			# Ajouter une délétion
			if pos != 0 :
				new_sequence = new_sequence[:pos]+new_sequence[pos:]
			else :
				new_sequence = new_sequence[1:]

	# Retourner la nouvelle séquence
	return new_sequence


def gnrErronedAmplicons(sequence: str) -> tp.List[str] :
	
	""" Génère une liste de séquences avec erreurs, à partir d'une séquence modèle."""

	# Initialiser les variables globales
	global n_sequences
	global forwards
	global reverses
	
	# Retourner la liste des séquences
	return [forwards+addError(sequence=sequence)+reverses for _ in range(0, n_sequences, 1)]


def generateSpeciesAmplicons(outFile: str) -> None :

	""" Crée les données de séquençage factices. """

	# Initialiser les variables globales
	global species
	global amplicon_pos
	global forwards
	global reverses
	global n_sample
	global max_length
	global path
	global random_state

	# Initialiser les variables
	dataset = pd.DataFrame(columns=[i for i in range(0, max_length, 1)])
	i = 0
	specie = ''
	start = 0
	end = 0
	sequence = ''
	sequences = []
	
	# Initialiser une graine aléatoire pour la reproductibilité des résultats
	rd.seed(random_state)

	# Pour chaque espèce
	for i, specie in enumerate(species) :
		
		# Récupérer les position pour générer l'amplicon de PCR
		start, end = amplicon_pos[i][0], amplicon_pos[i][1]

		# Récupérer l'amplicon
		sequence = readFasta(filename=f"{path}/{specie}")[0][start:end]
		
		# Générer les séquences factices
		sequences = gnrErronedAmplicons(sequence=sequence)

		# Ajouter les séquences au dataframe
		sequences = init_dataFrame(sequences=sequences, n_features=max_length).assign(species=specie.split(sep='.')[0])
		sequences.columns = [i for i in range(0, max_length+1, 1)]
		dataset = pd.concat(objs=[dataset, sequences], axis=0).reset_index(drop=True)

	# Sauver le dataset de séquences factices
	dataset.to_csv(path_or_buf=outFile, index=False)

################################### SIMULATED SPECIES #####################################


if __name__ == "__main__" :

	# Définir la graine aléatoire
	rd.seed(random_state)

################################## UNIDENTIFIED SPECIES ###################################

	generateUnidentifiedAmplicons(outFile=f"{res_directory}/fakes_unidentified_amplicons.csv")

################################## UNIDENTIFIED SPECIES ###################################


################################### SIMULATED SPECIES #####################################

	generateSpeciesAmplicons(outFile=f"{res_directory}/fakes_species_amplicons.csv")

################################### SIMULATED SPECIES #####################################