#!/bin/python3
#-*- coding:utf-8 -*-


################################# IMPORTATION DES MODULES #################################

import pandas as pd
import subprocess as sp
import typing as tp

################################# IMPORTATION DES MODULES #################################


################################### VARIABLES GLOBALES ####################################

cherif: tp.Dict[str, list] = {"forwards" : ["GCGGTAATTCCAGCTCCA"],
							  "reverses" : ["CACTCAATTAAGCGCACACG"]}
bijection: tp.Dict[str, int] = {'n' : 0,
								'a' : 1,
								't' : 2,
								'c' : 3,
								'g' : 4}
invBijection: tp.Dict[int, str] = {0 : 'n',
								   1 : 'a',
								   2 : 't',
								   3 : 'c',
								   4 : 'g'}

################################### VARIABLES GLOBALES ####################################


####################################### FONCTIONS #########################################

# Exécute des commandes dans le shell
exec_cmd: tp.Callable[[str], None] = lambda cmd : sp.run(args=[cmd], shell=True)


def readFasta(filename: str) -> tp.List[str] :

	""" Récupère les séquences d'un fichier fasta. """

	# Initialiser les variables 
	sequences: tp.List[str] = []
	fasta: tp.TextIO[str] = ''
	current_line: str = ''
	current_seq: str = ''

	# Ouvrir le fichier
	with open(file=filename, mode='r') as fasta :

		# Récupérer la ligne courrante
		current_line = fasta.readline().rstrip()

		# Initialiser la séquence courrante
		current_seq = ''

		# Tant que l'on n'a pas parcouru tout le fichier
		while len(current_line) != 0 :

			# Si la première lettre de la ligne n'est pas un nucléotide
			if current_line[0] not in ('>', '\n', '\t', '\r', '', ' ') :

				# Ajouter la ligne à la séquence courrante
				current_seq+=current_line
			
			# Sinon
			else :

				# Si on n'est pas à la fin du fichier
				if len(current_seq) != 0 :

					# Ajouter la séquence courrante à la liste des séquences
					sequences.append(current_seq)

					# Réinitialiser la séquence courrante
					current_seq = ''
			
			# Récupérer la nouvelle ligne courrante
			current_line = fasta.readline().rstrip()

			# S'il n'y a plus de lignes dans le fichier, mais que la séquence courrante contient une séquence
			if len(current_line) == 0 and len(current_seq) != 0 :

				# Ajouter la séquence à la liste des séquences
				sequences.append(current_seq)

	# Retourner les séquences sous forme d'une liste triée
	return sorted(sequences)


def convert_sequence(sequence: str) -> tp.List[int] :
	
	""" Convertie une séquence nucléotidique en liste d'entiers. """

	# Récupérer la variable globale "bijection"
	global bijection

	# Récupérer la liste d'entier
	return [bijection[nuc.lower()] for nuc in list(sequence)]


def init_dataFrame(sequences: tp.List[str], n_features: int) -> pd.DataFrame :

	""" Initialiser un dataframe, à partir d'une liste de séquences. """  

	# Retourner le dataframe
	return pd.DataFrame(data=[convert_sequence(sequence=sequence+"n"*(n_features-len(sequence))) for sequence in sequences])


def reconstructSequence(data: pd.DataFrame, index: int) -> str :

	""" Reconstruit la séquence d'une ligne du dataset. """

	# Récupérer la variable globale "invBijection"
	global invBijection

	# Retourner la séquence reconstituée
	return ''.join([invBijection[value] for value in data.iloc[index, :] if value in (1, 2, 3, 4)])

####################################### FONCTIONS #########################################