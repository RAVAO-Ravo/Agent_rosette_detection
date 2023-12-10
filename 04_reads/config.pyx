#-*- coding:utf-8 -*-
#cython:language_level=3

# Importation des modules
import numpy as np


cdef class GlobalValue(object):
	"""
	Classe représentant des valeurs globales.

	Attributes:
		PRIMERS(dict): Les primers disponibles.
		RANGES(dict): Les ranges pour chaque primer disponible.
		ENCODER (dict): L'encodage des espèces.
		COLUMNS (dict): Les éléments composants un read.
		CLF_FILE (str): Fichier contenant le classifier.
		NONE (int): Valeur par défaut pour 'NONE'.

	Methods:
		get_values(): Retourne une valeur globale, en fonction du nom passé en paramètre.
	"""
	# Initialisation des variables
	cdef dict PRIMERS
	cdef dict RANGES
	cdef dict ENCODER
	cdef dict COLUMNS
	cdef str CLF_FILE
	cdef int NONE

	def __init__(self) -> None:
		"""
		Construit une instance de la classe GlobalValue.
		"""
		# Les primers disponibles
		self.PRIMERS = {
			"CHERIF" : {
				"forwards": ["GCGGTAATTCCAGCTCCA"],
				"reverses": ["CACTCAATTAAGCGCACACG"]
			},

			"CARPIO" : {
				"forwards" : ["GGTGGGTTTCAGTAGACAATGC"],
				"reverses" : ["GGCGGCAATAACAAATGGTAGT"]
			},

			"GOZLAN" : {
				"forwards" : ["AATCGTATGACATTTTGTCGAC", "CAGGGCTTTTTAAGTCTT"],
				"reverses" : ["GAAGTCACAGGCGATTCGG", "TGATGGAGTCATAGAATTAACATCC"]
			},

			"PARVA" : {
				"forwards" : ["CGAGCCAAATAACAGAGGGT"],
				"reverses" : ["CAGGCGAGGCTTATGTTTGC"]
			}
		}

		# Les ranges pour chaque primer disponible
		self.RANGES = {
			"CHERIF" : (150, 160),
			"CARPIO" : (70, 85),
			"GOZLAN" : (580, 1120),
			"PARVA"  : (140, 155)
		}

		# Encodage des différentes espèces
		self.ENCODER = {
						"Anurofeca_richardsi": 0, "Dermocystidium_salmonis": 1, "Ichthyophonus_hoferi": 2, 
						"Pseudoperkinsus_tapetis": 3, "Psorospermium_haeckelii": 4, "Rhinosporidium_cygnus": 5,
						"Rhinosporidium_seeberi": 6, "Sphaeroforma_spB7": 7, "Sphaeroforma_spCRG3": 8,
						"Sphaerothecum_destruens": 9, "unidentified": 10
						}

		# Les éléments composants un read
		self.COLUMNS = {
			key: value for value, key in enumerate([
				"id", "sequence", "quality",
				"prefix", "suffix", "length",
				'A', 'T', 'C', 'G', 'N',
				"GC", "AT", "Q_Score",
				"std_QScore", "q25", "q50",
				"q75", "q90", "ScoreDist"
			])
		}

		# Fichier contenant le classifier
		self.CLF_FILE = "./clf_ra.sav"

		# Valeur par défaut pour 'none'
		self.NONE = (-99)

	cpdef object get_values(self, str name):
		"""
		Récupère les valeurs associées à un nom spécifique.

		Args:
			name (str): Le nom de la valeur à récupérer ("primers", "ranges", "encoder", "columns", "clf_file", "none").

		Returns:
			object: La valeur associée au nom spécifié.

		Raises:
			ValueError: Si le nom spécifié n'est pas valide.
		"""
		if name == "primers":
			return self.PRIMERS
		elif name == "ranges":
			return self.RANGES
		elif name == "encoder":
			return self.ENCODER
		elif name == "columns":
			return self.COLUMNS
		elif name == "clf_file":
			return self.CLF_FILE
		elif name == "none":
			return self.NONE
		else:
			raise ValueError(f"Nom de valeur non valide : {name}")


cpdef str getComplementarySeq(str sequence):
	"""
	Retourne la séquence complémentaire d'une séquence nucléotidique.

	Args:
		sequence (str): La séquence nucléotidique d'origine.

	Returns:
		str: La séquence complémentaire.
	"""
	# Initialiser les variables
	cdef str complementary = ''
	cdef str nuc = ''

	# Pour chaque base de la séquence nucléotique (itérée dans le sens inverse)
	for nuc in sequence[len(sequence)::-1]:
		# Récupérer le complémentaire de la base
		if nuc == 'T':
			complementary += 'A'
		elif nuc == 't':
			complementary += 'a'
		elif nuc == 'A':
			complementary += 'T'
		elif nuc == 'a':
			complementary += 't'
		elif nuc == 'G':
			complementary += 'C'
		elif nuc == 'g':
			complementary += 'c'
		elif nuc == 'C':
			complementary += 'G'
		else:
			complementary += 'g'

	# Retourner le complémentaire
	return complementary

cpdef int encode(str nuc):
	"""
	Récupère l'encodage d'un nucléotide.

	Args:
		nuc (str): Le nucléotide à encoder (A, T, C ou G, en majuscules ou minuscules).

	Returns:
		int: L'encodage correspondant au nucléotide (1 pour 'A', 2 pour 'T', 3 pour 'C', 4 pour 'G', 0 sinon).
	"""
	# Retourner la valeur
	return 1 if nuc.lower() == 'a' else \
 		   2 if nuc.lower() == 't' else \
 		   3 if nuc.lower() == 'c' else \
 		   4 if nuc.lower() == 'g' else \
 		   0

cpdef list transformer(str sequence, int n_features):
	"""
	Transforme une séquence nucléotidique en un vecteur.

	Args:
		sequence (str): La séquence nucléotidique d'origine.
		n_features (int): Le nombre de caractéristiques (dimensions) du vecteur.

	Returns:
		list: Un vecteur basé sur la séquence nucléotidique.
	"""
	# Retourner le vecteur
	return [encode(nuc=nuc) for nuc in sequence + 'n' * (n_features - len(sequence))]