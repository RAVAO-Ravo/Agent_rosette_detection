#-*- coding:utf-8 -*-
#cython:language_level=3

# Importation des variables
import numpy as np
import pandas as pd
import random as rd
from tqdm import tqdm

# Module config
from config import (
	GlobalValue,
	getComplementarySeq
)

# Module AhoCorasick
from ahocorasick import (
	AhoCorasick,
	getPosPrimer
)

# Module read
from read import (
	readFastq,
	createRead,
	headcroping,
	tailcroping,
	makeKmers
)

# Variables globales
GVALUES: GlobalValue = GlobalValue()
cdef dict COLUMNS = GVALUES.get_values("columns")
cdef str CLF_FILE = GVALUES.get_values("clf_file")
cdef int NONE = GVALUES.get_values("none")


cdef class Reads(object):
	"""
	Classe Reads modélisant un ensemble de reads, et permettant de manipuler l'ensemble.

	Attributes:
		filename (str): Le nom du fichier contenant les reads.
		reads (list): La liste des reads extraits du fichier.
		n_reads (int): Le nombre de reads dans l'ensemble.

	Methods:
		Reads(): Instancie un objet Reads.
		getFileName(): Renvoie le nom du fichier.
		getReads(): Renvoie la liste de reads.
		getN_reads(): Renvoie le nombre de reads.
		setReads(reads): Met à jour la liste de reads.
		setN_reads(n_reads): Met à jour le nombre de reads.
		getColumn(): Renvoie la liste des valeurs contenus dans une colonne.
		to_numpy(): Retourne les statistiques des reads sous forme d'un tableau numpy.
		to_DataFrame(): Retourne les statistiques des reads sous forme d'un dataframe.
		processing_reads(): Retourne les reads avec une barre de chargement.
		filter_Duplicates(): Filtre les reads dupliqués.
		filter_size(): Filtre les reads en fonction de la taille.
		filter_Q(): Filtre les reads en fonction de la valeur moyenne du Q_score.
		filter_quality(): Filtre les reads en fonction de la qualité.
		filter_N(): Filtre les reads en fonction du nombre de bases indeterminées (Q_score < 15).
		filter_GC(): Filtre les reads en fonction du taux de GC.
		HTC(): Excise les reads en début et/ou en fin de séquence.
		cutPrimers(): Excise les reads en fonction de primers.
		multiFilter(): Filtre les reads en fonction de plusieurs paramètres.
		readsToKmers(): Infère tous les K-mers possibles à partir de tous les reads.
		write(): Permet de sauvegarder les reads.
		writeWithIdx(): Permet de sauvegarder les reads, en fonction d'une liste d'index.
		makeSample(): Permet d'échantillonner des reads à partir des reads contenus dans l'instance Reads courrante.
	"""
	# Initialiser les variables de classes
	cdef str filename
	cdef list reads
	cdef int n_reads

	def __init__(self, str filename):
		"""
		Construit une instance de la classe Reads à partir d'un fichier fastq.

		Args:
			filename (str): Le chemin du fichier contenant les reads.
		"""
		# Nom du fichier
		self.filename = filename

		# Ensemble des reads
		self.reads = list(readFastq(filename))

		# Nombre de reads
		self.n_reads = len(self.reads)

	cpdef str getFileName(self):
		"""
		Renvoie le nom du fichier.

		Returns:
			str: Le nom du fichier.
		"""
		return self.filename

	cpdef list getReads(self):
		"""
		Renvoie la liste de reads.

		Returns:
			List[Tuple[Any]]: La liste des reads.
		"""
		return self.reads

	cpdef int getN_reads(self):
		"""
		Renvoie le nombre de reads.

		Returns:
			int: Le nombre de reads.
		"""
		return self.n_reads

	cpdef void setReads(self, list reads):
		"""
		Met à jour la liste de reads.

		Args:
			List[Tuple[Any]]: La nouvelle liste de reads.
		"""
		self.reads = reads

	cpdef void setN_reads(self, int n_reads):
		"""
		Met à jour le nombre de reads.

		Args:
			n_reads (int): Le nouveau nombre de reads.
		"""
		self.n_reads = n_reads

	cpdef list getColumn(self, str column):
		"""
		Renvoie la liste des valeurs contenus dans une colonne.

		Returns:
			List[str]: La liste de la colonne demandée.
		"""
		# Vérifier si la colonne demandées est dans la liste des colonnes
		if column in COLUMNS.keys():
			return [read[COLUMNS[column]] for read in self.reads]
		else:
			# Retourner une erreur
			raise ValueError(f"'{column}' ne fait pas partir des valeurs autorisées, voici la liste :\n{COLUMNS.keys()}")
	
	def to_numpy(self) -> np.ndarray:
		"""
		Retourne les statistiques des reads sous forme d'un tableau numpy.

		Returns:
			np.ndarray: Un tableau numpy contenant les statistiques des reads.
		"""
		return np.array(self.reads)

	def to_DataFrame(self) -> pd.DataFrame:
		"""
		Retourne les statistiques des reads sous forme d'un dataframe.

		Returns:
			pd.DataFrame: Un dataframe contenant les statistiques des reads.
		"""
		return pd.DataFrame(data=self.reads, columns=COLUMNS.keys())

	cdef object processing_reads(self, str step_name = "Processing"):
		"""
		Retourne les reads pour un traitement, en y ajoutant une barre de chargement.

		Args:
			step_name (str): Nom de l'étape.

		Returns:
			tqdm: Les reads avec une barre de chargement.
		"""
		return tqdm(iterable=self.reads, desc=step_name, unit=" Reads")

	cpdef Reads filter_Duplicates(self):
		"""
		Supprime les éléments en double (ou plus) de l'ensemble des reads.

		Returns:
			Reads: L'instance Reads sans les reads dupliqués.
		"""
		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []
		cdef str prefix = ''
		cdef set unique_ = set()
		
		# Pour chaque read
		for read in self.processing_reads("Filter Duplicates"):
			# Récupérer son préfixe
			prefix = read[COLUMNS["prefix"]]

			# Si le préfixe n'est pas dans la liste unique_
			if prefix not in unique_:
				# Ajouter le préfixe à unique_
				unique_.add(prefix)

				# Ajouter le read à nouvelle liste des reads
				new_reads.append(read)
			
			# Sinon, ne pas y ajouter le read 
			else:
				self.n_reads -= 1

		# Mettre à jour les reads
		self.reads = new_reads
		
		# Retourner l'objet
		return self

	cpdef Reads filter_size(self, float inf = 0.0, float sup = 1.0):
		"""
		Filtrer les reads en fonction de leur taille.

		Args:
			inf (float): La limite inférieure de la taille de read.
			sup (float): La limite supérieure de la taille de read.

		Returns:
			Reads: L'instance Reads avec les tailles filtrées.
		"""
		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []

		# Récupérer les valeurs des quantiles
		if (0.0 <= inf and inf < 1.0) or (0.0 < sup and sup <= 1.0):
			lengths: np.ndarray = np.array([read[COLUMNS["length"]] for read in self.reads], dtype=float)
			if 0.0 <= inf and inf < 1.0:
				inf = np.quantile(a=lengths, q=inf)
			if 0.0 < sup and sup <= 1.0:
				sup = np.quantile(a=lengths, q=sup)

		# Récupérer les reads avec une taille comprise entre les bornes
		for read in self.processing_reads("Filter Size"):
			if inf <= read[COLUMNS["length"]] and read[COLUMNS["length"]] <= sup:
				new_reads.append(read)
			else:
				self.n_reads -= 1
	
		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self

	cpdef Reads filter_Q(self, float seuil = 20.0):
		"""
		Filtrer les reads en fonction de leur score Q.

		Args:
			seuil (float): Le score Q minimum pour conserver les reads.

		Returns:
			Reads: L'instance Reads avec les reads filtrés par la valeur de Q.
		"""
		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []

		# Récupérer les reads dont le score est supérieur à la borne minorante
		for read in self.processing_reads("Filter Q"):
			if seuil <= read[COLUMNS["Q_Score"]]:
				new_reads.append(read)
			else:
				self.n_reads -= 1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self

	cpdef Reads filter_quality(self, float seuil = 1.0):
		"""
		Filtrer les reads en fonction de leur qualité.

		Args:
			seuil (float): La limite supérieure du score de qualité des reads à conserver.

		Returns:
			Reads: L'instance Reads avec les reads filtrés par la qualité.
		"""
		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []

		# Récupérer la valeur du quantile
		scores: np.ndarray = np.array([read[COLUMNS["ScoreDist"]] for read in self.reads], dtype=float)
		seuil = np.quantile(a=scores, q=seuil)

		# Récupérer les reads dont le score est inférieur à la borne majorante
		for read in self.processing_reads("Filter Quality"):
			if read[COLUMNS["ScoreDist"]] <= seuil:
				new_reads.append(read)
			else:
				self.n_reads -= 1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self

	cpdef Reads filter_N(self, float seuil = 1.0):
		"""
		Filtrer les reads en fonction de la quantité de bases indéterminées (Q_Score(bases) < 15).

		Args:
			seuil (float): La limite supérieure de la quantité de bases indéterminées pour les reads à conserver.

		Returns:
			Reads: L'instance Reads avec les reads filtrés par la quantité de bases indéterminées.
		"""
		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []

		# Récupérer la valeur du quantile
		N_content: np.ndarray = np.array([read[COLUMNS['N']] for read in self.reads], dtype=float)
		seuil = np.quantile(a=N_content, q=seuil)

		# Récupérer les reads ayant les plus faibles quantités de bases indécises
		for read in self.processing_reads("Filter N"):
			if read[COLUMNS['N']] <= seuil:
				new_reads.append(read)
			else:
				self.n_reads -= 1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self

	cpdef Reads filter_GC(self, float minGC = 0.0, float maxGC = 1.0):
		"""
		Filtrer les reads en fonction de leur contenu en GC.

		Args:
			minGC (float): La teneur minimale en GC pour les reads à conserver.
			maxGC (float): La teneur maximale en GC pour les reads à conserver.

		Returns:
			Reads: L'instance Reads avec les reads filtrés par la quantitée de GC.
		"""
		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []

		# Récupérer les reads ayant un taux de GC contenu dans les bornes
		for read in self.processing_reads("Filter GC"):
			if minGC <= read[COLUMNS["GC"]] and read[COLUMNS["GC"]] <= maxGC:
				new_reads.append(read)
			else:
				self.n_reads -= 1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self

	cpdef Reads HTC(self, float head = 0.0, float tail = 0.0):
		"""
		Rogner les extrémités de chaque read d'un certain nombre de bases.

		Args:
			head (float): Le nombre de bases à rogner depuis le début de chaque read.
			tail (float): Le nombre de bases à rogner depuis la fin de chaque read.

		Returns:
			Reads: L'instance Reads avec les reads excisés.
		"""
		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []
		cdef str id_ = ''
		cdef str sequence = ''
		cdef str quality = ''
		cdef float len_seq = 0.0

		# Amputer chaque read, selon les paramètres définis
		for read in self.processing_reads("Head/Tail Croping"):
			# Récupérer les éléments nécéssaire du read
			id_ = read[COLUMNS["id"]]
			sequence = read[COLUMNS["sequence"]]
			quality = read[COLUMNS["quality"]]
			len_seq = float(read[COLUMNS["length"]])

			# Exciser au début du read
			if 0.0 < head and head < 1.0:
				sequence, quality = headcroping(sequence=sequence, quality=quality, value=(len_seq*head))
			elif 1.0 <= head:
				sequence, quality = headcroping(sequence=sequence, quality=quality, value=head)

			# Exciser à la fin du read
			if 0.0 < tail and tail < 1.0:
				sequence, quality = tailcroping(sequence=sequence, quality=quality, value=(len_seq*tail))
			elif 1.0 <= tail:
				sequence, quality = tailcroping(sequence=sequence, quality=quality, value=tail)

			# Créer le nouveau read
			if len(sequence) != 0:
				new_reads.append(createRead((id_, sequence, quality)))
			else:
				self.n_reads -= 1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self

	cpdef Reads cutPrimers(self, list forwards, list reverses, int seuil_length = 1):
		"""
		Filtrer les reads en fonction des séquences d'amorces.

		Args:
			forwards (list): Liste des amorces forwards.
			reverses (list): Liste des amorces reverses.
			seuil_length (int): La longueur minimale d'une read à conserver.

		Returns:
			Reads: L'instance Reads avec les reads excisés au niveau des primers.
		"""
		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []
		cdef str id_ = ''
		cdef str sequence = ''
		cdef str quality = ''
		
		# Créer les automates de recherche
		reverses = [getComplementarySeq(sequence=sequence) for sequence in reverses]
		forwardMachine: AhoCorasick = AhoCorasick(words=forwards)
		reverseMachine: AhoCorasick = AhoCorasick(words=reverses)

		# Récupérer les reads ayant les primers
		for read in self.processing_reads("CutPrimers"):
			# Récupérer la séquence du read
			id_ = read[COLUMNS["id"]]
			sequence =  read[COLUMNS["sequence"]]
			quality = read[COLUMNS["quality"]]

			# Récupérer la postion du forward et du reverse
			min_idx = getPosPrimer(sequence=sequence, searchMachine=forwardMachine, primers=forwards, forward=True)
			max_idx = getPosPrimer(sequence=sequence, searchMachine=reverseMachine, primers=reverses, forward=False)

			# Si les deux primers sont présents, dans le bon sens
			if (min_idx != NONE) and (max_idx != NONE) and (min_idx < max_idx):
				# Couper au niveau des primers
				sequence = sequence[min_idx:max_idx]
				quality = quality[min_idx:max_idx]

				# Si la taille est suffisante, récupérer le read
				if seuil_length <= len(sequence):
					new_reads.append(createRead((id_, sequence, quality)))

				else:
					self.n_reads -= 1
			else:
				self.n_reads -= 1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self

	cpdef Reads multiFilter(self,
							float size_inf = 0.0,
							float size_sup = 1.0,
							float seuil_qual = 1.0,
							float seuil_N = 1.0,
							float head = 0.0,
							float tail = 0.0,
							float minGC = 0.0,
							float maxGC = 1.0,
							bint  duplicated = False):
		"""
		Filtre les reads en fonction de plusieurs paramètres.

		Args:
			size_inf (float): La limite inférieure de la taille des reads.
			size_sup (float): La limite supérieure de la taille des reads.
			seuil_qual (float): La limite supérieure du score de qualité pour les reads.
			seuil_N (float): La limite supérieure de la quantité de bases indéterminées pour les reads.
			head (float): Le nombre de bases à couper du début de chaque read.
			tail (float): Le nombre de bases à couper de la fin de chaque read.
			minGC (float): La teneur minimale en GC pour les reads.
			maxGC (float): La teneur maximale en GC pour les reads.
			duplicated (bool): Si True, exclut les reads en double.

		Returns:
			Reads: L'instance Reads avec les reads filtrés en fonction des paramètres.
		"""
		# Initialiser les variables
		df: pd.DataFrame  = self.to_DataFrame().loc[:, ("length", "ScoreDist", "N", "GC")]
		cdef list new_reads = []
		cdef tuple read = ()
		cdef float len_seq = 0.0

		# Récupérer les valeurs de quantiles pour la taille des reads
		if (0.0 <= size_inf and size_inf < 1.0) or (0.0 < size_sup and size_sup <= 1.0):
			if 0.0 <= size_inf and size_inf < 1.0:
				size_inf = df["length"].quantile(q=size_inf)
			if 0.0 < size_sup and size_sup <= 1.0:
				size_sup = df["length"].quantile(q=size_sup)

		# Récupérer les valeurs de quantiles pour la qualité
		seuil_qual = df["ScoreDist"].quantile(q=seuil_qual)
		
		# Récupérer les valeurs de quantiles pour le pourcentage de bases indécises
		seuil_N = df['N'].quantile(q=seuil_N)

		# Récupérer les reads selon les filtres
		for read in self.processing_reads("MultiFilter"):
			if (size_inf <= read[COLUMNS["length"]]) and (read[COLUMNS["length"]] <= size_sup) and \
			   (read[COLUMNS["ScoreDist"]] <= seuil_qual) and (read[COLUMNS['N']] <= seuil_N) and \
			   (minGC <= read[COLUMNS["GC"]]) and (read[COLUMNS["GC"]] <= maxGC):
				
				id_ = read[COLUMNS["id"]]
				sequence = read[COLUMNS["sequence"]]
				quality = read[COLUMNS["quality"]]
				len_seq = read[COLUMNS["length"]]

				if 0.0 < head and head < 1.0:
					sequence, quality = headcroping(sequence=sequence, quality=quality, value=(len_seq*head))
				elif 1.0 <= head:
					sequence, quality = headcroping(sequence=sequence, quality=quality, value=head)

				if 0.0 < tail and tail < 1.0:
					sequence, quality = tailcroping(sequence=sequence, quality=quality, value=(len_seq*tail))
				elif 1.0 <= tail:
					sequence, quality = tailcroping(sequence=sequence, quality=quality, value=tail)

				if len(sequence) != 0:
					new_reads.append(createRead((id_, sequence, quality)))

				else:
					self.n_reads -= 1
			else:
				self.n_reads -= 1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Filtrer les duplications
		if duplicated == True:
			self.filter_Duplicates()

		# Retourner l'objet
		return self

	cpdef list readsToKmers(self, int k = 21):
		"""
		Crée tout les k-mers possibles.

		Args:
			k (int): La longueur des k-mers. (Par défaut = 21).

		Returns:
			List[str]: Une liste de tout les k-mers.
		"""
		# Retourner la liste triée des kmers
		return sorted(set([kmer for read in self.reads for kmer in makeKmers(sequence=read[COLUMNS["sequence"]], k=k)]))

	cpdef void write(self, str filename, str outputFormat = "fq"):
		"""
		Enregistre les reads dans un fichier dans un format biologique spécifié.

		Args:
			filename (str): Le nom du fichier de sortie.
			outputFormat (str): Le format de sortie, 'fq' pour fastq, 'fa' pour fasta.

		Returns:
			None
		"""
		# Initialiser les variables
		cdef tuple read = ()

		# Créer le fichier de sortie
		with open(file=filename, mode='w') as outputFile:
			# Si l'on souhaite une sortie fastq
			if outputFormat == "fq":
				# Enregistrer selon le format fastq
				for read in self.reads:
					outputFile.write(f"""{read[COLUMNS["id"]]}\n{read[COLUMNS["sequence"]]}\n+\n{read[COLUMNS["quality"]]}\n""")

			# Si l'on souhaite une sortie fastq
			elif outputFormat == "fa":
				# Enregistrer selon le format fasta
				for read in self.reads:
					outputFile.write(f""">{read[COLUMNS["id"]][1:]}\n{read[COLUMNS["sequence"]]}\n""")

			# Sinon
			else:
				# Ne pas enregistrer le fichier
				print("Format inexistant, ou pas pris en compte. L'enregistrement est abordé.")

	cpdef void writeWithIdx(self, str filename, list idx_, str outputFormat = "fq"):
		"""
		Enregistre les reads sélectionnées dans un fichier en fonction de leurs identifiants.

		Args:
			filename (str): Le nom du fichier de sortie.
			idx_ (list): Liste d'identifiants de reads à enregistrer.
			outputFormat (str): Le format de sortie, 'fq' pour fastq, 'fa' pour fasta.

		Returns:
			None
		"""
		# Initialiser les variables
		cdef tuple read = ()
		cdef str id_ = ''

		# Créer le fichier de sortie
		with open(file=filename, mode='w') as outputFile:
			# Enregistrer les reads sélectionnés, selon le format fastq
			if outputFormat == "fq":
				for read in self.reads:
					id_ = read[COLUMNS["id"]]
					if id_ in idx_:
						outputFile.write(f"""{id_}\n{read[COLUMNS["sequence"]]}\n+\n{read[COLUMNS["quality"]]}\n""")

			# Enregistrer les reads sélectionnés, selon le format fasta
			elif outputFormat == "fa":
				for read in self.reads:
					id_ = read[COLUMNS["id"]]
					if id_ in idx_:
						outputFile.write(f""">{id_[1:]}\n{read[COLUMNS["sequence"]]}\n""")

			# Sinon
			else:
				# Ne pas enregistrer le fichier
				print("Format inexistant, ou pas pris en compte. L'enregistrement est abordé.")

	cpdef void makeSample(self, int n_sample, str filename = '', str outputFormat = "fq", object random_state = None):
		"""
		Sous-échantillonne des reads à partir de l'ensemble complet de reads.

		Args:
			n_sample (int): Nombre de reads à sous-échantillonner.
			filename (str): Le nom du fichier de sortie. Si non fourni, le nom de fichier est généré automatiquement.
			outputFormat (str): Le format de sortie, 'fq' pour fastq, 'fa' pour fasta.
			random_state (object): Graine aléatoire pour la reproductibilité.

		Returns:
			None
		"""
		# Définir un graine aléatoire
		rd.seed(random_state)

		# Sélectionner les reads
		cdef list selected_reads = self.getColumn("id")
		rd.shuffle(selected_reads)

		# Enregistrer le sous-échantillon
		if filename == '':
			self.writeWithIdx(filename=self.getFileName()+"_sample", idx_=selected_reads[:n_sample], outputFormat=outputFormat)
		else:
			self.writeWithIdx(filename=filename, idx_=selected_reads[:n_sample], outputFormat=outputFormat)