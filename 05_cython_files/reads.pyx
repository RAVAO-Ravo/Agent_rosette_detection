#-*- coding:utf-8 -*-
#cython: language_level=3


############################### IMPORTATIONS ####################################

import numpy as np
import pandas as pd
import pickle as pk
import random as rd
import typing as tp
from my_ahocorapy import ahocorasick
from sklearn.ensemble import IsolationForest
from xgboost import XGBClassifier

############################### IMPORTATIONS ####################################


################################ AHOCORASICK ####################################

cdef int none = (-99)


cdef str getComplementarySeq(str sequence) :

	""" Retourne le complémentaire d'un séquence nucléotidique. """

	# Initialiser les variables
	cdef str complementary = ''
	cdef str nuc = ''

	# Pour chaque base de la séquence nucléotique (itérée dans le sens inverse)
	for nuc in sequence[len(sequence)::-1] :

		# Récupérer le complémentaire de la base
		if nuc == 'T' :
			complementary+='A'
		elif nuc == 't' :
			complementary+='a'
		elif nuc == 'A' :
			complementary+='T'
		elif nuc == 'a' :
			complementary+='t'
		elif nuc == 'G' :
			complementary+='C'
		elif nuc == 'g' :
			complementary+='c'
		elif nuc == 'C' :
			complementary+='G'
		else :
			complementary+='g'

	# Retourner le complémentaire
	return complementary


cdef dict search(str text, object searchMachine, dict dict_results) :

	""" Effectue la recherche de motifs, à partir d'un automate de recherche ahocorasick. """

	# Initialiser les variables
	cdef list results = searchMachine.search_all(text)
	cdef str motif = ''
	cdef int position = 0

	# Pour les chaque résultat, ajouter le au dictionnaire
	for motif, position in results :
		dict_results[motif].append(position)

	# Retourner le dictionnaire
	return dict_results


cdef int getPosPrimer(str sequence, object searchMachine, list primers, bint forward=True) :

	""" Récupère le positionnement d'un primer. """

	# Initialiser les variables
	cdef str primer = ''
	cdef dict results = {}

	# Récupérer les positions des mots contenus dans la machine
	results = search(text=sequence,
					 searchMachine=searchMachine, 
					 dict_results={primer.lower() : [] for primer in primers})

	# Récupérer les primers présents
	results = {primer : positions for (primer, positions) in results.items() if len(positions) != 0}

	# Récupérer les résultats
	if len(results) == 0 :
		return none
	else:
		if forward == True :
			return min([position for positions in results.values() for position in positions])
		else :
			return max([position+len(primer) for (primer, positions) in results.items() for position in positions])

################################ AHOCORASICK ####################################


################################ CLASSIFIEUR ####################################

cdef dict cherif = {"forwards" : ["GCGGTAATTCCAGCTCCA"],
	   	  			"reverses" : ["CACTCAATTAAGCGCACACG"]}


cdef int encode(str nuc) :
	
	""" Récupère l'encodge d'un nucléotide. """

	# Retourner la valeur
	return 1 if nuc.lower() == 'a' else \
 		   2 if nuc.lower() == 't' else \
 		   3 if nuc.lower() == 'c' else \
 		   4 if nuc.lower() == 'g' else \
 		   0


cdef object transformer(str sequence, int n_features) :

	""" Transforme une séquence nucléotidique en vecteur. """

	# Retourner le vecteur
	return np.array([[encode(nuc=letter) for letter in sequence+'n'*(n_features-len(sequence))]])

################################ CLASSIFIEUR ####################################


################################ STRUCTURE READ #################################

cdef dict columns = {key : value for value, key in enumerate(["id",
					 										  "sequence",
					 										  "quality",
					 										  "prefix",
					 										  "suffix",
					 										  "length",
					 										  'A',
					 										  'T',
					 										  'C',
					 										  'G',
					 										  'N',
					 										  "GC",
					 										  "AT",
					 										  "Q_Score",
					 										  "std_QScore",
					 										  "q25",
					 										  "q50",
					 										  "q75",
					 										  "q90",
					 										  "ScoreDist"])}


cdef tuple computeQuantValues(list q_scoreRead, int len_seq) :

	""" Retourne un tuple de valeurs de plusieurs quantiles, pour un ensemble de données triées. """

	# Récupérer les valeurs
	cdef float q25 = q_scoreRead[int(0.25*len_seq)]
	cdef float q50 = q_scoreRead[int(0.50*len_seq)]
	cdef float q75 = q_scoreRead[int(0.75*len_seq)]
	cdef float q90 = q_scoreRead[int(0.90*len_seq)]

	# Retourner le tuple des valeurs
	return (q25, q50, q75, q90)


cdef float dist(list p1, list p2) :

	""" Calcule la distance entre deux points. """
	
	# Initialiser les variables
	cdef float pc1 = 0.0
	cdef float pc2 = 0.0
	cdef float square_sum = 0.0
	
	# Calculer la somme des différences au carrés de chaque composante dimensionnelle
	for pc1, pc2 in zip(p1, p2) :
		square_sum+=(pc1-pc2)**2

	# Retourner la distance
	return square_sum**0.5


cdef float computeDistScore(list values) :

	""" Calcule le 'score de distance qualité'. """

	# Initialiser les variables
	cdef float value = 0.0

	# Retourner le score
	return dist(p1=[0]*len(values), p2=[1/value for value in values])


cpdef tuple createRead(tuple read) :

	""" Créer une sorte d'objet read. """

	# Initialiser les variables
	cdef str id_ = read[0]
	cdef str sequence = read[1]
	cdef str quality = read[2]
	cdef int len_seq = len(sequence)
	cdef list q_scoreRead = [0]*len_seq
	cdef int i = 0
	cdef float a = 0.0
	cdef float t = 0.0
	cdef float c = 0.0
	cdef float g = 0.0
	cdef float n = 0.0
	cdef float gc_content = 0.0
	cdef float at_content = 0.0
	cdef float sum_ = 0.0
	cdef float sum_square = 0.0
	cdef float moy = 0.0
	cdef float std_ = 0.0
	cdef float q25 = 0.0
	cdef float q50 = 0.0
	cdef float q75 = 0.0
	cdef float q90 = 0.0
	cdef float score = 0.0

	# Pour chaque base et lettre de score qualité
	for i in range(0, len_seq, 1) :

		# Compter chaque base du read
		if sequence[i] == 'A' :
			a+=1
		elif sequence[i] == 'T' :
			t+=1
		elif sequence[i] == 'C' :
			c+=1
		else :
			g+=1

		# Récupérer le score qualité de chaque base, du read, dans une liste
		q_scoreRead[i] = ord(quality[i])-33
		sum_+=q_scoreRead[i]
		sum_square+=q_scoreRead[i]**2
		n += 1 if q_scoreRead[i] < 15 else 0

	# Calculer les pourcentages de chaque base
	a/=len_seq
	t/=len_seq
	c/=len_seq
	g/=len_seq
	n/=len_seq

	# Calculer les proportions de GC et AT
	gc_content = (c+g)
	at_content = 1-gc_content

	# Calculer la moyenne et l'écart-type
	moy = sum_/len_seq
	std_ = ((sum_square/len_seq)-moy**2)**0.5

	# Calculer les valeurs de quantiles, et le score de distance qualité
	q_scoreRead.sort()
	q25, q50, q75, q90 = computeQuantValues(q_scoreRead=q_scoreRead, len_seq=len_seq)
	score = computeDistScore([moy, q25, q50, q75, q90])

	# Retourner la 'structure' read
	return (id_,
			sequence,
			quality,
			sequence[:50],
			sequence[-50:],
			len_seq,
			a,
			t,
			c,
			g,
			n,
			gc_content,
			at_content,
			moy,
			std_,
			q25,
			q50,
			q75,
			q90,
			score)


cdef list makeKmers(str sequence, int k=21) :

	""" Crée les k-mers possibles pour une séquence. """

	# Initialiser les variables
	cdef int len_seq = len(sequence)
	cdef int i = 0

	# Retourner la liste des k-mers
	return [sequence[i:i+k] for i in range(0, len_seq-k+1, 1)]


cdef tuple headcroping(str sequence, str quality, float value) :

	""" Ampute plusieurs bases au début de la séquence. """

	# Initialiser les variables
	cdef int n_bases = len(sequence)
	cdef int finalValue = 0

	# Si la valeur passer est entre 0 et 1, récupérer le nombre de bases associé à cette valeur
	if 0.0 <= value and value < 1.0 :
		value = n_bases*value

	# Valeur final pour l'amputation
	finalValue = int(round(value))

	# Retourner la séquence nouvelle séquence
	return (sequence[finalValue:], quality[finalValue:])


cdef tuple tailcroping(str sequence, str quality, float value) :

	""" Amputer plusieurs bases à la fin de la séquence. """

	# Initialiser les variables
	cdef int n_bases = len(sequence)
	cdef int finalValue = 0

	# Si la valeur passer est entre 0 et 1, récupérer le nombre de bases associé à cette valeur
	if 0.0 <= value and value < 1.0 :
		value = n_bases*value

	# Valeur final pour l'amputation
	finalValue = int(round(value))

	# Retourner la séquence nouvelle séquence
	return (sequence[:n_bases-finalValue], quality[:n_bases-finalValue])

################################ STRUCTURE READ #################################


################################## CLASS READS ##################################


def readFastq(filename: str) -> tp.Generator :

	""" Récupère les reads d'un fichier fastq. """

	# Initialiser les variables
	cdef str id_ = ''
	cdef str seq = ''
	cdef str qual = ''
	fastq: tp.TextIO[str] = None

	# Ouvrir le fastq
	with open(file=filename, mode='r') as fastq :

		# Tant qu'il y a des lignes
		while True :

			# Récupérer la séquence, et les infos associées
			id_ = fastq.readline().rstrip()
			seq = fastq.readline().rstrip()
			fastq.readline()
			qual = fastq.readline().rstrip()

			# Stopper la lecture quand il n'y a plus de lignes
			if len(seq) == 0 :
				break

			yield createRead((id_, seq, qual))


cdef class Reads:

	""" Classe Reads modélisant un ensemble de reads, et permet de manipuler l'ensemble. """

	### Attributs de la classe ###

	# Nom du fichier
	cdef str filename
	
	# Ensemble des reads
	cdef list reads
	
	# Nombre de reads
	cdef int n_reads


	### Constructeur de la classe ###

	def __init__(self, str filename) :

		""" Construit un objet Reads, à partir d'un fichier de reads. """

		self.filename = filename
		self.reads = list(readFastq(filename))
		self.n_reads = len(self.reads)

	### Getters de la classe ###

	cpdef str getFileName(self) :
		return self.filename

	cpdef list getReads(self) :
		return self.reads

	cpdef int getN_reads(self) :
		return self.n_reads


	### Setters et updaters de la classe ###

	cpdef void setReads(self, list reads) :
		self.reads = reads

	cpdef void setN_reads(self, int n_reads) :
		self.n_reads = n_reads


	### Getters spéciaux de la classe ###

	cpdef list getIDs(self) :
		cdef tuple read = ()
		return [read[columns["id"]] for read in self.reads]

	cpdef list getSequences(self) :
		cdef tuple read = ()
		return [read[columns["sequence"]] for read in self.reads]

	cpdef list getQualities(self) :
		cdef tuple read = ()
		return [read[columns["Q_Score"]] for read in self.reads]

	cpdef list getPrefixes(self) :
		cdef tuple read = ()
		return [read[columns["prefix"]] for read in self.reads]

	cpdef list getSuffixes(self) :
		cdef tuple read = ()
		return [read[columns["suffix"]] for read in self.reads]


	### Méthodes de la classe ###
	
	def to_numpy(self) -> np.ndarray :

		""" Retourne les stats des reads sous forme d'un tableau numpy. """

		cdef tuple read = ()
		return np.array([read for read in self.reads])


	def to_DataFrame(self) -> pd.DataFrame :

		""" Retourne les stats des reads sous forme d'un dataframe. """

		cdef tuple read = ()
		return pd.DataFrame(data=[read for read in self.reads], columns=columns.keys())


	cpdef Reads filter_Duplicated(self) :

		""" Retire les éléments dupliqués de l'ensemble des reads. """

		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []
		cdef str prefix = ''
		cdef set unique_ = set()
		
		# Pour chaque read
		for read in self.reads :

			# Récupérer son préfixe
			prefix = read[columns["prefix"]]

			# Si le préfixe n'est pas dans la liste unique_
			if prefix not in unique_ :

				# Ajouter le préfixe à unique_
				unique_.add(prefix)

				# Ajouter le read à nouvelle liste des reads
				new_reads.append(read)
			
			# Sinon, ne pas y ajouter le read 
			else :
				self.n_reads-=1

		# Mettre à jour les reads
		self.reads = new_reads
		
		# Retourner l'objet
		return self


	cpdef Reads filter_size(self, float inf=0.0, float sup=1.0) :

		""" Filtre les reads selon leur taille. """

		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []

		# Récupérer les valeurs des quantiles
		if (0.0 <= inf and inf < 1.0) or (0.0 < sup and sup <= 1.0) :
			lengths: np.ndarray = np.array([read[columns["length"]] for read in self.reads], dtype=float)
			if 0.0 <= inf and inf < 1.0 :
				inf = np.quantile(a=lengths, q=inf)
			if 0.0 < sup and sup <= 1.0 :
				sup = np.quantile(a=lengths, q=sup)

		# Récupérer les reads avec une taille comprise entre les bornes
		for read in self.reads :
			if inf <= read[columns["length"]] and read[columns["length"]] <= sup :
				new_reads.append(read)
			else :
				self.n_reads-=1
	
		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self


	cpdef Reads filter_Q(self, float seuil=20.0) :

		""" Filtre les reads par leur Q score. """

		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []

		# Récupérer les reads dont le score est supérieur à la borne minorante
		for read in self.reads :
			if seuil <= read[columns["Q_Score"]] :
				new_reads.append(read)
			else :
				self.n_reads-=1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self


	cpdef Reads filter_quality(self, float seuil=1.0) :

		""" Filtre les reads par leur qualité. """

		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []

		# Récupérer la valeur du quantile
		scores: np.ndarray = np.array([read[columns["ScoreDist"]] for read in self.reads], dtype=float)
		seuil = np.quantile(a=scores, q=seuil)

		# Récupérer les reads dont le score est inférieur à la borne majorante
		for read in self.reads :
			if read[columns["ScoreDist"]] <= seuil :
				new_reads.append(read)
			else :
				self.n_reads-=1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self


	cpdef Reads filter_N(self, float seuil=1.0) :

		""" Filtre les reads selon leur quantité de bases indécises (Q_Score(bases) < 15). """

		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []

		# Récupérer la valeur du quantile
		N_content: np.ndarray = np.array([read[columns['N']] for read in self.reads], dtype=float)
		seuil = np.quantile(a=N_content, q=seuil)

		# Récupérer les reads ayant les plus faibles quantités de bases indécises
		for read in self.reads :
			if read[columns['N']] <= seuil :
				new_reads.append(read)
			else :
				self.n_reads-=1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self


	cpdef Reads filter_GC(self, float minGC=0.0, float maxGC=1.0) :

		""" Filtre les reads selon leur quantité de GC. """

		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []

		# Récupérer les reads ayant un taux de GC contenu dans les bornes
		for read in self.reads :
			if minGC <= read[columns["GC"]] and read[columns["GC"]] <= maxGC :
				new_reads.append(read)
			else :
				self.n_reads-=1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self


	cpdef Reads HTC(self, float head=0.0, float tail=0.0) :

		""" Ampute les extrémités de chaque read, d'un certain nombre de bases. """

		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []
		cdef str id_ = ''
		cdef str sequence = ''
		cdef str quality = ''
		cdef float len_seq = 0.0

		# Amputer chaque read, selon les paramètres définis
		for read in self.reads :
			
			# Récupérer les éléments nécéssaire du read
			id_ = read[columns["id"]]
			sequence = read[columns["sequence"]]
			quality = read[columns["quality"]]
			len_seq = float(read[columns["length"]])

			# Exciser au début du read
			if 0.0 < head and head < 1.0 :
				sequence, quality = headcroping(sequence=sequence, quality=quality, value=(len_seq*head))
			elif 1.0 <= head :
				sequence, quality = headcroping(sequence=sequence, quality=quality, value=head)

			# Exciser à la fin du read
			if 0.0 < tail and tail < 1.0 :
				sequence, quality = tailcroping(sequence=sequence, quality=quality, value=(len_seq*tail))
			elif 1.0 <= tail :
				sequence, quality = tailcroping(sequence=sequence, quality=quality, value=tail)

			# Créer le nouveau read
			if len(sequence) != 0 :
				new_reads.append(createRead((id_, sequence, quality)))
			else:
				self.n_reads-=1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self


	cpdef Reads cutPrimers(self, list forwards, list reverses, int seuil_length=1) :

		""" Filtre les reads selon des primers. """

		# Initialiser les variables
		cdef tuple read = ()
		cdef list new_reads = []
		cdef str id_ = ''
		cdef str sequence = ''
		cdef str quality = ''
		
		# Créer les automates de recherche
		reverses = [getComplementarySeq(sequence=sequence) for sequence in reverses]
		forwardMachine: ahocorasick = ahocorasick(words=forwards)
		reverseMachine: ahocorasick = ahocorasick(words=reverses)

		# Récupérer les reads ayant les primers
		for read in self.reads :

			# Récupérer la séquence du read
			id_ = read[columns["id"]]
			sequence =  read[columns["sequence"]]
			quality = read[columns["quality"]]

			# Récupérer la postion du forward et du reverse
			min_idx = getPosPrimer(sequence=sequence, searchMachine=forwardMachine, primers=forwards, forward=True)
			max_idx = getPosPrimer(sequence=sequence, searchMachine=reverseMachine, primers=reverses, forward=False)

			# Si les deux primers sont présents, dans le bon sens
			if (min_idx != none) and (max_idx != none) and (min_idx < max_idx) :
				
				# Couper au niveau des primers
				sequence = sequence[min_idx:max_idx]
				quality = quality[min_idx:max_idx]

				# Si la taille est suffisante, récupérer le read
				if seuil_length <= len(sequence) :
					new_reads.append(createRead((id_, sequence, quality)))

				else :
					self.n_reads-=1
			else :
				self.n_reads-=1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner l'objet
		return self


	cpdef tuple clf_RA(self, float seuil=0.90) :

		""" Filtre les reads à partir du classifieur 'clf_RA'. """

		# Initialiser les variables
		cdef list reverses = [getComplementarySeq(sequence=sequence) for sequence in cherif["reverses"]]
		cdef list forwards = cherif["forwards"]
		cdef tuple read = ()
		cdef list new_reads = []
		cdef str id_ = ''
		cdef str sequence = ''
		cdef str quality = ''
		cdef int initial = 0
		cdef int primer_cutting = 0
		cdef int n_features = 160
		cdef int size_filter = 0
		pred: np.ndarray = None

		# Créer les automates de recherche
		forwardMachine: ahocorasick = ahocorasick(words=forwards)
		reverseMachine: ahocorasick = ahocorasick(words=reverses)

		# Récupérer le classifieur
		with open(file="./clf_RA.sav", mode="rb") as model_clf :
			clf: XGBClassifier = pk.load(file=model_clf)

		# Récupérer le nombre de reads initial
		initial = self.n_reads

		# Classifier les reads
		for read in self.reads :

			# Récupérer la séquence
			id_ = read[columns["id"]]
			sequence = read[columns["sequence"]]
			quality = read[columns["quality"]]

			# Récupérer la postion du forward et du reverse
			min_idx = getPosPrimer(sequence=sequence, searchMachine=forwardMachine, primers=forwards, forward=True)
			max_idx = getPosPrimer(sequence=sequence, searchMachine=reverseMachine, primers=reverses, forward=False)

			# Si les deux primers sont présents, dans le bon sens
			if (min_idx != none) and (max_idx != none) and (min_idx < max_idx) :
				
				# Couper au niveau des primers
				sequence = sequence[min_idx:max_idx]
				quality = quality[min_idx:max_idx]

				# Incrémenter 'primer_cutting'
				primer_cutting+=1

				# Si la taille est correcte
				if n_features-10 <= len(sequence) <= n_features :

					# Prédire la probalité d'appartenance à l'agent rosette
					pred = clf.predict_proba(transformer(sequence=sequence, n_features=n_features))

					# Incrémenter size_filtrer
					size_filter+=1

					# Si la probabilité est supérieur au seuil, récupérer la séquence
					if seuil <= pred[0][9] :
						new_reads.append(createRead((id_, sequence, quality)))
					
					else :
						self.n_reads-=1
				else :
					self.n_reads-=1
			else :
				self.n_reads-=1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Retourner les résultats
		return (initial, primer_cutting, size_filter, self.n_reads)


	cpdef Reads multiFilter(self,
							float size_inf=0.0,
							float size_sup=1.0,
							float seuil_qual=1.0,
							float seuil_N=1.0,
							float head=0.0,
							float tail=0.0,
							float minGC=0.0,
							float maxGC=1.0,
							bint duplicated=False) :

		""" Filtre les reads selon plusieurs paramètres. """
		
		# Initialiser les variables
		df: pd.DataFrame  = self.to_DataFrame().loc[:, ("length", "ScoreDist", "N", "GC")]
		cdef list new_reads = []
		cdef tuple read = ()
		cdef float len_seq = 0.0

		# Récupérer les valeurs de quantiles pour la taille des reads
		if (0.0 <= size_inf and size_inf < 1.0) or (0.0 < size_sup and size_sup <= 1.0) :
			if 0.0 <= size_inf and size_inf < 1.0 :
				size_inf = df["length"].quantile(q=size_inf)
			if 0.0 < size_sup and size_sup <= 1.0 :
				size_sup = df["length"].quantile(q=size_sup)

		# Récupérer les valeurs de quantiles pour la qualité
		seuil_qual = df["ScoreDist"].quantile(q=seuil_qual)
		
		# Récupérer les valeurs de quantiles pour le pourcentage de bases indécises
		seuil_N = df['N'].quantile(q=seuil_N)

		# Récupérer les reads selon les filtres
		for read in self.reads :

			if (size_inf <= read[columns["length"]]) and (read[columns["length"]] <= size_sup) and \
			   (read[columns["ScoreDist"]] <= seuil_qual) and (read[columns['N']] <= seuil_N) and \
			   (minGC <= read[columns["GC"]]) and (read[columns["GC"]] <= maxGC) :
				
				id_ = read[columns["id"]]
				sequence = read[columns["sequence"]]
				quality = read[columns["quality"]]
				len_seq = read[columns["length"]]

				if 0.0 < head and head < 1.0 :
					sequence, quality = headcroping(sequence=sequence, quality=quality, value=(len_seq*head))
				elif 1.0 <= head :
					sequence, quality = headcroping(sequence=sequence, quality=quality, value=head)

				if 0.0 < tail and tail < 1.0 :
					sequence, quality = tailcroping(sequence=sequence, quality=quality, value=(len_seq*tail))
				elif 1.0 <= tail :
					sequence, quality = tailcroping(sequence=sequence, quality=quality, value=tail)

				if len(sequence) != 0 :
					new_reads.append(createRead((id_, sequence, quality)))

				else:
					self.n_reads-=1
			else :
				self.n_reads-=1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Filtrer les duplications
		if duplicated == True :
			self.filter_Duplicated()

		# Retourner l'objet
		return self


	cpdef Reads multiFilterIF(self, 
							  object seuil="auto",
							  tuple dims=("length", "ScoreDist", "N", "GC"),
							  float head=0.0,
							  float tail=0.0,
							  bint duplicated=False,
							  object random_state=None) :

		""" Filtre les reads selon l'algorithme de l'isolation forest. """

		# Initialiser les variables
		cdef str dim = ''
		cdef tuple read = ()
		cdef list outliers = [[read[columns[dim]] for dim in dims] for read in self.reads]
		cdef list new_reads = []
		cdef int outlier = 0
		cdef str id_ = ''
		cdef str sequence = ''
		cdef str quality = ''
		cdef float len_seq = 0.0

		# Récupérer les outliers
		if isinstance(seuil, str) :
			outliers = list(IsolationForest(random_state=random_state, contamination=seuil).fit_predict(X=outliers))
		else :
			outliers = list(IsolationForest(random_state=random_state, contamination=(1-seuil)).fit_predict(X=outliers))

		# Récupérer les reads non 'outliers'
		for read, outlier in zip(self.reads, outliers) :

			if outlier == 1 :

				id_ = read[columns["id"]]
				sequence = read[columns["sequence"]]
				quality = read[columns["quality"]]
				len_seq = read[columns["length"]]

				if 0.0 < head and head < 1.0 :
					sequence, quality = headcroping(sequence=sequence, quality=quality, value=(len_seq*head))
				elif 1.0 <= head :
					sequence, quality = headcroping(sequence=sequence, quality=quality, value=head)

				if 0.0 < tail and tail < 1.0 :
					sequence, quality = tailcroping(sequence=sequence, quality=quality, value=(len_seq*tail))
				elif 1.0 <= tail :
					sequence, quality = tailcroping(sequence=sequence, quality=quality, value=tail)

				if len(sequence) != 0 :
					new_reads.append(createRead((id_, sequence, quality)))

				else:
					self.n_reads-=1
			else :
				self.n_reads-=1

		# Mettre à jour les attributs
		self.reads = new_reads

		# Filtrer les duplications
		if duplicated == True :
			self.filter_Duplicated()

		# Retourner l'objet
		return self


	cpdef list readsToKmers(self, int k=21) :

		""" Retourne l'ensemble des k-mers uniques possibles pour les reads. """

		# Initialiser les variables
		cdef tuple read = ()
		cdef set kmers = set()
		cdef str kmer = ''

		# Retourner la liste triée des kmers
		return sorted(set([kmer for read in self.reads for kmer in makeKmers(sequence=read[columns["sequence"]], k=k)]))


	cpdef void write(self, str filename, str outputFormat="fq") :

		""" Enregistre les reads, selon un format biologique défini. """

		# Initialiser les variables
		cdef tuple read = ()

		# Créer le fichier de sortie
		with open(file=filename, mode='w') as outputFile :

			# Si l'on souhaite une sortie fastq
			if outputFormat == "fq" :
				# Enregistrer selon le format fastq
				for read in self.reads :
					outputFile.write(f"""{read[columns["id"]]}\n{read[columns["sequence"]]}\n+\n{read[columns["quality"]]}\n""")

			# Si l'on souhaite une sortie fastq
			elif outputFormat == "fa" :
				# Enregistrer selon le format fasta
				for read in self.reads :
					outputFile.write(f""">{read[columns["id"]][1:]}\n{read[columns["sequence"]]}\n""")

			# Sinon
			else :
				# Ne pas enregistrer le fichier
				print("Format inexistant, ou pas pris en compte. L'enregistrement est abordé.")


	cpdef void writeWithIdx(self, str filename, list idx_, str outputFormat="fq") :

		""" Enregistre les reads, selon un format défini, à partir d'identifiants. """

		# Initialiser les variables
		cdef tuple read = ()
		cdef str id_ = ''

		# Créer le fichier de sortie
		with open(file=filename, mode='w') as outputFile :

			# Enregistrer les reads sélectionnés, selon le format fastq
			if outputFormat == "fq" :
				for read in self.reads :
					id_ = read[columns["id"]]
					if id_ in idx_ :
						outputFile.write(f"""{id_}\n{read[columns["sequence"]]}\n+\n{read[columns["quality"]]}\n""")

			# Enregistrer les reads sélectionnés, selon le format fasta
			elif outputFormat == "fa" :
				for read in self.reads :
					id_ = read[columns["id"]]
					if id_ in idx_ :
						outputFile.write(f""">{id_[1:]}\n{read[columns["sequence"]]}\n""")

			# Sinon
			else :
				# Ne pas enregistrer le fichier
				print("Format inexistant, ou pas pris en compte. L'enregistrement est abordé.")


	cpdef void makeSample(self, int n_sample, str filename='', str outputFormat="fq", object random_state=None) :

		""" Sous-échantillonne des reads depuis l'ensemble des reads. """

		# Définir un graine aléatoire
		rd.seed(random_state)

		# Sélectionner les reads
		cdef list selected_reads = self.getIDs()
		rd.shuffle(selected_reads)

		# Enregistrer le sous-échantillon
		if filename == '' :
			self.writeWithIdx(filename=self.getFileName()+"_sample", idx_=selected_reads[:n_sample], outputFormat=outputFormat)
		else :
			self.writeWithIdx(filename=filename, idx_=selected_reads[:n_sample], outputFormat=outputFormat)

################################## CLASS READS ##################################