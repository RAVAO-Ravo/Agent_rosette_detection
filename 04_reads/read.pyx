#-*- coding:utf-8 -*-
#cython:language_level=3

# Importation des modules
from typing import Generator, TextIO

cdef tuple computeQuantValues(list q_scoreRead, int len_seq):
	"""
	Retourne un tuple de valeurs de plusieurs quantiles, pour un ensemble de données triées.

	Args:
		q_scoreRead (List[int]): Liste triée des scores de qualité.
		len_seq (int): La longueur de la séquence.

	Returns:
		Tuple[float, float, float, float]: Un tuple de valeurs de quantiles (q25, q50, q75, q90).
	"""
	# Récupérer les valeurs
	cdef float q25 = q_scoreRead[int(0.25*len_seq)]
	cdef float q50 = q_scoreRead[int(0.50*len_seq)]
	cdef float q75 = q_scoreRead[int(0.75*len_seq)]
	cdef float q90 = q_scoreRead[int(0.90*len_seq)]

	# Retourner le tuple des valeurs
	return (q25, q50, q75, q90)

cdef float dist(list p1, list p2):
	"""
	Calcule la distance entre deux points.

	Args:
		p1 (List[float]): Les coordonnées du premier point.
		p2 (List[float]): Les coordonnées du deuxième point.

	Returns:
		float: La distance entre les deux points.
	"""
	# Initialiser les variables
	cdef float pc1 = 0.0
	cdef float pc2 = 0.0
	cdef float square_sum = 0.0
	
	# Calculer la somme des différences au carrés de chaque composante dimensionnelle
	for pc1, pc2 in zip(p1, p2):
		square_sum += (pc1 - pc2) ** 2

	# Retourner la distance
	return square_sum ** 0.5

cdef float computeDistScore(list values):
	"""
	Calcule le 'score de distance qualité'.

	Args:
		values (List[float]): Une liste de valeurs à partir desquelles calculer le score.

	Returns:
		float: Le score de distance qualité.
	"""
	# Retourner le score
	return dist(p1=[0]*len(values), p2=[1/value for value in values])


cpdef tuple createRead(tuple read):
	"""
	Crée une structure de données représentant un read.

	Args:
		read (Tuple[str, str, str]): Un tuple contenant l'identifiant, la séquence et la qualité du read.

	Returns:
		Tuple[Any]: Une structure de données représentant le read avec diverses caractéristiques calculées.
	"""
	# Déclarer les types pour les variables
	cdef str id_, sequence, quality
	cdef int len_seq
	cdef list q_scoreRead

	# Initialiser les variables
	id_, sequence, quality = read
	len_seq = len(sequence)
	q_scoreRead = [0] * len_seq

	# Déclarer les types pour les variables utilisées dans les boucles
	cdef int i
	cdef float a = 0.0
	cdef float t = 0.0
	cdef float c = 0.0
	cdef float g = 0.0
	cdef float n = 0.0
	cdef float gc_content, at_content
	cdef float sum_ = 0.0
	cdef float sum_square = 0.0
	cdef float moy, std_
	cdef float q25, q50, q75, q90
	cdef float score

	# Pour chaque base et lettre de score qualité
	for i in range(0, len_seq, 1):
		# Compter chaque base du read
		if sequence[i] == 'A':
			a += 1
		elif sequence[i] == 'T':
			t += 1
		elif sequence[i] == 'C':
			c += 1
		else:
			g += 1

		# Récupérer le score qualité de chaque base, du read, dans une liste
		q_scoreRead[i] = ord(quality[i]) - 33
		sum_ += q_scoreRead[i]
		sum_square += q_scoreRead[i] ** 2
		n += 1 if q_scoreRead[i] < 15 else 0

	# Calculer les pourcentages de chaque base
	a /= len_seq
	t /= len_seq
	c /= len_seq
	g /= len_seq
	n /= len_seq

	# Calculer les proportions de GC et AT
	gc_content = (c + g)
	at_content = 1 - gc_content

	# Calculer la moyenne et l'écart-type
	moy = sum_ / len_seq
	std_ = ((sum_square / len_seq) - moy ** 2) ** 0.5

	# Calculer les valeurs de quantiles, et le score de distance qualité
	q_scoreRead.sort()
	q25, q50, q75, q90 = computeQuantValues(q_scoreRead=q_scoreRead, len_seq=len_seq)
	score = computeDistScore([moy, q25, q50, q75, q90])

	# Retourner la structure de données du read
	return (
		id_, sequence, quality, sequence[:50], sequence[-50:], len_seq, a, t, c, g, n,
		gc_content, at_content, moy, std_, q25, q50, q75, q90, score
	)

def readFastq(filename: str) -> Generator:
	"""
	Récupère les reads depuis un fichier fastq.

	Args:
		filename (str): Le chemin du fichier Fastq.

	Yields:
		Tuple[Any]: Un générateur de récupération des reads.
	"""
	# Initialiser les variables
	cdef str id_ = ''
	cdef str seq = ''
	cdef str qual = ''
	fastq: TextIO[str] = None

	# Ouvrir le fastq
	with open(file=filename, mode='r') as fastq:
		# Tant qu'il y a des lignes
		while True:
			# Récupérer la séquence, et les infos associées
			id_ = fastq.readline().rstrip()
			seq = fastq.readline().rstrip()
			fastq.readline()
			qual = fastq.readline().rstrip()

			# Stopper la lecture quand il n'y a plus de lignes
			if len(seq) == 0:
				break

			yield createRead((id_, seq, qual))

cpdef tuple headcroping(str sequence, str quality, float value):
	"""
	Ampute plusieurs bases au début de la séquence.

	Args:
		sequence (str): La séquence à modifier.
		quality (str): La qualité associée à la séquence.
		value (float): La valeur pour l'amputation.

	Returns:
		Tuple[str, str]: La séquence et la qualité modifiées.
	"""
	# Initialiser les variables
	cdef int n_bases = len(sequence)
	cdef int finalValue = 0

	# Si la valeur passer est entre 0 et 1, récupérer le nombre de bases associé à cette valeur
	if 0.0 <= value and value < 1.0:
		value = n_bases * value

	# Valeur final pour l'amputation
	finalValue = int(round(value))

	# Retourner la séquence nouvelle séquence
	return (sequence[finalValue:], quality[finalValue:])

cpdef tuple tailcroping(str sequence, str quality, float value):
	"""
	Ampute plusieurs bases à la fin de la séquence.

	Args:
		sequence (str): La séquence à modifier.
		quality (str): La qualité associée à la séquence.
		value (float): La valeur pour l'amputation.

	Returns:
		Tuple[str, str]: La séquence et la qualité modifiées.
	"""
	# Initialiser les variables
	cdef int n_bases = len(sequence)
	cdef int finalValue = 0

	# Si la valeur passer est entre 0 et 1, récupérer le nombre de bases associé à cette valeur
	if 0.0 <= value and value < 1.0:
		value = n_bases * value

	# Valeur final pour l'amputation
	finalValue = int(round(value))

	# Retourner la séquence nouvelle séquence
	return (sequence[:n_bases-finalValue], quality[:n_bases-finalValue])

cpdef list makeKmers(str sequence, int k = 21):
	"""
	Crée les k-mers possibles pour une séquence.

	Args:
		sequence (str): La séquence à partir de laquelle générer les k-mers.
		k (int): La longueur des k-mers. (Par défaut = 21).

	Returns:
		List[str]: Une liste de k-mers générés à partir de la séquence.
	"""
	# Initialiser les variables
	cdef int len_seq = len(sequence)

	# Retourner la liste des k-mers
	return [sequence[i:i+k] for i in range(0, len_seq-k+1, 1)]