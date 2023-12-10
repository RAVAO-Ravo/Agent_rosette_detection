#-*- coding:utf-8 -*-
#cython:language_level=3

# Importation des modules
from collections import deque
from config import GlobalValue

# Initialisation d'un objet "GlobalValue" et récupération de la valeur du none
GVALUES: GlobalValue = GlobalValue()
cdef int NONE = GVALUES.get_values("none")

cdef class AhoCorasick(object):
	"""
	Classe AhoCorasick permettant de créer un automate de recherche de motifs, utilisant l'algorithme 'Aho-Corasick'.

	Attributes:
		automaton (list): La structure de données représentant l'automate de recherche.

	Methods:
		next_node(): Récupère l'index du nœud suivant.
		add(): Ajoute un mot dans l'automate.
		adds(): Ajoute un liste de mots dans l'automate
		set_fail_transitions(): Ajoute les liens echecs dans l'automate
		search_all(): Permet la recherche des mots, à partir d'un texte.
	"""
	# L'automate de recherche, sous forme d'une liste de dictionnaires, chaque dictionnaire étant un nœud de Trie
	cdef list automaton

	def __init__(self, list words):
		"""
		Construit une instance de l'automate Aho-Corasick pour la recherche de motifs.

		Args:
			words (List[str]): Une liste de mots (chaînes de caractères) à rechercher.

		Returns:
			AhoCorasick: Une instance de l'automate Aho-Corasick prête à être utilisée.
		"""
		# Initialiser la racine du Trie
		self.automaton = [{"value": '', "next_nodes": [], "fail_state": 0, "outputs": []}]

		# Ajouter les mots au Trie
		self.adds(words=words)
		
		# Créer les liens echecs
		self.set_fail_transitions()

	cdef int next_node(self, int current_node, str value):
		"""
		Récupère l'index du nœud suivant contenant la valeur demandée ('value') à partir d'un point de départ du Trie ('current_node').

		Args:
			current_node (int): L'index du nœud de départ dans le Trie.
			value (str): La valeur recherchée dans les nœuds fils.

		Returns:
			int: L'index du nœud suivant contenant la valeur demandée ('value') s'il est trouvé, ou 'NONE' si la valeur n'est pas trouvée.
		"""
		# Initialiser la variable 'id_node'
		cdef int id_node = 0
		
		# Pour chaque nœud fils du nœud courrant
		for id_node in self.automaton[current_node]["next_nodes"]:
			# Si un nœud fils contient la valeur demandée, retourner le
			if self.automaton[id_node]["value"] == value:
				return id_node
		
		# Retourner 'NONE' si la valeur n'a pas été trouvée
		return NONE

	cdef void add(self, str word):
		"""
		Ajoute un mot à l'automate de recherche.

		Args:
			word (str): Le mot à ajouter à l'automate.

		Returns:
			None
		"""
		# Initialiser les variables
		word = word.lower()
		cdef int current_node = 0
		cdef int i = 0
		cdef int pos = 0
		cdef int len_word = len(word)
		cdef int next_node = self.next_node(current_node, word[pos])

		# Tant qu'une partie du mot (préfixe) est présent dans le Trie
		while next_node != NONE:
			# Passer au nœud, et à la lettre, suivant(e)
			current_node = next_node
			pos+=1

			# Sortir de la boucle si on a parcouru toute les lettres du mot
			if pos < len_word:
				next_node = self.next_node(current_node, word[pos])
			else:
				break

		# Pour chaque lettre, du mot, non présent dans le Trie
		for i in range(pos, len_word, 1):
			# Ajouter la lettre
			self.automaton.append({"value": word[i], "next_nodes": [], "fail_state": 0, "outputs": []})
			self.automaton[current_node]["next_nodes"].append(len(self.automaton)-1)
			current_node = len(self.automaton)-1

		# Ajouter le mot à l'outputs du dernier nœud
		self.automaton[current_node]["outputs"].append(word)

	cdef void adds(self, list words):
		"""
		Ajoute une liste de mots à l'automate de recherche.

		Args:
			words (List[str]): Une liste de mots (chaînes de caractères) à ajouter à l'automate.

		Returns:
			None
		"""
		# Initialiser la variables 'word'
		cdef str word = ''

		# Ajouter chaque mot de la liste au Trie
		for word in words:
			self.add(word=word)

	cdef void set_fail_transitions(self):
		"""
		Met en place les liens d'échec pour les nœuds de l'automate Aho-Corasick.

		Returns:
			None
		"""
		# Initialiser les variables
		q: deque = deque()
		cdef int r = 0
		cdef int id_node = 0
		cdef int current_node = 0

		# Pour chaque nœud fils du premier nœud
		for id_node in self.automaton[current_node]["next_nodes"]:
			# Initialiser les liens echecs à zéro
			q.append(id_node)
			self.automaton[id_node]["fail_state"] = 0

		# Tant que l'on a des index de nœuds dans le dèque
		while q:
			# Récupérer le premier index de nœuds du dèque
			r = q.popleft()

			# Pour chaque nœud fils de ce nœud
			for id_node in self.automaton[r]["next_nodes"]:
				# Ajouter les nœuds fils de ce nœud au dèque
				q.append(id_node)
				current_node = self.automaton[r]["fail_state"]
				
				# Parcourir le Trie selon les liens echecs
				while self.next_node(current_node, self.automaton[id_node]["value"]) == NONE and current_node != 0:
					current_node = self.automaton[current_node]["fail_state"]
				
				# Ajouter les liens echecs
				self.automaton[id_node]["fail_state"] = self.next_node(current_node, self.automaton[id_node]["value"])
				if self.automaton[id_node]["fail_state"] == NONE:
					self.automaton[id_node]["fail_state"] = 0

				# Ajouter les mots contenus dans d'autres mots, dans l'outputs des mots les contenant
				self.automaton[id_node]["outputs"].extend(self.automaton[self.automaton[id_node]["fail_state"]]["outputs"])

	cpdef list search_all(self, str text):
		"""
		Recherche tous les mots contenus dans l'automate dans un texte donné.

		Args:
			text (str): Le texte dans lequel effectuer la recherche.

		Returns:
			List[tuple]: Une liste de tuples (mot, position) représentant les mots trouvés dans le texte.
		"""
		# Initialiser les variables
		text = text.lower()
		cdef list results = []
		cdef int current_node = 0
		cdef int pos = 0
		cdef str word = ''

		# Pour chaque lettre du texte
		for pos in range(0, len(text), 1):
			# Aller aux liens echecs quand il n'y a plus matchs avec les mots contenus dans l'automate
			while self.next_node(current_node, text[pos]) == NONE and current_node != 0:
				current_node = self.automaton[current_node]["fail_state"]
				
			# Aller au nœud suivant s'il y a match
			current_node = self.next_node(current_node, text[pos])
			
			# Ajouter la position du mot aux résultats, si le mots est présent dans le texte
			if current_node == NONE:
				current_node = 0
			else:
				for word in self.automaton[current_node]["outputs"]:
					results.append((word, pos-len(word)+1))

		# Retourner le résultats sous forme d'une liste de tuples (mot, position)
		return results


cdef dict search(str text, AhoCorasick searchMachine, dict dict_results):
	"""
	Effectue la recherche de motifs à partir d'un automate de recherche Ahocorasick.

	Args:
		text (str): Le texte dans lequel effectuer la recherche.
		searchMachine (AhoCorasick): L'automate de recherche Ahocorasick préparé pour la recherche.
		dict_results (Dict[str, List[int]]): Un dictionnaire pour stocker les résultats de la recherche.

	Returns:
		Dict[str, List[int]]: Le dictionnaire mis à jour avec les résultats de la recherche.
	"""
	# Initialiser les variables
	cdef list results = searchMachine.search_all(text)
	cdef str motif = ''
	cdef int position = 0

	# Pour chaque résultat, ajouter le au dictionnaire
	for motif, position in results:
		dict_results[motif].append(position)

	# Retourner le dictionnaire
	return dict_results

cpdef int getPosPrimer(str sequence, AhoCorasick searchMachine, list primers, bint forward = True):
	"""
	Récupère le positionnement d'un primer.

	Args:
		sequence (str): La séquence dans laquelle rechercher les primers.
		searchMachine (AhoCorasick): L'automate de recherche préparé pour la recherche de primers.
		primers (List[str]): Une liste de primers à rechercher.
		forward (bool, optional): Indique si la recherche doit être effectuée en avant. (Par défaut = True).

	Returns:
		int: La position du premier primer trouvé.
	"""
	# Initialiser le dictionnaires
	cdef dict results = {}

	# Récupérer les positions des mots contenus dans la machine
	results = search(text=sequence,
					 searchMachine=searchMachine, 
					 dict_results={primer.lower(): [] for primer in primers})

	# Récupérer les primers présents
	results = {primer: positions for (primer, positions) in results.items() if len(positions) != 0}

	# Récupérer les résultats
	if len(results) == 0:
		return NONE
	else:
		if forward == True:
			return min([position for positions in results.values() for position in positions])
		else:
			return max([position+len(primer) for (primer, positions) in results.items() for position in positions])