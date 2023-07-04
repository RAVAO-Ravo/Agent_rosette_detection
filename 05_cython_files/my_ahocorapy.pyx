#-*- coding:utf-8 -*-
#cython: language_level=3


# Importation des modules
from collections import deque


# Valeur du none
cdef int none = (-99)


cdef class ahocorasick:

	""" Classe ahocorasick permettant de créer un automate de recherche de motifs, selon l'algorithme 'aho-corasick'. """

	# L'automate de recherche, sous forme d'une liste de dictionnaires, chaque dictionnaire étant un noeud de Trie
	cdef list automaton


	def __init__(self, list words) :
		
		""" Constructeur de la machine de recherche. """

		# Initialiser la racine du Trie
		self.automaton = [{"value" : '', "next_nodes" : [], "fail_state" : 0, "outputs" : []}]

		# Ajouter les mots au Trie
		self.adds(words=words)
		
		# Créer les liens echecs
		self.set_fail_transitions()


	cdef int next_node(self, int current_node, str value) :

		""" Récupère l'index du noeud suivant contenant la valeur demandée ('value'), depuis un point de départ du Trie ('current_node'). """

		# Initialiser la variable id_node
		cdef int id_node = 0
		
		# Pour chaque noeud fils du noeud courrant
		for id_node in self.automaton[current_node]["next_nodes"] :

			# Si un noeud fils contient la valeur demandée, retourner le
			if self.automaton[id_node]["value"] == value :
				return id_node
		
		# Retourner None si la valeur n'a pas été trouvée
		return none


	cdef void add(self, str word) :
		
		""" Ajoute un mot à l'automate de recherche. """
		
		# Initialiser les variables
		word = word.lower()
		cdef int current_node = 0
		cdef int i = 0
		cdef int pos = 0
		cdef int len_word = len(word)
		cdef int next_node = self.next_node(current_node, word[pos])

		# Tant qu'une partie du mot (préfixe) est présent dans le Trie
		while next_node != none :
			
			# Passer au noeud, et à la lettre, suivant(e)
			current_node = next_node
			pos+=1

			# Sortir de la boucle si on a parcouru toute les lettres du mot
			if pos < len_word :
				next_node = self.next_node(current_node, word[pos])
			else:
				break

		# Pour chaque lettre, du mot, non présent dans le Trie
		for i in range(pos, len_word, 1) :

			# Ajouter la lettre
			self.automaton.append({"value" : word[i], "next_nodes" : [], "fail_state" : 0, "outputs" : []})
			self.automaton[current_node]["next_nodes"].append(len(self.automaton)-1)
			current_node = len(self.automaton)-1

		# Ajouter le mot à l'outputs du dernier noeud
		self.automaton[current_node]["outputs"].append(word)


	cdef void adds(self, list words) :

		""" Ajoute une liste de mots à l'automate de recherche """

		# Initialiser la variables 'word'
		cdef str word = ''

		# Ajouter chaque mot de la liste au Trie
		for word in words :
			self.add(word=word)


	cdef void set_fail_transitions(self) :

		""" Met en place les liens echecs. """

		# Initialiser les variables
		q: deque = deque()
		cdef int r = 0
		cdef int id_node = 0
		cdef int current_node = 0

		# Pour chaque noeud fils du premier noeud
		for id_node in self.automaton[0]["next_nodes"] :
			
			# Initialiser les liens echecs à zéro
			q.append(id_node)
			self.automaton[id_node]["fail_state"] = 0

		# Tant que l'on n'a des index de noeuds dans le dèque
		while q :
			
			# Récupérer le premier index de noeuds du dèque
			r = q.popleft()

			# Pour chaque noeud fils de ce noeud
			for id_node in self.automaton[r]["next_nodes"] :
				
				# Ajouter les noeuds fils de ce noeud au dèque
				q.append(id_node)
				current_node = self.automaton[r]["fail_state"]
				
				# Parcourir le Trie selon les liens echecs
				while self.next_node(current_node, self.automaton[id_node]["value"]) == none and current_node != 0 :
					current_node = self.automaton[current_node]["fail_state"]
				
				# Ajouter les liens echecs
				self.automaton[id_node]["fail_state"] = self.next_node(current_node, self.automaton[id_node]["value"])
				if self.automaton[id_node]["fail_state"] == none :
					self.automaton[id_node]["fail_state"] = 0

				# Ajouter les mots contenus dans d'autres mots, dans l'outputs des mots les contenant
				self.automaton[id_node]["outputs"].extend(self.automaton[self.automaton[id_node]["fail_state"]]["outputs"])


	cpdef list search_all(self, str text) :
		
		""" Recherche les mots contenus dans l'automate dans un texte. """
		
		# Initialiser les variables
		text = text.lower()
		cdef list results = []
		cdef int current_node = 0
		cdef int pos = 0
		cdef str word = ''

		# Pour chaque lettre du texte
		for pos in range(0, len(text), 1) :
			
			# Aller aux liens echecs quand il n'y a plus matchs avec les mots contenus dans l'automate
			while self.next_node(current_node, text[pos]) == none and current_node != 0 :
				current_node = self.automaton[current_node]["fail_state"]
				
			# Aller au noeud suivant s'il y a match
			current_node = self.next_node(current_node, text[pos])
			
			# Ajouter la position du mot aux résultats, si le mots est présent dans le texte
			if current_node == none :
				current_node = 0
			else :
				for word in self.automaton[current_node]["outputs"] :
					results.append((word, pos-len(word)+1))

		# Retourner le résultats sous forme d'une liste de tuples (mot, position)
		return results