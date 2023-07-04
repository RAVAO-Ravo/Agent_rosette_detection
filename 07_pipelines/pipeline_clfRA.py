#!/bin/python3
#-*- coding:utf-8 -*-


############################### IMPORTATIONS ####################################

import matplotlib.patches
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import subprocess as sp
import time as tm
import typing as tp
from argparse import ArgumentParser, HelpFormatter
from datetime import date
from math import log10
from reads import Reads

############################### IMPORTATIONS ####################################


########################## ARGUMENTS DU PROGRAMME ###############################

formatter: tp.Callable[[str], HelpFormatter] = lambda prog : HelpFormatter(prog=prog, max_help_position=100, width=200)
parser: ArgumentParser = ArgumentParser(description="Ce programme permet la detection de l'agent rosette dans un run de séquençage.", formatter_class=formatter)
parser.add_argument("-r", "--run", type=str, help="Nom du dossier à utiliser, les sous-dossiers doivent être du style 'barcode01'")
parser.add_argument("-o", "--outputDir", type=str, help="Nom du dossier de sortie")
parser.add_argument("-s", "--seuil", type=float, default=0.90, help="Le seuil pour retenir un read")
parser.add_argument("-b", "--begin", type=int, default=1, help="Premier barcode à utiliser")
parser.add_argument("-e", "--end", type=int, help="Dernier barcode à utiliser")
parser.add_argument("-t", "--threads", type=int, default=18, help="Nombre de threads pour la génération des fastqc")
parser.add_argument("-f", "--figsize", type=int, nargs=2, default=(30, 20), help="Taille de la figure")
args = parser.parse_args()

########################## ARGUMENTS DU PROGRAMME ###############################


########################## PARAMÈTRES PRÉDICTIONS ###############################

run: str = args.run
outputDir: str = args.outputDir
seuil: float = args.seuil
start: int = args.begin
end: int = args.end
barcodes: tp.List[str] = [f"0{i}" if i < 10 else f"{i}" for i in range(start, end+1, 1)]
threads: int = args.threads

########################## PARAMÈTRES PRÉDICTIONS ###############################


################################# SUMMARIES #####################################

col_names: tp.List[str] = ["barcode",
				   		   "n_reads_initial",
				   		   "n_reads_withPrimers",
				   		   "n_reads_withGoodSize",
				   		   f"RA_amplicons P(>={seuil})",
				   		   "pourcentage"]
summary: pd.DataFrame = pd.DataFrame(columns=col_names)
figsize: tp.Tuple[int] = tuple(args.figsize)
colors: tp.Tuple[str, str, str, str] = ("black", "blue", "green", "red")

################################# SUMMARIES #####################################


########################## FONCTIONS UTILITAIRES ################################

# Exécute une commande dans le shell
exec_cmd: tp.Callable[[str], None] = lambda cmd : sp.run(args=[cmd], shell=True)


# Retourne le log de 10 d'un nombre ou 0 sinon
log10Transfomer : tp.Callable[[int], float] = lambda x : log10(x) if x != 0 else 0


def autolabel(rects: tp.List[matplotlib.patches.Rectangle]) -> None:

	""" Annote les rectangles du barplot avec les pourcentages associés. """

	# Pour chaque rectangle
	for rect in rects :
		
		# Annoter le rectangle
		height = rect.get_height()
		rect.axes.annotate(text=round(height, 2),
						   xy=((rect.get_x()+rect.get_width()/2), height),
						   xytext=(0, 3), 
						   textcoords="offset points",
						   ha="center", 
						   va="bottom")


def makePlotSummary(summaryFile: str, 
					figsize: tp.Tuple[int, int]=figsize,
					colors: tp.Tuple[str, str, str, str]=colors,
					sep: str='\t', 
					alpha: float=0.75,
					style: str="dark",
					ylabel: str="Log10 du nombre de reads",
					xlabel: str="Barcodes",
					title: str="Log10 du nombre de reads, selon l'espèce, et selon le barcode",
					edgeColor: str="black") -> None :
	
	""" Créer le barplot du summary final. """

	# Récupérer le nom des colonnes
	global seuil, col_names

	# Ouvrir le fichier summary
	data = pd.read_csv(filepath_or_buffer=summaryFile, sep=sep)
	data['A'] = data["n_reads_initial"].apply(log10Transfomer)
	data['B'] = data["n_reads_withPrimers"].apply(log10Transfomer)
	data['C'] = data["n_reads_withGoodSize"].apply(log10Transfomer)
	data['D'] = data[f"RA_amplicons P(>={seuil})"].apply(log10Transfomer)

	# Créer une figure
	plt.figure(figsize=figsize)
	
	# Créer le barplot
	sns.set_style(style=style)

	# Pour chaque élément de la pile
	for i, y in enumerate(['A', 'B', 'C', 'D']) :

		# Construire le barplot
		fig = sns.barplot(data=data, x=col_names[0], y=y, color=colors[i], alpha=alpha, label=col_names[i+1])

		# Annoter le barplot
		autolabel(rects=[rect for rect in fig.patches if rect.get_height() != 0.0])

		# Changer le contour du barplot
		[child.set_edgecolor(color=edgeColor) for child in fig.get_children() if isinstance(child, matplotlib.patches.Rectangle)]
	
	# Paramètrer la figure
	plt.legend(loc="best")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	
	# Sauvegarder la figure
	plt.savefig(summaryFile.split(sep='.')[0]+".png")
	
########################## FONCTIONS UTILITAIRES ################################


################################# PIPELINE ######################################

def makeWorkDir(workDir: str, barcode: str, all_reads: str) -> None :
	
	""" Crée le répertoire de travail pour l'exécution du pipeline. """

	# Initialiser les variables globales
	global run, exec_cmd

	# Si le répertoire de travaille n'existe pas, créer le répertoire
	if workDir not in os.listdir() :
		exec_cmd(f"mkdir {workDir}")
	
	# Copier les '*.gz' contenant les reads dans le répertoire de travail 
	exec_cmd(f"cp {run}/barcode{barcode}/*.gz {workDir}")

	# Décompresser les '*.gz'
	exec_cmd(f"gzip -d {workDir}/*.gz")

	# Concatener tout les fastqs dans un seuil fichier
	exec_cmd(f"cat {workDir}/*.fastq > {all_reads}.protected")

	# Supprimer tout les fastqs excepté le fichier concaténé
	exec_cmd(f"rm {workDir}/*.fastq")
	exec_cmd(f"mv {all_reads}.protected {all_reads}")
	

def qualityCheck(barcode: str, preFastqc: str, all_reads: str, reads_cleaning: str, postFastqc: str) -> None :

	""" Effectue le contrôle qualité, et filtre les reads. """

	# Définir les variables globales
	global seuil, threads, barcodes, col_names, summary

	# Générer les fastqc avant prétraitements
	exec_cmd(f"mkdir {preFastqc} ; fastqc -t {threads} {all_reads} -o {preFastqc}")

	# Filtrer les reads
	# ============================================================================================== #

	fastq = Reads(filename=all_reads)
	print(f"Nombre de reads : {fastq.getN_reads()}")

	initial, primer_cutting, size_filter, final = fastq.clf_RA(seuil=seuil)
	print(f"Nombre de reads avec primers: {primer_cutting}")
	print(f"Nombre de reads avec la bonne taille: {size_filter}")
	print(f"Nombre d'amplicons d'agent rosette: {final}")

	# Ajouter une nouvelle ligne au summary
	new_line = pd.DataFrame(data={col_names[0] : [barcode],
				   				  col_names[1] : [initial],
								  col_names[2] : [primer_cutting],
								  col_names[3] : [size_filter],
								  col_names[4] : [final],
								  col_names[5] : [(final/size_filter)*100 if size_filter != 0 else 0]})
	summary = pd.concat(objs=[summary, new_line], axis=0)

	# Enregistrer les reads restants
	fastq.write(filename=reads_cleaning, outputFormat="fq")

	# ============================================================================================== #

	# Générer les fastqc après prétraitements
	exec_cmd(f"mkdir {postFastqc} ; fastqc -t {threads} {reads_cleaning} -o {postFastqc} ; gzip {all_reads} ; gzip {reads_cleaning}")


def pullRes(day: str, map_res: str) -> None :

	""" Regroupe les résultats de prédictions de chaque barcode. """

	# Définir les variables globales
	global col_names, summary, figsize, colors

	# Créer le répertoire de regroupement
	if map_res not in os.listdir() :
		exec_cmd(f"mkdir {map_res}")

	# Pour chaque élément du répertoire courrant
	for element in os.listdir() :

		# Regrouper dans le répertoire dédié les résultats de prédictions de chaque barcode
		if "barcode" in element and "_predictions" in element :
			exec_cmd(f"cp -r {element}/ {map_res}/")
			exec_cmd(f"rm -r {element}/")

	# Enregistrer le summary global sous le format '.csv'
	summary.reset_index(drop=True).to_csv(path_or_buf=f"{map_res}/summary_predictions_{day}.csv", sep='\t', index=False)
	
	# Créer le barplot du summary global
	makePlotSummary(summaryFile=f"{map_res}/summary_predictions_{day}.csv")

################################# PIPELINE ######################################


if __name__ == "__main__" :

	# Définir le jour actuel
	day = date.today().strftime("%Y_%m_%d")

	# Définir le répertoire de regroupement
	map_res = f"{outputDir}_{day}"

	# Définir la balise start
	d = tm.perf_counter()

	# Pour chaque barcode
	for barcode in barcodes :

		# Définir les informations nécéssaires pour les prédictions
		workDir = f"barcode{barcode}_predictions/"
		all_reads = f"{workDir}/all_reads_BC{barcode}.fastq"
		postFastqc = f"{workDir}/postCleaningGraphs"
		preFastqc = f"{workDir}/preCleaningGraphs"
		reads_cleaning = f"{workDir}/reads_BC{barcode}_PC_SF.fastq"

		# Créer le répertoire de travail pour le barcode
		makeWorkDir(workDir=workDir, barcode=barcode, all_reads=all_reads)
		
		# Effectuer le contrôle qualité, ainsi que les différents filtrages
		qualityCheck(barcode=barcode, preFastqc=preFastqc, all_reads=all_reads, reads_cleaning=reads_cleaning, postFastqc=postFastqc)
	
	# Regrouper tout les résultats de prédictions, de chaque barcode, dans le répertoire de regroupement
	pullRes(day=day, map_res=map_res)

	# Définir la balise stop
	f = tm.perf_counter()

	# Afficher le temps d'exécution
	print(f"\nTime : {round(number=f-d, ndigits=0)} s\n")