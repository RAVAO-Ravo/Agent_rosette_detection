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


################################## PRIMERS ######################################

carpio: tp.Dict[str, list] = {"forwards" : ["GGTGGGTTTCAGTAGACAATGC"],
	   	  					  "reverses" : ["GGCGGCAATAACAAATGGTAGT"]}
carpioSize_Range: tp.Tuple[int, int] = (70, 85)


cherif: tp.Dict[str, list] = {"forwards" : ["GCGGTAATTCCAGCTCCA"],
	   	  					  "reverses" : ["CACTCAATTAAGCGCACACG"]}
cherifSize_range: tp.Tuple[int, int] = (150, 160)


gozlan: tp.Dict[str, list] = {"forwards" : ["AATCGTATGACATTTTGTCGAC", "CAGGGCTTTTTAAGTCTT"],
	   	  					  "reverses" : ["GAAGTCACAGGCGATTCGG", "TGATGGAGTCATAGAATTAACATCC"]}
gozlanSize_range: tp.Tuple[int, int] = (580, 1120)


parva: tp.Dict[str, list] = {"forwards" : ["CGAGCCAAATAACAGAGGGT"],
	   	 					 "reverses" : ["CAGGCGAGGCTTATGTTTGC"]}
parvaSize_range: tp.Tuple[int, int] = (140, 155)

################################## PRIMERS ######################################


########################## ARGUMENTS DU PROGRAMME ###############################

formatter: tp.Callable[[str], HelpFormatter] = lambda prog : HelpFormatter(prog=prog, max_help_position=100, width=200)
parser: ArgumentParser = ArgumentParser(description="Description: Ce programme permet le mapping d'un run sur des espèces", formatter_class=formatter)
parser.add_argument("-r", "--run", type=str, help="Nom du dossier à utiliser, les sous-dossiers doivent être du style 'barcode01'")
parser.add_argument("-s", "--speciesDir", type=str, default="species/", help="Le dossier contenant les espèces à mapper")
parser.add_argument("-o", "--outputDir", type=str, help="Nom du dossier de sortie")
parser.add_argument("-p", "--primers", type=str, help="Nom des couples d'amorces à utiliser")
parser.add_argument("-b", "--begin", type=int, default=1, help="Premier barcode à utiliser")
parser.add_argument("-e", "--end", type=int, help="Dernier barcode à utiliser")
parser.add_argument("-k", "--kmer_size", type=int, default=28, help="Taille des kmers")
parser.add_argument("-t", "--threads", type=int, default=18, help="Nombre de threads pour la génération des fastqc")
parser.add_argument("-f", "--figsize", type=int, nargs=2, default=(30, 20), help="Taille de la figure")
parser.add_argument("-l", "--spToPlot", type=str, nargs='*', default=["Sphaerothecum_destruens01"], 
					help="Les espèces à plotter, 3 espèces max, noms des fichiers sans les extensions")
args = parser.parse_args()

########################## ARGUMENTS DU PROGRAMME ###############################


############################ PARAMÈTRES MAPPING #################################

run: str = args.run
speciesDir: str = args.speciesDir
outputDir: str = args.outputDir
primers: str = args.primers
K: int = args.kmer_size
threads: int = args.threads
start: int = args.begin
end: int = args.end
barcodes: tp.List[str] = [f"0{i}" if i < 10 else f"{i}" for i in range(start, end+1, 1)]
allSpecies: tp.List[str] = sorted(os.listdir(path=speciesDir))

############################ PARAMÈTRES MAPPING #################################


################################# SUMMARIES #####################################

columns: tp.List[str] = ["barcodes",
						 "species",
						 "QC_passedReads",
						 "QC_failedReads",
						 "mappedReads(count)",
						 "mappedReads(%)",
						 "startPos",
						 "endPos",
						 "covBases",
						 "coverage",
						 "meanDepth",
						 "meanBaseQ",
						 "meanMapQ" ]
summary: pd.DataFrame = pd.DataFrame(columns=columns)
figsize: tp.Tuple[int, int] = tuple(args.figsize)
species: tp.List[str] = args.spToPlot
colors: tp.Tuple[str, str, str] = ("blue", "red", "green")

################################# SUMMARIES #####################################


########################## FONCTIONS UTILITAIRES ################################

# Exécute des commandes dans le shell
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
					species: tp.List[str]=species,
					figsize: tp.Tuple[int]=figsize,
					colors: tp.Tuple[str]=colors,
					sep: str='\t',
					alpha: float=0.75,
					style: str="dark",
					ylabel: str="Log10 du nombre de reads mappés",
					xlabel: str="Barcodes",
					title: str="Log10 du nombre de reads mappés, selon l'espèce, et selon le barcode",
					edgeColor: str="black") -> None :

	""" Crée le barplot du summary file global. """

	# Récupérer le fichier csv sous forme de dataframe
	df = pd.read_csv(filepath_or_buffer=str(summaryFile), sep=sep)

	# Filtrer le dataframe, en fonction des espèces cibles
	df = df[df["species"].isin(values=species)]

	# Créer la colonne du log10 du nombre de reads mappés
	df['A'] = df["mappedReads(count)"].apply(log10Transfomer)

	# Initialiser une figure
	plt.figure(figsize=figsize)

	# Définir le style de la figure
	sns.set_style(style=style)

	# Apposer une grille
	plt.grid()

	# Créer le barplot
	fig = sns.barplot(data=df, x="barcodes", y='A', hue="species", palette=colors, alpha=alpha)
	
	# Définir la couleur du contour des rectangles
	[child.set_edgecolor(color=edgeColor) for child in fig.get_children() if isinstance(child, matplotlib.patches.Rectangle)]

	# Annoter les axes
	plt.ylabel(ylabel=ylabel)
	plt.xlabel(xlabel=xlabel)

	# Titrer la figure
	plt.title(label=title)

	# Annoter les rectangles
	autolabel(rects=[rect for rect in fig.patches if rect.get_height() != 0.0])

	# Sauvegarder les figures
	plt.savefig(summaryFile.split(sep='.')[0]+".png")


def PMR_Summary(summary_file: str, outputFile: str, sep: str='\t') -> None :

	""" Crée le summary du pourcentage de reads mappés, selon l'espèce, et selon le barcode. """
	
	# Initialiser les variables globales 
	global start, end

	# Ouvrir le fichier summary global 
	df  = pd.read_csv(filepath_or_buffer=summary_file, sep=sep)
	
	# Récupérer le nom des espèces, et les barcodes 
	species = sorted(list(set(df["species"])))
	bc = pd.DataFrame(data={"barcode" : [i for i in range(start, end+1, 1)]}).reset_index(drop=True)
	
	# Initialiser le PMR summary 
	pmr_summary = pd.DataFrame(columns=species)

	# Pour chaque barcode 
	for i in range(start, end+1, 1) :

		# Récupérer le pourcentage de reads mappés
		tmp = df.query(f"barcodes == {i}").transpose().iloc[1:, :].reset_index(drop=False)
		tmp = tmp.query("index == 'mappedReads(%)'").iloc[:, 1:]
		
		# Attribuer chaque colonne à l'espèce associé
		tmp.columns = species
		
		# Ajouter la nouvelle ligne au PMR summary
		pmr_summary = pd.concat(objs=[pmr_summary, tmp], axis=0)

	# Enregister le PMR summary
	pmr_summary = pmr_summary.reset_index(drop=True)
	pmr_summary = pd.concat(objs=[bc, pmr_summary], axis=1).reset_index(drop=True).to_csv(path_or_buf=outputFile, sep='\t', index=False)


def getAccesNumbers(outputFile: str) -> None :
	
	""" Récupère les espèces utilisées pour le mapping. """

	# Définir la variable global specieDir
	global speciesDir
	
	# Initialiser le dictionnaire
	dic = {"Accession" : [], "Species" : [], "Description" : []}

	# Pour chaque espèce, ajouter ses informations au dictionnaire
	for file in os.listdir(path=speciesDir) :
		if file.endswith(".fasta") :
			with open(file=speciesDir+file, mode='r') as inputFile :
				firstLine = inputFile.readline().split()
				dic["Accession"].append(firstLine[0][1:])
				dic["Species"].append('_'.join(firstLine[1:3]))
				dic["Description"].append('_'.join(firstLine[3:]))

	# Enregistrer le résultat
	pd.DataFrame(data=dic).to_csv(path_or_buf=outputFile, sep='\t', index=False)


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

	# Concatener tout les fastqs dans un seul fichier
	exec_cmd(f"cat {workDir}/*.fastq > {all_reads}.protected")

	# Supprimer tout les fastqs excepté le fichier concaténé
	exec_cmd(f"rm {workDir}/*.fastq")
	exec_cmd(f"mv {all_reads}.protected {all_reads}")


def filter_reads(fastq: Reads, forwards: tp.List[str], reverses: tp.List[str], inf: int, sup: int) -> None :

	""" Filtre les reads selon des primers, et selon leur taille. """

	# Exciser les reads, au niveau des primers, et filtrer selon la taille des reads
	fastq.cutPrimers(forwards=forwards, reverses=reverses, seuil_length=inf).filter_size(inf=inf, sup=sup)

	# Afficher le nombre de reads ayant passés les filtres
	print(f"Nombre de reads : {fastq.getN_reads()}")


def qualityCheck(preFastqc: str, all_reads: str, reads_cleaning: str, postFastqc: str, primers: str) -> None :

	""" Effectue le contrôle qualité, et filtre les reads. """

	# Définir les variable globales
	global gozlan, gozlanSize_range
	global cherif, cherifSize_range
	global carpio, carpioSize_Range
	global parva, parvaSize_range
	global barcodes, threads, exec_cmd

	# Générer les fastqc avant prétraitements
	exec_cmd(f"mkdir {preFastqc} ; fastqc -t {threads} {all_reads} -o {preFastqc}")

	# Filtrer les reads
	# ==================================================================================== #

	fastq = Reads(filename=all_reads)
	print(f"Nombre de reads : {fastq.getN_reads()}")

	if primers == "carpio" :
		forwards, reverses, inf, sup = carpio["forwards"], carpio["reverses"], carpioSize_Range[0], carpioSize_Range[1]
		filter_reads(fastq=fastq, forwards=forwards, reverses=reverses, inf=inf, sup=sup)

	elif primers == "cherif" :
		forwards, reverses, inf, sup = cherif["forwards"], cherif["reverses"], cherifSize_range[0], cherifSize_range[1]
		filter_reads(fastq=fastq, forwards=forwards, reverses=reverses, inf=inf, sup=sup)

	elif primers == "gozlan" :
		forwards, reverses, inf, sup = gozlan["forwards"], gozlan["reverses"], gozlanSize_range[0], gozlanSize_range[1]
		filter_reads(fastq=fastq, forwards=forwards, reverses=reverses, inf=inf, sup=sup)

	elif primers == "parva" :
		forwards, reverses, inf, sup = parva["forwards"], parva["reverses"], parvaSize_range[0], parvaSize_range[1]
		filter_reads(fastq=fastq, forwards=forwards, reverses=reverses, inf=inf, sup=sup)

	fastq.write(filename=reads_cleaning, outputFormat="fq")

	# ==================================================================================== #

	# Générer les fastqc après prétraitements
	exec_cmd(f"mkdir {postFastqc} ; fastqc -t {threads} {reads_cleaning} -o {postFastqc} ; gzip {all_reads}")


def mapping(workDir: str, reads_cleaning: str, barcode: str) -> None :

	""" Mappe les reads sur un ensemble d'espèces. """

	# Définir les variables globales
	global allSpecies, columns, summary, exec_cmd

	# Créer le répertoire contenant les résultats de mapping
	exec_cmd(f"mkdir {workDir}/mapping_res")

	# Initialiser le summary local (celui du barcode traité actuellement)
	summaryBC = pd.DataFrame(columns=columns)

	# Pour chaque espèce de l'ensemble d'espèces
	for specie in allSpecies :

		# Créer son répertoire résultat
		specieFileRes = f"{workDir}/mapping_{specie.split(sep='.')[0]}"
		exec_cmd(f"mkdir {specieFileRes}")

		# Mapper les reads sur l'espèce
		exec_cmd(f"minimap2 -a -k {K} -x map-ont {speciesDir}/{specie} {reads_cleaning} \
				| samtools view -S -b \
				| samtools sort -o {specieFileRes}.bam")

		# Enregistrer les résultats de mapping
		exec_cmd(f"samtools coverage {specieFileRes}.bam > {specieFileRes}.cov")
		exec_cmd(f"samtools coverage -m {specieFileRes}.bam > {specieFileRes}.covHist")
		exec_cmd(f"samtools depth {specieFileRes}.bam > {specieFileRes}.depth")
		exec_cmd(f"samtools flagstat -O tsv {specieFileRes}.bam > {specieFileRes}.stat")

		# Récupérer les informations contenues dans le fichier '.stat'
		df = pd.read_csv(filepath_or_buffer=f"{specieFileRes}.stat", names=["col1", "col2", "descriptor"], sep='\t').fillna(value=0)
		passed = int(df.at[0, "col1"])
		failed = int(df.at[0, "col2"])
		mapped = int(df.at[4, "col1"])
		mappedPerc = float(df.at[5, "col1"].strip('%')) if isinstance(df.at[5, "col1"], str) else float(df.at[5, "col1"])

		# Récupérer les informations contenues dans le fichier '.cov'
		df = pd.read_csv(filepath_or_buffer=f"{specieFileRes}.cov", sep='\t')

		# Récupérer la range de mapping
		startPos = int(df.loc[:, ("startpos")].min())
		endPos = int(df.loc[:, ("endpos")].max())

		# Récupérer que les régions non-nulles
		df = df.query("covbases != 0 & coverage != 0 & meandepth != 0 & meanbaseq != 0 & meanbaseq != 0 & meanmapq != 0")
		if df.empty == False :
			covBases = float(df.loc[:, ("covbases")].sum())
			coverage = float(df.loc[:, ("coverage")].sum())
			meanDepth = float(df.loc[:, ("meandepth")].mean())
			meanBaseQ = float(df.loc[:, ("meanbaseq")].mean())
			meanMapQ = float(df.loc[:, ("meanmapq")].mean())
		else :
			covBases = 0
			coverage = 0
			meanDepth = 0
			meanBaseQ = 0
			meanMapQ = 0

		# Créer une nouvelle ligne pour le summary local
		df = pd.DataFrame(data={"barcodes" : [barcode],
								"species" : [specie.split(sep='.')[0]],
								"QC_passedReads" : [passed],
								"QC_failedReads" : [failed],
								"mappedReads(count)" : [mapped],
								"mappedReads(%)" : [mappedPerc],
								"startPos" : [startPos],
								"endPos" : [endPos],
								"covBases" : [covBases],
								"coverage" : [coverage],
								"meanDepth" : [meanDepth],
								"meanBaseQ" : [meanBaseQ],
								"meanMapQ" : [meanMapQ]})

		# Ajouter les informations au summary local, et summary global (tout les barcodes pour chaque espèce)
		summaryBC = pd.concat([summaryBC, df], axis=0)
		summary = pd.concat([summary, df], axis=0)

		# Ajouter les résultats de mapping au répertoire résultat dédié
		exec_cmd(f"mv 	{specieFileRes}.bam \
						{specieFileRes}.cov \
						{specieFileRes}.covHist \
						{specieFileRes}.depth \
						{specieFileRes}.stat \
						{specieFileRes}")
		exec_cmd(f"mv {specieFileRes} {workDir}/mapping_res")
	
	# Enregistrer le summary local sous le format '.csv', et gzipper les reads filtrés
	summaryBC.to_csv(f"{workDir}/summary_BC{barcode}.csv", sep='\t', index=False)
	exec_cmd(f"gzip {reads_cleaning}")


def pullRes(day: str, map_res: str) -> None :

	""" Regroupe les résultats de mapping de chaque barcode. """

	# Définir les variables globales
	global summary, exec_cmd

	# Créer le répertoire de regroupement
	if map_res not in os.listdir() :
		exec_cmd(f"mkdir {map_res}")

	# Pour chaque élément du répertoire courrant
	for element in os.listdir() :

		# Regrouper dans le répertoire dédié les résultats de mapping de chaque barcode
		if "barcode" in element and "_mapping" in element :
			exec_cmd(f"cp -r {element}/ {map_res}/")
			exec_cmd(f"rm -r {element}/")

	# Enregistrer le summary global sous le format '.csv'
	summary.reset_index(drop=True).to_csv(path_or_buf=f"{map_res}/summary_mapping_{day}.csv", sep='\t', index=False)
	
	# Créer le barplot du summary global
	makePlotSummary(summaryFile=f"{map_res}/summary_mapping_{day}.csv")

	# Créer le PMR summary
	PMR_Summary(summary_file=f"{map_res}/summary_mapping_{day}.csv", outputFile=f"{map_res}/summary_PMR_{day}.csv")

	# Ajouter le fichier d'epsèces
	getAccesNumbers(outputFile=f"{map_res}/species.csv")

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

		# Définir les informations nécéssaires pour le mapping
		workDir = f"barcode{barcode}_mapping/"
		all_reads = f"{workDir}/all_reads_BC{barcode}.fastq"
		postFastqc = f"{workDir}/postCleaningGraphs"
		preFastqc = f"{workDir}/preCleaningGraphs"
		reads_cleaning = f"{workDir}/reads_BC{barcode}_PC_SF.fastq"

		# Créer le répertoire de travaille pour le barcode
		makeWorkDir(workDir=workDir, barcode=barcode, all_reads=all_reads)
		
		# Effectuer le contrôle qualité, ainsi que les différents filtrages
		qualityCheck(preFastqc=preFastqc, all_reads=all_reads, reads_cleaning=reads_cleaning, postFastqc=postFastqc, primers=primers)

		# Mapper les reads nettoyés sur l'ensemble des espèces d'intérêts
		mapping(workDir=workDir, reads_cleaning=reads_cleaning, barcode=barcode)
	
	# Regrouper tout les résultats de mapping, de chaque barcode, dans le répertoire de regroupement
	pullRes(day=day, map_res=map_res)

	# Définir la balise stop
	f = tm.perf_counter()

	# Afficher le temps d'exécution
	print(f"\nTime : {round(number=f-d, ndigits=0)} s\n")