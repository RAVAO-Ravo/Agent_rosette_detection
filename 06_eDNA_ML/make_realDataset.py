#!/bin/python3
#-*- coding:utf-8 -*-


################################# IMPORTATION DES MODULES #################################

import os
import typing as tp
from argparse import ArgumentParser, HelpFormatter
from reads import Reads
from utilitaire import cherif, exec_cmd, init_dataFrame

################################# IMPORTATION DES MODULES #################################


################################# ARGUMENTS DU PROGRAMME ##################################

formatter: tp.Callable[[str], HelpFormatter] = lambda prog : HelpFormatter(prog=prog, max_help_position=100, width=200)
parser = ArgumentParser(description="Description: Ce programme permet de transformer les fastq des barcodes 6 et du 12 du run de 28/02/2023 en datasets.",
						formatter_class=formatter)
parser.add_argument("-K", "--kmer_size", type=int, default=28, help="k-mers size for minimap2")
parser.add_argument("-Q", "--mapQ", type=float, default=60, help="minimum mapQ to keep sequences")
parser.add_argument("-@", "--n_threads", type=int, default=0, help="number of threads to use for the samtools step")
parser.add_argument("-m", "--min_length", type=int, default=150, help="minimum length to keep sequences")
parser.add_argument("--max_length", type=int, default=160, help="maximum length to keep sequences")
args = parser.parse_args()

################################# ARGUMENTS DU PROGRAMME ##################################


################################### VARIABLES GLOBALES ####################################

K: int = args.kmer_size
mapQ: float = args.mapQ
n_threads: int = args.n_threads
min_length: int = args.min_length
max_length: int = args.max_length
res_directory: str = "dataset"

################################### VARIABLES GLOBALES ####################################


####################################### FONCTIONS #########################################

def inferDataset(outputFile: str) -> None :

	""" Permet la création d'un dataset. """

	# Initialiser les variables globales
	global cherif
	global min_length
	global max_length

	# Récupérer les primers
	forwards = cherif["forwards"]
	reverses = cherif["reverses"]

	# Ouvrir le fastq temporaire
	fastq = Reads(filename="tmp.fastq")
	print(f"Nombre de reads : {fastq.getN_reads()}")

	# Filtrer par rapport aux primers, et filtrer par rapport à la taille des reads
	fastq.cutPrimers(forwards=forwards, reverses=reverses, seuil_length=min_length).filter_size(inf=min_length, sup=max_length)
	print(f"Nombre de reads : {fastq.getN_reads()}")

	# Créer le dataset
	dataset = init_dataFrame(sequences=fastq.getSequences(), n_features=max_length).assign(species="Sphaerothecum_destruens")
	dataset.columns = [i for i in range(0, max_length+1, 1)]
	dataset.to_csv(path_or_buf=outputFile, index=False)

####################################### FONCTIONS #########################################


if __name__ == "__main__" :

	# Créer le répertoire résultats
	if "dataset" not in os.listdir(path='.') :
		exec_cmd(cmd=f"mkdir {res_directory}")

	########## BARCODE 06 ##########

	# Mapper les reads sur la séquence de référence
	exec_cmd(cmd=f"minimap2 -a -x map-ont -k {K} species/Sphaerothecum_destruens.fasta all_reads_BC06.fastq.gz \
				   | samtools view -S -b -F 4 -q {mapQ} \
				   | samtools sort > tmp.bam ; \
				   \
				   samtools fastq -@ {n_threads} tmp.bam > tmp.fastq")
	
	# Générer le testSet
	inferDataset(outputFile=f"{res_directory}/Sd_testSet_BC06.csv")
	
	########## BARCODE 06 ##########


	########## BARCODE 12 ##########

	# Mapper les reads sur la séquence de référence
	exec_cmd(cmd=f"minimap2 -a -x map-ont -k {K} species/Sphaerothecum_destruens.fasta all_reads_BC12.fastq.gz \
				   | samtools view -S -b -F 4 -q {mapQ} \
				   | samtools sort > tmp.bam ; \
				   \
				   samtools fastq -@ {n_threads} tmp.bam > tmp.fastq")

	# Générer le validationSet
	inferDataset(outputFile=f"{res_directory}/Sd_testSet_BC12.csv")

	########## BARCODE 12 ##########

	# Effacer les fichiers temporaires
	exec_cmd(cmd="rm *.bam ; rm *.fastq")