#!/bin/python3
#-*- coding:utf-8 -*-


############################### IMPORTATIONS ####################################

import os
import subprocess as sp
import typing as tp
from argparse import ArgumentParser, HelpFormatter

############################### IMPORTATIONS ####################################


################################ ARGUMENTS ######################################

formatter: tp.Callable[[str], HelpFormatter] = lambda prog : HelpFormatter(prog=prog, max_help_position=100, width=200)
parser: ArgumentParser = ArgumentParser(description="Description: Ce programme permet de concaténer des dossiers de runs.", formatter_class=formatter)
parser.add_argument("-r", "--runs", type=str, nargs='*', help="Liste des dossiers de runs à concaténer")
parser.add_argument("-b", "--begin", type=int, default=1, help="Index du premier barcode à utiliser")
parser.add_argument("-e", "--end", type=int, help="Index du dernier barcode à utiliser")
parser.add_argument("-o", "--outputDir", type=str, help="Nom du dossier de parent")
args = parser.parse_args()

################################ ARGUMENTS ######################################


############################ VARIABLES GLOBALES #################################

runs: tp.List[str] = args.runs
begin: int = args.begin
end: int = args.end
outputDir: str = args.outputDir
barcodes: tp.List[str] = [f"0{i}" if i < 10 else f"{i}" for i in range(begin, end+1, 1)]

############################ VARIABLES GLOBALES #################################


########################## FONCTIONS UTILITAIRES ################################

exec_cmd: tp.Callable[[str], None] = lambda cmd : sp.run(args=[cmd], shell=True)

########################## FONCTIONS UTILITAIRES ################################


if __name__ == "__main__" :

	# Créer le dossier de pulling
	if outputDir not in os.listdir() :
		exec_cmd(f"mkdir {outputDir}")
	
	# Créer les sous-dossiers
	for barcode in barcodes :
		exec_cmd(f"mkdir {outputDir}/barcode{barcode}")

	# Puller les runs
	for run in  runs :
		for barcode in barcodes :
			if len(os.listdir(f"{run}/barcode{barcode}")) != 0 :
				exec_cmd(f"cp {run}/barcode{barcode}/* {outputDir}/barcode{barcode}")