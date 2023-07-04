#!/bin/python3
#-*- coding : utf-8 -*-


# Importation
import subprocess as sp
import typing as tp


# Fonctions utilitaires
exec_cmd: tp.Callable[[str], None] = lambda cmd : sp.run(args=[cmd], shell=True)


# Commandes
cmds: tp.List[str] = ["./pipeline_mapping.py -r run_2023_02_28/ -s species_destruens_01/ -o run_2023_02_28_gozlan_map -p gozlan -b 1 -e 6 -l Sphaerothecum_destruens01 Pimephales_promelas",
					  "./pipeline_mapping.py -r run_2023_02_28/ -s species_destruens_01/ -o run_2023_02_28_cherif_map -p cherif -b 7 -e 12 -l Sphaerothecum_destruens01 Pimephales_promelas",
					  "./pipeline_clfRA.py -r run_2023_02_28/ -o run_2023_02_28_pred -s 0.9 -b 1 -e 12",
					  "./pipeline_mapping.py -r run_2023_04_12/ -s species_destruens_02/ -o run_2023_04_12_cherif_map -p cherif -b 1 -e 12 -l Sphaerothecum_destruens01",
					  "./pipeline_clfRA.py -r run_2023_04_12/ -o run_2023_04_12_pred -s 0.9 -b 1 -e 12",
					  "./pipeline_mapping.py -r run_2023_04_25/ -s species_destruens_02/ -o run_2023_04_25_cherif_map -p cherif -b 1 -e 23 -l Sphaerothecum_destruens01",
					  "./pipeline_mapping.py -r run_2023_04_25/ -s species_carpio/ -o run_2023_04_25_carpio_map -p carpio -b 24 -e 24 -l Cyprinus_carpioCytB",
					  "./pipeline_clfRA.py -r run_2023_04_25/ -o run_2023_04_25_pred -s 0.9 -b 1 -e 23",
					  "./pipeline_mapping.py -r run_2023_05_11/ -s species_destruens_02/ -o run_2023_05_11_cherif_map -p cherif -b 1 -e 24 -l Sphaerothecum_destruens01",
					  "./pipeline_clfRA.py -r run_2023_05_11/ -o run_2023_05_11_pred -s 0.9 -b 1 -e 24",
					  "./pipeline_mapping.py -r run_2023_04_27/ -s species_parva/ -o run_2023_04_27_parva_map -p parva -b 1 -e 24 -l Pseudorasbora_parva16S",
					  "./pipeline_mapping.py -r run_2023_05_09/ -s species_parva/ -o run_2023_05_09_parva_map -p parva -b 1 -e 24 -l Pseudorasbora_parva16S",
					  "./pipeline_mapping.py -r run_2023_06_20/ -s species_parva/ -o run_2023_06_20_parva_map01 -p parva -b 1 -e 12 -l Pseudorasbora_parva16S",
					  "./pipeline_mapping.py -r run_2023_06_20/ -s species_parva/ -o run_2023_06_20_parva_map02 -p parva -b 14 -e 24 -l Pseudorasbora_parva16S"]


if __name__ == "__main__" :

	# Exécution des commandes
	[exec_cmd(cmd) for cmd in cmds]