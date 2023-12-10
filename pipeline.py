#!/bin/python3
#-*- coding:utf-8 -*-

# Importation des modules
import os
import gzip
import pandas as pd
import pickle as pk
import subprocess as sp
import matplotlib.pyplot as plt
from argparse import ArgumentParser, MetavarTypeHelpFormatter
from typing import Dict, List, Tuple
from math import log10
from pathlib import Path
from tqdm import tqdm
from reads import Reads
from config import GlobalValue, transformer

# Gestion des erreurs pour le package "xgboost"
import warnings
warnings.filterwarnings(action="ignore", module="xgboost")
import xgboost as xgb_

parser = ArgumentParser(description="Exécute, au choix, le pipeline de mapping ou le pipeline de classification sur des runs.",
                        formatter_class=lambda prog : MetavarTypeHelpFormatter(prog=prog, max_help_position=50, width=500))
parser.add_argument("-c", "--classif", action="store_true", help="Choix entre le pipeline de mapping et le pipeline de classification")
parser.add_argument("-s", "--species_folder", type=str, default="", help="Chemin du dossier d'espèces à utiliser (défaut = '') (NOTE : Utile que pour le mapping)")
parser.add_argument("-r", "--run_folder", type=str, help="Chemin du run à analyser")
parser.add_argument("-o", "--result_folder", type=str, help="Nom du dossier qui contiendra les résultats")
parser.add_argument("-a", "--primers", type=str, default="cherif", choices=["cherif", "gozlan", "carpio", "parva"], help="Choix des primers [cherif, gozlan, carpio, parva]")
parser.add_argument("-t", "--n_threads", type=int, default=64, help="Nombre de threads pour la création de graphique fastQC (défaut = 64)")
parser.add_argument("-q", "--mapq", type=int, default=40, help="Valeur de MapQ pour garder un read, pour le mapping (défaut = 40)")
parser.add_argument("-k", "--kmers", type=int, default=28, help="Taille des k-mers pour le pipeline de mapping (défaut = 28)")
parser.add_argument("-p", "--seuil", type=float, default=0.90, help="Seuil de probabilité pour garder un read, pour la classification (défaut = 0.90)")
args = parser.parse_args()

# Initialisation d'un objet "GlobalValue"
GVALUES: GlobalValue = GlobalValue()

# Choix des primers
CHOICE: str = str(args.primers).upper()
PRIMERS: Dict[str, List[str]] = GVALUES.get_values("primers")[CHOICE]
FORWARDS: str = PRIMERS["forwards"]
REVERSES: str = PRIMERS["reverses"]

# Récupération des ranges de taille de séquences
RANGES: Dict[str, Tuple[int, int]] = GVALUES.get_values("ranges")[CHOICE]
MIN_LENGTH: int = RANGES[0]
MAX_LENGTH: int = RANGES[1]

# Récupération l'encodeur des noms d'espèces, pour la classification
ENCODER: Dict[str, int] = GVALUES.get_values("encoder")

# Récupération du nom du classifier
CLF_FILE: str = GVALUES.get_values("clf_file")

# Nom du dossier qui contiendra les FastQC et nombre de threads à utiliser
REPORT_FOLDER: str = "fastqc_reports"
N_THREADS: int = args.n_threads

# Nom du dossier qui contiendra les résultats d'analyse
RESULT_FOLDER: str = args.result_folder

# Nom du dossier de run à utiliser
RUN_FOLDER: str = args.run_folder

# Choix entre les deux pipelines
CLASSIF: bool = args.classif

# Éléments pour le pipeline de mapping
SPECIES_FOLDER: str = args.species_folder
K: int = args.kmers
MAPQ: int = args.mapq

# Éléments pour le pipeline de classification
SEUIL: float = args.seuil
with open(file=CLF_FILE, mode="rb") as xgb_object:
    CLF: xgb_.XGBClassifier = pk.load(file=xgb_object)


def make_result_folder(run_folder: str, result_folder: str) -> None:
    """
    Créer le dossier résultat, ainsi que ses sous-dossiers barcode, et également les fichiers fastq concaténés.

    Args:
    - run_folder (str): Chemin vers le dossier de run, contenant des sous-dossiers barcode, contenant des fichiers 'fastq.gz'.
    - result_folder (str): Chemin du dossier qui contiendra les résultats d'analyses.

    Returns:
    - None
    """
    # Créer le dossier des résultats
    Path(result_folder).mkdir(parents=True, exist_ok=True)

    # Itérer les éléments du dossier de run
    for barcode_folder in os.listdir(path=run_folder):
        print(f"\n{barcode_folder}")

        # Vérifier si c'est un dossier barcode
        if barcode_folder.startswith("barcode"):
            # Récupérer le chemin du dossier barcode du run
            run_barcode_folder = os.path.join(run_folder, barcode_folder)
            
            # Créer le fichier barcode résultat
            result_barcode_folder = os.path.join(result_folder, barcode_folder)
            Path(result_barcode_folder).mkdir(parents=True, exist_ok=True)
            
            # Créer le chemin du fichier fastq concaténé
            concatenated_fastq_file = os.path.join(result_barcode_folder, f"{barcode_folder}.fastq")

            # Créer et ouvrir le fichier fastq combiné en écriture
            with open(file=concatenated_fastq_file, mode="wt") as text_writer:
                # Itérer les éléments du dossier barcode du run
                for fastq_file in tqdm(os.listdir(path=run_barcode_folder), desc="Concaténation", unit=" fichiers"):
                    # Si c'est un fichier fastq.gz
                    if fastq_file.endswith(".fastq.gz"):
                        # Récupérer le chemin du fichier
                        fastq_file = os.path.join(run_barcode_folder, fastq_file)

                        # Ouvrir le fichier en lecture binaire
                        with gzip.open(filename=fastq_file, mode="rb") as binary_reader:
                            # Copier le contenu dans le fichier concaténé
                            text_writer.write(binary_reader.read().decode(encoding="utf-8"))

    print("\nOpération terminée. Les fichiers ont été concaténé avec succès.")

def run_fastqc(result_folder: str, filtered_file: bool = False, n_threads: int = 64) -> None:
    """
    Lance FastQC sur les fichiers fastq du dossier d'entrée et enregistre les résultats dans le dossier de sortie.

    Args:
    - result_folder (str): Chemin du dossier résultat, contenant les sous-dossiers barcode, contenant des fichiers fastq.
    - filtered_file (bool): Indique s'il faut traiter le fichier filtré ou pas.
    - n_threads (int): Nombre de threads à utiliser pour la commande.

    Returns:
    - None
    """
    # Itérer les éléments du dossier résultat
    for barcode_folder in os.listdir(path=result_folder):
        print(f"\n{barcode_folder}")

        # Créer le dossier report
        fastqc_folder = os.path.join(result_folder, barcode_folder, REPORT_FOLDER)
        Path(fastqc_folder).mkdir(parents=True, exist_ok=True)

        # S'il faut traiter le fichier filtré ou le fichier d'avant filtrage
        if filtered_file == False:
            fastq_file = os.path.join(result_folder, barcode_folder, f"{barcode_folder}.fastq")
        else:
            fastq_file = os.path.join(result_folder, barcode_folder, f"{barcode_folder}_filtered.fastq")

        # Construire la commande FastQC avec le chemin du fichier fastq à traiter et du dossier report
        fastqc_cmd = [f"fastqc {fastq_file} -t {n_threads} -o {fastqc_folder}"]

        # Exécutez la commande FastQC en utilisant subprocess
        sp.run(args=fastqc_cmd, shell=True)

    print("\nOpération terminée. Les fastqc ont été crée avec succès.")

def filter_reads(result_folder: str) -> None:
    """
    Filtre les reads à partir des fichiers FastQ dans le dossier résultat.

    Args:
    - result_folder (str): Chemin du dossier résultat, contenant les sous-dossiers barcode, contenant des fichiers fastq.

    Returns:
    - None
    """
    # Itérer les éléments du dossier résultat
    for barcode_folder in os.listdir(path=result_folder):
        print(f"\n{barcode_folder}")
        
        # Récupérer le chemin du fichier fastq concaténé pour le barcode actuel
        fastq_file = os.path.join(result_folder, barcode_folder, f"{barcode_folder}.fastq")

        # Créer un objet Reads pour manipuler les reads
        reads_object = Reads(filename=fastq_file)

        # Filtrer les reads selon des amorces
        reads_object.cutPrimers(forwards=FORWARDS, reverses=REVERSES, seuil_length=MIN_LENGTH)
        
        # Filtrer les reads selon leur taille
        reads_object.filter_size(inf=MIN_LENGTH, sup=MAX_LENGTH)

        # Enregistrez les lectures filtrées dans un nouveau fichier fastq
        fastq_filtered_file = os.path.join(result_folder, barcode_folder, f"{barcode_folder}_filtered.fastq")
        reads_object.write(filename=fastq_filtered_file)

        # Compressez le fichier fastq original
        gzip_cmd = [f"gzip {fastq_file}"]
        sp.run(args=gzip_cmd, shell=True)

    print("\nOpération terminée. Les reads ont été filtré avec succès.")

def mapping_specie(result_folder: str, specie_file: str, k: int = 28, mapQ: int = 40) -> pd.Series:
    """
    Effectue le mappage des lectures filtrées d'une espèce spécifique et retourne une série Pandas.

    Args:
    - result_folder (str): Chemin du dossier résultat, contenant les sous-dossiers barcode, contenant des fichiers fastq.
    - specie_file (str): Chemin vers le fichier FASTA représentant l'espèce de référence pour le mapping.
    - k (int): Longueur des k-mers utilisés par minimap2 (par défaut = 28).
    - mapQ (int): Valeur de qualité minimale du mapping (par défaut = 40).

    Returns:
    - pd.Series: Série Pandas contenant le nombre de reads mappés pour chaque barcode, en échelle logarithmique.
    """
    # Créer un DataFrame Pandas
    df = pd.DataFrame()

    # Extraire le nom de l'espèce à partir du chemin du fichier FASTA
    specie_name = specie_file.replace(".fasta", '')
    specie_name = specie_name.split('/')[-1]
    specie_name = specie_name.replace('_', ' ')

    # Créer le nom de la colonne en fonction de l'espèce et du seuil de qualité de mapping
    column = f"{specie_name} ({mapQ} <= MapQ)"
    df[column] = 0

    # Récupérer la liste des dossiers de barcode dans le dossier résultats
    barcodes = [barcode_folder for barcode_folder in os.listdir(path=result_folder) if Path(os.path.join(result_folder, barcode_folder)).is_dir()]

    # Itérer les éléments du dossier résultat
    for i, barcode_folder in tqdm(enumerate(barcodes), desc="Mapping", unit=" Barcodes"):
        # Récupérer le fichier contenant les reads à mapper
        fastq_file = os.path.join(result_folder, barcode_folder, f"{barcode_folder}_filtered.fastq")

        # Utiliser minimap2 et samtools pour effectuer le mapping et obtenir le nombre de reads mappés
        minimap2_cmd = [f"minimap2 -a -x map-ont -k {k} {specie_file} {fastq_file} | samtools view -S -F 4 -q {mapQ} | wc -l"]
        n_reads_mapped = sp.run(args=minimap2_cmd, shell=True, capture_output=True, text=True).stdout.strip()
        df.at[i, column] = int(n_reads_mapped)

    # Définir les index du DataFrame à partir des noms des barcodes
    df.index = [barcode_folder for barcode_folder in barcodes]

    # Appliquez une transformation logarithmique aux valeurs de la colonne (sauf zéro)
    return df[column].apply(lambda x: 0 if x == 0 else log10(x))

def mapping_pipeline(result_folder: str, species_folder: str, k: int = 28, mapQ: int = 40) -> pd.DataFrame:
    """
    Applique le pipeline de mapping sur plusieurs espèces et enregistre les résultats dans un fichier CSV.

    Args:
    - result_folder (str): Chemin du dossier résultat, contenant les sous-dossiers barcode, contenant des fichiers fastq.
    - species_folder (str): Chemin vers le dossier contenant les fichiers FASTA représentant les espèces de référence.
    - k (int): Longueur des k-mers utilisés par minimap2 (par défaut = 28).
    - mapQ (int): Valeur de qualité minimale du mappage (par défaut = 40).

    Returns:
    - pd.DataFrame: Dataframe dont les colonnes sont les espèces et les index sont les barcodes.
    """
    # Initialiser un DataFrame
    df = pd.DataFrame()

    # Itérer les fichiers de chaque espèce traitée
    for specie_file_name in os.listdir(path=species_folder):
        print(f"\nSpecie: {specie_file_name}")

        # Construire le chemin complet du fichier de l'espèce
        specie_file_path = os.path.join(species_folder, specie_file_name)

        # Exécuter la fonction de mapping pour une espèce et concaténer les résultats au DataFrame principal
        df = pd.concat([df, mapping_specie(result_folder=result_folder, specie_file=specie_file_path, k=k, mapQ=mapQ)], axis=1)

    # Construire le nom du fichier résultat
    result_name = "result_mapping.csv"
    result_path = os.path.join(result_folder, result_name)

    # Enregistrer le DataFrame dans un fichier CSV
    df = df.astype(dtype=float)
    df.to_csv(path_or_buf=result_path)
    
    print("\nOpération terminée. Les reads ont été mappé avec succès.")

    return df

def inv_trans(value: int) -> str:
    """
    Convertit une valeur encodée en espèce correspondante.

    Args:
    - value (int): Valeur encodée de l'espèce.

    Returns:
    - str: Nom de l'espèce.
    """
    return next((k for k, v in ENCODER.items() if v == value), None)

def get_sequences(fastq_file: str) -> List[str]:
    """
    Récupère les séquences à partir d'un fichier FastQ.

    Args:
    - fastq_file (str): Chemin vers le fichier FastQ.

    Returns:
    - List[str]: Liste des séquences.
    """
    # Ouvrir le fichier en lecture
    with open(file=fastq_file, mode="r") as reader:
        # Récupérer l'ensemble des lignes
        lines = reader.readlines()
        
        # Récupérer les séquences nucléotidiques
        return [lines[i + 1].strip() for i in range(0, len(lines), 4)]

def predict_barcode(result_folder: str, barcode_folder: str, seuil: float) -> pd.DataFrame:
    """
    Effectue une prédiction avec le classifier pour un dossier barcode.

    Args:
    - result_folder (str): Chemin du dossier résultat, contenant les sous-dossiers barcode, contenant des fichiers fastq.
    - barcode_folder (str): Nom du dossier barcode.
    - seuil (float): Seuil de probabilité pour la prédiction.

    Returns:
    - pd.DataFrame: DataFrame des occurrences d'espèces pour le code-barres spécifique.
    """    
    # Initialiser la variable de comptage pour chaque espèce
    df_species = {specie : 0 for specie in ENCODER.keys()}
    
    # Récupérer le chemin du fichier fastq à traité
    fastq_file = os.path.join(result_folder, barcode_folder, f"{barcode_folder}_filtered.fastq")
    
    # Récupérer les séquences du fichier fastq, sous forme vectorisée
    sequences = [transformer(sequence=sequence, n_features=MAX_LENGTH) for sequence in get_sequences(fastq_file=fastq_file)]
    
    # Transformer le tout en dataframe
    tmp = pd.DataFrame(data=sequences, columns=list(range(MAX_LENGTH)))

    # Si le dataframe est non vide
    if not tmp.empty:
        # Faire la prédiction avec le classifieur
        y_pred = CLF.predict(X=tmp)
        
        # Récupérer les espèces associées à chaque séquence si la probabilité d'appartenance (à l'espèce) de la séquence dépasse le seuil
        y_pred = [inv_trans(value=value) for probas, value in zip(CLF.predict_proba(X=tmp), y_pred) if seuil <= probas[value]]

        # Incrémenter les variables de comptage si y_pred est non vide
        if len(y_pred) != 0:
            for specie in y_pred:
                df_species[specie] += 1

    # Retourner le dataframe de comptage pour ce barcode 
    return pd.DataFrame({specie : [log10(count)] if count != 0 else [0] for specie, count in df_species.items()})

def classif_pipeline(result_folder: str, seuil: float = 0.90) -> pd.DataFrame:
    """
    Applique le pipeline d'analyse sur tout les barcodes et enregistre le résultat dans un fichier CSV.

    Args:
    - result_folder (str): Chemin du dossier résultat, contenant les sous-dossiers barcode, contenant des fichiers fastq.
    - seuil (float): Seuil de probabilité pour la prédiction.

    Returns:
    - pd.DataFrame: Dataframe dont les colonnes sont les espèces et les index sont les barcodes.
    """
    # Initialiser la liste des lignes du dataframe final
    df_lines = []

    # Initialiser la liste des index à utiliser pour le dataframe final
    new_indexes = []

    # Récupérer les noms des dossiers barcode
    barcodes = [barcode_folder for barcode_folder in os.listdir(path=result_folder) if Path(os.path.join(result_folder, barcode_folder)).is_dir()]

    # Itérer les barcodes
    for barcode_folder in barcodes:
        print(f"\n Predictions : {barcode_folder}")

        # Récupérer une ligne
        df_line = predict_barcode(result_folder=result_folder, barcode_folder=barcode_folder, seuil=seuil)
        
        # Ajouter la ligne à liste des lignes
        df_lines.append(df_line)
        
        # Ajouter le nom du barcode à 'new_indexes'
        new_indexes.append(barcode_folder)

    # Concaténer les lignes dans le dataframe final
    df = pd.concat(df_lines, axis=0, ignore_index=True)

    # Renommer les colonnes en fonction des spécifications
    columns = [column.replace('_', ' ') for column in df.columns]
    df.columns = [f"{column} ({seuil} <= p)" for column in columns]

    # Renommer les index selon les spécifications
    df.index = new_indexes

    # Construire le nom du fichier résultat
    result_name = "result_classif.csv"
    result_path = os.path.join(result_folder, result_name)

    # Enregistrer le DataFrame dans un fichier CSV
    df = df.astype(dtype=float)
    df.to_csv(path_or_buf=result_path)

    print("\nOpération terminée. Les reads ont été classifié avec succès.")

    return df

def plot_bar(df: pd.DataFrame, fig_name: str, result_folder: str) -> None:
    """
    Trace un diagramme à barres pour chaque espèce en utilisant des sous-graphiques.

    Args:
    - df (pd.DataFrame): DataFrame avec les espèces en tant que colonnes, les codes-barres en tant qu'index et les valeurs en log10.
    - fig_name (str): Nom du fichier pour sauvegarder la figure.
    - result_folder (str): Répertoire de résultats pour sauvegarder la figure.

    Returns:
    - None
    """
    # Obtenir la liste des espèces
    species_list = df.columns.tolist()

    # Définir le nombre de lignes et de colonnes pour les sous-graphiques
    n_rows = len(species_list) // 3 + len(species_list) % 3
    n_cols = 3

    # Créer les sous-graphiques
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(30, 10 * n_rows))

    # Ajuster les espacements vertical et horizontal
    fig.subplots_adjust(hspace=0.25, wspace=0.25)

    # Aplatir les axes s'il y a plus d'une ligne
    if n_rows > 1:
        axes = axes.flatten()

    # Itérer sur les espèces et créer un sous-graphique pour chacune
    for i, species in enumerate(iterable=species_list, start=0):
        # Sélectionner le sous-graphique actuel
        ax = axes[i] if n_rows > 1 else axes

        # Tracer le diagramme à barres empilées pour l'espèce actuelle
        df_species = df[[species]]
        df_species.plot(kind="bar", ax=ax, legend=None, color="skyblue", edgecolor="black")

        # Définir les labels et le titre du sous-graphique
        ax.set_xlabel("Barcodes")
        ax.set_ylabel("Nombre de reads (en log10)")
        ax.set_title(species)

        # Ajouter une ligne en pointillés à y = 1.0
        ax.axhline(y=1.0, color="black", linestyle="dashed", linewidth=1)

        # Ajouter une grille
        ax.grid(True, axis='y', linestyle='-', alpha=0.7)

    # Ajuster la disposition et afficher le graphique
    plt.savefig(os.path.join(result_folder, fig_name))
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

if __name__ == "__main__":
    # Crée le dossier de résultats, et crée les fichiers concaténés
    make_result_folder(run_folder=RUN_FOLDER, result_folder=RESULT_FOLDER)

    # Exécute FastQC sur les fichiers concaténés avant le filtrage
    run_fastqc(result_folder=RESULT_FOLDER, filtered_file=False, n_threads=N_THREADS)

    # Filtre les reads en fonction des critères définis dans la fonction
    filter_reads(result_folder=RESULT_FOLDER)

    # Exécute FastQC sur les fichiers concaténés après le filtrage
    run_fastqc(result_folder=RESULT_FOLDER, filtered_file=True, n_threads=N_THREADS)

    # Vérifie si le pipeline de classification est choisie
    if CLASSIF == True:
        # Exécute le pipeline de classification avec le seuil spécifié
        df = classif_pipeline(result_folder=RESULT_FOLDER, seuil=SEUIL)
        # Définit le nom du fichier de la figure pour le pipeline de classification
        fig_name = "classif_fig.png"
    else:
        # Exécute le pipeline de mapping avec les paramètres spécifiés
        df = mapping_pipeline(result_folder=RESULT_FOLDER, species_folder=SPECIES_FOLDER, k=K, mapQ=MAPQ)
        # Définit le nom du fichier de la figure pour le pipeline de mapping
        fig_name = "mapping_fig.png"

    # Génère le diagramme à barres en utilisant les résultats du pipeline
    plot_bar(df=df, fig_name=fig_name, result_folder=RESULT_FOLDER)