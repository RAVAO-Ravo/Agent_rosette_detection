# Test du classifieur

Ce dossier est dédié à l'analyse du comportement du classifieur "clf_RA" pour la classification.

## Structure du Projet

1. **`datasets/`**
   - Dossier contenant les données de barcodes (du 29 février 2023) vectorisées (barcodes 01 à 12).

2. **`results/`**
   - Dossier qui contiendra des sous-dossiers, sous-dossiers contenant des résultats d'analyses du comportement du classifieur (1 sous-dossier par barcode).

3. **`make_datas.py`**
   - Fichier permettant la vectorisation des données de barcodes issues du dossier `../01_clf_ra/fastqs`.

4. **`test_clf.ipynb`**
   - Notebook permettant d'inférer les résultats d'analyse du comportement du classifieur.

## Déroulement de l'Analyse

1. **Vectorisation des Données de Barcodes :**
   - Utilisation du fichier `make_datas.py` pour effectuer la vectorisation des données de barcodes provenant du dossier `../01_clf_ra/fastqs`.

2. **Inférence des Résultats d'Analyse :**
   - Utilisation du notebook `test_clf.ipynb` pour inférer les résultats d'analyse du comportement du classifieur pour la prédiction de chaque barcode. Les résultats sont enregistrés dans le dossier `results`.

## Prérequis

Assurez-vous d'inclure les fichiers `.so` issus de la compilation des `.pyx` du dossier `00_reads` pour garantir le bon fonctionnement du projet. De plus, veillez à avoir le fichier du classifieur "clf_RA.sav".