# Agent_rosette_detection

Le script `Pipeline.py` a été développé pour la détection d'espèces : Sphaerothecum destruens (Agent Rosette) ; Cyprinus carpio ; et Pseudorasbora parva. Le pipeline utilise le classifieur "clf_RA", ainsi que le module "reads", ainsi que les outils "fastqc", "samtools" et "minimap2". Voici des informations sur la configuration et l'utilisation du script :

## Téléchargement du dépôt

Pour récupérer le dépôt "Agent_rosette_detection" depuis la ligne de commande, vous pouvez utiliser la commande `git clone`. Exécutez la commande suivante dans votre terminal :

```bash
git clone https://github.com/RAVAO-Ravo/Agent_rosette_detection.git
```

Assurez-vous d'avoir Git installé sur votre machine pour utiliser cette commande. Une fois la commande exécutée, le dépôt sera cloné dans un nouveau dossier appelé "Agent_rosette_detection" dans votre répertoire actuel.

## Dépendances

### Dépendances Python

Pour installer les dépendances Python nécessaires au projet, utilisez le fichier `requirements.txt`. Exécutez la commande suivante dans l'emplacement du fichier `requirements.txt` :

```bash
pip3 install -r requirements.txt
```

Assurez-vous d'avoir `Python3` et `pip3` installés sur votre machine. Cette commande installera les bibliothèques Python requises pour le bon fonctionnement du projet, en se basant sur les spécifications du fichier `requirements.txt`.

### Autres dépendances

Pour les autres outils, vous pouvez exécuter ces commandes dans un terminal de commandes :

```bash
sudo apt-get install samtools=1.13
```

```bash
sudo apt-get install fastqc=0.11.9
```

```bash
sudo apt-get install minimap2=2.24-r1122
```

## Utilisation

Le script prend en charge deux pipelines : le pipeline de mapping et le pipeline de classification. Vous pouvez choisir entre ces deux options en utilisant l'argument `-c` ou `--classif`.

### Options Principales

- **`-c` ou `--classif`**
  - Choix entre le pipeline de mapping et le pipeline de classification.

- **`-s [SPECIES_FOLDER]` ou `--species_folder [SPECIES_FOLDER]`**
  - Chemin du dossier d'espèces à utiliser (utile uniquement pour le pipeline de mapping).

- **`-r [RUN_FOLDER]` ou `--run_folder [RUN_FOLDER]`**
  - Chemin du run à analyser.

- **`-o [RESULT_FOLDER]` ou `--result_folder [RESULT_FOLDER]`**
  - Nom du dossier qui contiendra les résultats.

- **`-a [PRIMERS]` ou `--primers [PRIMERS]`**
  - Choix des primers parmi [cherif, gozlan, carpio, parva]. (par défaut, "cherif") (NOTE : ne pas utilisé si `-c` est spécifié, laissé sur la valeur par défaut)

- **`-t [N_THREADS]` ou `--n_threads [N_THREADS]`**
  - Nombre de threads pour la création de graphiques FastQC (par défaut, 64).

- **`-q [MAPQ]` ou `--mapq [MAPQ]`**
  - Valeur de MapQ pour garder un read lors du mapping (par défaut, 40).

- **`-k [KMERS]` ou `--kmers [KMERS]`**
  - Taille des k-mers pour le pipeline de mapping (par défaut, 28).

- **`-p [SEUIL]` ou `--seuil [SEUIL]`**
  - Seuil de probabilité pour garder un read lors de la classification (par défaut, 0.90).

### Exemples d'utilisation

Pour l'utilisation du classifieur :

```bash
python3 Pipeline.py -c -r /chemin/vers/run -o resultats_classif -p 0.95
```
Pour l'utilisation du mapping :

```bash
python3 Pipeline.py -r /chemin/vers/run -s /chemin/vers/species -o resultats_mapping -k 28 -q 60 -t 64 -a cherif
```

## Licence

Ce projet est sous licence [Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)](https://creativecommons.org/licenses/by-sa/4.0/).

Vous êtes libre de :

- **Partager** : copier et redistribuer le matériel sous n'importe quel format ou médium.
- **Adapter** : remixer, transformer et construire à partir du matériel.

Selon les conditions suivantes :

- **Attribution** : Vous devez donner le crédit approprié, fournir un lien vers la licence et indiquer si des modifications ont été faites. Vous devez le faire d'une manière raisonnable, mais d'une manière qui n'implique pas que l'auteur vous approuve ou approuve votre utilisation du matériel.
- **ShareAlike** : Si vous remixez, transformez ou construisez à partir du matériel, vous devez distribuer vos contributions sous la même licence que l'original.

Veuillez consulter le fichier [LICENSE](LICENSE) pour plus de détails.
