# Project eDNA

Ce dossier contient les éléments du projet eDNA. Cela comprend la bibliographie utilisée, la description des amorces de PCR, les propriétés du matériel, les logiciels de traitement, les runs sous format fastq compressé (*.gz), les séquences ADN des espèces d'intérêts sous format fasta, les scripts/pipelines en langage python3, et les résultats obtenus.

Pour la reproduction des résultats. Copier dans un dossier (A) tous les dossiers présents dans le dossier 03_runs/, ainsi que ceux du dossier 04_species/.

Une fois les éléments précédents copiés, ajouter (dans A) les éléments du dossier 07_pipelines/.

Ouvrir un terminal de commandes dans A, et lancer le script "process_run.py" via la commande : 

```bash
    python3 process_runs.py
```