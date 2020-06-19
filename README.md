Développement (from scratch) d’un pipeline d’analyse de données de séquençage whole exome (WES).

Le but est d’automatiser un pipeline d’analyse de données de Séquençage Whole Exome (WES) en utilisant python. Le pipeline doit être exécutable dans un environnement de calcul haute performance (HPC) utilisant Slurm comme gestionnaire de tâches. Il devra lancer des scripts bash et gérer les dépendances et les entrées/sorties des outils tiers utilisés.

Les principales étapes de l’analyse sont :
• QC (fastqc).
• Le trimming (cutadapt).
• L’alignement (bwa, bowtie).
• Le traitement des fichiers bams (Picard, Samtools).
• Variant calling (GATK,Varscan).
• L’annotations de variants (Annovar).
• CNV calling (CNVkit).
