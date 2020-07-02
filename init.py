"""
melek soumia
khelil yasmine
script python pour le lancment d'analyses d'un pipeline WES
gestion des fichiers d'entrée et de sortie
"""

import logging, os, threading, subprocess


class Run:
    def __init__(self):
        self.__process = None

