import copy as c
import math
import pandas as p
import time

p.set_option('display.max_columns', 500)
p.options.display.float_format = '{:,.3f}'.format

class Sequence:
    """ 
    Classe qui représente un objet Séquence d'acides aminées 
    """

    def __init__(self, title = "", seq = ""):
        """ Crée un objet Séquence """

        self.__seq = seq # Séquence : string
        self.__title = title #Titre de la séquence

    def get_acids(self):
        """ Renvoie les lettres de la séquence """
        return self.__seq 

    def get_acid(self, i):
        """ Renvoie l'acide aminé à la position i """
        return self.__seq[i]

    def length(self):
        """ Renvoie la taille de la séquence """
        return len(self.__seq)

    def get_title(self):
        """ Renvoie le titre de la séquence """
        return self.__title 

    def set_title(self, title):
        """ Donne un nom a la séquence """
        self.__title = title 
    
    def set_seq(self, seq):
        """ met à jour la séquence """
        self.__seq = seq 

    def display(self):
        """ Affiche la séquence au format FASTA """
        print(self.__title)
        print(self.__seq, end="\n\n")


class ParserSequence:
    """ 
    Classe qui représente un Parser qui va lire un fichier avec des séquences
    """

    def __init__(self, file_name):
        """ Crée un objet Parser """

        self.lines = open(file_name, "r").readlines()  #Lignes du fichier
        self.nb_lines = len(self.lines)
        self.sequences = [] #Liste d'objets Séquence

    def create_seq(self, seq, title):
        """ Ajoute une séquence a la liste de séquences """

        new_S = Sequence()   # Crée un nouvel objet Séquence
        new_S.set_seq(seq)   # lui donne la séquence correspondante
        new_S.set_title(title) # titre de la séquence

        self.sequences.append(new_S) #ajout d'un nouvel objet séquence

    def get_seq(self, i):
        """ Renvoie une séquence de la liste """
        return self.sequences[i]
    
    def get_all_seq(self):
        """ Renvoie toute la liste de séquences parsée """
        return self.sequences
        
    def parse(self):        
        """ Récupère les séquences du fichier """

        tmp = ""
        title = ""
        for line in self.lines: # lignes du fichier
            if line[0] == ">":
                if tmp != "":    
                    self.create_seq(tmp, title)    # on crée un objet séquence          
                    tmp = ""
                title = line
            else:
                tmp += line.strip("\n")

        self.create_seq(tmp, title)

class Matrix:
    """ 
    Classe qui représente un objet Matrice 
    """

    def __init__(self):
        """ Crée un objet Matrice """

        self.mat = []
        self.n = 0 # colonnes (seq1)
        self.m = 0 # lignes (seq2)
        self.letters_seq1 = {} # dictionnaire clé = lettre, valeur = indice dans matrice
        self.letters_seq2 = {}  #     "
        self.seq1 = Sequence() # Première séquence
        self.seq2 = Sequence() # Deuxième séquence

    def get_acid_score(self, i, j):
        """ Renvoie le score d'une cellule représentée par 2 lettres """
        if i == "-" or j == "-":
            return -1
        else:
            return self.mat[self.letters_seq1[i]][self.letters_seq2[j]]

    def get_score(self, i, j):
        """ Renvoie le score d'une cellule en mat[i][j] """
        return self.mat[i][j]

    def set_score(self, i, j, score):
        """ donne un score a une cellule """
        self.mat[i][j] = score

    def addline(self):
        """ ajoute une ligne dans la matrice """
        self.mat.append([])
        self.m += 1

    def add_cell(self, i, score):
        """ ajoute un élément dans une ligne (score) """
        if isinstance(score, str):
            score = int(score)

        self.mat[i].append(score)

    def set_lign(self, scores):
        """ ajoute toute une ligne avec des scores """
        self.mat.append(scores)

    def set_nb_col(self):
        """ détermine le nombre de colonnes """
        self.n = len(self.mat[0]) 

    def get_letters(self, i, j):
        """ Renvoie les lettres correspondantes a la cellule i,j """

        return self.seq2.get_acid(i) + self.seq1.get_acid(j)

    def setup_seq(self, seq1, seq2):
        """ remplit les dictionnaires par
            clé : acide aminée, valeur : indice dans la matrice
        """

        self.seq1.set_seq(seq1)
        self.seq2.set_seq(seq2)
        self.letters_seq1 = dict(zip(seq1, (i for i in range(len(seq1)))))
        self.letters_seq2 = dict(zip(seq2, (i for i in range(len(seq2)))))

    def panda(self, len_seq1=None, len_seq2=None):
        """ Moyen d'afficher la matrice de manière plus compacte """
        s1 = list(self.seq1.get_acids())
        s2 = list(self.seq2.get_acids())
        if len_seq1 != None:
            s1 = [i for i in range(len_seq1)]
        if len_seq2 != None:
            s2 = [j for j in range(len_seq2)]
        print(p.DataFrame(self.mat, s2, p.Index(s1, name="*")), end="\n\n")


    def get_max(self):
        """ Renvoie la position de l'élément maximal de la matrice """

        maxi = -float("inf") # - inf 
        i, j = 0, 0
        for line in range(self.m):
            current_max = max(self.mat[line])
            if current_max > maxi:
                maxi = current_max
                i = line 
                j = self.mat[line].index(current_max)

        return (i,j)

    
    def set_zero(self, positions):
        """ Met à 0 les éléments donnés en paramètre dans la matrice """

        for pos in positions:
            self.mat[pos[0]][pos[1]] = 0


class MatSubstitution(Matrix):
    """ 
    Classe qui représente un Parser qui va lire un fichier avec une matrice
    de substitution
    Cette classe hérite de Matrix
    """

    def __init__(self, file_name):
        """ Crée un objet Parser """

        Matrix.__init__(self) 
        self.lines = open(file_name, "r").readlines()  #Lignes du fichier
        self.nb_lines = len(self.lines)

    def get_mat_sub(self):
        """ Renvoie la matrice de subsitution parsée """

        return self.mat

    def parse(self):        
        """ Récupère la matrice du fichier """

        i = 0
        for line in self.lines: # lignes du fichier

            line = line.strip("\n")
            if len(line) > 0:
                if line[0] == " ":
                    self.setup_seq(line.replace(" ", ""), line.replace(" ", ""))
                       # lettres de seq 1 et 2
                else:
                    if line[0] != "#":
                        self.addline()
                        all_line = line.split(" ")
                        for l in all_line:
                            if l != "" and not l.isalpha() and l != "*":
                                self.add_cell(i, l)
                        i += 1

        self.set_nb_col()


class MatScoring(Matrix):
    """ Classe qui représente une matrice contenant les scores avec une pénalité
        de gap affine
        Cette classe hérite de Matrix
    """ 

    def __init__(self, I, E, n, m):
        Matrix.__init__(self)

        self.n = n 
        self.m = m 
        self.I = I
        self.E = E 

    def init_S_global(self):
        """ Initialise 1ere ligne et 1ere colonne de S pour la méthode globale """

        newline = [0]
        for j in range(self.n-1):
            newline.append(-self.I - j*self.E)
        self.set_lign(newline)
        for i in range(self.m-1):
            self.set_lign([-self.I - i*self.E] + (self.n-1)*[""])

    def init_S_local(self):
        """ Initialise 1ere ligne et 1ere colonne de S pour la méthode locale """

        for i in range(self.m):
            newline = []
            if i == 0:
                newline += [0]*self.n 
            else:
                newline += [0] + (self.n-1)*[""]

            self.set_lign(newline)


class MatV(Matrix):
    """ Classe qui représente une matrice permettant de sauvegarder des valeurs 
        liées aux lignes lors du calcul de la matrice de scoring 
        Cette classe hérite de Matrix
    """

    def __init__(self, n, m):
        Matrix.__init__(self)

        self.n = n 
        self.m = m 

    def init_V(self):
        """ Initialise 1ere ligne et 1ere colonne de V """

        for i in range(self.m):
            if i == 0:
                self.set_lign([-float("inf")]*self.n) # -inf -inf -inf ...
            else:
                self.set_lign([0] + (self.n-1)*[""]) # 0 ...


class MatW(Matrix):
    """ Classe qui représente une matrice permettant de sauvegarder des valeurs 
        liées aux colonnes lors du calcul de la matrice de scoring 
        Cette classe hérite de Matrix
    """

    def __init__(self, n, m):
        Matrix.__init__(self)

        self.n = n
        self.m = m 

    def init_W(self):
        """ Initialise 1ere ligne et 1ere colonne de W """

        for i in range(self.m):
            if i == 0:
                self.set_lign([-float("inf")] + [0]*(self.n-1)) # -inf 0 0 0 0 ...
            else:
                self.set_lign([-float("inf")] + (self.n-1)*[""]) # 0 ...