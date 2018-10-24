import copy as c
import math
import pandas as p
import time

p.set_option('display.max_columns', 500)
p.options.display.float_format = '{:,.0f}'.format


class Sequence:
    """ 
    Classe qui représente un objet Séquence d'acides aminées 
    """

    def __init__(self, title = "", seq = ""):
        """ Crée un objet Séquence """

        if title != "" and seq != "":
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

    def set_letters_seq(self, seq1, seq2):
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

        # if len_seq != None:
        #     print(p.DataFrame(self.mat, [i for i in range(len_seq)], \
        #     p.Index([j for j in range(len_seq)], name="*")), end="\n\n")

        # else:
        #     print(p.DataFrame(self.mat, list(self.seq2.get_acids()), \
        #     p.Index(list(self.seq1.get_acids()), name="*")), end="\n\n")

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
                    self.set_letters_seq(line.replace(" ", ""), line.replace(" ", ""))
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


class Alignment:
    """ Classe représentant un objet qui trouvera l'alignement de 2 séquences
        d'acides aminées
    """

    def __init__(self, I, E, mat_file, seq1, seq2, p):
        
        self.I = I 
        self.E = E
        self.p = p
        self.seq1 = seq1
        self.seq2 = seq2
        self.n = self.seq1.length()+1
        self.m = self.seq2.length()+1

        if self.p == 1 or p == 2: # Si on veut imprimer les informations
            print("Séquence 1 de longueur {0}: ".format(self.n))
            self.seq1.display() 
            print("Séquence 2 de longueur {0}: ".format(self.m))
            self.seq2.display()
            print("matrice de substitution utilisée: {0}".format(mat_file))
            print("Pénalité de gap affine: I = {0} | E = {1}".format(self.I, self.E))

        # =================== SUBSITUTION ==============================
        self.t = MatSubstitution(mat_file)
        self.t.parse()

        # =================== SCORING ===============================

        self.S = MatScoring(I, E, self.n, self.m) # Création de la matrice scoring

        # ===================== V ET W ===================================

        self.V = MatV(self.n, self.m)
        self.W = MatW(self.n, self.m)

        self.V.init_V() # Initialise V
        self.V.set_letters_seq("-"+self.seq1.get_acids(), "-"+self.seq2.get_acids())
        self.W.init_W() # Initialise W
        self.W.set_letters_seq("-"+self.seq1.get_acids(), "-"+self.seq2.get_acids())

        self.current_sol = [] # Pour le backtracking
        self.all_solutions = []

    def get_v(self, i, j):
        return max(self.S.get_score(i-1, j)-self.I,self.V.get_score(i-1,j)-self.E)

    def get_w(self, i, j):
        return max(self.S.get_score(i, j-1)-self.I,self.W.get_score(i, j-1)-self.E)

    def is_previous(self, mat, i, j):
        """ Regarde si l'élement ij résulte de l'élément en diagonale,
            en haut ou a gauche 
        """

        res = False
        if mat == "v" and self.V.get_score(i,j) == self.S.get_score(i,j):
                res = True 

        elif mat == "w" and self.W.get_score(i,j) == self.S.get_score(i,j):
                res = True 

        elif mat == "s":
            letters_ij = self.S.get_letters(i,j) # 'AB' par exemple
            t_ij = self.t.get_acid_score(letters_ij[0], letters_ij[1])  # t(i,j)

            if self.S.get_score(i-1, j-1) + t_ij == self.S.get_score(i,j):
                res = True 

        return res 

    def add_sol(self, i, j, pos):
        """ Ajoute une solution à la liste des solutions
            De la forme: ('AB', score)
        """

        letters_ij = self.S.get_letters(i,j) 
        if pos == 0:
            letters_ij = letters_ij[0]+"-" # A- par exemple => Si gap vient de V
        elif pos == 1:
            letters_ij = "-"+letters_ij[1] # -A par exemple => Si gap vient de W
        self.current_sol.append((letters_ij, \
                            self.t.get_acid_score(letters_ij[0], letters_ij[1])))


class GlobalAlignment(Alignment):
    """
    Classe qui va trouver un aligement de séquences d'acides aminées avec
    une pénalité affine et une méthode globale
    """

    def __init__(self, k, I, E, mat_file, seq1, seq2, p):
        Alignment.__init__(self, I, E, mat_file, seq1, seq2, p)
        self.k = k

        self.S.init_S_global()  # On initialise la matrice de scoring
        self.S.set_letters_seq("-"+self.seq1.get_acids(), "-"+self.seq2.get_acids())

    def get_s_global(self, i, j, v_ij, w_ij):
        letters_ij = self.S.get_letters(i,j) # 'AB' par exemple
        t_ij = self.t.get_acid_score(letters_ij[0], letters_ij[1])

        s_ij = max( self.S.get_score(i-1,j-1) + t_ij, v_ij, w_ij )
        # self.S.set_score(i,j, s_ij)

        return s_ij

    def Needleman_Wunsch(self):
        """ 
        Algorithme qui calcule la matrice de scoring pour l’alignement global 
        en utilisant la pénalité affine puis fait un backtracking pour récupérer
        tous les alignements optimaux possibles
        """

        # ==================== CREATION SCORING ===========================

        for i in range(1, self.m):
            for j in range(1, self.n):
                v_ij = self.get_v(i, j)   # V(i,j)
                self.V.set_score(i,j, v_ij)
                w_ij = self.get_w(i, j)
                self.W.set_score(i,j, w_ij) # W(i,j)

                # self.update_s_global(i, j, v_ij, w_ij)
                s_ij = self.get_s_global(i, j, v_ij, w_ij)
                self.S.set_score(i, j, s_ij) # S(i,j)

        if self.p == 2: # Si on veut afficher les matrices
            self.V.panda()
            self.W.panda()
            self.S.panda()

        self.backtracking_global(self.m-1, self.n-1) #appel sur element en derniere
                                                     #ligne, derniere colonne

        R = Result(self.all_solutions, self.p)
        return R.bind()

    def backtracking_global(self, i, j):
        """ Remonte la matrice de scoring a partir du dernier élément jusqu'à [0][0]
            pour avoir les k alignements
        """
        #print("current: ", self.current_sol)

        if i == 0 or j == 0:
            if not (i==0 and j==0):
                # self.current_sol.append((i,j))
                self.add_sol(i,j,2)
            if len(self.all_solutions) == self.k: # Si on a deja trouvé k alignements
                return 
            if self.current_sol not in self.all_solutions:
                #print("1 solution trouvé: ", self.current_sol)
                self.all_solutions.append(c.deepcopy(self.current_sol))           

        else:
            for pos in range(3):

                new_i = i 
                new_j = j 
                valid = False 
                if pos == 0 and self.is_previous("v", i, j):  # haut
                    new_i -= 1 # i-1
                    valid = True 
                elif pos == 1 and self.is_previous("w", i, j): # gauche
                    new_j -= 1 # j-1  
                    valid = True 
                elif pos == 2 and self.is_previous("s", i, j): # diagonale
                    new_i -= 1 # i - 1
                    new_j -= 1  # j - 1
                    valid = True 

                if valid: 
                    # self.current_sol.append((i, j))
                    self.add_sol(i, j, pos)
                    self.backtracking_global(new_i, new_j) # appel sur cellule suivante
                    # if len(self.current_sol) > 1:
                    self.current_sol.pop() # destruction sol partielle


class LocalAlignment(Alignment):
    """
    Classe qui qui va trouver un aligement de séquences d'acides aminées avec
    une pénalité affine et une méthode locale
    """

    def __init__(self, l, I, E, mat_file, seq1, seq2, p):
        Alignment.__init__(self, I, E, mat_file, seq1, seq2, p)
        self.l = l

        self.S.init_S_local()  # On initialise la matrice de scoring
        self.S.set_letters_seq("-"+self.seq1.get_acids(), "-"+self.seq2.get_acids())

        self.zeros = []
        self.found = False

    def get_s_local(self, i, j, v_ij, w_ij):
        """ détermine la valeur de S en fonction de V et W et de t"""

        letters_ij = self.S.get_letters(i,j) # 'AB' par exemple
        t_ij = self.t.get_acid_score(letters_ij[0], letters_ij[1])
        s_ij = max( self.S.get_score(i-1,j-1) + t_ij, v_ij, w_ij, 0 )
        # self.S.set_score(i,j, s_ij)    
    
        return s_ij

    def compute_scoring(self, start_i, start_j):
        """ ReCalcule la matrice de
        scoring pour les alignements locaux  avec pénalité affine 
        après avoir trouvé un alignement local
        """

        # ==================== CREATION SCORING ===========================

        for i in range(start_i, self.m):
            for j in range(start_j, self.n):
                if (i,j) not in self.zeros:
                    v_ij = self.get_v(i, j)     # V(i,j)
                    self.V.set_score(i,j, v_ij)
                    w_ij = self.get_w(i, j)     # W(i,j)
                    self.W.set_score(i,j, w_ij)

                    # self.update_s_local(i, j, v_ij, w_ij)
                    s_ij = self.get_s_local(i, j, v_ij, w_ij)
                    self.S.set_score(i, j, s_ij)

        if self.p == 2:
            print("la matrice de Scoring: ")
            self.S.panda()

    def is_modified(self, i, j):
        res = False 
        if (i,j) in self.zeros:
            res = True
        else: 
            v_ij = self.get_v(i,j)
            w_ij = self.get_w(i,j)
            s_ij = self.get_s_local(i,j, v_ij, w_ij)
            # print("score de {0}, {1} : {2}".format(i, j, s_ij))
            # if v_ij != self.V.get_score(i,j) and \
            #     w_ij != self.W.get_score(i,j) and \
            if s_ij != self.S.get_score(i,j):

                self.V.set_score(i,j,v_ij) # si c'est modifié, on modifie
                self.W.set_score(i,j,w_ij)
                self.S.set_score(i,j,s_ij)
                res = True 
        # print(" {0} {1} : changé ? {2}".format(i, j, res))
        return res 

    def recompute_scoring(self, start_i, start_j):
        """ test de recalcul de la matrice de scoring"""
        
        # start_i = start_i 
        # start_j = start_j 
        i = start_i
        j = start_j 
        while i != self.m and j != self.n:
            # print("start ij : ({0},{1})".format(i, j))
            while i != self.m and j != self.n and self.is_modified(i, j):
                # print("ij 1 : ({0},{1})".format(i, j))
                j += 1 # colonne suivante 
            j = start_j 
            i += 1
            while i != self.m and j != self.n and self.is_modified(i, j):
                # print("ij 2: ({0},{1})".format(i, j))
                i += 1 # ligne suivante 
            start_i += 1
            start_j += 1
            i = start_i # sous carré suivant
            j = start_j
        if self.p == 2:
            # print("la matrice de Scoring: ")
            self.S.panda()

    
    def sol_found(self, i,j):
        """ Détermine si on a fini un alignement local """

        return (i == 0 or j == 0 or self.S.get_score(i,j) == 0)

    def bottom_up(self, i,j):
        """ Remonte la matrice de scoring a partir du max de la matrice
            de scoring jusqu'à un élément de la 1ere ligne ou 1ere colonne
            ou un élément 0
        """
        while not self.sol_found(i,j):
            prev_i = i
            prev_j = j
            pos = 2
            if self.is_previous("v", i, j):    # haut
                i-=1
                pos = 0
            elif self.is_previous("w", i, j):  # gauche
                j-=1
                pos = 1
            elif self.is_previous("s", i, j):  # diagonale
                i-=1
                j-=1
            if self.S.get_score(i,j) != 0:
                self.add_sol(prev_i,prev_j,pos)
                self.zeros.append((prev_i,prev_j))

        self.all_solutions.append(self.current_sol)
        self.add_sol(prev_i, prev_j, 2)
        self.zeros.append((prev_i, prev_j))

    
    def Smith_Waterman(self):
        self.compute_scoring(1,1)

        for i in range(self.l):
            current_max = self.S.get_max()
            self.bottom_up(current_max[0], current_max[1])
            if i == self.l-1: # Si on a trouvé l solutions
                break
            # une fois le backtrack fini, on met a 0 l'alignement 
            sols = len(self.current_sol)
            self.S.set_zero(self.zeros[len(self.zeros)-sols:])
            self.V.set_zero(self.zeros[len(self.zeros)-sols:])
            self.W.set_zero(self.zeros[len(self.zeros)-sols:])

            self.current_sol = []
            if self.p == 2:
                print("mise a 0: ")
                self.S.panda()

                print("recalcul de la matrice: ")
            # self.compute_scoring(self.zeros[-1][0], self.zeros[-1][1])
            self.recompute_scoring(self.zeros[-1][0], self.zeros[-1][1])
        R = Result(self.all_solutions, self.p)
        return R.bind()


class Result:
    """ Classe représentant un objet Résultat dans laquelle on va créer le
        résultat de l'alignement global ou local
    """

    def __init__(self, all_sol, p):
        self.all_sol = all_sol 
        self.p = p 
        # self.alignments = [] # contiendra les alignements # ex: 'AB-:. CDE'
        self.scores = [] # score de similarité de tous les alignements

    def bind(self):
        """ Crée les liaisons entre 2 séquences + leurs scores de similarité """

        for sol in self.all_sol:
            seq1, seq2, links = "", "", ""
            for i in range(len(sol)-1, -1, -1):
                seq1 += sol[i][0][1]  # 2eme lettre
                seq2 += sol[i][0][0]  # 1ere lettre

                if sol[i][0][1] == "-" or sol[i][0][0] == "-": # SI gap
                    links += " "   # Pas de correspondance
                else:
                    if sol[i][0][1] == sol[i][0][0]:
                        links += ":"       # Identiques
                    elif sol[i][1] >= 0:   # Si le score est positif => similaires
                        links += "."
                    else:
                        links += " "

            # pourcentage de similarité
            similarity = ((links.count(".")+links.count(":"))/len(links)) * 100
            identity = (links.count(":")/len(links)) * 100 # pourcentage d'identité
            self.scores.append(similarity)
            
            if self.p == 1 or self.p == 2: # Si on veut afficher les résultats
                self.print_result(seq1, seq2, links, similarity, identity)

        return self.scores

    def print_result(self, seq1, seq2, links, similarity, identity):
        """ Imprime de manière jolie les 2 séquences alignées 
            seq est de la forme: [('AB', scoreAB), .. ]
        """

        print("\n>>> Alignement: ")

        i, j = 0, 0
        while i < (len(seq1) // 60):
            print(seq1[j:j+60]+"\n"+links[j:j+60]+"\n"+seq2[j:j+60]+"\n\n")
            i += 1
            j += 60
        end = len(seq1) - len(seq2)%60
        print(seq1[end:]+"\n"+links[end:]+"\n"+seq2[end:]+"\n\n")

        print("==> similarité: {0} %".format(similarity))
        print("==> Identité: {0} %".format(identity))

def get_best_globalalignments(S, I, E):

    # sims = Matrix()
    best_score = [0, 0, 0] # les 2 derniers 0 accueilleront les numéros de séquences dont
                             # le score est maximal
    
    for i in range(len(S)):
        # sims.addline()
        for j in range(len(S)):
             # on ne compare pas une séquence avec elle-même
            if i != j:
                #print(i, "  ", j)
                G_tmp = GlobalAlignment(1, I, E, "blosum62.txt", S[i], S[j], 0)
                tmp_score = G_tmp.Needleman_Wunsch()
                # sims.add_cell(i, round(tmp_score[0], 2))
                if tmp_score[0] > best_score[0]:
                    best_score[0] = tmp_score[0]
                    best_score[1] = i
                    best_score[2] = j
            # else:
            #     sims.add_cell(i, 100)
    # sims.set_nb_col()
    # sims.panda(len(S))
    
    return best_score

def is_very_similar(S, I, E, seq):
    """
    Regarde si une séquence est similaire a toutes les autres séquences
    """

    for s in S:
        if s.get_acids() == seq.get_acids():
            return True 
        G_tmp = GlobalAlignment(1, I, E, "blosum62.txt", s, seq, 0)
        tmp_score = G_tmp.Needleman_Wunsch()
        if tmp_score[0] > 60:
            return True 
    return False

def remove_high_scores(S, I, E):
    # S = liste de séquences
    reduced = [S[0]]
    for i in range(1, len(S)):
        print("i = ", i)
        if not is_very_similar(reduced, I, E, S[i]):
            print("ajouté")
            reduced.append(S[i])
        print("len reduced: ", len(reduced))

    reduced_file = open("to-be-aligned-reduced.fasta", "w") 
    for s in reduced:
        reduced_file.write(s.get_title())
        reduced_file.write(s.get_acids()+"\n")
    reduced_file.close()

class Profile:
    """ Classe représentant un PSSM, un profil, qui représente de manière
        compacte un alignement de séquence multiple
    """

    def __init__(self, MSA):
        self.prof = Matrix()
        self.MSA = MSA # Liste des séquences alignées de manière multiple
        self.acids = "RHKDESTNQCGPAILMFWYV"
        self.prof.set_letters_seq("", self.acids)
        self.Nseq = self.MSA[0].length() # longueur des séquences
        self.Nacids = 20 # nombre d'acides aminées

        # pseudocounts
        self.alpha = self.Nseq - 1  # facteur de cadrage pour les données observées
        self.beta = math.sqrt(self.Nseq)  # facteur de cadrage
        # Ala (A) 8.25   Gln (Q) 3.93   Leu (L) 9.65   Ser (S) 6.62
        # Arg (R) 5.53   Glu (E) 6.73   Lys (K) 5.81   Thr (T) 5.35
        # Asn (N) 4.05   Gly (G) 7.07   Met (M) 2.41   Trp (W) 1.09
        # Asp (D) 5.46   His (H) 2.27   Phe (F) 3.86   Tyr (Y) 2.91
        # Cys (C) 1.38   Ile (I) 5.92   Pro (P) 4.73   Val (V) 6.86

        self.pa = {"R": 5.53, "H": 2.27,"K": 5.81,"D":5.46,"E":6.73,"S":6.62,\
        "T":5.35,"N":4.05,"Q":3.93,"C":1.38,"G":7.07,"P":4.73,"A":8.25,"I":5.92,\
        "L":9.65,"M":2.41,"F":3.86,"W":1.09,"Y":2.91,"V":6.86} # probabilité d'apparition des acides aminées

    def number_acids(self, u, a):
        """
        Calcule le nombre d'acide a dans la colonne u
        """
        Nua = 0 
        for i in range(self.Nacids):
            if self.MSA[i].get_acid(u) == a: # acide aminé a la séquence i, col u
                Nua += 1
        return Nua
    
    def frequency(self, u, b):
        """ 
        Calcule la fréquence d'apparition d'un acide aminé b dans la colonne b
        parmi toutes les séquences de l'alignement multiple
        """

        Nub = 0 # nombre d'acide b dans la colonne u
        for i in range(self.Nacids):
            if self.MSA[i].get_acid(u) == b: # acide aminé a la séquence i, col u
                Nub += 1
        
        return Nub / self.Nseq

    def pseudocount_q(self, u, a):
        """
        Calcule le pseudocount associé à l'acide aminé a
        """
        return (self.number_acids(u, a) + self.beta*(self.pa[a]/100)) \
                / (self.Nseq + self.beta )
        
    def compute_profile(self):
        print("length seq: ", self.Nseq)
        # self.prof.panda(self.Nseq)  # on affiche le profil
        self.MSA[0].display()         # on affiche la première séquence

        for a in range(self.Nacids):   # lignes
            for u in range(self.Nseq): # colonnes
                Qua = self.pseudocount_q(u, self.acids[a])
                self.prof[a][u] = math.log(Qua/(self.pa[a]/100))

        

def main():

    # ================================= GLOBAL ==============================
    P = ParserSequence("BRD-sequence.fasta")
    P.parse()

    P1 = ParserSequence("to-be-aligned.fasta") 
    P1.parse()

    # print("========================= TEST BEST SCORE to be aligned ========\n\n")
    # best= get_best_globalalignments(P1.get_all_seq(), 4, 4)
    # print("le meilleur score est celui-ci: ", best[0])

    # reduce
    # begin = time.time()
    # remove_high_scores(P1.get_all_seq(), 4, 4)
    # print("temps écoulé: ", time.time() - begin)

    # test si le reduce a bien max 60%
    # begin = time.time()
    # Ptest = ParserSequence("to-be-aligned-reduced.fasta")
    # Ptest.parse()
    # best= get_best_globalalignments(Ptest.get_all_seq(), 4, 4)
    # print("le meilleur score est celui-ci: ", best[0])
    # print("temps écoulé: ", time.time() - begin)

    # print(" ==============================================TEST AZAP et AI: \n\n")
    # seq1 = Sequence("Séquence 1 de l'exemple","AZAP")
    # seq2 = Sequence("Séquence 2 de l'exemple", "AI")
    # G0 = GlobalAlignment(3, 12, 2, "blosum62.txt", seq1, seq2, 2)
    # G0.Needleman_Wunsch()

    # print(" =======================================TEST BEST SCORE: \n\n")
    # best= get_best_globalalignments(P.get_all_seq(), 12, 2, None)
    # print("le meilleur score est celui-ci: ", best[0])
    # G1 = GlobalAlignment(1, 12, 2, "blosum62.txt", P.get_seq(best[1]), P.get_seq(best[2]), 1)
    # G1 = G1.Needleman_Wunsch()

    # print(" =======================================TEST MGGETFA et GGVTTF: \n\n")
    # seq1 = Sequence("Séquence 1 de l'exemple","MGGETFA")
    # seq2 = Sequence("Séquence 2 de l'exemple", "GGVTTF")
    # G0 = GlobalAlignment(3, 12, 2, "blosum62.txt", seq1, seq2, 2)
    # G0.Needleman_Wunsch()

    # print(" ==============================================TEST seq0 et seq0: \n\n")
    # G2 = GlobalAlignment(1, 12, 2, "blosum62.txt", P.get_seq(0), P.get_seq(0), 1)
    # G2.Needleman_Wunsch()

    # print(" ==============================================TEST seq0 et seq1: \n\n")
    # G3 = GlobalAlignment(1, 12, 2, "blosum62.txt", P.get_seq(0), P.get_seq(1), 1)
    # G3.Needleman_Wunsch()

    # print(" ==============================================TEST seq0 et seq2: \n\n")
    # G4 = GlobalAlignment(1, 12, 2, "blosum62.txt", P.get_seq(0), P.get_seq(2), 1)
    # G4.Needleman_Wunsch()

    # # =============================== LOCAL ===============================
    # P2 = ParserSequence("protein-sequences.fasta")
    # P2.parse()

    # seq1 = Sequence("Séquence 1 de l'exemple","ISALIGNED")
    # seq2 = Sequence("Séquence 2 de l'exemple", "THISLINE")
    # L0 = LocalAlignment(2, 4, 4, "blosum62.txt", seq1, seq2, 2)
    # L0.Smith_Waterman()

    # L1 = LocalAlignment(1, 12, 2, "blosum62.txt", P2.get_seq(0), P2.get_seq(0), 1)
    # L1.Smith_Waterman()

    # L3 = LocalAlignment(3, 12, 2, "blosum62.txt", P2.get_seq(0), P2.get_seq(1), 1)
    # L3.Smith_Waterman()


    # ================================== PROFIL =============================

    # P3 = ParserSequence("msaresults-MUSCLE.fasta") 
    # P3.parse() # alignements multiples de séquences 

    # Prof = Profile(P3.get_all_seq())
    # Prof.compute_profile()

if __name__ == "__main__":
    main()