import copy as c
import pandas

pandas.set_option('display.max_columns', 500)
pandas.options.display.float_format = '{:,.0f}'.format


class Sequence:
    """ 
    Classe qui représente un objet Séquence d'acides aminées 
    """

    def __init__(self, title = None, seq = None):
        """ Crée un objet Séquence """

        if title != None and seq != None:
            self.__seq = seq # Séquence : string
            self.__title = title #Titre de la séquence

    def get_acids(self):
        """ Renvoie les lettres de la séquence """
        return self.__seq 

    def length(self):
        """ Renvoie la taille de la séquence """
        return len(self.__seq)

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
        
    def parse(self):        
        """ Récupère les séquences du fichier """

        tmp = ""
        title = ""
        for line in self.lines: # lignes du fichier
            if line[0] == ">":
                # print(tmp + "\n")
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
        self.letters_orders1 = "" #Les lettres des colonnes mis dans l'ordre 
        self.letters_orders2 = "" #Les lettres des lignes mis dans l'ordre


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

    def inc_nb_col(self):
        """ augmente le nb de colonnes de 1 """
        self.n += 1

    def set_nb_col(self):
        """ détermine le nombre de colonnes """
        self.n = len(self.mat[0])

    def inc_nb_lign(self):
        """ augmente le nb de lignes de 1 """
        self.m += 1  

    def get_letters(self, i, j):
        """ Renvoie les lettres correspondantes a la cellule i,j """
        return self.letters_orders2[i] + self.letters_orders1[j]

    def set_letters_seq(self, seq1, seq2):
        """ remplit les dictionnaires par
            clé : acide aminée, valeur : indice dans la matrice
        """

        seq1 = seq1.replace(" ", "")
        seq2 = seq2.replace(" ", "")
        self.letters_orders1 = seq1
        self.letters_orders2 = seq2

        self.letters_seq1 = dict(zip(seq1, (i for i in range(len(seq1)))))
        self.letters_seq2 = dict(zip(seq2, (i for i in range(len(seq2)))))

    def display(self):
        """ Affiche la matrice """

        for l in self.letters_orders1:
            if l == "-":
                print("\t-", end="\t")
            else:
                print(l, end= "\t")
        print("\n")

        for i in range(self.m): # lignes
            print(self.letters_orders2[i], end ="\t")
            if self.letters_orders2[i] == "-":
                #print("  ", end= "")
                pass
            for j in range(self.n): # colonne
                print(self.mat[i][j], end= "\t")
            print("\n")
        print("\n") 

    def panda(self):
        """ Autre moyen d'afficher la matrice de manière plus compacte """
        
        print(pandas.DataFrame(self.mat, list(self.letters_orders2), \
        pandas.Index(list(self.letters_orders1), name="*")), end="\n\n")

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
                    self.set_letters_seq(line, line)   # lettres de seq 1 et 2
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

        first_col = [0]
        for i in range(self.m):

            newline = []
            if i == 0:
                newline += [0, -self.I]    # 1ere ligne [0, -I, -I-E, ...]
                for j in range(2, self.n):
                    newline.append(newline[j-1] - self.E)
                

            elif i == 1:
                newline += [-self.I] + (self.n-1)*[""]
                first_col += [-self.I]
            else:
                first_col.append(first_col[i-1] - self.E)
                newline += [first_col[i]] + (self.n-1)*[""]
                
            self.set_lign(newline)

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

        self.S = MatScoring(I, E, self.n, self.m) 
        # on crée un objet matrice de scoring

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
        return max( self.S.get_score(i-1, j) - self.I, 
                    self.V.get_score(i-1,j) - self.E )

    def get_w(self, i, j):
        return max( self.S.get_score(i, j-1) - self.I,
                    self.W.get_score(i, j-1) - self.E )

    def is_previous(self, mat, i, j):
        """ Regarde si l'élement ij résulte de l'élément en diagonale,
            en haut ou a gauche 
        """

        res = False
        if mat == "v":
            if self.V.get_score(i,j) == self.S.get_score(i,j):
                res = True 

        elif mat == "w":
            if self.W.get_score(i,j) == self.S.get_score(i,j):
                res = True 

        elif mat == "s":
            letters_ij = self.S.get_letters(i,j) # 'AB' par exemple
            t_ij = self.t.get_acid_score(letters_ij[0], letters_ij[1])  # t(i,j)

            if self.S.get_score(i-1, j-1) + t_ij == self.S.get_score(i,j):
                res = True 

        return res 



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

    def update_s_global(self, i, j, v_ij, w_ij):
        letters_ij = self.S.get_letters(i,j) # 'AB' par exemple
        t_ij = self.t.get_acid_score(letters_ij[0], letters_ij[1])

        s_ij = max( self.S.get_score(i-1,j-1) + t_ij, v_ij, w_ij )
        self.S.set_score(i,j, s_ij)

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

                self.update_s_global(i, j, v_ij, w_ij)

        if self.p == 2: # Si on veut afficher les matrices
            self.V.display()
            self.W.panda()
            self.S.panda()

        self.current_sol.append((self.m-1, self.n-1)) 
        self.backtracking_global(self.m-1, self.n-1) #appel sur element en derniere
                                                     #ligne, derniere colonne

        
        R = Result(self.S, self.t, self.all_solutions, self.p)
        res = R.compute_result()
        return res

    def backtracking_global(self, i, j):
        """ Remonte la matrice de scoring a partir du dernier élément jusqu'à [0][0]
            pour avoir les k alignements
        """

        if i == 0 or j == 0:
            if i == 0 and j == 0:
                if len(self.current_sol) > 0:
                    self.current_sol.pop() # Si on est en (0,0)
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
                    self.current_sol.append((new_i, new_j))
                    self.backtracking_global(new_i, new_j) # appel sur cellule suivante
                    if len(self.current_sol) != 0:
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

    def update_s_local(self, i, j, v_ij, w_ij):
        """ détermine la valeur de S en fonction de V et W et de t"""

        letters_ij = self.S.get_letters(i,j) # 'AB' par exemple
        t_ij = self.t.get_acid_score(letters_ij[0], letters_ij[1])
        s_ij = max( self.S.get_score(i-1,j-1) + t_ij, v_ij, w_ij, 0 )
        self.S.set_score(i,j, s_ij)    

    
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

                    self.update_s_local(i, j, v_ij, w_ij)

        if self.p == 2:
            print("la matrice de Scoring: ")
            self.S.panda()

    
    def sol_found(self, i,j):
        """ Détermine si on a fini un alignement local """

        return (i == 0 or j == 0 or self.S.get_score(i,j) == 0)

    def bottom_up(self, i,j):
        """ Remonte la matrice de scoring a partir du max de la matrice
            de scoring jusqu'à un élément de la 1ere ligne ou 1ere colonne
            ou un élément 0
        """
        self.current_sol.append((i,j))
        while not self.sol_found(i,j):
            if self.is_previous("v", i, j):    # haut
                i-=1
            elif self.is_previous("w", i, j):  # gauche
                j-=1
            elif self.is_previous("s", i, j):  # diagonale
                i-=1
                j-=1
            if self.S.get_score(i,j) != 0:
                self.current_sol.append((i,j))

        self.all_solutions.append(self.current_sol)
        self.zeros += self.current_sol
        self.current_sol = []

    
    def Smith_Waterman(self):
        self.compute_scoring(1,1)

        for i in range(self.l):
            current_max = self.S.get_max()
            self.bottom_up(current_max[0], current_max[1])
            if i == self.l-1: # Si on a trouvé l solutions
                break
            # une fois le backtrack fini, on met a 0 l'alignement 
            self.S.set_zero(self.all_solutions[-1])
            self.V.set_zero(self.all_solutions[-1])
            self.W.set_zero(self.all_solutions[-1])
            if self.p == 2:
                print("mise a 0: ")
                self.S.panda()
                self.V.panda()
                self.W.panda()

                print("recalcul de la matrice: ")
            

            self.compute_scoring(self.all_solutions[-1][-1][0], \
                                self.all_solutions[-1][-1][1])

        R = Result(self.S, self.t, self.all_solutions, self.p)
        return R.compute_result()



class Result:
    """ Classe représentant un objet Résultat dans laquelle on va aligner
        2 séquences selon la matrice de scoring obtenue
    """

    def __init__(self, S, t, all_sol, p):
        self.S = S  # matrice de scoring 
        self.t = t # matrice de substitution
        self.all_solutions = all_sol
        self.gap = "-"
        self.p = p # Pour savoir si on veut print le résultat ou non
        
    def compute_result(self):
        """ trouve les correspondances entre lettres de l'alignement """
        
        scores_sim = [] # stocke le score de similarité des résultats

        for sol in self.all_solutions:
            used = {} # dictionnaire avec indices deja utilisés 
                    # clé = colonnes, valeur = lignes
            res = []
            for pos in range(len(sol)-1, -1, -1):
                
                letters = self.S.get_letters(sol[pos][0], sol[pos][1])
                score = self.t.get_acid_score(letters[0], letters[1])
                if sol[pos][0] in used.values():
                    letters = self.gap + letters[1] # '-B' par ex
                if sol[pos][1] in used.keys():
                    letters = letters[0] + self.gap # 'B-' par ex

                res.append((letters, score)) # On ajoute son score aussi

                used[sol[pos][1]] = sol[pos][0]

            #self.bind(res)
            scores_sim.append(self.bind(res))
        
        return scores_sim


    def bind(self, seq):
        """ crée la liaison entre les 2 séquences (similarité, identité, ..)
            ainsi que les scores et pourcentages
        """

        # print(seq)
        seq1, seq2, links = "", "", ""
        for i in range(len(seq)):
            # print("seq: ", seq[i][0], " score: ", seq[i][1])
            seq1 += seq[i][0][1]  # 2eme lettre
            seq2 += seq[i][0][0]  # 1ere lettre

            if seq1[i] == self.gap or seq2[i] == self.gap:
                links += " "   # Pas de correspondance
            else:
                if seq1[i] == seq2[i]:
                    links += ":"       # Identiques
                elif seq[i][1] >= 0:   # Si le score est positif => similaires
                    links += "."
                else:
                    links += " "
                    
        similarity = ((links.count(".")+links.count(":"))/len(links)) * 100 # pourcentage de similarité
        identity = (links.count(":")/len(links)) * 100 # pourcentage d'identité
        
        if self.p == 1 or self.p == 2:
            self.print_result(seq1, seq2, links, similarity, identity)
            
        return similarity

    
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



def main():
    # seq1 = Sequence("Séquence 1 de l'exemple","AZAP")
    # seq2 = Sequence("Séquence 2 de l'exemple", "AI")
    # G0 = GlobalAlignment(3, 12, 2, "blosum62.txt", seq1, seq2, 2)
    # G0.Needleman_Wunsch()

    P2 = ParserSequence("protein-sequences.fasta")
    P2.parse()

    # seq1 = Sequence("Séquence 1 de l'exemple","ISALIGNED")
    # seq2 = Sequence("Séquence 2 de l'exemple", "THISLINE")
    # L0 = LocalAlignment(2, 4, 4, "blosum62.txt", seq1, seq2, 2)
    # L0.Smith_Waterman()

    L3 = LocalAlignment(3, 12, 2, "blosum62.txt", P2.get_seq(0), P2.get_seq(1), 1)
    L3.Smith_Waterman()

if __name__ == "__main__":
    main()