import copy as c
import math
import pandas as p
import time
import other_ADT as o

p.set_option('display.max_columns', 500)
p.options.display.float_format = '{:,.3f}'.format

class Alignment:
    """ Classe représentant un objet qui trouvera l'alignement de 2 séquences
        d'acides aminées
    """

    def __init__(self, I, E, mat_file, seq1, seq2, p, prof=None):
        
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
            print("Pénalité de gap : I = {0} | E = {1}".format(self.I, self.E))

        # =================== SUBSITUTION ==============================
        self.t = o.MatSubstitution(mat_file)
        self.t.parse()

        # =================== SCORING ===============================

        self.S = o.MatScoring(I, E, self.n, self.m) # Création matrice scoring

        # ===================== V ET W ===================================

        self.V = o.MatV(self.n, self.m)
        self.W = o.MatW(self.n, self.m)

        self.V.init_V() # Initialise V
        self.W.init_W() # Initialise W
        self.V.setup_seq("-"+self.seq1.get_acids(), "-"+self.seq2.get_acids())
        self.W.setup_seq("-"+self.seq1.get_acids(), "-"+self.seq2.get_acids())

        self.current_sol = [] # Pour le backtracking
        self.all_solutions = []
        
        self.prof = prof # Si on veut aligner avec une matrice PROFIL

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
        
        elif mat == "vp" and self.S.get_score(i-1,j)-self.I == self.S.get_score(i,j):
            res = True 
        elif mat == "wp" and self.S.get_score(i,j-1)-self.I == self.S.get_score(i,j):
            res = True
        elif mat == "sp": # Si on a aligné à un profil
            pssm_ij = self.prof.get_cell(\
                self.prof.acid_pos(self.seq2.get_acid(i-1)),j-1)
            if self.S.get_score(i-1, j-1) + pssm_ij == self.S.get_score(i,j):
                res = True 

        return res 

    def add_sol(self, i, j, pos, p):
        """ Ajoute une solution à la liste des solutions
            De la forme: ('AB', score)
        """
        if p: # Si on aligne au profil
            # self.current_sol.append((self.seq2.get_acid(i), i))
            self.current_sol.append(i)
        else:
            letters_ij = self.S.get_letters(i,j) 
            if pos == 0:
                letters_ij = letters_ij[0]+"-" # A- par exemple => Si gap vient de V
            elif pos == 1:
                letters_ij = "-"+letters_ij[1] # -A par exemple => Si gap vient de W
            # if p: # Si on aligne à un profil
            #     self.current_sol.append((letters_ij, i)) # (A-, 0) par ex
            else:
                self.current_sol.append((letters_ij, \
                                self.t.get_acid_score(letters_ij[0], letters_ij[1])))


class GlobalAlignment(Alignment):
    """
    Classe qui va trouver un aligement de séquences d'acides aminées avec
    une pénalité affine et une méthode globale
    """

    def __init__(self, k, I, E, mat_file, seq1, seq2, p, prof=None):
        Alignment.__init__(self, I, E, mat_file, seq1, seq2, p, prof)
        self.k = k

        self.S.init_S_global()  # On initialise la matrice de scoring
        self.S.setup_seq("-"+self.seq1.get_acids(), "-"+self.seq2.get_acids())

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
                s_ij = self.get_s_global(i, j, v_ij, w_ij)
                self.S.set_score(i, j, s_ij) # S(i,j)

        if self.p == 2: # Si on veut afficher les matrices
            self.V.panda()
            self.W.panda()
            self.S.panda()

        # self.backtracking_global(self.m-1, self.n-1) #appel sur dernier element
        self.bottom_up(self.m-1, self.n-1)

        R = Result(self.all_solutions, self.p)
        return R.bind()

    def backtracking_global(self, i, j):
        """ Remonte la matrice de scoring a partir du dernier élément jusqu'à [0][0]
            pour avoir les k alignements
        """

        if i == 0 or j == 0:
            if not (i==0 and j==0):
                # self.current_sol.append((i,j))
                self.add_sol(i,j,2, self.prof != None)
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
                    self.add_sol(i, j, pos, self.prof != None)
                    self.backtracking_global(new_i, new_j) # appel sur cellule suivante
                    # if len(self.current_sol) > 1:
                    self.current_sol.pop() # destruction sol partielle

    def bottom_up(self, i,j):
        """ 
        Remonte la matrice de scoring du dernier élément vers l'élement (0,0) 
        pour construire l'alignement 
        """
        while not (i == 0 or j == 0):
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
            self.add_sol(prev_i,prev_j,pos, self.prof != None)

        if not (i == 0 and j == 0):
            self.add_sol(i, j, 2, self.prof != None)
        self.all_solutions.append(self.current_sol)

class LocalAlignment(Alignment):
    """
    Classe qui qui va trouver un aligement de séquences d'acides aminées avec
    une méthode locale
    """

    def __init__(self, l, I, E, mat_file, seq1, seq2, p, prof=None):
        Alignment.__init__(self, I, E, mat_file, seq1, seq2, p, prof)
        self.l = l

        self.S.init_S_local()  # On initialise la matrice de scoring
        self.S.setup_seq("-"+self.seq1.get_acids(), "-"+self.seq2.get_acids())

        self.zeros = []
        self.found = False

    def get_s_local(self, i, j, v_ij, w_ij):
        """ détermine la valeur de S en fonction de V et W et de t"""

        if self.prof == None: # Si on n'aligne pas à un profil
            letters_ij = self.S.get_letters(i,j) # 'AB' par exemple
            t_ij = self.t.get_acid_score(letters_ij[0], letters_ij[1])
            return max( self.S.get_score(i-1,j-1) + t_ij, v_ij, w_ij, 0 )
    
        else: 

            pssm_ij = self.prof.get_cell(\
                self.prof.acid_pos(self.seq2.get_acid(i-1)),j-1)

            return max(self.S.get_score(i-1,j-1) + pssm_ij, v_ij, w_ij, 0)

    def compute_scoring(self, start_i, start_j):
        """ 
        Calcule la matrice de scoring pour l'alignement local
        """

        # ==================== CREATION SCORING ===========================

        for i in range(start_i, self.m):
            for j in range(start_j, self.n):
                if (i,j) not in self.zeros:
                    if self.prof == None:
                        v_ij = self.get_v(i, j)     # V(i,j)
                        self.V.set_score(i,j, v_ij)
                        w_ij = self.get_w(i, j)     # W(i,j)
                        self.W.set_score(i,j, w_ij)
                    elif self.prof != None:
                        v_ij = self.S.get_score(i-1,j) - self.I
                        w_ij = self.S.get_score(i,j-1) - self.I
                    s_ij = self.get_s_local(i, j, v_ij, w_ij)
                    self.S.set_score(i, j, s_ij)

        if self.p == 2:
            print("la matrice de Scoring: ")
            self.S.panda()

    # def is_modified(self, i, j):
    #     res = False 
    #     if (i,j) in self.zeros:
    #         res = True
    #     else: 
    #         if self.prof == None:
    #             v_ij = self.get_v(i,j)
    #             w_ij = self.get_w(i,j)
    #         elif self.prof != None:
    #             v_ij = self.S.get_score(i-1,j) - self.I
    #             w_ij = self.S.get_score(i,j-1) - self.I
    #         s_ij = self.get_s_local(i,j, v_ij, w_ij)
    #         if s_ij != self.S.get_score(i,j):
    #             if self.prof == None:
    #                 self.V.set_score(i,j,v_ij) # si c'est modifié, on modifie
    #                 self.W.set_score(i,j,w_ij)
    #             self.S.set_score(i,j,s_ij)
    #             res = True 
    #     return res 

    # def recompute_scoring(self, start_i, start_j):
    #     """ test de recalcul de la matrice de scoring"""

    #     i = start_i
    #     j = start_j 
    #     while i != self.m and j != self.n:
    #         while i != self.m and j != self.n and self.is_modified(i, j):
    #             j += 1 # colonne suivante 
    #         j = start_j 
    #         i += 1
    #         while i != self.m and j != self.n and self.is_modified(i, j):
    #             i += 1 # ligne suivante 
    #         start_i += 1
    #         start_j += 1
    #         i = start_i # sous carré suivant
    #         j = start_j
    #     if self.p == 2:
    #         self.S.panda()

    
    def sol_found(self, i,j):
        """ Détermine si on a fini un alignement local """

        return (i == 0 or j == 0 or self.S.get_score(i,j) == 0)

    def bottom_up(self, i,j):
        """ Remonte la matrice de scoring a partir du max de la matrice
            de scoring jusqu'à un élément de la 1ere ligne ou 1ere colonne
            ou un élément 0
        """
        while not self.sol_found(i,j):
            # print(i)
            prev_i = i
            prev_j = j
            pos = 2

            if self.prof == None and self.is_previous("v", i, j) or \
                self.prof != None and self.is_previous("vp", i, j):   # haut
                # print("v")
                i-=1
                pos = 0
            elif self.prof == None and self.is_previous("w", i, j) or \
                self.prof != None and self.is_previous("wp", i, j): # gauche
                # print("w")
                j-=1
                pos = 1
            elif self.prof == None and self.is_previous("s", i, j) or \
                self.prof != None and self.is_previous("sp", i, j):  # diagonale
                # print("sp")
                i-=1
                j-=1
            if self.S.get_score(i,j) != 0:
                self.add_sol(prev_i,prev_j,pos, self.prof != None)
                self.zeros.append((prev_i,prev_j))

        self.all_solutions.append(self.current_sol)
        self.add_sol(prev_i, prev_j, 2, self.prof != None)
        self.zeros.append((prev_i, prev_j))

    
    def Smith_Waterman(self):
        self.compute_scoring(1,1)
        # self.S.panda()
        if self.prof != None:
            self.all_solutions.append(self.seq2)

        for i in range(self.l):
            current_max = self.S.get_max()
            print("max : ", current_max[0])
            self.bottom_up(current_max[0], current_max[1])
            if i == self.l-1: # Si on a trouvé l solutions
                break
            # une fois le backtrack fini, on met a 0 l'alignement 
            sols = len(self.current_sol)
            self.S.set_zero(self.zeros[len(self.zeros)-sols:])
            if self.prof == None:
                self.V.set_zero(self.zeros[len(self.zeros)-sols:])
                self.W.set_zero(self.zeros[len(self.zeros)-sols:])

            self.current_sol = []
            if self.p == 2:
                print("mise a 0: ")
                self.S.panda()

                print("recalcul de la matrice: ")
            # self.recompute_scoring(self.zeros[-1][0], self.zeros[-1][1])
            self.compute_scoring(self.zeros[-1][0], self.zeros[-1][1])

        R = Result(self.all_solutions, self.p)
        if self.prof == None:
            return R.bind()
        else:
            #print(self.all_solutions)
            return R.display_local_prof()


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

    def display_local_prof(self):
        """
        Imprime les sous-chaînes issues de l'alignement local entre une séquence
        de prot et un Profil
        """
        for sol in range(1, len(self.all_sol)):
            tmp_res = ""
            # print(self.all_sol[sol][-1], self.all_sol[sol][0])
            for i in range(self.all_sol[sol][-1]-1, self.all_sol[sol][0]-1):
                tmp_res += self.all_sol[0].get_acid(i) 
            # for s in sol:
                # tmp_res = s[0][0] + tmp_res # : le A de A- par ex

                # tmp_res = s[0] + tmp_res
            print("Sous séquence obtenue: ", tmp_res, len(tmp_res))
            # print("=> intervalle: [{0} -> {1}]".format(sol[-1][1]-1, sol[0][1]-1))
            print("=> intervalle: [{0} -> {1}]".format(
                self.all_sol[sol][-1], self.all_sol[sol][0]-1))

def get_best_globalalignments(S, I, E):

    # sims = o.Matrix()
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

    print("j:", end="")
    for j in range(len(S)):
        if S[j].get_acids() == seq.get_acids():
            print("\n")
            return True 
        G_tmp = GlobalAlignment(1, I, E, "pam120.txt", S[j], seq, 0)
        tmp_score = G_tmp.Needleman_Wunsch()
        if tmp_score[0] > 60:
            print("\n")
            return True 
        print(j, end=", ", flush=True)
    return False

def remove_high_scores(S, I, E):
    # S = liste de séquences
    reduced = [S[0]]
    for i in range(1, len(S)):
        print("i = ", i)
        if not is_very_similar(reduced, I, E, S[i]):
            print(" ==> ajouté")
            reduced.append(S[i])
        print("len reduced: ", len(reduced))

    reduced_file = open("reduced_pam120.fasta", "w") 
    for s in reduced:
        reduced_file.write(s.get_title())
        reduced_file.write(s.get_acids()+"\n")
    reduced_file.close()

class Profile:
    """ Classe représentant un PSSM, un profil, qui représente de manière
        compacte un alignement de séquence multiple et donc un bromodomain
    """

    def __init__(self, MSA, pen):
        self.prof = o.Matrix()
        self.pen = pen # pénalité de gap
        self.MSA = MSA # Liste des séquences alignées de manière multiple
        self.acids = "RHKDESTNQCGPAILMFWYV-"
        self.prof.setup_seq("", self.acids)
        self.Nseq = len(self.MSA) # nombre de séquences
        self.Nacids = 20 # nombre d'acides aminées
        self.len_seq = self.MSA[0].length() # longueur des séquences du MSA
        y = self.MSA[0].get_acids()

        # pseudocounts
        self.alpha = self.Nseq - 1  # facteur de cadrage pour les données observées
        self.beta = math.sqrt(self.Nseq)  # facteur de cadrage
        # self.beta = 1

        self.pa = {"R": 5.53, "H": 2.27,"K": 5.81,"D":5.46,"E":6.73,"S":6.62,\
        "T":5.35,"N":4.05,"Q":3.93,"C":1.38,"G":7.07,"P":4.73,"A":8.25,"I":5.92,\
        "L":9.65,"M":2.41,"F":3.86,"W":1.09,"Y":2.91,"V":6.86} # probabilité d'apparition des acides aminées

        self.consencus = ""

    def get_cell(self, i, j):
        """
        Renvoie l'élément en position i,j de la matrice Profil
        """
        return self.prof.get_score(i,j)

    def acid_pos(self, acid):
        """
        Renvoie la position de l'acide aminé (la ligne dans la matrice)
        """
        return self.acids.index(acid)

    def number_acids(self, u, a):
        """
        Calcule le nombre d'acide a dans la colonne u
        """
        Nua = 0 
        for i in range(self.Nseq):
            if self.MSA[i].get_acid(u) == a: # acide aminé a la séquence i, col u
                Nua += 1
        return Nua
    
    def frequency(self, u, b):
        """ 
        Calcule la fréquence d'apparition d'un acide aminé b dans la colonne b
        parmi toutes les séquences de l'alignement multiple
        """

        Nub = 0 # nombre d'acide b dans la colonne u
        for i in range(self.Nseq):
            if self.MSA[i].get_acid(u) == b: # acide aminé a la séquence i, col u
                Nub += 1
        
        return Nub / self.Nseq


    def q1(self, u, a):
        """
        Calcule les pseudocounts associés à l'acide aminé a
        """
        return (self.number_acids(u, a) + self.beta*(self.pa[a]/100)) \
                / (self.Nseq + self.beta )

    def q2(self, u, a):
        """
        Calcule les pseudocounts associés à l'acide aminé a
        """
        # print(self.frequency(u,a))
        return (self.alpha*self.frequency(u, a) + self.beta*(self.pa[a]/100)) \
                / (self.alpha + self.beta)
        
    def compute_profile(self, q):
        print("length seq: ", self.Nseq)
        print("beta: ", self.beta)
        print("alpha: ", self.alpha)

        for a in range(self.Nacids):   # lignes
            self.prof.addline()
            for u in range(self.len_seq): # colonnes
                if q == 1:
                    Qua = self.q1(u, self.acids[a])
                elif q == 2:    
                    Qua = self.q2(u, self.acids[a])
                Mua = math.log10(Qua/(self.pa[self.acids[a]]/100))
                # print("u = {0}, a = {1}, Qua = {2}, Mua = {3}".format(\
                    # u, self.acids[a], Qua, Mua))
                self.prof.add_cell(a, Mua)
        
        self.prof.set_lign([self.pen]*self.len_seq)
        # print("1ere ligne: ", self.prof[0])
        # self.prof.panda(self.len_seq)

    def compute_consencus(self):
        """ Renvoie la séquence dont chaque acide aminé possède le score maximal
            dans sa position dans le profil
        """

        for u in range(self.len_seq):
            # print("u: ", u)
            tmp_max = ["-", -float("inf"), 0]
            for a in range(self.Nseq):
                # print(self.MSA[a].get_acid(u), end=", ")
                if self.MSA[a].get_acid(u) != "-":
                    Mua = self.prof.get_score(self.acids.index(self.MSA[a].get_acid(u)),u)
                    # print(self.MSA[a].get_acid(u), " Mua: ", Mua)
                    if Mua > tmp_max[1]:
                        # tmp_max = (self.MSA[a].get_acid(u), Mua)
                        tmp_max[0] = self.MSA[a].get_acid(u)
                        tmp_max[1] = Mua
                    tmp_max[2] += 1
                # else:
                #     tmp_max[2] += 1
                #     break
            # print("\n")
            if tmp_max[2] <= 10 and self.Nseq > 3:
                tmp_max[0] = "-"
            # print("u: {0}, nb: {1}".format(u, tmp_max[2]))
            self.consencus += tmp_max[0]
        
        # print("consencus: ", self.consencus, len(self.consencus))
        self.display_consencus()

    def display_consencus(self):
        """
        Affiche de manière jolie le consencus 
        """
        print("consencus: ")
        # i = 
        step = 0
        i = 0
        while i < self.len_seq // 40:
            for j in range(step, step+40):
                print(self.consencus[j], end= " ")
            print("\n")
            # # print(self.consencus[step:step+40])
            # for j in range(step, step+40, 5):
            #     print(j, end="          ")
            # print("\n")
            step += 40
            i += 1
        for j in range(step, self.len_seq):
            print(self.consencus[j], end= " ")
        print("\n")
        # for j in range(step, self.len_seq, 5):
        #     print(j, end="     ")
        # print("\n")



def main():
    # =======================================================================
    # ================================= GLOBAL ==============================
    # =======================================================================

    # P = o.ParserSequence("BRD-sequence.fasta")
    # P.parse()

    # P1 = o.ParserSequence("to-be-aligned.fasta") 
    # P1.parse()

    # print("========================= TEST BEST SCORE to be aligned ========\n\n")
    # best= get_best_globalalignments(P.get_all_seq(), 4, 4)
    # print("le meilleur score est celui-ci: ", best[0])

    # REDUCE
    # begin = time.time()
    # remove_high_scores(P1.get_all_seq(), 2, 2)
    # print("temps écoulé: {0} secondes".format(time.time() - begin))

    ## TEST si le reduce a bien max 60%
    # begin = time.time()
    # Ptest = o.ParserSequence("reduced_penalite2.fasta")
    # Ptest.parse()
    # best= get_best_globalalignments(Ptest.get_all_seq(), 2, 2)
    # print("le meilleur score est celui-ci: ", best[0])
    # print("temps écoulé: ", time.time() - begin)

    # print(" ==============================================TEST AZAP et AI: \n\n")
    # seq1 = o.Sequence("Séquence 1 de l'exemple","AZAP")
    # seq2 = o.Sequence("Séquence 2 de l'exemple", "AI")
    # G0 = GlobalAlignment(3, 12, 2, "blosum62.txt", seq1, seq2, 2)
    # G0.Needleman_Wunsch()

    # print(" =======================================TEST BEST SCORE: \n\n")
    # best= get_best_globalalignments(P.get_all_seq(), 12, 2, None)
    # print("le meilleur score est celui-ci: ", best[0])
    # G1 = GlobalAlignment(1, 12, 2, "blosum62.txt", P.get_seq(best[1]), P.get_seq(best[2]), 1)
    # G1 = G1.Needleman_Wunsch()

    # print(" =======================================TEST MGGETFA et GGVTTF: \n\n")
    # seq1 = o.Sequence("Séquence 1 de l'exemple","MGGETFA")
    # seq2 = o.Sequence("Séquence 2 de l'exemple", "GGVTTF")
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




    # =======================================================================
    # # =============================== LOCAL ===============================
    # =======================================================================
    # P2 = o.ParserSequence("protein-sequences_part1.fasta")
    # P2.parse()

    # seq1 = o.Sequence("Séquence 1 de l'exemple","ISALIGNED")
    # seq2 = o.Sequence("Séquence 2 de l'exemple", "THISLINE")
    # L0 = LocalAlignment(2, 4, 4, "blosum62.txt", seq1, seq2, 2)
    # L0.Smith_Waterman()

    # Identiques
    # L1 = LocalAlignment(1, 12, 2, "blosum62.txt", P2.get_seq(0), P2.get_seq(0), 1)
    # L1.Smith_Waterman()

    # séquence 1 avec séquence 2 du fichier
    # L3 = LocalAlignment(3, 12, 2, "blosum62.txt", P2.get_seq(0), P2.get_seq(1), 1)
    # L3.Smith_Waterman()


    # =======================================================================
    # ================================== PROFIL =============================
    # =======================================================================

    MSA1 = o.ParserSequence("msaresults-MUSCLE.fasta") 
    MSA1.parse() # alignements multiples de séquences 
    # MSA1 = o.ParserSequence("clustal.fasta")
    # MSA1.parse()

    MSA2 = o.ParserSequence("msaresults-reduced-MUSCLE.fasta") 
    MSA2.parse() # alignements multiples de séquences réduits

    # EXEMPLE slides
    # S1 = o.Sequence("1", "TGVEAENLLL")
    # S2 = o.Sequence("2", "PRAKAEESLS")
    # S3 = o.Sequence("3", "GRKDAERQLL")
    # Prof0 = Profile([S1, S2, S3], 4)
    # Prof0.compute_profile(2)
    # Prof0.compute_consencus()

    # PROFIL MSA
    # Prof = Profile(MSA1.get_all_seq(), 4)
    # Prof.compute_profile(2)
    # Prof.compute_consencus()

    # #PROFIL MSA REDUCED
    Prof1 = Profile(MSA2.get_all_seq(), 4)
    Prof1.compute_profile(2)
    Prof1.compute_consencus()



    # =========================================================================
    # ================================ LOCAL + PROFIL =========================
    # =========================================================================

    P3 = o.ParserSequence("protein-sequences.fasta") # protéines à aligner avec profil
    P3.parse()

    # print("Alignement Séquence 1 avec profil")
    # Lp0 = LocalAlignment(2, 4, 4, "blosum62.txt", o.Sequence("", "-"*(MSA1.get_seq(0).length())),
    #                         P3.get_seq(0), 1, Prof)
    # Lp0.Smith_Waterman()

    # print("Alignement Séquence 2 avec profil")
    # Lp1 = LocalAlignment(6, 2, 2, "blosum62.txt", o.Sequence("", "-"*(MSA1.get_seq(0).length())),
    #                         P3.get_seq(1), 1, Prof)
    # Lp1.Smith_Waterman()

    # print("Alignement Séquence 1 avec profil de reduced")
    # Lp2 = LocalAlignment(2, 4, 4, "blosum62.txt", o.Sequence("", "-"*(MSA2.get_seq(0).length())),
    #                         P3.get_seq(0), 1, Prof1)
    # Lp2.Smith_Waterman()

    # print("Alignement Séquence 2 avec profil de reduced")
    # Lp3 = LocalAlignment(1, 4, 4, "blosum62.txt", o.Sequence("", "-"*(MSA2.get_seq(0).length())),
    #                         P3.get_seq(1), 1, Prof1)
    # Lp3.Smith_Waterman()


    # print("Alignement Séquence EXEMPLE avec profil EXEMPLE")
    # Lptest = LocalAlignment(2, 4, 4, "blosum62.txt", o.Sequence("", 10*"-"),
    #                         o.Sequence("", "SRNAAEYLLS"), 1, Prof0)
    # Lptest.Smith_Waterman()

if __name__ == "__main__":
    main()

