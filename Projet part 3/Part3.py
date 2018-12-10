import math as m
import time as t
import statistics


class ParserDSSP:
    """
    Classe représentant un parser de fichier .dssp
    """

    def __init__(self, cath_file, dssp):
        try:
            self.cath = open(cath_file, "r") # fichier avec les identifiants dssp
        except:
            print("Entrez un fichier valide s.v.p")

        self.secondary_struct = {"H": "H", "G": "H", "I": "H",\
                                 "E": "E", "B": "E"\
                                ,"T": "C", "C": "C", "S": "C", " ": "C"}
        
        self.dssp_dir = dssp

    def get_title(self, protein_name, organism):
        """
        Renvoie le titre de la prot sous forme:
        > Identifiant|Nom Prot|Organisme
        """

        res = "> "
        res += protein_name[-2] + "|"      # identifiant
        for n in protein_name[1:-3]:
            res += n + " "
        res = res[:-1] + "|" + " ".join(organism) # organisme
        return res

    def parse(self, dssp_file, chain):
        """
        Parse un fichier dssp et renvoie les informations
        """
        try:
            res = ["> ", "", ""] # contient l'entête de la prot, la prot et structure sec
            d = open(dssp_file, "r")
        except:
            print("fichier dssp invalide")
        
        # ====== TITLE =======
        lines = d.readlines()
        res[0] = self.get_title(lines[2].split(), (lines[4].split(":")[1]).split()[:-1])

        # ====== PROT + STRUCTURE SECONDAIRE ======
        for i in range(28, len(lines)):
            # Si la chaine est celle qu'on donne et le résidu n'est pas X,Z, B
            if lines[i][11] == chain and lines[i][13] not in "XZB":
                if lines[i][13].islower():
                    res[1] += "C"
                else:
                    res[1] += lines[i][13] # acide aminé
                res[2] += self.secondary_struct[lines[i][16]] # structure secondaire

        return res[0].strip() + "\n" + res[1] + "\n" + res[2] + "\n"

    def create_proteins(self, prot_file):
        """
        Crée les protéines contenues dans les fichier dssp parsées et les met
        dans un fichier
        """
        try:
            p = open(prot_file, "w")
            for line in self.cath:
                l = line.split()[0]
                dssp_file = "dataset/dataset/{}/{}.dssp".format(self.dssp_dir, l[:-1])
                p.write(self.parse(dssp_file, l[-1]))
        except:
            print("Fichier proteins invalide")

class Protein:
    """
    Classe représentant une protéine en format:
    > identifier|Protein name|Organism
    Protein
    Secondary struct
    """

    def __init__(self, title="", seq="", struct="", prediction = ""):
        self.title = title
        self.seq = seq
        self.struct = struct 
        self.prediction = prediction
    
    def get_seq(self):
        """Renvoie la séquence"""
        return self.seq
    
    def set_seq(self, seq):
        """Met une séquence"""
        self.seq = seq

    def get_aa(self, i):
        """Renvoie l'acide aminé à la position i dans la séquence"""
        return self.seq[i]

    def set_prediction(self, p):
        """ Met une prédiction """
        self.prediction = p

    def get_prediction(self):
        """Renvoie la prédiction"""
        return self.prediction 
    
    def get_pred(self, i):
        """Renvoie la structure i de la prédiction"""
        return self.prediction[i]

    def get_structures(self):
        """Renvoie la structure secondaire """
        return self.struct

    def get_struct(self, i):
        """Renvoie la structure à la position i"""
        return self.struct[i]

    def set_struct(self, s):
        """Met une nouvelle structure"""
        self.struct = s

    def __str__(self):
        """Affiche la protéine"""
        return self.title + "\n" + self.seq + "\n" + self.struct

class ParserProteins:
    """
    Classe représentant un parser de protéines issues des fichiers dssp de la forme
    > identifier|Protein name|Organism
    Protein
    Secondary struct
    """

    def __init__(self, prot_file):
        try:
            self.file = open(prot_file, "r")
        except:
            print("fichier invalide")

        self.proteins = [] # liste de protéines
        
    def get_prot(self, i):
        """
        Renvoie la protéine i
        """
        return self.proteins[i]

    def get_proteins(self):
        """
        Renvoie la liste des protéines
        """
        return self.proteins

    def parse(self):
        """
        Parse le fichier et met les protéines dans une liste
        """
        lines = self.file.readlines()
        for i in range(0, len(lines), 3):
            self.proteins.append(Protein(lines[i].strip(), lines[i+1].strip(), \
                                    lines[i+2].strip()))

    def display(self):
        for i in self.proteins:
            print(i, end="\n\n")
            
class Counter:
    """
    Classe représentant un compteur, celui-ci va servir à comptabiliser les 
    fréquences 
    """

    def __init__(self, proteins):
        self.aa = {"R","H","K","D","E","S","T","N","Q","C","G","P","A","I","L","M",\
              "F","W","Y","V"}
        self.o = {"H": "EC", "E": "HC", "C": "HE"} # les n - S possibles

        self.F_s = {"H": 0, "E": 0, "C": 0}
        self.F_sr = {"H": {}, "E": {}, "C": {}}
        self.F_srr = {"H": {}, "E": {}, "C": {}}
        for a in self.aa:
            for s in "HEC":
                self.F_sr[s][a] = 0
                self.F_srr[s][a] = {}

        self.proteins = proteins # liste de protéines (+/- 3000)

    def Freq_s(self, s, o = False):
        """
        Renvoie la fréquence d'apparition de la structure s
        """
        if o: # Si n-s
            return self.F_s[self.o[s][0]] + self.F_s[self.o[s][1]]
        else:
            return self.F_s[s]

    def Freq_sr(self, s, r, o = False):
        """
        Renvoie la fréquence d'apparition de la structure s avec un résidu r
        """
        if o: # Si n-s
            return self.F_sr[self.o[s][0]][r] + self.F_sr[self.o[s][1]][r]
        else:
            return self.F_sr[s][r]

    def Freq_srr(self, s, r, rm, o = False):
        """
        Renvoie la fréquence d'apparition de la structure s avec un résidu r
        et un autre résidu r dans la position m dans son voisinnage
        """
        if o: # Si n-s
            return self.F_srr[self.o[s][0]][r].get(rm, 1) + \
                   self.F_srr[self.o[s][1]][r].get(rm, 1)
        else:
            return self.F_srr[s][r].get(rm, 1)

    def compute_frequencies(self):
        """
        Calcule les fréquences (F_sr, F_s et F_s,rj+m,rj)
        """

        for p in self.proteins:                                     # protéines
            for i in range(len(p.get_structures())):
                self.F_s[p.get_struct(i)] += 1                      # F_s
                self.F_sr[p.get_struct(i)][p.get_aa(i)] += 1        # F_sr

                for j in range(-8, 9):
                    # Si pas 0 et que le voisin existe
                    if j != 0 and i+j >= 0 and i+j < len(p.get_structures()): 
                        if (j, p.get_aa(i+j)) not in \
                            self.F_srr[p.get_struct(i)][p.get_aa(i)]:

                            self.F_srr[p.get_struct(i)][p.get_aa(i)]\
                                [(j, p.get_aa(i+j))] = 1
                        
                        else:
                            # Incrémentation du compteur F_srr
                            self.F_srr[p.get_struct(i)][p.get_aa(i)]\
                                [(j, p.get_aa(i+j))] += 1

        # print("moi: ", self.F_s)
        # print("ric: ", {'E': 151521, 'C': 249950, 'H': 223419})
        # print("nico: ", {'C': 247262, 'E': 151521, 'H': 223301})
        # print("sim: ", {'H': 228538, 'E': 154221, 'C': 253836})
        # self.print_fsr()

    def print_fsr(self):
        l = []
        for i in self.F_srr:
            for a in self.F_srr[i]:
                l += self.F_srr[i][a].values()
            # l += self.F_sr[i].values()
        # print(sorted(l))
        l = sorted(l)
        t = open("compteur_fsrr.fasta", "w")
        for j in l:
            t.write("{}, ".format(j))


class gorIII:
    """
    Classe représentant un prédicteur de structure utilisant la méthode GOR III
    """
    def __init__(self, counter, prot):
        self.c = counter
        self.prot = prot

    def info_individuelle(self, s, r):
        """
        Calcule l'information individuelle:
        log(F_sr/Fn-s,r) + log(Fn-s/Fs)
        """

        return (m.log10(self.c.Freq_sr(s, r)) - m.log10(self.c.Freq_sr(s,r,True)))\
             + (m.log10(self.c.Freq_s(s,True)) - m.log10(self.c.Freq_s(s)))

    def info_directionnelle(self, s, r, rj):
        """
        Renvoie l'information directionnelle (prenant en compte le voisinnage)
        I(Sj, Rj+m, Rj) qui vaut
        log(F_srr/Fn-s,rr) + log(F_n-s,r / F_sr)
        """

        return (m.log10(self.c.Freq_srr(s,r,rj)) - m.log10(self.c.Freq_srr(s,r,rj,True)))\
             + (m.log10(self.c.Freq_sr(s, r, True)) - m.log10(self.c.Freq_sr(s, r)))

    def probability(self, s, r, j):
        """
        calcule la probabilité qu'un résidu adopte une certaine structure
        I(Sj,Rj) + sum(I(Sj,Rj+m, Rj))
        """

        res = self.info_individuelle(s, r)
        for m in range(-8, 9):
            if m != 0 and j+m >= 0 and j+m < len(self.prot.get_seq()):
                res += self.info_directionnelle(s, r, (m, self.prot.get_aa(j+m)))

        return res

    def predict(self):
        """
        prédit la structure secondaire d'une séquence d'acide aminé sur base
        des informations individuelles et directionnelles
        """

        predicted = ""
        for j in range(len(self.prot.get_seq())):
            scores = [0, 0, 0] # [0] = "H", [1] = "E", [2]: "C"
            scores[0] = self.probability("H", self.prot.get_aa(j), j)
            scores[1] = self.probability("E", self.prot.get_aa(j), j)
            scores[2] = self.probability("C", self.prot.get_aa(j), j)

            # print(scores)
            predicted += "HEC"[scores.index(max(scores))]
        
        return predicted

class PredictionQuality:
    """ 
    Classe représentant un évaluateur de la qualité de plusieurs prédictions
    """

    def __init__(self, prot, counter):
        self.prot = prot
        self.counter = counter
        self.scores_q3 = []
        self.scores_mcc = {"H": [], "E": [], "C": []}

    def combine_prot(self):
        """
        Combine toutes les protéines pour former qu'une protéine
        """
        struct = pred = ""
        for p in self.prot:
            struct += p.get_structures()
            pred += p.get_prediction()
        return Protein("", "", struct, pred)

    def Q3(self, prot = None):
        """
        Mesure de qualité d'une prédiction
        => nombre de résidus correctement prédits / nombre total de résidus
        """
        if prot == None: # Si on veut Q3 sur l'ensemble des protéines
            prot = self.combine_prot()

        corrects = 0
        for i in range(len(prot.get_prediction())):
            if prot.get_pred(i) == prot.get_struct(i):
                corrects += 1

        return corrects/len(prot.get_prediction())

    def MCC(self, p = None):
        """
        Renvoie le calcul du mcc sur base d'une matrice de confusion
        """
        mcc = {"H": 0, "E": 0, "C": 0}
        if p == None:
            p = self.combine_prot()

        for x in "HEC":
            TP = FP = FN = TN = 0
            for i in range(len(p.get_prediction())):
                # Si prédit x alors que x
                if p.get_pred(i) == x and p.get_pred(i) == p.get_struct(i): 
                    TP += 1
                # Si prédit x alors que pas x
                elif p.get_pred(i) == x and p.get_struct(i) != x: 
                    FP += 1
                # Si prédit pas x alors que x
                elif p.get_pred(i) != x and p.get_struct(i) == x: 
                    FN += 1
                # Si prédit pas x alors que pas x
                elif p.get_pred(i) != x and p.get_struct(i) != x: 
                    TN += 1
            den = m.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
            if den != 0:
                mcc[x] += (TP*TN-FP*FN) / den
        return mcc

    def test_prediction(self):
        mcc = {"H": [0, 0, 0, 0], "E": [0, 0, 0, 0], "C": [0, 0, 0, 0]}
        conf = [0, 0, 0, 0]
        for p in self.prot:
            tmp_g = gorIII(self.counter, p)
            p.set_prediction(tmp_g.predict())
            self.scores_q3.append(self.Q3(p))
            mcc = self.MCC(p)
            for x in "HEC":
                self.scores_mcc[x].append(mcc[x])

        mcc = self.MCC()
        q3 = self.Q3()
        
        print("Q3 global :", q3)
        print("mcc global : ", mcc)
        

    def ecart_type(self, option):
        """
        Renvoie l'écart type d'un ensemble de données
        """        
        if option == 1:
            return statistics.stdev(self.scores_q3)
        else:
            return [statistics.stdev(self.scores_mcc["H"]),\
                    statistics.stdev(self.scores_mcc["E"]),\
                    statistics.stdev(self.scores_mcc["C"])]
    

def main():

    # ==================== compteurs 3000 prots =======================

    # p = ParserDSSP("dataset/dataset/CATH_info.txt", "dssp")
    # p.create_proteins("proteins.fasta")
    # NICO
    # prot3 = ParserProteins("dataset_nico.fasta")
    # prot3.parse()
    # RIC
    # prot3 = ParserProteins("proteins_ric.fasta")
    # prot3.parse()
    # MOI
    prot3 = ParserProteins("proteins.fasta") # les 3000 prots
    prot3.parse()
    c1 = Counter(prot3.get_proteins()[:3000])
    d = t.time()
    c1.compute_frequencies()
    print("Temps écoulé pour compteurs des 3000 prots: {} sec".format(t.time()-d))
    
    # ======================== PARTIE 2 ==============================

    d2 = t.time()
    # test_prediction(prot3.get_proteins()[3000:], c1) # les 713 autres prots
    qua1 = PredictionQuality(prot3.get_proteins()[3000:], c1)
    qua1.test_prediction()
    print(qua1.ecart_type(1))
    print(qua1.ecart_type(2))
    print("Temps écoulé pour prédire les 713 prots: {} sec".format(t.time()-d2))

    # ======================= PARTIE 3 ==============================

    # p1 = ParserDSSP("dataset/dataset/CATH_info_test.txt", "dssp_test")
    # p1.create_proteins("proteins_test.fasta")
    # prot2 = ParserProteins("proteins_test.fasta") # les 5 prots
    # prot2.parse()
    # d = t.time()
    # test_prediction(prot2.get_proteins(), c1, True)
    # print("Temps écoulé pour prédire les 5 prots: {} sec".format(t.time()-d))

if __name__ == "__main__":
    main()