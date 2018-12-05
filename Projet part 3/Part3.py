import math as m
import time as t


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
        res = res[:-1] + "|" + organism[0][1:] # organisme
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
        res[0] = self.get_title(lines[2].split(), (lines[4].split(":")[1]).split(';'))

        # ====== PROT + STRUCTURE SECONDAIRE ======
        for i in range(28, len(lines)):
            # Si la chaine est celle qu'on donne et le résidé n'est pas X,Z ou B
            if lines[i][11] == chain and lines[i][13] not in "XZ":
                if lines[i][13].islower():
                    res[1] += "C"
                else:
                    res[1] += lines[i][13] # acide aminé
                res[2] += self.secondary_struct[lines[i][16]] # structure secondaire

        # print(res[0] + "\n" + res[1] + "\n" + res[2] + "\n")
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
                # dssp_file = "dataset/dataset/{}/".format(dssp) + l[:-1] + ".dssp"
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

    def __init__(self, title, seq, struct):
        self.title = title
        self.seq = seq
        self.struct = struct 
    
    def get_seq(self):
        """
        Renvoie la séquence
        """
        return self.seq

    def get_aa(self, i):
        """
        Renvoie l'acide aminé à la position i dans la séquence
        """
        return self.seq[i]

    def get_structures(self):
        """
        Renvoie la structure secondaire 
        """
        return self.struct

    def get_struct(self, i):
        """
        Renvoie la structure à la position i
        """
        return self.struct[i]

    def __str__(self):
        """
        Affiche la protéine
        """
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
              "F","W","Y","V","B"}
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

        # print("")
        # print("F_s = {}".format(self.F_s))
        # print("F_sr = {}".format(self.F_sr))
        # print("F_srr = {}".format(self.F_srr))
        # print(self.F_srr["H"]["A"][(-5, "F")])
        # print(self.Freq_srr("H", "A", (-5, "F"), True))
        # print(self.F_srr["C"]["T"][(-4, "V")])
        # print(self.F_srr["E"]["I"][(3, "N")])
        # print(self.F_srr["E"])
        # print(self.F_srr["C"]["G"])
        # print(self.Freq_s("C"))

    def __str__(self):
        """
        Affiche le tableau
        """

        res = ""
        # F_s
        res += "F_s: \nH : {}   |   E : {}   |   C : {}\n\n   "\
            .format(self.F_s["H"], self.F_s["E"], self.F_s["C"])

        # F_sr
        for a in self.aa:
            res += a + "  "
        res += "\nF_sr: \n"
        for s in "HEC":
            res += s + "| "
            for a in self.aa:
                res += "{}  ".format(self.F_sr[s][a])
            res += "\n"
        # F_srr
        # res += "F_srr: \n"
        # for s in "HEC":
        #     res += s + ": \n"
        #     for a in self.aa:
        #         res += "   {} : ".format(a)
        #         for f in self.F_srr[s][a]:
        #             res += "{}({}):{}, ".format(f[1], f[0], self.F_srr[s][a][f])
        #         res += "\n"
        #     res += "\n"
        return res


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

        # return m.log10(self.c.Freq_sr(s,r)/self.c.Freq_sr(s,r,True)) + m.log10(\
        #                self.c.Freq_s(s, True)/self.c.Freq_s(s))

        return (m.log10(self.c.Freq_sr(s, r)) - m.log10(self.c.Freq_sr(s,r,True)))\
             + (m.log10(self.c.Freq_s(s,True)) - m.log10(self.c.Freq_s(s)))

    def info_directionnelle(self, s, r, rj):
        """
        Renvoie l'information directionnelle (prenant en compte le voisinnage)
        I(Sj, Rj+m, Rj) qui vaut
        log(F_srr/Fn-s,rr) + log(F_n-s,r / F_sr)
        """

        # return m.log10(self.c.Freq_srr(s, r, rj)/ self.c.Freq_srr(s, r, rj, True))\
        #      + m.log10(self.c.Freq_sr(s, r, True) / self.c.Freq_sr(s, r))

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
    Classe représentant un évaluateur de la qualité d'une prédiction
    """

    def __init__(self, prot, predicted):
        self.prot = prot
        self.predicted = predicted

    def Q3(self):
        """
        Mesure de qualité d'une prédiction
        => nombre de résidus correctement prédits / nombre total de résidus
        """

        corrects = 0
        for i in range(len(self.predicted)):
            if self.predicted[i] == self.prot.get_struct(i):
                corrects += 1
        
        return corrects/len(self.predicted)

    def MCC(self):
        """
        Mesure de qualité d'une prédiction prenant en compte 
        TP : prédit x alors que x
        TN : prédit pas x alors que pas x
        FP : prédit x alors que pas x
        FN : prédit pas x alors que pas x
        """

        mcc = {"H": 0, "E": 0, "C": 0}

        for x in "HEC":
            conf = {"TP": 0, "FP": 0, "FN": 0, "TN": 0} # mat de confusion
            for i in range(len(self.predicted)):

                # Si prédit x alors que x
                if self.predicted[i] == x and self.predicted[i] == \
                    self.prot.get_struct(i):
                    conf["TP"] += 1
                # Si prédit x alors que pas x
                elif self.predicted[i] == x and self.prot.get_struct(i) != x:
                    conf["FP"] += 1
                # Si prédit pas x alors que x
                elif self.predicted[i] != x and self.prot.get_struct(i) == x:
                    conf["FN"] += 1
                # Si prédit pas x alors que pas x
                elif self.predicted[i] != x and self.prot.get_struct(i) != x \
                    and self.predicted[i] == self.prot.get_struct(i):
                    conf["TN"] += 1

            # Calcul du mcc
            # (TPxTN-FPxFN) / sqrt( (TP+FP)(TP+FN)(TN+FP)(TN+FN) ) 

            # dénominateur
            den = m.sqrt((conf["TP"]+conf["FP"])*(conf["TP"]+conf["FN"])*\
                         (conf["TN"]+conf["FP"])*(conf["TN"]+conf["FN"]))
            if den != 0:
                mcc[x] += (conf["TP"]*conf["TN"]-conf["FP"]*conf["FN"]) / den
        
        return mcc


def test_prediction(prot, counter, display=False):
    """
    Teste la prédiction GOR3 sur les protéines données et les évalue
    avec des scores Q3 et MCC
    """

    average_scores = {"Q3": 0, "MCC": {"H": 0, "E": 0, "C": 0}}
    for p in prot:
        tmp_g = gorIII(counter, p)
        prediction = tmp_g.predict()
        qua = PredictionQuality(p, prediction)
        q3 = qua.Q3()
        average_scores["Q3"] += q3
        mcc = qua.MCC()
        average_scores["MCC"]["H"] += mcc["H"]
        average_scores["MCC"]["E"] += mcc["E"]
        average_scores["MCC"]["C"] += mcc["C"]

        if display:
            print(p)
            print("predicted: ")
            print(prediction)
            print("Q3: {}%, MCC: [H]: {}  [E]: {}  [C]: {}".format(\
            q3*100, mcc["H"], mcc["E"], mcc["C"]))
    
    average_scores["Q3"] /= len(prot)
    for s in average_scores["MCC"]:
        average_scores["MCC"][s] /= len(prot) # pour la moyenne

    print("Scores moyens: ", average_scores)
    

def main():
    # d = t.time()
    # p = ParserDSSP("dataset/dataset/CATH_info.txt", "dssp")
    # p.create_proteins("proteins.fasta")
    # print("Temps écoulé: {} sec".format(t.time()-d))
    # print(p.parse("dataset/dataset/dssp_test/1AVA.dssp", "C"))


    # prot1 = ParserProteins("proteins2.fasta")
    # prot1.parse()


    # ======================== TESTs pour partie 3 ============================
    # p1 = ParserDSSP("dataset/dataset/CATH_info_test.txt", "dssp_test")
    # p1.create_proteins("proteins_test.fasta")
    prot2 = ParserProteins("proteins_test.fasta") # les 5 prots
    prot2.parse()
    # c = Counter(prot2.get_proteins())
    # c.compute_frequencies()
    
    # =============================== 1ere prot des 5 de test ================
    # prot_test = prot2.get_prot(0)
    # g = gorIII(c, prot_test)
    # predicted = g.predict()
    # print(prot2.get_prot(0))
    # print("predicted: ")
    # print(predicted)
    # qua = PredictionQuality(prot2.get_prot(0), predicted)
    # print("Q3: {}%".format(qua.Q3()*100))
    # print("MCC: ", qua.MCC())

    # ========================== 3000 prots =======================
    prot3 = ParserProteins("proteins.fasta") # les 3000 prots
    prot3.parse()
    c1 = Counter(prot3.get_proteins()[:3000])
    d = t.time()
    c1.compute_frequencies()
    print("Temps écoulé pour compteurs des 3000 prots: {} sec".format(t.time()-d))
    
    # d2 = t.time()
    # test_prediction(prot3.get_proteins()[3000:], c1) # les 713 autres prots
    # print("Temps écoulé pour prédire les 713 prots: {} sec".format(t.time()-d2))

    # ======================= PARTIE 3 ==============================

    d = t.time()
    test_prediction(prot2.get_proteins(), c1, True)
    print("Temps écoulé pour prédire les 5 prots: {} sec".format(t.time()-d))

if __name__ == "__main__":
    main()