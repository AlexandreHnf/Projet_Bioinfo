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
            if lines[i][11] == chain and lines[i][13] not in "XZB":
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

        self.proteins = []
        
    def get_prot(self, i):
        """
        Renvoie la protéine i
        """
        return self.proteins[i]

    def parse(self):
        """
        Parse le fichier et met les protéines dans une liste
        """
        lines = self.file.readlines()
        print(len(lines))
        for i in range(1, len(lines), 3):
            self.proteins.append(lines[i] + lines[i+1].strip())
        
        self.display()

    def display(self):
        print(len(self.proteins))
        for i in self.proteins:
            print(i, end=" <<<<<<<<<<<<<<<<\n\n")
            

class gorIII:
    """
    Classe représentant un prédicteur de structure utilisant la méthode GOR III
    """
    def __init__(self, seq, struct):
        self.seq = seq
        self.struct = struct
        self.predicted = ""

    def frequency(self, s, use_s, r = None):
        """
        Calcule la fréquence d'apparition:
        - d'une structure
        - d'un acide aminé dans une structure
        - des autres structures (use_s = False)
        """

        res = 0
        if r == None: # f_s
            for i in self.struct:
                if use_s == True and i == s or use_s == False and i != s:
                    res += 1
        else:
            pass


    def info_individuelle(self, j):
        pass

    def info_directionnelle(self, j):
        pass

    def predict(self):
        pass


class PredictionQuality:
    """ 
    Classe représentant un évaluateur de la qualité d'une prédiction
    """

    def __init__(self, seq, struct, predicted):
        self.seq = seq
        self.struct = struct 
        self.predicted = predicted

    def Q3(self):
        """
        Mesure de qualité d'une prédiction
        => nombre de résidus correctement prédits / nombre total de résidus
        """

        incorrects = 0
        for i in range(len(self.predicted)):
            if self.predicted[i] != self.struct[i]:
                incorrects += 1
        
        return incorrects/len(self.predicted)

    def MCC(self, x):
        """
        Mesure de qualité d'une prédiction prenant en compte 
        TP : prédit x alors que x
        TN : prédit pas x alors que pas x
        FP : prédit x alors que pas x
        FN : prédit pas x alors que pas x
        """

        mcc = {"H": 0, "E": 0, "C": 0}

        for x in "HEC":
            conf = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
            for i in range(len(self.predicted)):

                # Si prédit x alors que x
                if self.predicted[i] == x and self.predicted[i] == self.struct[i]:
                    conf["TP"] += 1
                # Si prédit x alors que pas x
                elif self.predicted[i] == x and self.struct[i] != x:
                    conf["FP"] += 1
                # Si prédit pas x alors que x
                elif self.predicted[i] != x and self.struct[i] == x:
                    conf["FN"] += 1
                # Si prédit pas x alors que pas x
                elif self.predicted[i] != x and self.struct[i] != x \
                    and self.predicted[i] == self.struct[i]:
                    conf["TN"] += 1

            # Calcul du mcc
            # (TPxTN-FPxFN) / sqrt( (TP+FP)(TP+FN)(TN+FP)(TN+FN) ) 

            m_c_c = (conf["TP"]*conf["TN"]-conf["FP"]*conf["FN"]) / \
            m.sqrt((conf["TP"]+conf["FP"])*(conf["TP"]+conf["FN"])*\
            (conf["TN"]+conf["FP"])*(conf["TN"]+conf["FN"]))

            mcc[x] += m_c_c

        print(mcc)

        
def main():
    # d = t.time()
    # p = ParserDSSP("dataset/dataset/CATH_info.txt", "dssp")
    # p.create_proteins("proteins3.txt")
    # print(t.time()-d)
    # print(p.parse("dataset/dataset/dssp_test/1AVA.dssp", "C"))


    # prot1 = ParserProteins("proteins2.txt")
    # prot1.parse()





    # ======================== TESTs pour partie 3 ============================
    # p1 = ParserDSSP("dataset/dataset/CATH_info_test.txt", "dssp_test")
    # p1.create_proteins("proteins_test.txt")
    # prot2 = ParserProteins("proteins_test.txt")
    # prot2.parse()
    # g2 = gorIII()


if __name__ == "__main__":
    main()