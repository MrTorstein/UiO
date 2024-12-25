class DFK():
    """
    Tar imot liste over konstanter, leser dem fra en textfil og definerer dem som klassevariabler
    """
    
    def __init__(self, liste_konstanter):
        verdinavn = []
        verdier = []
        i = 0

        with open("Konstanter.txt", "r") as innfil:
            for linje in innfil:
                if len(linje.split()) < 2:
                    pass
                elif linje.split()[0] in liste_konstanter:
                    exec("self.%s = %s"%(linje.split()[0], linje.split()[1]))