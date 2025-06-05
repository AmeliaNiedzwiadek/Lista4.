class BioSequence:
    
    """
    Bazowa klasa reprezentująca sekwencję biologiczną (DNA, RNA, białko).
    
    @param identifier: Identyfikator sekwencji.
    @param data: Sekwencja biologiczna (np. DNA).
    @param valid_chars: Zbiór dozwolonych znaków.
    
    @raises ValueError: Jeśli sekwencja zawiera nieprawidłowe znaki.
    """

    def __init__(self, identifier: str, data: str, valid_chars: set):
        
        self.identifier = identifier
        self.data = data.upper()
        self.valid_chars = valid_chars
        if not all(char in valid_chars for char in self.data):
            raise ValueError(f"Nieprawidłowe znaki w sekwencji: {self.data}")

    def length(self):
        """
        Zwraca długość sekwencji.
        
        @return: Liczba znaków w sekwencji.
        """
        return len(self.data)

    def __str__(self):
        
        """
        Zwraca reprezentację sekwencji w formacie FASTA.
        
        @return: Sekwencja w formacie FASTA.
        
        """
        
        return f">{self.identifier}\n{self.data}"

    def mutate(self, position: int, value: str):
        
        """
        Modyfikuje znak sekwencji na podanej pozycji.
        
        @param position: Indeks znaku do modyfikacji.
        @param value: Nowy znak (musi być dozwolony).
        
        @raises IndexError: Jeśli pozycja jest poza zakresem.
        @raises ValueError: Jeśli znak jest nieprawidłowy.
        
        """
        if position < 0 or position >= len(self.data):
            raise IndexError("Pozycja poza zakresem.")
        if value not in self.valid_chars:
            raise ValueError(f"Znak {value} nie jest dozwolony.")
        self.data = self.data[:position] + value + self.data[position+1:]

    def find_motif(self, motif: str) -> int:
        """
        Szuka motywu w sekwencji.
        
        @param motif: Motyw (podciąg) do wyszukania.
        @return: Indeks pierwszego wystąpienia lub -1.
        
        """
        return self.data.find(motif)


class DNASequence(BioSequence):
    
    """
    Reprezentuje sekwencję DNA.
    
    @param identifier: Identyfikator sekwencji.
    @param data: Sekwencja DNA.
    
    """

    def __init__(self, identifier: str, data: str):
        
        super().__init__(identifier, data, {'A', 'T', 'G', 'C'})

    def complement(self) -> str:
        
        """
        Zwraca komplementarną, odwróconą sekwencję DNA.
        
        @return: Komplementarna sekwencja DNA.
        """
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        comp = ''.join(complement_map[base] for base in reversed(self.data))
        return comp

    def transcribe(self):
        
        """
        Transkrybuje DNA do RNA.
        
        @return: Obiekt RNASequence.
        
        """
        complement_map = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
        rna_data = ''.join(complement_map[base] for base in reversed(self.data))
        return RNASequence(self.identifier, rna_data)


class RNASequence(BioSequence):
    
    """
    Reprezentuje sekwencję RNA.
    
    @param identifier: Identyfikator sekwencji.
    @param data: Sekwencja RNA.
    
    """

    def __init__(self, identifier: str, data: str):
        
        super().__init__(identifier, data, {'A', 'U', 'G', 'C'})

    def complement(self) -> str:
        
        """
        Zwraca komplementarną, odwróconą sekwencję RNA.
        
        @return: Komplementarna sekwencja RNA.
        
        """
        complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        comp = ''.join(complement_map[base] for base in reversed(self.data))
        return comp

    def transcribe(self):
        """
        Tłumaczy RNA na sekwencję białkową (do kodonu STOP).
        
        @return: Obiekt klasy ProteinSequence.
        
        """
        codon_table = {
            'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'UAU': 'Y', 'UAC': 'Y', 'UGU': 'C', 'UGC': 'C',
            'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
            'UAA': '*', 'UAG': '*', 'UGA': '*'
        }

        protein = ''
        for i in range(0, len(self.data) - 2, 3):
            codon = self.data[i:i+3]
            amino_acid = codon_table.get(codon, 'X')
            if amino_acid == '*':
                break
            protein += amino_acid
        return ProteinSequence(self.identifier, protein)


class ProteinSequence(BioSequence):
    
    """
    Reprezentuje sekwencję białkową.
    
    @param identifier: Identyfikator sekwencji.
    @param data: Sekwencja aminokwasów.
    """

    def __init__(self, identifier: str, data: str):
        
        allowed = set("ACDEFGHIKLMNPQRSTVWY")
        super().__init__(identifier, data, allowed)


# ===================== TESTY ======================

def assert_equals(expected, actual, test_name):
    result = "OK" if expected == actual else f"FAIL — Expected: {expected}, Actual: {actual}"
    print(f"{test_name}: {result}")

def main():
    dna = DNASequence("DNA1", "ATGC")
    test_basic(dna)
    test_mutations(dna)
    test_motifs(dna)
    test_transcription_and_complement(dna)

def test_basic(dna):
    print("\n=== Test podstawowy ===")
    print(dna)
    print("Długość:", dna.length)

def test_mutations(dna):
    print("\n=== Testowanie mutacji ===")

    print("-> Mutacja pozycji 2 na 'T'")
    try:
        dna.mutate(2, 'T')
        print("Po mutacji:", dna)
        assert_equals("ATTC", dna.data, "Mutacja poprawna")
    except Exception as e:
        print(" Nieoczekiwany błąd przy mutacji pozycji 2:", e)

    print("-> Próba mutacji na pozycji 10 ")
    try:
        dna.mutate(10, 'A')
        print(" Błąd! Oczekiwano wyjątku, ale mutacja się powiodła.")
    
    except Exception as e:
        print(f" Niespodziewany typ wyjątku: {e}")

    print("-> Próba mutacji pozycji 0 na niedozwolony znak 'X'")
    try:
        dna.mutate(0, 'X')
        print(" Błąd! Oczekiwano wyjątku, ale mutacja się powiodła.")
    except ValueError as e:
        print(f" {e}")
    except Exception as e:
        print(f" Niespodziewany typ wyjątku: {e}")


def test_motifs(dna):
    print("\n=== Testowanie motywów ===")
    pos = dna.find_motif("TT")
    print("Motyw 'TT' na pozycji:", pos)
    assert_equals(1, pos, "Motyw TT")

    no_pos = dna.find_motif("GGG")
    print("Motyw 'GGG' na pozycji:", no_pos)
    assert_equals(-1, no_pos, "Brak motywu")

def test_transcription_and_complement(dna):
    print("\n=== Transkrypcja i komplement ===")
    comp = dna.complement()
    print("Komplement DNA:", comp)
    assert_equals("GAAT", comp, "Komplement DNA (odwrócony)")

    rna = dna.transcribe()
    print(f"RNA (z {dna.data}): {rna.data}")
    assert_equals("GAAU", rna.data, "Transkrypcja RNA")

    rna_comp = rna.complement()
    print("Komplement RNA:", rna_comp)
    assert_equals("AUUC", rna_comp, "Komplement RNA")

    protein = rna.transcribe()
    print("Białko:", protein)
    print("Sekwencja białka:", protein.data)

if __name__ == "__main__":
    main()

#Podczas tworzenia kodu korzystałam z następujących źródeł:
#https://docs.python.org/3/tutorial/classes.html
#https://docs.python.org/3/library/functions.html#super
#https://docs.python.org/3/library/functions.html#all
#https://docs.python.org/3/library/stdtypes.html#str.find
#https://docs.python.org/3/library/stdtypes.html#text-sequence-type-str
#https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#https://realpython.com/python3-object-oriented-programming/
#ChatGPT
#https://docs.python.org/3/tutorial/inputoutput.html#formatted-string-literals
#https://docs.python.org/3/tutorial/errors.html
#https://docs.python.org/3/reference/compound_stmts.html#the-try-statement
#https://realpython.com/python-testing/#writing-basic-tests-with-assert
#https://realpython.com/python-main-function/
#https://docs.python.org/3/reference/expressions.html#comparisons
#https://github.com/AmeliaNiedzwiadek/lista-2
