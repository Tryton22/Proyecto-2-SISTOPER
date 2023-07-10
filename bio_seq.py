from bio_structs import DNA_Codons, RNA_Codons, NUCLEOTIDE_BASE
from collections import Counter


class bio_seq:

    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        """Inicializacion de secuencia, validacion."""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Data entregada no parece ser una secuencia tipo {self.seq_type} "

    def __validate(self):
        """Corrobora si secuencia esta escrita en lenguaje nucleotidico"""
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    def get_seq_biotype(self):
        """Retorna el tipo de secuencia"""
        return self.seq_type

    def get_seq_info(self):
        '''Retorna 4 strings asociados a informacion de secuencia'''
        return f"ID: {self.label}\nSecuencia: {self.seq}\nTipo: {self.seq_type}\nLonguitud: {len(self.seq)}"

    def nucleotide_frequency(self):
        '''Cuenta los nucleotidos en una secuenia dada. Retorna un diccionario'''
        return dict(Counter(self.seq))

    def transcription(self):
        """Transcripcion de ADN a ARN. Reemplaza Timina por Uracilo"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "No es una secuencia de ADN"

    def reverse_complement(self):
        """
        Intercambia adenina con timina y guanina con citosina.
        Invierte la cadena recien generada.
        """
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        """Contenido de GC en una secuencia de ADN/ARN"""
        return round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100)

    def gc_content_subsec(self, k=20):
        """Contenido de GC en una subsecuencia de ADN/ARN de longitud k. k=20 por defecto"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(
                round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100))
        return res

    def translate_seq(self, init_pos=0):
        """Traduce una secuencia de ADN en una secuencia de aminoacidos"""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    def codon_usage(self, aminoacid):
        """Proporciona la frecuencia de cada codon que codifica un aminoacido 
        en particular en una secuencia de ADN"""
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        freqDict = dict(Counter(tmpList))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict

    def gen_reading_frames(self):
        """Genere los seis marcos de lectura de una secuencia de ADN, incluido el complemento inverso"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames

    def proteins_from_rf(self, aa_seq):
        """Calcule todas las proteinas posibles en una secuencia de aminoacidos. 
        Devuelva una lista de posibles proteinas"""
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                # DETENER la acumulacion de aminoacidos si se encontro codon de termino
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                # INICIAR la acumulacion de aminoacidos si se encontro codon de inicio
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        """Calcular todas las proteinas posibles para todos los marcos de lectura abiertos"""
        """Base de datos de busqueda de proteina: https://www.ncbi.nlm.nih.gov/nucore/NM_001185097.2"""
        """API se puede usar para extraer informacion de proteinas"""
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(
                self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res

