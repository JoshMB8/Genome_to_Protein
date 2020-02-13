import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Enter a valid file name including a valid format. e.g .fasta or .txt. Input file must include '>' if first line contains file info e.g. I.claudius genome sequence must be structed as >I.claudius genome sequence.")
parser.add_argument("-o", "--output", help="Enter a name for your new file and include a valid format for that file. e.g. .fasta or .txt", default="Proteins.fasta")
parser.add_argument("-s", "--species", help="Enter the name of your species or an abbreviation of your choosing as an identifier. e.g. I_claudius or IC", type=str, default="species")
parser.add_argument("-f", "--frame", help="Enter the frame you wish to retrieve. Valid options are 1, 2, and 3 for the forward 3 frames, -1, -2, and -3 for the reverse 3 frames, and 0 for all frame", type=int, default=0)
parser.add_argument("-m", "--morfl", help="Enter a value for minimum ORF length. Minimum ORF length is the number of amino-acids e.g. 100", type=int, default=50)
args = parser.parse_args()

fileIn = args.input              # Assigns the argument (input file name) input by user in command line to 'fileIn'. e.g. 'genome.fasta'
fileOut = args.output            # Assigns the argument (output file name) input by the user in command line to 'fileOut'. e.g. 'species_proteins.fasta'
spec = args.species              # Assigns the argument (species name/ abbreviation) input by the user in command line to 'spec'. e.g 'IC'
frame = args.frame               # Assigns the argument (frame) input by the user in command line to 'frame'. e.g '1' for forward frame 1 only
min_orf_length = args.morfl      # Assigns the argument (value of min orf length) input by the user in command line to min_orf_length. 'e.g. 100'

    
def readFile(fileName):                        # Function that reads in the file.
    Genome = open(fileName) # Read in Genome Sequence file.
    gen_seq = "" # Create string for full genome sequence to be added to.

    for seq in Genome:           
        seq_cmb = seq.rstrip()   # Lines are then joined using .rstrip() to remove the 'new line' auto-formatting.
        if seq[0] != ">":        # If sequence starts with '>' it is ignored and not added to the growing string.
            gen_seq += seq_cmb.upper()    # The sequence is then added to 'gen_seq' line-by-line until all are incorporated in uppercase.
    return gen_seq

genetic_code = {
    "AAA" : "K", "AAC" : "N", "AAG" : "K", "AAT" : "N",
    "ACA" : "T", "ACC" : "T", "ACG" : "T", "ACT" : "T",
    "AGA" : "R", "AGC" : "S", "AGG" : "R", "AGT" : "S",
    "ATA" : "I", "ATC" : "I", "ATG" : "M", "ATT" : "I",
    "CAA" : "Q", "CAC" : "H", "CAG" : "Q", "CAT" : "H",
    "CCA" : "P", "CCC" : "P", "CCG" : "P", "CCT" : "P",
    "CGA" : "R", "CGC" : "R", "CGG" : "R", "CGT" : "R",
    "CTA" : "L", "CTC" : "L", "CTG" : "L", "CTT" : "L",
    "GAA" : "E", "GAC" : "D", "GAG" : "E", "GAT" : "D",
    "GCA" : "A", "GCC" : "A", "GCG" : "A", "GCT" : "A",
    "GGA" : "G", "GGC" : "G", "GGG" : "G", "GGT" : "G",
    "GTA" : "V", "GTC" : "V", "GTG" : "V", "GTT" : "V",
    "TAA" : "*", "TAC" : "Y", "TAG" : "*", "TAT" : "Y",
    "TCA" : "S", "TCC" : "S", "TCG" : "S", "TCT" : "S",
    "TGA" : "*", "TGC" : "C", "TGG" : "W", "TGT" : "C",
    "TTA" : "L", "TTC" : "F", "TTG" : "L", "TTT" : "F" }     # Dictionary of genetic code, codon to amino-acid translations.
def revComp(seq):                              # Function to generate the reverse compliment of the input sequence.
    transTab = str.maketrans('ATGC', 'TACG')
    comp = seq.translate(transTab)
    revcomp = comp[::1]
    return revcomp

def dna_to_prot(seq, revseq, spec, min_orf_length, frame_start):  # Function that will read through the file, separate open reading frames and translate to protein sequences. revseq argument tells function if the input sequence is the reverse compliment. spec argument allows user to choose name of spcies in output. 
    frame = 0
    orf = ""
    n_orf = 0
    prot_seq = ""
    for pos in range(frame_start, len(seq), 3):
        codon = seq[pos:pos+3]
        if pos+3 > len(seq):
            break
        if(codon == "ATG" and not frame):     # Reads through the sequence codon by codon to find a start codon 'ATG'.
            orf = genetic_code[codon]         # If found ORF if started, translated to amino-acid code via dictionary, and recorded.
            frame += 1
        elif(codon == "TGA" or codon == "TAA" or codon == "TAG"): # If a stop codon is found in the sequence then the ORF stops recording and moves to the next ORF.
            if frame:
                orf += genetic_code[codon]    # Translates codons to amino-acid code, and adds to the growing ORF.
                if len(orf) > min_orf_length:
                    n_orf += 1
                    if revseq == False:       # If the sequence is the reverse compliment (== True) then print with "-" frame_start +1 to give frames -1, -2, -3.
                        prot_seq += ">%sGenome|Frame:%s|ORF:%s|Length:%s|\n%s\n" % (spec, frame_start+1, n_orf, len(orf), orf)
                    elif revseq == True:
                        prot_seq += ">%sGenome|Frame:-%s|ORF:%s|Length:%s|\n%s\n" % (spec, frame_start+1, n_orf, len(orf), orf)
                frame = 0                     # Stops recording codons to sequence.
        elif frame:                           # Any other codon combination is added to the growing sequence, including if another 'ATG' is found withing an ORF.
            try:
                orf += genetic_code[codon]    # Translates codons to amino-acid code.
            except KeyError:
                frame = 0                     # If genetic code contains anything not in dictionary e.g 'NGC' doesn't add and instead moves to next sequence.
    return prot_seq

def dna2prot_allframes(seq, spec, frame, min_orf_length=50):   # Function that retrieves the translated sequences for all 6 open reading frames, of a certain length (the default is 50 amino acids), and concatinates them to 'Frames'. IF provided a valid input for argument 'frame' will only return frame associated with that input.
    Frames = ""
    if frame == 1:
        Frames += dna_to_prot(seq, False, spec, min_orf_length, 0)
        return Frames
    elif frame == 2:
        Frames += dna_to_prot(seq, False, spec, min_orf_length, 1)
        return Frames
    elif frame == 3:
        Frames += dna_to_prot(seq, False, spec, min_orf_length, 2)
        return Frames
    elif frame == -1:
        Frames += dna_to_prot(revComp(seq), True, spec, min_orf_length, 0)
        return Frames
    elif frame == -2:
        Frames += dna_to_prot(revComp(seq), True, spec, min_orf_length, 1)
        return Frames
    elif frame == -3:
        Frames += dna_to_prot(revComp(seq), True, spec, min_orf_length, 2)
        return Frames
    elif frame == 0:
        Frames += dna_to_prot(seq, False, spec, min_orf_length, 0)
        Frames += dna_to_prot(seq, False, spec, min_orf_length, 1)
        Frames += dna_to_prot(seq, False, spec, min_orf_length, 2)
        Frames += dna_to_prot(revComp(seq), True, spec, min_orf_length, 0)
        Frames += dna_to_prot(revComp(seq), True, spec, min_orf_length, 1)
        Frames += dna_to_prot(revComp(seq), True, spec, min_orf_length, 2)
        return Frames

def main(fileIn, fileOut, spec, frame, min_orf_length):             # Main Function which takes in file name, e.g. genome.fasta, output file name, species name, frame number, and a minimum orf length value as arguments and outputs the protein sequences, for all orf's, with a orf length defined by the user into a fasta file.
    gen_seq = readFile(fileIn)
    Frames = dna2prot_allframes(gen_seq, spec, frame, min_orf_length)
    proteome = open(fileOut, "w")          # Opens a file, with name defined by user, to write the output to.
    proteome.write(Frames)                 # Writes AllFrames into the Proteins file.
    proteome.close                         # Closes the Proteins file so no more is written to it.

main(fileIn, fileOut, spec, frame, min_orf_length)         # Calls the main function on 'fileIn', 'fileOut', 'spec', 'frame' and 'min_orf_length' as defined earlier.