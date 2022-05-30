/*--------------------------------------------------------------------------
 Program name     : SatuRation.cpp
 
 Version          : 1.00
 
 Author           : Lars S Jermiin
 
 Institutions     : Australian National University
                    Research School of Biology
                    Acton, ACT 2601, Australia
 
                    Univerity College Dublin
                    School of Biology & Environmental Science
                    Belfield, Dublin 4, Ireland
 
 Emails           : lars.jermiin [at] anu.edu.au
                    lars.jermiin [at] ucd.ie

 URL              : https://github.com/lsjermiin/SatuRation_v1.0
 
 Date begun       : 18 April, 2018
 
 Date modified    : 28 May, 2022
 
 Copyright        : Copyright © 2019-2022 Lars Sommer Jermiin.
                    All rights reserved.
 
 Responsibility   : The copyright holders take no legal responsibility for
                    the correctness of results obtained using this program.
 
 Summary          : SatuRation estimates the degree of satuation observed
                    between two sequences of nucleotides or amino acids.
 
                    Sequences must be stored in the FASTA format.
 
                    Characters are converted to integers to speed up the
                    program.

 Data types       : Sequences can be read a strings of one, two, or three
                    nucleotides, strings of 10- or 14-state genotypes, or 
                    as strings of amino acids.

 Nucleotides______
         Singlets : Alphabet: [A,C.G,T/U,-] = [0,1,2,3,4].
 
                    Ambiguous characters (i.e., ?, N, B, D, H, K, M, R, S,
                    V, W and Y) are treated as if they were alignment gaps
                    (-) (i.e., as missing data).
 
          Duplets : Alphabet: [AA,AC,AG,AT/AU,CA,CC,CG,CT/CU,GA,GC,GG,
                    GT/GU,TA/UA,TC/UC,TG/UG,TT/TU/UT/UU,--] = [0,1,2,3,4,5,
                    6,7,8,9,10,11,12,16,14,15,16]
 
                    Ambiguous characters (i.e., ?, N, B, D, H, K, M, R, S,
                    V, W and Y) are treated as if they were alignment gaps
                    (-) (i.e., as missing data).
 
         Triplets : Alphabet: All codons
 
 Genotypes________
        10 states : Alphabet: [A,C,G,K,M,R,S,T/U,W,Y,-] =
                    [0,1,2,3,4,5,6,7,8,9,10].
 
                    Ambiguous characters (i.e., ?, N, B, D, G and V) are
                    treated as if they were alignment gaps (-) (i.e., as
                    missing data).
 
        14 states : Alphabet: [A,C,G,T/U,K,M,R,S,W,Y,B,D,H,V,-] =
                    [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14].
 
                    Ambiguous characters (i.e., ? and N) are treated as if
                    they were alignment gaps (-) (i.e., as missing data).
 
 Amino acids      : Alphabet: [A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,-] =
                    [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
  
                    Ambiguous characters (i.e., ?, X and Z) are treated as
                    if they were alignment gaps (-) (i.e., as missing data).
 
 Constant sites   : A site is defined as a constant site if it contains the
                    same character in every sequence at a given site, or if
                    that character in some of the sequences is replaced by
                    an alignment gap (i.e., '-').
 
 Manuscript       : Jermiin et al. (2019)
 
                    Zhang...
 
 ----------------------------------------------------------------------------*/

#include <cctype>
#include <cmath>
#include <limits>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#define SQR(a) ((a) * (a))

// The following variables are declared externally because they
// are needed in different functions
const unsigned TWO(2);          // for 2-state alphabet (recoded DNA)
const unsigned THREE(3);        // for 3-state alphabet (recoded DNA)
const unsigned FOUR(4);         // for 4-state alphabet (DNA)
const unsigned SIX(6);          // for 6-state alphabet (recoded amino acids)
const unsigned TEN(10);         // for 10-state alphabet (genotype data)
const unsigned FOURTEEN(14);    // for 14-state alphabet (genotype data)
const unsigned SIXTEEN(16);     // for 16-state alphabet (subset of codons)
const unsigned TWENTY(20);      // for 20-state alphabet (amino acids)
const unsigned SIXTYFOUR(64);   // for 64-state alphabet (codons)
const unsigned max_array(65);   // corresponding to amino acids & gaps
std::string sites;                  // string controlling whether a site is used or not
std::vector<std::string> taxon;           // 2D container for sequence names
std::vector<std::vector<int> > alignment; // 2D container for sequence data



// This function translates a string of characters into a vector of integers
std::vector<int> Translator(unsigned datatype, std::string seq) {
    int unit; // integer for singlet, duplet or triplet (codon)
    std::string duplet(""), triplet(""); // strings for dinucleotides and codons
    std::vector<int> seq_data;
    
    switch (datatype) {
        case 1: // Nucleotides (A|C|G|T)
             for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 2: // Nucleotides (C|T|R)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(2); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'R': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 3: // Nucleotides (A|G|Y)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(2); break;
                    case 'U': seq_data.push_back(2); break;
                    case 'Y': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 4: // Nucleotides (A|T|S)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'C': seq_data.push_back(2); break;
                    case 'S': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 5: // Nucleotides (C|G|W)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(2); break;
                    case 'U': seq_data.push_back(2); break;
                    case 'W': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 6: // Nucleotides (A|C|K)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(2); break;
                    case 'U': seq_data.push_back(2); break;
                    case 'K': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 7: // Nucleotides (G|T|M)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'G': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(2); break;
                    case 'C': seq_data.push_back(2); break;
                    case 'M': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 8: // Nucleotides (K|M)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'G': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(0); break;
                    case 'U': seq_data.push_back(0); break;
                    case 'K': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'M': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 9: // Nucleotides (R|Y)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(0); break;
                    case 'R': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'Y': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 10: // Nucleotides (S|W)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(0); break;
                    case 'R': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'Y': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 11: // Nucleotides (A|B)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'B': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 12: // Nucleotides (C|D)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'D': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 13: // Nucleotides (G|H)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'G': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'H': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 14: // Nucleotides (T|V)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'T': seq_data.push_back(0); break;
                    case 'U': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'V': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 15: // Di-nucleotides (AA|AC|..|TG|TT)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 2) {
                for (std::string::size_type j = i; j != i + 2; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': duplet.push_back('0'); break;
                        case 'C': duplet.push_back('1'); break;
                        case 'G': duplet.push_back('2'); break;
                        case 'T': duplet.push_back('3'); break;
                        case 'U': duplet.push_back('3'); break;
                        default : duplet.push_back('4'); break;
                    }
                }
                unit = stoi(duplet);
                duplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 16: // Di-nucleotides 1st position (A|C|G|T)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 2) {
                switch (toupper(seq[i])) {
                    case 'A' : seq_data.push_back(0); break; // 1st A
                    case 'C' : seq_data.push_back(1); break; // 1st C
                    case 'G' : seq_data.push_back(2); break; // 1st G
                    case 'T' : seq_data.push_back(3); break; // 1st T
                    default  : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 17: // Di-nucleotides 2nd position (A|C|G|T)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 2) {
                switch (toupper(seq[i+1])) {
                    case 'A' : seq_data.push_back(0); break; // 2nd A
                    case 'C' : seq_data.push_back(1); break; // 2nd C
                    case 'G' : seq_data.push_back(2); break; // 2nd G
                    case 'T' : seq_data.push_back(3); break; // 2nd T
                    default  : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 18: // Codons (AAA|AAC|...|TTG|TTT)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 000: seq_data.push_back(0);  break; // AAA
                    case 001: seq_data.push_back(1);  break; // AAC
                    case 002: seq_data.push_back(2);  break; // AAG
                    case 003: seq_data.push_back(3);  break; // AAT
                    case 010: seq_data.push_back(4);  break; // ACA
                    case 011: seq_data.push_back(5);  break; // ACC
                    case 012: seq_data.push_back(6);  break; // ACG
                    case 013: seq_data.push_back(7);  break; // ACT
                    case 020: seq_data.push_back(8);  break; // AGA
                    case 021: seq_data.push_back(9);  break; // AGC
                    case 022: seq_data.push_back(10); break; // AGG
                    case 023: seq_data.push_back(11); break; // AGT
                    case 030: seq_data.push_back(12); break; // ATA
                    case 031: seq_data.push_back(13); break; // ATC
                    case 032: seq_data.push_back(14); break; // ATG
                    case 033: seq_data.push_back(15); break; // ATT
                    case 100: seq_data.push_back(16); break; // CAA
                    case 101: seq_data.push_back(17); break; // CAC
                    case 102: seq_data.push_back(18); break; // CAG
                    case 103: seq_data.push_back(19); break; // CAT
                    case 110: seq_data.push_back(20); break; // CCA
                    case 111: seq_data.push_back(21); break; // CCC
                    case 112: seq_data.push_back(22); break; // CCG
                    case 113: seq_data.push_back(23); break; // CCT
                    case 120: seq_data.push_back(24); break; // CGA
                    case 121: seq_data.push_back(25); break; // CGC
                    case 122: seq_data.push_back(26); break; // CGG
                    case 123: seq_data.push_back(27); break; // CGT
                    case 130: seq_data.push_back(28); break; // CTA
                    case 131: seq_data.push_back(29); break; // CTC
                    case 132: seq_data.push_back(30); break; // CTG
                    case 133: seq_data.push_back(31); break; // CTT
                    case 200: seq_data.push_back(32); break; // GAA
                    case 201: seq_data.push_back(33); break; // GAC
                    case 202: seq_data.push_back(34); break; // GAG
                    case 203: seq_data.push_back(35); break; // GAT
                    case 210: seq_data.push_back(36); break; // GCA
                    case 211: seq_data.push_back(37); break; // GCC
                    case 212: seq_data.push_back(38); break; // GCG
                    case 213: seq_data.push_back(39); break; // GCT
                    case 220: seq_data.push_back(40); break; // GGA
                    case 221: seq_data.push_back(41); break; // GGC
                    case 222: seq_data.push_back(42); break; // GGG
                    case 223: seq_data.push_back(43); break; // GGT
                    case 230: seq_data.push_back(44); break; // GTA
                    case 231: seq_data.push_back(45); break; // GTC
                    case 232: seq_data.push_back(46); break; // GTG
                    case 233: seq_data.push_back(47); break; // GTT
                    case 300: seq_data.push_back(48); break; // TAA
                    case 301: seq_data.push_back(49); break; // TAC
                    case 302: seq_data.push_back(50); break; // TAG
                    case 303: seq_data.push_back(51); break; // TAT
                    case 310: seq_data.push_back(52); break; // TCA
                    case 311: seq_data.push_back(53); break; // TCC
                    case 312: seq_data.push_back(54); break; // TCG
                    case 313: seq_data.push_back(55); break; // TCT
                    case 320: seq_data.push_back(56); break; // TGA
                    case 321: seq_data.push_back(57); break; // TGC
                    case 322: seq_data.push_back(58); break; // TGG
                    case 323: seq_data.push_back(59); break; // TGT
                    case 330: seq_data.push_back(60); break; // TTA
                    case 331: seq_data.push_back(61); break; // TTC
                    case 332: seq_data.push_back(62); break; // TTG
                    case 333: seq_data.push_back(63); break; // TTT
                    default:  seq_data.push_back(64); break; // In case of other characters
                }
            }
            break;
        case 19: // Codons 1st + 2nd positions (AA|AC|...|TG|TT)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(2,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 20: // Codons 1st + 3rd positions (AA|AC|...|TG|TT)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(1,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 21: // Codons 2nd + 3rd positions (AA|AC|...|TG|TT)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(0,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 22: // Codons 1st + 2nd positions (A|C|G|T)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    if (j == i || j == i + 1) {
                        switch (toupper(seq[j])) {
                            case 'A': seq_data.push_back(0); break;
                            case 'C': seq_data.push_back(1); break;
                            case 'G': seq_data.push_back(2); break;
                            case 'T': seq_data.push_back(3); break;
                            case 'U': seq_data.push_back(3); break;
                            default : seq_data.push_back(4); break;
                        }
                    }
                }
            }
            break;
        case 23: // Codons 1st + 3rd positions (A|C|G|T)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    if (j == i || j == i + 2) {
                        switch (toupper(seq[j])) {
                            case 'A': seq_data.push_back(0); break;
                            case 'C': seq_data.push_back(1); break;
                            case 'G': seq_data.push_back(2); break;
                            case 'T': seq_data.push_back(3); break;
                            case 'U': seq_data.push_back(3); break;
                            default : seq_data.push_back(4); break;
                        }
                    }
                }
            }
            break;
        case 24: // Codons 2nd + 3rd positions (A|C|G|T)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    if (j == i + 1 || j == i + 2) {
                        switch (toupper(seq[j])) {
                            case 'A': seq_data.push_back(0); break;
                            case 'C': seq_data.push_back(1); break;
                            case 'G': seq_data.push_back(2); break;
                            case 'T': seq_data.push_back(3); break;
                            case 'U': seq_data.push_back(3); break;
                            default : seq_data.push_back(4); break;
                        }
                    }
                }
            }
            break;
        case 25: // Codons 1st position (A|C|G|T)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(2,1);
                triplet.erase(1,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case   0: seq_data.push_back(0); break;
                    case   1: seq_data.push_back(1); break;
                    case   2: seq_data.push_back(2); break;
                    case   3: seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 26: // Codons 2nd position (A|C|G|T)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(2,1);
                triplet.erase(0,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case   0: seq_data.push_back(0); break;
                    case   1: seq_data.push_back(1); break;
                    case   2: seq_data.push_back(2); break;
                    case   3: seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 27: // Codons 3rd position (A|C|G|T)
            for (std::string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (std::string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(1,1);
                triplet.erase(0,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case   0: seq_data.push_back(0); break;
                    case   1: seq_data.push_back(1); break;
                    case   2: seq_data.push_back(2); break;
                    case   3: seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 28: // 10-state genotype data
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    case 'K': seq_data.push_back(4); break;
                    case 'M': seq_data.push_back(5); break;
                    case 'R': seq_data.push_back(6); break;
                    case 'Y': seq_data.push_back(7); break;
                    case 'S': seq_data.push_back(8); break;
                    case 'W': seq_data.push_back(9); break;
                    default : seq_data.push_back(10);break; // In case of other characters
                }
            }
            break;
        case 29: // 14-state genotype data
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    case 'K': seq_data.push_back(4); break;
                    case 'M': seq_data.push_back(5); break;
                    case 'R': seq_data.push_back(6); break;
                    case 'Y': seq_data.push_back(7); break;
                    case 'S': seq_data.push_back(8); break;
                    case 'W': seq_data.push_back(9); break;
                    case 'B': seq_data.push_back(10);break;
                    case 'D': seq_data.push_back(11);break;
                    case 'H': seq_data.push_back(12);break;
                    case 'V': seq_data.push_back(13);break;
                    default : seq_data.push_back(14);break; // In case of other characters
                }
            }
            break;
        case 30: // amino acids (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'D': seq_data.push_back(2); break;
                    case 'E': seq_data.push_back(3); break;
                    case 'F': seq_data.push_back(4); break;
                    case 'G': seq_data.push_back(5); break;
                    case 'H': seq_data.push_back(6); break;
                    case 'I': seq_data.push_back(7); break;
                    case 'K': seq_data.push_back(8); break;
                    case 'L': seq_data.push_back(9); break;
                    case 'M': seq_data.push_back(10);break;
                    case 'N': seq_data.push_back(11);break;
                    case 'P': seq_data.push_back(12);break;
                    case 'Q': seq_data.push_back(13);break;
                    case 'R': seq_data.push_back(14);break;
                    case 'S': seq_data.push_back(15);break;
                    case 'T': seq_data.push_back(16);break;
                    case 'V': seq_data.push_back(17);break;
                    case 'W': seq_data.push_back(18);break;
                    case 'Y': seq_data.push_back(19);break;
                    default : seq_data.push_back(20);break; // In case of other characters
                }
            }
            break;
        default: // Dayhoff-6 (AGPST|DENQ|HKR|MIVL|WFY|C)
            for (std::string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(0); break;
                    case 'P': seq_data.push_back(0); break;
                    case 'S': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(0); break;
                    case 'D': seq_data.push_back(1); break;
                    case 'E': seq_data.push_back(1); break;
                    case 'N': seq_data.push_back(1); break;
                    case 'Q': seq_data.push_back(1); break;
                    case 'H': seq_data.push_back(2); break;
                    case 'K': seq_data.push_back(2); break;
                    case 'R': seq_data.push_back(2); break;
                    case 'M': seq_data.push_back(3); break;
                    case 'I': seq_data.push_back(3); break;
                    case 'L': seq_data.push_back(3); break;
                    case 'V': seq_data.push_back(3); break;
                    case 'F': seq_data.push_back(4); break;
                    case 'W': seq_data.push_back(4); break;
                    case 'Y': seq_data.push_back(4); break;
                    case 'C': seq_data.push_back(5); break;
                    default : seq_data.push_back(6); break; // In case of other characters
                }
            }
            break;
    }
    return(seq_data);
}


// This function translates a vector of intergers into a string of characters
std::string Back_translator(unsigned datatype, std::vector<int> seq_data) {
    unsigned number(0);
    std::string str("");
    
    switch (datatype) {
        case 1:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('G'); break;
                        case 3:  str.push_back('T'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 2:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('C'); break;
                        case 1:  str.push_back('T'); break;
                        case 2:  str.push_back('R'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 3:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('G'); break;
                        case 2:  str.push_back('Y'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 4:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('T'); break;
                        case 2:  str.push_back('S'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 5:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('C'); break;
                        case 1:  str.push_back('G'); break;
                        case 2:  str.push_back('W'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 6:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('K'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 7:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('G'); break;
                        case 1:  str.push_back('T'); break;
                        case 2:  str.push_back('M'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 8:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('K'); break;
                        case 1:  str.push_back('M'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 9:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('R'); break;
                        case 1:  str.push_back('Y'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 10:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('S'); break;
                        case 1:  str.push_back('W'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 11:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('B'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 12:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('C'); break;
                        case 1:  str.push_back('D'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 13:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('G'); break;
                        case 1:  str.push_back('H'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 14:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('T'); break;
                        case 1:  str.push_back('V'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 15:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 00: str = str + "AA"; break;
                        case 01: str = str + "AC"; break;
                        case 02: str = str + "AG"; break;
                        case 03: str = str + "AT"; break;
                        case 10: str = str + "CA"; break;
                        case 11: str = str + "CC"; break;
                        case 12: str = str + "CG"; break;
                        case 13: str = str + "CT"; break;
                        case 20: str = str + "GA"; break;
                        case 21: str = str + "GC"; break;
                        case 22: str = str + "GG"; break;
                        case 23: str = str + "GT"; break;
                        case 30: str = str + "TA"; break;
                        case 31: str = str + "TC"; break;
                        case 32: str = str + "TG"; break;
                        case 33: str = str + "TT"; break;
                        default: str = str + "--"; break;
                    }
                }
            }
            break;
        case 16:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('G'); break;
                        case 3:  str.push_back('T'); break;
                        default: str.push_back('-'); break;
                    }
                }
                number = seq_data[i];
            }
            break;
        case 17:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('G'); break;
                        case 3:  str.push_back('T'); break;
                        default: str.push_back('-'); break;
                    }

                }
            }
            break;
        case 18:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 000: str = str + "AAA"; break;
                        case 001: str = str + "AAC"; break;
                        case 002: str = str + "AAG"; break;
                        case 003: str = str + "AAT"; break;
                        case 010: str = str + "ACA"; break;
                        case 011: str = str + "ACC"; break;
                        case 012: str = str + "ACG"; break;
                        case 013: str = str + "ACT"; break;
                        case 020: str = str + "AGA"; break;
                        case 021: str = str + "AGC"; break;
                        case 022: str = str + "AGG"; break;
                        case 023: str = str + "AGT"; break;
                        case 030: str = str + "ATA"; break;
                        case 031: str = str + "ATC"; break;
                        case 032: str = str + "ATG"; break;
                        case 033: str = str + "ATT"; break;
                        case 100: str = str + "CAA"; break;
                        case 101: str = str + "CAC"; break;
                        case 102: str = str + "CAG"; break;
                        case 103: str = str + "CAT"; break;
                        case 110: str = str + "CCA"; break;
                        case 111: str = str + "CCC"; break;
                        case 112: str = str + "CCG"; break;
                        case 113: str = str + "CCT"; break;
                        case 120: str = str + "CGA"; break;
                        case 121: str = str + "CGC"; break;
                        case 122: str = str + "CGG"; break;
                        case 123: str = str + "CGT"; break;
                        case 130: str = str + "CTA"; break;
                        case 131: str = str + "CTC"; break;
                        case 132: str = str + "CTG"; break;
                        case 133: str = str + "CTT"; break;
                        case 200: str = str + "GAA"; break;
                        case 201: str = str + "GAC"; break;
                        case 202: str = str + "GAG"; break;
                        case 203: str = str + "GAT"; break;
                        case 210: str = str + "GCA"; break;
                        case 211: str = str + "GCC"; break;
                        case 212: str = str + "GCG"; break;
                        case 213: str = str + "GCT"; break;
                        case 220: str = str + "GGA"; break;
                        case 221: str = str + "GGC"; break;
                        case 222: str = str + "GGG"; break;
                        case 223: str = str + "GGT"; break;
                        case 230: str = str + "GTA"; break;
                        case 231: str = str + "GTC"; break;
                        case 232: str = str + "GTG"; break;
                        case 233: str = str + "GTT"; break;
                        case 300: str = str + "TAA"; break;
                        case 301: str = str + "TAC"; break;
                        case 302: str = str + "TAG"; break;
                        case 303: str = str + "TAT"; break;
                        case 310: str = str + "TCA"; break;
                        case 311: str = str + "TCC"; break;
                        case 312: str = str + "TCG"; break;
                        case 313: str = str + "TCT"; break;
                        case 320: str = str + "TGA"; break;
                        case 321: str = str + "TGC"; break;
                        case 322: str = str + "TGG"; break;
                        case 323: str = str + "TGT"; break;
                        case 330: str = str + "TTA"; break;
                        case 331: str = str + "TTC"; break;
                        case 332: str = str + "TTG"; break;
                        case 333: str = str + "TTT"; break;
                        default:  str = str + "---"; break;
                    }
                }
            }
            break;
        case 19:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 00: str = str + "AA"; break;
                        case 01: str = str + "AC"; break;
                        case 02: str = str + "AG"; break;
                        case 03: str = str + "AT"; break;
                        case 10: str = str + "CA"; break;
                        case 11: str = str + "CC"; break;
                        case 12: str = str + "CG"; break;
                        case 13: str = str + "CT"; break;
                        case 20: str = str + "GA"; break;
                        case 21: str = str + "GC"; break;
                        case 22: str = str + "GG"; break;
                        case 23: str = str + "GT"; break;
                        case 30: str = str + "TA"; break;
                        case 31: str = str + "TC"; break;
                        case 32: str = str + "TG"; break;
                        case 33: str = str + "TT"; break;
                        default: str = str + "--"; break;
                    }
                }
            }
            break;
        case 20:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 00: str = str + "AA"; break;
                        case 01: str = str + "AC"; break;
                        case 02: str = str + "AG"; break;
                        case 03: str = str + "AT"; break;
                        case 10: str = str + "CA"; break;
                        case 11: str = str + "CC"; break;
                        case 12: str = str + "CG"; break;
                        case 13: str = str + "CT"; break;
                        case 20: str = str + "GA"; break;
                        case 21: str = str + "GC"; break;
                        case 22: str = str + "GG"; break;
                        case 23: str = str + "GT"; break;
                        case 30: str = str + "TA"; break;
                        case 31: str = str + "TC"; break;
                        case 32: str = str + "TG"; break;
                        case 33: str = str + "TT"; break;
                        default: str = str + "--"; break;
                    }
                }
            }
            break;
        case 21:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 00: str = str + "AA"; break;
                        case 01: str = str + "AC"; break;
                        case 02: str = str + "AG"; break;
                        case 03: str = str + "AT"; break;
                        case 10: str = str + "CA"; break;
                        case 11: str = str + "CC"; break;
                        case 12: str = str + "CG"; break;
                        case 13: str = str + "CT"; break;
                        case 20: str = str + "GA"; break;
                        case 21: str = str + "GC"; break;
                        case 22: str = str + "GG"; break;
                        case 23: str = str + "GT"; break;
                        case 30: str = str + "TA"; break;
                        case 31: str = str + "TC"; break;
                        case 32: str = str + "TG"; break;
                        case 33: str = str + "TT"; break;
                        default: str = str + "--"; break;
                    }
                }
            }
            break;
        case 22:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('G'); break;
                        case 3:  str.push_back('T'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 23:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('G'); break;
                        case 3:  str.push_back('T'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 24:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('G'); break;
                        case 3:  str.push_back('T'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 25:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('G'); break;
                        case 3:  str.push_back('T'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 26:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('G'); break;
                        case 3:  str.push_back('T'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 27:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('G'); break;
                        case 3:  str.push_back('T'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 28:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case 0:  str.push_back('A'); break;
                        case 1:  str.push_back('C'); break;
                        case 2:  str.push_back('G'); break;
                        case 3:  str.push_back('T'); break;
                        case 4:  str.push_back('K'); break;
                        case 5:  str.push_back('M'); break;
                        case 6:  str.push_back('R'); break;
                        case 7:  str.push_back('Y'); break;
                        case 8:  str.push_back('S'); break;
                        case 9:  str.push_back('W'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 29:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case  0: str.push_back('A'); break;
                        case  1: str.push_back('C'); break;
                        case  2: str.push_back('G'); break;
                        case  3: str.push_back('T'); break;
                        case  4: str.push_back('K'); break;
                        case  5: str.push_back('M'); break;
                        case  6: str.push_back('R'); break;
                        case  7: str.push_back('Y'); break;
                        case  8: str.push_back('S'); break;
                        case  9: str.push_back('W'); break;
                        case 10: str.push_back('B'); break;
                        case 11: str.push_back('D'); break;
                        case 12: str.push_back('H'); break;
                        case 13: str.push_back('V'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        case 30:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case  0: str.push_back('A'); break;
                        case  1: str.push_back('C'); break;
                        case  2: str.push_back('D'); break;
                        case  3: str.push_back('E'); break;
                        case  4: str.push_back('F'); break;
                        case  5: str.push_back('G'); break;
                        case  6: str.push_back('H'); break;
                        case  7: str.push_back('I'); break;
                        case  8: str.push_back('K'); break;
                        case  9: str.push_back('L'); break;
                        case 10: str.push_back('M'); break;
                        case 11: str.push_back('N'); break;
                        case 12: str.push_back('P'); break;
                        case 13: str.push_back('Q'); break;
                        case 14: str.push_back('R'); break;
                        case 15: str.push_back('S'); break;
                        case 16: str.push_back('T'); break;
                        case 17: str.push_back('V'); break;
                        case 18: str.push_back('W'); break;
                        case 19: str.push_back('Y'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
        default:
            for (std::vector<int>::size_type i = 0; i != seq_data.size(); ++i) {
                if (sites[i] == '1') {
                    number = seq_data[i];
                    switch (number) {
                        case  0: str.push_back('1'); break;
                        case  1: str.push_back('2'); break;
                        case  2: str.push_back('3'); break;
                        case  3: str.push_back('4'); break;
                        case  4: str.push_back('5'); break;
                        case  5: str.push_back('6'); break;
                        default: str.push_back('-'); break;
                    }
                }
            }
            break;
    }
    return(str);
}


// Function that reads input file and stores data in two 2D containers
unsigned long Read_Input(std::string inname, unsigned datatype){
    unsigned long alignment_length(0);
    unsigned long counter(0);
    std::string seq(""), str(""), tmp(""); // temporary string used to store input
    std::vector<int> sequence;    // temporary vector used to store input
    std::ifstream infile;
    
    infile.open(inname.c_str());
    while (getline(infile, str)) {
        if (!str.empty()) {
            // remove blank space in string
            tmp.clear();
            for (std::string::size_type i = 0; i != str.size(); ++i) {
                if (!isblank(str[i])) {
                    tmp.push_back(str[i]);
                }
            }
            if (tmp[0] == '>') {
                if (seq.size() > 0) {
                    if (datatype > 14 && datatype < 18) {
                        if (seq.size() % 2 != 0) {
                            std::cerr << "\nERROR: expected sequence of di-nucleotides" << "\n" << std::endl;
                            exit(1);
                        }
                    }
                    if (datatype > 17 && datatype < 28) {
                        if (seq.size() % 3 != 0) {
                            std::cerr << "\nERROR: expected sequence of codons" << "\n" << std::endl;
                            exit(1);
                        }
                    }
                    sequence = Translator(datatype, seq);
                    alignment.push_back(sequence); // stores sequence in vector
                    if (alignment_length == 0)
                        alignment_length = sequence.size();
                    sequence.clear();
                    seq.clear();
                }
                tmp.erase(tmp.begin()); // removes first character from name
                taxon.push_back(tmp); // stores sequence name in vector
            } else {
                seq += tmp;
            }
            str.clear();
        }
    }
    // Store last sequence in vector
    if (seq.size() > 0) {
        sequence = Translator(datatype, seq);
        alignment.push_back(sequence);
    } else {
        std::cerr << "\nERROR: last sequence empty" << "\n" << std::endl;
        exit(1);
    }
    //Check whether the sequence names are unique
    for (std::vector<std::string>::const_iterator iter1 = taxon.begin(); iter1 != taxon.end(); ++iter1) {
        for (std::vector<std::string>::const_iterator iter2 = iter1 + 1; iter2 != taxon.end(); ++iter2) {
            if (*iter1 == *iter2) {
                std::cerr << "\nERROR: sequence names not unique -- look for " << *iter1 << "\n" << std::endl;
                exit(1);
            }
        }
    }
    // Check whether the sequences have the same length
    for (std::vector<std::vector<int> >::const_iterator iter = alignment.begin()+1; iter != alignment.end(); ++iter) {
        ++counter;
        sequence = *iter;
        if (sequence.size() != alignment_length) {
            std::cerr << "\nERROR: sequences 1 and " << counter + 1 << " differ in length!\n" << std::endl;
            exit(1);
        }
    }
    return(alignment_length);
}


// Function used to identify variant sites
void Identify_Variant_Sites(std::string choice_of_sites, unsigned states, unsigned long alignment_length) {
    std::vector<int> column;

    for (std::string::size_type i = 0; i != alignment_length; ++i) {
        sites.push_back('1');
    }
    if (toupper(choice_of_sites[0]) == 'V') {
        for (std::string::size_type j = 0; j != sites.size(); ++j) {
            for (std::string::size_type i = 0; i != taxon.size(); ++i) {
                column.push_back(alignment[i][j]);
            }
            sort(column.begin(),column.end());
            std::vector<int>::iterator end_unique = unique(column.begin(),column.end());
            column.erase(end_unique,column.end());
            if (column.size() == 1 || (column.size() == 2 && column[1] == states)) {
               sites[j] = '0'; // exclude constant sites from the analysis
            }
            column.clear();
        }
    }
}


// Main program that calls the various functions above
int main(int argc, char** argv){
    unsigned alphabet(0);
    unsigned dataType(0);
    unsigned long alignment_length(0), sum_dm(0), sum_diag_dm(0);
    unsigned long varSites(0), total;
    unsigned long dm[max_array][max_array];            // 2D divergence matrix
    unsigned long row_sum[max_array], col_sum[max_array];
    double d_obs(0), d_ran(0), lambda(0);
    double min_lambda(std::numeric_limits<double>::max());
    double max_lambda(-std::numeric_limits<double>::max());
    std::vector<int> sequence;
    std::vector<int> column;
    std::vector<double> row_of_double;
    std::vector<std::vector<double> > mat_dobs, mat_lambda, mat_dB;
    std::string name_1, name_2; // temporary variables holding names of two sequences
    std::string choice_of_sites, nature_of_data, brevity, characters;
    std::string inName, outName1, outName2, outName3, outName4, outName5;
    std::ifstream infile;
    std::ofstream outfile1, outfile2, outfile3, outfile4, outfile5;
   
    if(argc != 5) {
        std::cerr << "\nSatuRation v1.0 Copyright 2019-22, Lars Jermiin" << std::endl;
        std::cerr << "\nERROR -- use command: saturation <infile> <a|v> <b|f> <1|..|31>\n" << std::endl;
        std::cerr << "  infile   Fasta-formatted alignment" << std::endl;
        std::cerr << "     a|v   All or variant sites" << std::endl;
        std::cerr << "     b|f   Brief or full report of results" << std::endl;
        std::cerr << "       1   Nucleotides; 4 states (A|C|G|T)" << std::endl;
        std::cerr << "       2   Nucleotides; 3 states (C|T|AG)" << std::endl;
        std::cerr << "       3   Nucleotides; 3 states (A|G|CT)" << std::endl;
        std::cerr << "       4   Nucleotides; 3 states (A|T|CG)" << std::endl;
        std::cerr << "       5   Nucleotides; 3 states (C|G|AT)" << std::endl;
        std::cerr << "       6   Nucleotides; 3 states (A|C|GT)" << std::endl;
        std::cerr << "       7   Nucleotides; 3 states (G|T|AC)" << std::endl;
        std::cerr << "       8   Nucleotides; 2 states (GT|AC)" << std::endl;
        std::cerr << "       9   Nucleotides; 2 states (AG|CT)" << std::endl;
        std::cerr << "      10   Nucleotides; 2 states (GC|AT)" << std::endl;
        std::cerr << "      11   Nucleotides; 2 states (A|CGT)" << std::endl;
        std::cerr << "      12   Nucleotides; 2 states (C|AGT)" << std::endl;
        std::cerr << "      13   Nucleotides; 2 states (G|ACT)" << std::endl;
        std::cerr << "      14   Nucleotides; 2 states (T|ACG)" << std::endl;
        std::cerr << "      15   Di-nucleotides (pos 12); 16 states (AA|AC|...|TG|TT)" << std::endl;
        std::cerr << "      16   Di-nucleotides (pos  1);  4 states (A|C|G|T)" << std::endl;
        std::cerr << "      17   Di-nucleotides (pos  2);  4 states (A|C|G|T)" << std::endl;
        std::cerr << "      18   Codons (pos 123); 64 states (AAA|AAC|...|TTG|TTT)" << std::endl;
        std::cerr << "      19   Codons (pos  12); 16 states (AA|AC|...|TG|TT)" << std::endl;
        std::cerr << "      20   Codons (pos  13); 16 states (AA|AC|...|TG|TT)" << std::endl;
        std::cerr << "      21   Codons (pos  23); 16 states (AA|AC|...|TG|TT)" << std::endl;
        std::cerr << "      22   Codons (pos  12);  4 states (A|C|G|T)" << std::endl;
        std::cerr << "      23   Codons (pos  13);  4 states (A|C|G|T)" << std::endl;
        std::cerr << "      24   Codons (pos  23);  4 states (A|C|G|T)" << std::endl;
        std::cerr << "      25   Codons (pos   1);  4 states (A|C|G|T)" << std::endl;
        std::cerr << "      26   Codons (pos   2);  4 states (A|C|G|T)" << std::endl;
        std::cerr << "      27   Codons (pos   3);  4 states (A|C|G|T)" << std::endl;
        std::cerr << "      28   Genotypes; 10 states (A|C|G|T|K|M|R|Y|S|W)" << std::endl;
        std::cerr << "      29   Genotypes; 14 states (A|C|G|T|K|M|R|Y|S|W|B|D|H|V)" << std::endl;
        std::cerr << "      30   Amino acids; 20 states (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)" << std::endl;
        std::cerr << "      31   Amino acids;  6 states (AGPST|DENQ|HKR|MIVL|WFY|C) [D6]" << std::endl;
        std::cerr << std::endl;
        exit(1);
    }
    inName = argv[1];
    choice_of_sites = argv[2];
    brevity = argv[3];
    nature_of_data = argv[4];
    // Check availability of input file
    infile.open(inName.c_str());
    if (!infile) {
        std::cerr << "\nERROR: input file not found...\n" << std::endl;
        exit(1);
    }
    // Check choice of sites
    if (toupper(choice_of_sites[0]) != 'A' && toupper(choice_of_sites[0]) != 'V') {
        std::cerr << "\nERROR: incorrect choice of sites: [a|v]\n" << std::endl;
        exit(1);
    }
    // Check choice of output
    if (toupper(brevity[0]) != 'F' && toupper(brevity[0]) != 'B') {
        std::cerr << "\nERROR: incorrect choice of output: [b|f]\n" << std::endl;
        exit(1);
    }
    // Check choice of data and alphabet
    dataType = stoi(nature_of_data);
    if (dataType < 1 || dataType > 31) {
        std::cerr << "\nERROR: incorrect choice of data: [1|...|31]\n" << std::endl;
        exit(1);
    }
    if (toupper(brevity[0]) == 'F') {
        outName1.clear();
        for (std::string::size_type i = 0; i != inName.size() && inName[i] != '.'; ++i) {
            outName1 += inName[i];
        }
        outName5 = outName1 + "_sites_used.fst";
        outName4 = outName1 + "_lambda.csv";
        outName3 = outName1 + "_d_obs.dis";
        outName2 = outName1 + "_d_obs.csv";
        outName1 = outName1 + "_table.csv";
    }
    switch (dataType) {
        case  1: alphabet = FOUR; break;
        case  2: alphabet = THREE; break;
        case  3: alphabet = THREE; break;
        case  4: alphabet = THREE; break;
        case  5: alphabet = THREE; break;
        case  6: alphabet = THREE; break;
        case  7: alphabet = THREE; break;
        case  8: alphabet = TWO; break;
        case  9: alphabet = TWO; break;
        case 10: alphabet = TWO; break;
        case 11: alphabet = TWO; break;
        case 12: alphabet = TWO; break;
        case 13: alphabet = TWO; break;
        case 14: alphabet = TWO; break;
        case 15: alphabet = SIXTEEN; break;
        case 16: alphabet = FOUR; break;
        case 17: alphabet = FOUR; break;
        case 18: alphabet = SIXTYFOUR; break;
        case 19: alphabet = SIXTEEN; break;
        case 20: alphabet = SIXTEEN; break;
        case 21: alphabet = SIXTEEN; break;
        case 22: alphabet = FOUR; break;
        case 23: alphabet = FOUR; break;
        case 24: alphabet = FOUR; break;
        case 25: alphabet = FOUR; break;
        case 26: alphabet = FOUR; break;
        case 27: alphabet = FOUR; break;
        case 28: alphabet = TEN; break;
        case 29: alphabet = FOURTEEN; break;
        case 30: alphabet = TWENTY; break;
        default: alphabet = SIX; break;
    }
    if (toupper(brevity[0]) == 'F') {
        std::cout << "\nReading input file ..." << std::endl;
    }
    alignment_length = Read_Input(inName, dataType);
    if (toupper(brevity[0]) == 'F') {
        std::cout << "\nProcessing data ...\n" << std::endl;
    }
    Identify_Variant_Sites(choice_of_sites, alphabet, alignment_length);
    
    // Prime square matrices (mat_dobs, mat_lambda, mat_daitFS, and mat_daitMS) with 0.0
    for (std::vector<double>::size_type i = 0; i != taxon.size(); i++) {
        row_of_double.push_back(0.0);
    }
    for (std::vector<std::vector<double> >::size_type i = 0; i != taxon.size(); i++) {
        mat_dobs.push_back(row_of_double);
    }
    mat_lambda = mat_dobs;
    mat_dB = mat_dobs;
    if (toupper(brevity[0]) == 'F') {
        outfile1.open(outName1.c_str());
        outfile1 << "Taxon 1,Taxon 2,d_obs,d_ran,lambda" << std::endl;
    }
    // Start the generation of results
    total = taxon.size() * (taxon.size() - 1)/2;
    for (std::vector<std::vector<int> >::size_type iter1 = 0; iter1 != alignment.size(); ++iter1) {
        for (std::vector<std::vector<int> >::size_type iter2 = iter1 + 1; iter2 != alignment.size(); ++iter2) {
            // Prime divergence matrix with zero
            for (size_t i = 0; i != max_array; ++i) {
                for (size_t j = 0; j != max_array; ++j) {
                    dm[i][j] = 0;
                }
            }
            // Generate divergence matrix
            for (std::string::size_type i = 0; i != sites.size(); ++i) {
                if (sites[i] == '1') {
                    dm[alignment[iter1][i]][alignment[iter2][i]]++;
                }
            }
            // Generate the sum over all elements in the divergence matrix
            sum_dm = 0;
            for (size_t i = 0; i != alphabet; ++i) {
                for (size_t j = 0; j != alphabet; ++j) {
                    sum_dm += dm[i][j];
                }
            }
            // Generate the sum over diagonal elements in the divergence matrix
            sum_diag_dm = 0;
            for (size_t i = 0; i != alphabet; ++i) {
                sum_diag_dm += dm[i][i];
            }
            // Calculate the observed sequence divergence
            d_obs = (double)(sum_dm - sum_diag_dm)/sum_dm;
            mat_dobs[iter1][iter2] = d_obs;
            mat_dobs[iter2][iter1] = d_obs;
            
            // Prime marginal vectors of divergence matrix with zero
            for (size_t i = 0; i != alphabet; ++i) {
                row_sum[i] = 0;
                col_sum[i] = 0;
            }
            // Generate vectors of marginal frequencies
            for (size_t i = 0; i != alphabet; ++i) {
                for (size_t j = 0; j != alphabet; ++j) {
                    row_sum[i] += dm[i][j];
                }
            }
            for (size_t i = 0; i != alphabet; ++i) {
                for (size_t j = 0; j != alphabet; ++j) {
                    col_sum[i] += dm[j][i];
                }
            }
            // Calculate the divergence of randomized sequences
            d_ran = 0.0;
            for (size_t i = 0; i != alphabet; ++i) {
                d_ran += ((double)row_sum[i]/sum_dm)*((double)col_sum[i]/sum_dm);
            }
            d_ran = 1.0 - d_ran;
            // Calculate the level of saturation (lambda)
            lambda = d_obs/d_ran;
            if (lambda < min_lambda) {
                min_lambda = lambda;
            }
            if (lambda > max_lambda) {
                max_lambda = lambda;
                // Capture two numbers, pointing to the names of the relevant seequence pair
                name_1 = taxon[iter1];
                name_2 = taxon[iter2];
            }
            mat_lambda[iter1][iter2] = lambda;
            mat_lambda[iter2][iter1] = lambda;

            if (toupper(brevity[0]) == 'F') {
                outfile1 << taxon[iter1] << ",";
                outfile1 << taxon[iter2] << ",";
                outfile1 << d_obs << ",";
                outfile1 << d_ran << ",";
                outfile1 << lambda;
                outfile1 << std::endl;
            }
        }
    }
    if (toupper(brevity[0]) == 'F') {
        outfile1.close();
        // Printing tables to files
        outfile2.open(outName2.c_str());
        outfile2 << taxon.size() << std::endl;
        outfile3.open(outName3.c_str());
        outfile3 << taxon.size() << std::endl;
        outfile4.open(outName4.c_str());
        outfile4 << taxon.size() << std::endl;
        if (toupper(choice_of_sites[0]) == 'V') {
            outfile5.open(outName5.c_str());
        }
        for (std::vector<std::vector<int> >::size_type i = 0; i != taxon.size(); i++) {
            outfile2 << taxon[i];
            outfile3 << std::left << std::setw(10) << taxon[i];
            outfile4 << taxon[i];
            if (toupper(choice_of_sites[0]) == 'V') {
                outfile5 << ">" << taxon[i] << std::endl;
                sequence.clear();
                for (std::vector<std::vector<int> >::size_type j = 0; j != alignment_length; ++j) {
                    sequence.push_back(alignment[i][j]);
                }
                characters = Back_translator(dataType, sequence); // A|V sites
                outfile5 << characters << std::endl;
            }
            for (std::vector<std::vector<int> >::size_type j = 0; j != taxon.size(); j++) {
                outfile2 << "," << std::fixed << mat_dobs[i][j];
                outfile3 << "\t" << std::fixed << mat_dobs[i][j];
                outfile4 << "," << std::fixed << mat_lambda[i][j];
            }
            outfile2 << std::endl;
            outfile3 << std::endl;
            outfile4 << std::endl;
            if (toupper(choice_of_sites[0]) == 'V') {
                outfile5 << std::endl;
            }
        }
        outfile2.close();
        outfile3.close();
        outfile4.close();
        if (toupper(choice_of_sites[0]) == 'V') {
            outfile5.close();
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "   RESULTS FROM ANALYSIS OF SATURATION COMPLETE" << std::endl;
        std::cout << std::endl;
        std::cout << "   All estimates ............................. " << outName1 << std::endl;
        std::cout << "   Estimates of d_obs ........................ " << outName2 << std::endl;
        std::cout << "   Estimates of d_obs ........................ " << outName3 << std::endl;
        std::cout << "   Estimates of lambda ....................... " << outName4 << std::endl;
        if (toupper(choice_of_sites[0]) == 'V') {
            std::cout << "   Alignment of sites used in this analysis .. " << outName5 << std::endl;
        }
        std::cout << std::endl;
        if (min_lambda == std::numeric_limits<double>::max() || max_lambda == -std::numeric_limits<double>::max()) {
            std::cout << "   Unexpected value of lambda - check alignment" << std::endl;
        } else {
            std::cout << "   Min(lambda) ............................... " << min_lambda << std::endl;
            std::cout << "   Max(lambda) ............................... " << max_lambda << std::endl;
            std::cout << "   Range(lambda) ............................. " << max_lambda - min_lambda << std::endl;
        }
        std::cout << "   Sites considered during this analysis ..... ";
        if (toupper(choice_of_sites[0]) == 'V') {
            std::cout << "VARIANT" << std::endl;
        } else {
            std::cout << "ALL" << std::endl;
        }
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
    } else {
        std::cout << "File, #taxa, Sites, #sites, Min(lambda), Max(lambda), Pair (max(lambda))" << std::endl;
        for(std::string::size_type i = 0; i != sites.size(); ++i) {
            if(sites[i] == '1')
                ++varSites;;
        }
        std::cout << inName << "," << taxon.size() << ",";
        if (toupper(choice_of_sites[0]) == 'V') {
            std::cout << "VARIANT,";
        } else {
            std::cout << "All,";
        }
        std::cout << varSites << ",";
        if (min_lambda == std::numeric_limits<double>::max() || max_lambda == -std::numeric_limits<double>::max()) {
            std::cout << "Nan,Nan,Nan,ODD RESULT - check alignment" << std::endl;
        } else {
            std::cout << std::fixed << min_lambda << "," << std::fixed << max_lambda << ",";
            std::cout << name_1 << " vs " << name_2 << std::endl;
        }
    }
    return 0;
}
