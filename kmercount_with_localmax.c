/*  PROGRAM TO COUNT KMERS AND IDENTIFY LOCAL MAXIMA IN HUDDINGE SPACE        */
/*  ENCODES DNA AS 2 bit REPRESENTATION TO SPEED UP COUNTING                  */
/*  DIRECTLY USES KMER SEQUENCE VALUE AS INDEX, MEMORY REQUIREMENTS           */
/*  SCALE POORLY (k^2), GOOD UP TO 12-mers OR SO                              */
/*  BY J Taipale 2024, most code is from spacek40 & autoseed                  */
/*  example ./kmercount_with_localmax bak.seq sig.seq 40 8 8 100 0.35 n       */
/*  counts ungapped 8-mers from files containing 40 bp reads,                 */
/*  prints those with more than 100 counts                                    */
/*  compile with: cc -O3 -o kmercount_with_localmax kmercount_with_localmax.c */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

/* GLOBAL VARIABLES */
char *version = "kmercount_with_localmax, version 0.5 Oct 9, 2024";
char *usage = "\nUsage: ./kmercount_with_localmax [background filename] [signal filename] [read length (max 64)] [min_kmerlength] [max_kmerlength] [minimum count cutoff for printing] [base frequency cutoff for selecting a longer kmer as local maxima (default 0.35)] [gaps counted? counts center gap only by default, 'a'=count all gaps, 'n'=count only ungapped kmers]\nReads must be on individual lines with bases in all caps";
__uint128_t mask_ULL[42][42];   /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
__uint128_t lowmask_ULL[42];    /* LOW MASK FOR EACH KMER */
__uint128_t highmask_ULL[42];   /* HIGH MASK FOR EACH KMER */
short int Nlength;              /* LENGTH OF SEQUENCE READ FRAGMENTS */

/* SUBROUTINE THAT DETERMINES IF GAP IS AT EITHER OF THE CENTER POSITIONS, IF count_also is set to != 1 returns true */
short int Centergap (short int count_also_spaced_kmers, short int kmer_length, short int gap_position)
{
if (count_also_spaced_kmers != 1) return (1);
if (gap_position == kmer_length / 2) return (1);
if (kmer_length % 2 == 1 && gap_position == kmer_length / 2 - 1) return (1);
else return (0);
}

/* SUBROUTINE THAT DETERMINES IF A GAPPED KMER IS A LOCAL MAXIMUM WITHIN HUDDINGE DISTANCE OF 1 */
/* SEE NITTA ET AL. eLIFE 2015 Methods and Supplementary Figure 1 for algorithm description     */
/* https://doi.org/10.7554/eLife.04837.004                                                      */
short int Localmax(long int *****results, short int file_number, short int shortest_kmer, short int too_long_kmer, short int current_kmer_length,short int current_gap_position, short int current_gap_length, long int current_kmer, short int count_also_spaced_kmers, double kmer_length_difference_cutoff)
{
    long int kmer1_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
    long int kmer2_incidence;
    long int compared_kmer = current_kmer;
    signed short int position;
    short int counter;
    short int first_half;
    long int lowbit = 1;
    long int highbit = 2;
    short int shift;
    short int true_gap_position = current_gap_position;
    short int true_gap_length = current_gap_length;
    short int start;
    short int end;
    short int left = 0;
    short int right = 1;
    signed short int position2;
    
    /* Substitution; HAMMING OF 1, RETURNS 0 IF ANY KMER WITHIN HAMMING OF 1 HAS HIGHER COUNT */
    for(position=0; position < current_kmer_length; position++, lowbit <<= 2, highbit <<= 2)
    {
        compared_kmer = lowbit ^ current_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = highbit ^ current_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    /* Shift; FULL SHIFT BY ONE */
    shift = (current_gap_length != 0);
    current_gap_position += shift;
    
    /* ONLY LOOK AT FULL SHIFT FOR UNGAPPED KMERS, OR FOR GAPPED KMERS IF GAP CAN SHIFT ALSO (ALL GAP POSITIONS HAVE BEEN COUNTED) */
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 && current_gap_position < current_kmer_length) || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2 + 1))
    {
    /* COMPARED FULL SHIFT RIGHT (same shift in gap if any), RETURNS 0 IF ANY KMER WITHIN HAMMING OF 1 HAS HIGHER COUNT  */
    compared_kmer = (current_kmer >> 2) & lowmask_ULL[current_kmer_length-1];
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    lowbit >>= 2; highbit >>= 2;
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    current_gap_position -= shift;
    current_gap_position -= shift;
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 && current_gap_position > 0) || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2))
    {
    /* COMPARED FULL SHIFT LEFT (same shift in gap if any) */
    compared_kmer = (current_kmer << 2) & lowmask_ULL[current_kmer_length-1];
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    lowbit = 1; highbit = 2;
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    current_gap_position += shift;
    lowbit = 1; highbit = 2;
    if (count_also_spaced_kmers != 0)
    {
    if(current_gap_position == 0) current_gap_position = current_kmer_length / 2;
    
    /* Longer Gap; COMPARE TO KMER WITH LONGER GAP */
    current_gap_length++;
    if (current_gap_length < Nlength - current_kmer_length && true_gap_length != 0)
    {
    compared_kmer = (current_kmer & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    
    /* LOOP TO ANALYZE SHIFT OF EITHER HALF, STARTS WITH SECOND HALF */
    for(first_half = 0; first_half < 2; first_half++)
    {
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    
    /* SWITCHES TO FIRST HALF */
    compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
    }
    }
    
    /* Shorter Gap; COMPARE TO KMER WITH SHORTER GAP */
    current_gap_length--;
    current_gap_length--;

    lowbit = 1; highbit = 2;
    lowbit <<= ((current_kmer_length - true_gap_position - 1)*2); highbit <<= ((current_kmer_length - true_gap_position - 1)*2);
    if (current_gap_length == 0) current_gap_position = 0;
    if (current_gap_length >= 0)
    {
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer >> 2) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
        /* LOOP TO ANALYZE SHIFT OF EITHER HALF, STARTS WITH SECOND HALF */
        for(first_half = 0; first_half < 2; first_half++)
        {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            
            /* SWITCHES TO FIRST HALF */
            compared_kmer = lowmask_ULL[current_kmer_length-1] & ((current_kmer << 2) & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
            lowbit <<= 2; highbit <<= 2;
        }
    }
    
        /* COMPARES DIFFERENT GAP POSITIONS (SHIFTED BY ONE) */
        current_gap_length = true_gap_length;
        if ((count_also_spaced_kmers == 2 || (count_also_spaced_kmers == 1 && current_kmer_length % 2 == 1)) && true_gap_length > 0)
        {
        current_gap_position = true_gap_position + 1;
        lowbit = 1; highbit = 2;
        lowbit <<= ((current_kmer_length - current_gap_position)*2); highbit <<= ((current_kmer_length - current_gap_position)*2);
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));

        /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH RIGHT BY ONE */
        for(first_half = 0; first_half < 2; first_half++)
        {
        if (current_gap_position < current_kmer_length && current_gap_position > 0 && (count_also_spaced_kmers == 2 ||
        (count_also_spaced_kmers == 1 && ((current_gap_position == current_kmer_length / 2) || current_gap_position == current_kmer_length / 2 + 1))))
        {
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = highbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        }
            
        /* SWITCHES LEFT BY ONE */
        current_gap_position--;
        current_gap_position--;
        compared_kmer = ((current_kmer) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
        lowbit <<= 2; highbit <<= 2;
        }
        }
        
        start = 1;
        end = current_kmer_length;
        /* COMPARES UNGAPPED KMER TO ALL SINGLE GAPS IN ALL POSITIONS */
        if (count_also_spaced_kmers != 0 && true_gap_length == 0)
        {
            if(count_also_spaced_kmers == 1)
            {
                start = current_kmer_length / 2;
                end = start + 1 + (current_kmer_length % 2);
            }
        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++)
        {
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH BEGINNING */
            for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++)
            {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
           
        /* END */
        compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
        }
        }
        }
        
    }
/* Shorter; COMPARES KMER WITH ONE SHORTER */
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
current_kmer_length--;
end = current_kmer_length;
if (current_kmer_length >= shortest_kmer)
{
/* IF NO GAP, INSERTS GAP AT ALL ALLOWED POSITIONS */
if(current_gap_length == 0)
{
    if (count_also_spaced_kmers != 0)
    {
        if(count_also_spaced_kmers == 1)
        {
            start = current_kmer_length / 2;
            end = start + 1 + (current_kmer_length % 2);
        }
        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++)
        {
            compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                if(kmer2_incidence  * kmer_length_difference_cutoff > kmer1_incidence) return (0);
        }
    }
}

current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
if (count_also_spaced_kmers != 1) {left = 1; right = 1;}
if (current_gap_position == current_kmer_length / 2 && current_kmer_length % 2 == 0 && count_also_spaced_kmers == 1) {left = 1; right = 0;}
if (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1) left = 1;
    
/* LEFT PART */
if (current_gap_position < current_kmer_length)
{
    if (left == 1 || true_gap_length == 0)
    {
compared_kmer = (current_kmer >> 2);
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
    }
/* RIGHT PART */
if (current_gap_position != 1 && right == 1)
{
if (current_gap_position > 0) current_gap_position--;
    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
current_gap_position = true_gap_position;
/* Shorter with Longer Gap; LONGER GAP */
if(count_also_spaced_kmers != 0 && true_gap_position != 0)
{
current_gap_length++;
if (current_gap_position < current_kmer_length && left == 1)
{
/* RIGHT BASE GAPPED */
compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
/* LEFT BASE GAPPED */
if(current_gap_position > 1 && right == 1)
{
current_gap_position--;
compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
current_gap_length--;
current_gap_position++;
}
    /* COMPARES HANGING SINGLE BASE TO UNGAPPED KMER */
    current_gap_position = true_gap_position;
    current_gap_length = true_gap_length;
    if(current_gap_position == 1)
    {
    current_gap_length = 0;
    current_gap_position = 0;
    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
    else if(current_kmer_length-current_gap_position == 0)
    {
    current_gap_length = 0;
    current_gap_position = 0;
    compared_kmer = (current_kmer >> 2);
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
}
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
/* Longer; COMPARES KMER WITH ONE LONGER */
current_kmer_length++;
current_kmer_length++;
if (current_kmer_length < too_long_kmer)
{
compared_kmer = (current_kmer << 2);
/* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH BEGINNING */
for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++)
{
if(count_also_spaced_kmers != 1 || true_gap_length == 0 || (first_half == 1 || (count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2)))
{
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
compared_kmer = lowbit ^ compared_kmer;
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
compared_kmer = highbit ^ compared_kmer;
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
compared_kmer = lowbit ^ compared_kmer;
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
}
/* END */
compared_kmer = current_kmer;
lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
if (true_gap_length == 0) continue;
current_gap_position++;
if (current_gap_position > current_kmer_length || (count_also_spaced_kmers == 1 && current_gap_position != (current_kmer_length / 2 + (current_kmer_length % 2)))) break;
}
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
/* Longer with Shorter Gap; COMPARES TO ONE SHORTER GAP LENGTH */
if(count_also_spaced_kmers != 0 && true_gap_length >= 1)
{
current_gap_length--;
lowbit = 1; highbit = 2;
lowbit <<= ((current_kmer_length - current_gap_position-1)*2); highbit <<= ((current_kmer_length - current_gap_position-1)*2);
compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
/* NO GAP LEFT */
if (current_gap_length == 0) current_gap_position = 0;
    
    /* LOOP TO ANALYZE ADDED BASE ON EITHER SIDE OF GAP, STARTS WITH RIGHT */
    for(first_half = 0; first_half < 2; first_half++)
    {
        if (current_gap_position < current_kmer_length && (current_gap_position == 0 || Centergap (count_also_spaced_kmers, current_kmer_length, current_gap_position)))
        {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
        }
        
        /* SWITCHES ADDED BASE TO LEFT */
        current_gap_position++;
        if (current_gap_position > current_kmer_length || current_gap_position == 1) break;
        compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    }
current_gap_length++;
}
}
current_gap_position = true_gap_position;
return (1);
}

/* SUBROUTINE THAT GENERATES 128 bit BITMASKS FOR NUCLEOTIDES, KMER STRINGS AND DELETIONS */
void GenerateMask ()
{
    short int counter;
    short int start_position;
    short int position;
    short int current_kmer_length;
    /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
    for (mask_ULL[1][0] = 3, counter = 0; counter < Nlength; counter++) mask_ULL[1][counter+1] = mask_ULL[1][counter] << 2;
    /* GENERATES mask_ULLS FOR EXTRACTION OF EACH KMER STRING */
    for(current_kmer_length = 2; current_kmer_length <= Nlength; current_kmer_length++)
    {
        for (start_position = 0; start_position < Nlength-current_kmer_length; start_position++)
        {
            for (mask_ULL[current_kmer_length][start_position]=mask_ULL[1][start_position], position = start_position+1; position < current_kmer_length+start_position; position++) {mask_ULL[current_kmer_length][start_position] += mask_ULL[1][position];}
        }
    }
    /* GENERATES HIGH AND LOW mask_ULLS FOR DELETIONS */
    for (lowmask_ULL[0] = mask_ULL[1][0], position = 1; position < Nlength-2; position++) lowmask_ULL[position] = lowmask_ULL[position-1]+mask_ULL[1][position];
    for (highmask_ULL[Nlength-2] = mask_ULL[1][Nlength-2], position = Nlength-3; position > 0; position--) highmask_ULL[position] = highmask_ULL[position+1]+mask_ULL[1][position];
}

/* SUBROUTINE THAT CLEARS RESULTS TABLE */
void InitializeResults(long int *****results, short int shortest_kmer, short int too_long_kmer, short int count_also_spaced_kmers, short int number_of_files)
{
short int counter;
short int current_kmer_length;
short int current_gap_position;
short int current_gap_length;
short int kmer_length_size;
short int gap_position_size;
long long int variable_length;
long int *current_value;
long int current_kmer;
    
for(counter = 0; counter < number_of_files; counter++)
{
    results[counter] = malloc(sizeof(long int *) * too_long_kmer + 5);
    for(current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
    {
        kmer_length_size = ((current_kmer_length-2) * (count_also_spaced_kmers != 0)) + 2;
        results[counter][current_kmer_length] = malloc(sizeof(long int *) * kmer_length_size + 5);
        for(current_gap_position = 0; current_gap_position < kmer_length_size; current_gap_position++)
        {
            /* DO NOT ALLOCATE MEMORY FOR ALL GAPS IF ONLY CENTER GAPS COUNTED */
            if (count_also_spaced_kmers == 0 && current_gap_position > 1) continue;
            if (count_also_spaced_kmers == 1 && (current_gap_position != current_kmer_length / 2 && current_gap_position != current_kmer_length / 2 + current_kmer_length % 2 && current_gap_position > 1)) continue;
            
            gap_position_size = ((Nlength-current_kmer_length+1) * (count_also_spaced_kmers != 0)) + 1;
            results[counter][current_kmer_length][current_gap_position] = malloc(sizeof(long int *) * gap_position_size + 5);
            for(current_gap_length = 0; current_gap_length < gap_position_size; current_gap_length++)
            {
                /* DO NOT ALLOCATE MEMORY FOR ALL GAPS IF ONLY CENTER GAPS COUNTED */
                if (current_gap_position > 1 && current_gap_length == 0) continue;
                if (current_gap_position == 0 && current_gap_length > 1) continue;
                if (count_also_spaced_kmers != 2 && current_gap_position == 1 && current_gap_length == 1 && current_kmer_length > 3) break;
                /* ALLOCATES MEMORY AND ZEROES IT */
                variable_length = pow(4, current_kmer_length);
                results[counter][current_kmer_length][current_gap_position][current_gap_length] = malloc(variable_length * sizeof(long int) + 5);
                current_value = results[counter][current_kmer_length][current_gap_position][current_gap_length];
                for(current_kmer = 0; current_kmer < variable_length; current_kmer++) current_value[current_kmer] = 0;
            }
        }
    }
}
}
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
int main (int argc, char *argv[])
{
    long int counter;
    short int number_of_files = 2;
    FILE *open_file;
    char **file_name;
    file_name = malloc (sizeof (char *) * (number_of_files + 1) + 5);
    for (counter = 0; counter < number_of_files; counter++) {file_name[counter] = malloc(1000); strcpy(file_name[counter], "no file");}
    
    /* PARSES ARGUMENTS */
    if (argc < 7) {printf("%s%s\n", version, usage); exit(1);}
    strcpy(file_name[0], argv[1]);
    strcpy(file_name[1], argv[2]);
    Nlength = atoi(argv[3]) + 1;
    short int shortest_kmer = atoi(argv[4]);
    short int too_long_kmer = atoi(argv[5]) + 1;
    float kmer_length_difference_cutoff = 0.35;
    short int print_cutoff = atoi(argv[6]) + 1;
    if (argc > 7) kmer_length_difference_cutoff = atof(argv[7]);
    short int count_also_spaced_kmers = 1;
    if (argc > 8 && argv[8][0] == 'a') count_also_spaced_kmers = 2;
    if (argc > 8 && argv[8][0] == 'n') count_also_spaced_kmers = 0;

    /* INITIALIZES VARIABLES */
    short int gap_position_size;
    short int kmer_length_size;
    char text1;
    long int charcounter;
    short int nucleotide_value;
    short int file_number = 0;
    short int strand;
    long int read_index;
    short int eof_reached;
    long int current_kmer;
    short int current_kmer_length;
    short int current_gap_position;
    short int current_gap_length;
    signed short int position;
    short int start_position;
    short int end_position;
    long int kmer_incidence;
    short int deletion_size = 1;
    short int Nmer_position;
    short int too_long_nmer_position;
    char *current_sequence = malloc(10000);
    char *dnaforward = "ACGT";
    char *dnareverse = "TGCA";
    
    /* GENERATES BITMASKS AND DEFINES 128 bit ints FOR FULL READ SEQUENCE VALUES */
    GenerateMask ();
    __uint128_t current_sequence_value_ULL;
    __uint128_t forward_sequence_value_ULL;
    __uint128_t deleted_sequence_value_ULL;
    __uint128_t position_value;
    __uint128_t left_position_value = 1;
    left_position_value <<= ((Nlength-2) * 2);
    
    /* INITIALIZES RESULT TABLE */
    long int *****results;
    results = malloc(sizeof(long int *) * 3 + 5);
    InitializeResults(results, shortest_kmer, too_long_kmer, count_also_spaced_kmers, number_of_files);

    /* FILE MAIN LOOP */
    for (; file_number < 2; fclose(open_file), file_number++)
    {
        open_file = fopen(file_name[file_number], "r");
        if (open_file == (void *)0)
        {
            printf("File: %s not found\n\n", file_name[file_number]);
            exit(1);
        }
        
        /* SEQUENCE LINE LOADING LOOP */
        for(eof_reached = 0, read_index = 1; ; read_index++)
        {
        start_of_line_loop:
            /* TAKES ONE SEQUENCE FROM FILE */
            for(charcounter = 0; ; )
            {
                text1 = getc(open_file);
                if (text1 == EOF) {eof_reached = 1; break;}
                if (text1 == '\n') break;
                if (text1 == 'A' || text1 == 'C' || text1 == 'G' || text1 == 'T' || text1 == 'N') /* ONLY ACCEPTS NUCLEOTIDES IN CAPS, N ONLY ACCEPTED SO THAT ERROR WILL BE DIFFERENT */
                {
                    current_sequence[charcounter] = text1;
                    charcounter++;
                }
            }
            current_sequence[charcounter] = '\0';
            if(eof_reached == 0 && (strlen(current_sequence) != Nlength-1)) {printf("\nWrong sequence length on line %li", read_index); goto start_of_line_loop;}
            
            /* CHECKS IF AT END OF FILE */
            if (eof_reached == 1) {printf("\nEOF encountered in file %i on line %li\n", file_number, read_index); break;}
            
            /* STRAND LOOP */
            for(strand = 0; strand < 2; strand++)
            {
                /* CALCULATES INTEGER VALUE CORRESPONDING TO SEQUENCE N-mer */
                if (strand == 0)
                {
                    /* FORWARD STRAND */
                    for(current_sequence_value_ULL = 0, position = 0, position_value = left_position_value; position < Nlength-1; position++, position_value /= 4)
                    {
                        for (nucleotide_value = 0; nucleotide_value < 4 && current_sequence[position] != dnaforward[nucleotide_value]; nucleotide_value++);
                        if(nucleotide_value == 4) {printf("\nSEQUENCE ERROR AT POSITION %li, %i \n", read_index, position); goto start_of_line_loop;}
                        current_sequence_value_ULL += position_value * nucleotide_value;
                    }
                    forward_sequence_value_ULL = current_sequence_value_ULL;
                }
                else
                {
                    /* REVERSE STRAND */
                    for(current_sequence_value_ULL = 0, position = Nlength-2, position_value = left_position_value; position > -1; position--, position_value /= 4)
                    {
                        for (nucleotide_value = 0; nucleotide_value < 4 && current_sequence[position] != dnareverse[nucleotide_value]; nucleotide_value++);
                        if(nucleotide_value == 4) {printf("\nSEQUENCE ERROR AT POSITION %li, %i \n", read_index, position); goto start_of_line_loop;}
                        current_sequence_value_ULL += position_value * nucleotide_value;
                    }
                }
                                
                /* COUNTS KMERS */
                /* KMERS WITH NO GAPS */
                for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
                {
                    end_position = Nlength-current_kmer_length;
                    for (position = 0; position < end_position; position++)
                    {
                        results[file_number][current_kmer_length][0][0][(current_sequence_value_ULL & mask_ULL[current_kmer_length][position]) >> (position * 2)]++;
                    }
                }
            }
            if (count_also_spaced_kmers != 0)
            {
                /* KMERS WITH GAPS */
                /* GENERATES DELETIONS INTO CURRENT SEQUENCE */
                for (deletion_size = 1; deletion_size < Nlength - shortest_kmer; deletion_size++)
                {
                    Nmer_position = 1;
                    too_long_nmer_position = Nlength-1-deletion_size;
                    for (; Nmer_position < too_long_nmer_position; Nmer_position++)
                    {
                        deleted_sequence_value_ULL = (current_sequence_value_ULL & lowmask_ULL[Nmer_position-1]) ^ ((current_sequence_value_ULL & highmask_ULL[Nmer_position+deletion_size]) >> (deletion_size * 2));
                        
                        /* FINDS KMERS WITH GAPS FROM DELETED SEQUENCE */
                        for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
                        {
                            if (count_also_spaced_kmers == 1)
                            {
                                position = Nmer_position - current_kmer_length + current_kmer_length / 2;
                                end_position = position + 1 + current_kmer_length % 2;
                            }
                            else
                            {
                                position = Nmer_position - current_kmer_length + 1;
                                end_position = Nmer_position;
                            }
                            if (position < 0) position = 0;
                            if (end_position > Nlength - current_kmer_length - deletion_size) end_position = Nlength - current_kmer_length - deletion_size;
                            for(; position < end_position; position++)
                            {
                                results[file_number][current_kmer_length][current_kmer_length-Nmer_position+position][deletion_size][((deleted_sequence_value_ULL & mask_ULL[current_kmer_length][position])) >> (position * 2)]++;
                            }
                        }
                    }
                }
            }
        }
    }

    /* PRINTS KMER TABLE */
    printf("\nKmer_length\tGap_position\tGap_length\tBackground\tSignal\tLocal_max");
    for(current_gap_length = 0, current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++)
    {
        current_gap_position = 0;
        kmer_length_size = current_kmer_length;
        gap_position_size = (Nlength-current_kmer_length-1) * (count_also_spaced_kmers != 0) + 1;
        for(; current_gap_position < kmer_length_size; current_gap_position++)
        {
            for(current_gap_length = 0; current_gap_length < gap_position_size; current_gap_length++)
            {
                if (current_gap_length == 0 && current_gap_position == 0);
                else
                {
                    if(count_also_spaced_kmers == 1 && current_gap_position != current_kmer_length / 2 && current_gap_position != current_kmer_length / 2 + current_kmer_length % 2) continue;
                    if(current_gap_length == 0 && current_gap_position > 0) continue;
                    if(current_gap_position == 0 && current_gap_length > 0) break;
                }
                for(current_kmer = 0; current_kmer < pow(4, current_kmer_length); current_kmer++)
                {
                    /* PRINTS KMER */
                    kmer_incidence = results[1][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
                    if (kmer_incidence > print_cutoff)
                    {
                        printf("\n%i\t%i\t%i\t", current_kmer_length, current_gap_position, current_gap_length);
                        for(position = current_kmer_length-1; position > -1 ; position--)
                        {
                            if(current_kmer_length - position - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");
                            printf("%c", dnaforward[(current_kmer & mask_ULL[1][position]) >> (position * 2)]);
                        }
                        printf("\t%li\t%li", results[0][current_kmer_length][current_gap_position][current_gap_length][current_kmer], kmer_incidence);
                        if (Localmax(results, 1, shortest_kmer, too_long_kmer, current_kmer_length, current_gap_position, current_gap_length, current_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff) == 1) printf("\tlocal_max");
                    }
                }
            }
        }
    }
    printf("\n");
}
