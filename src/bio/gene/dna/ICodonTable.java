/*
    File:
        ICodonTable.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Interface of a codons table.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package bio.gene.dna;


public interface ICodonTable
{
    public static enum TYPE
    {
        FullName,
        OneLetterCode,
        ThreeLetterCode
    };
    
    /**
     *  Returns the name of the codons table.
     * 
     *  @return
     */
    public String getName();
    
    /**
     *  Checks whether or not the specified codon is a terminal codon.
     * 
     *  @param  strCodon
     *  @return true, if the codon is terminal and false otherwise.
     * 
     *  Note:
     *      if the codon contains invalid symbols false is returned.
     */
    public boolean isTerminal(String strCodon);
    
    /**
     *  Returns the amino acid in three-letter code, one-letter code or
     *  the full amino acid name, depending on the type value. 
     * 
     *  Remarks:
     *      Return value for terminal codons.
     *      
     *      type            | Value
     *      FullName        | Terminal
     *      OneLetterCode   | Ter
     *      ThreeLetterCode | *
     * 
     *  @param strCodon
     *  @param type
     *  @return
     */
    public String getAminoAcid(String strCodon, TYPE type);
    
    /**
     *  Checks whether or not the codons encode the same amino acid.
     *  If one of the codons is invalid false is returned.
     * 
     *  @param strCodon1
     *  @param strCodon2
     *  @return
     */
    public boolean areSynonymous(String strCodon1, String strCodon2);
    
    /**
     *  Returns the fold family of the amino acid encoded by the specified
     *  codon.
     * 
     *  @param  strCodon
     *  @return fold of the codon and 0 if the codon is invalid.
     */
    public int getFoldFamily(String strCodon);

    /**
     *  Returns true if the specified codon is a start codon and false otherwise.
     *
     *  @param strCodon
     *  @return
     */
    public boolean isStartCodon(String strCodon);
}
