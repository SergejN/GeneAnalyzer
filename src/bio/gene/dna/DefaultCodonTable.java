/*
    File:
        DefaultCodonTable.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Implements the genetic code table.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package bio.gene.dna;

import java.util.HashMap;

public class DefaultCodonTable implements ICodonTable
{
    /**
     *  Amino acid. 
     */
    private class Codon
    {
        private String  strName;
        private String  strTLC;
        private String  strOLC;
        private int     iFold;
        private boolean bTerminal;
        
        public Codon(String strName, String strTLC, 
                         String strOLC, int iFold, boolean bTerminal)
        {
            this.strName = strName;
            this.strTLC = strTLC;
            this.strOLC = strOLC;
            this.iFold = iFold;
            this.bTerminal = bTerminal;
        }
    };    
     
    /**
     *  Maps the codons to the corresponding amino acids.
     */
    private HashMap<String, Codon> codons = null;
    
    
    /**
     *  Initialize the codons table with codons and
     *  three-letter codes and the folds table with
     *  fold family values.
     *  NOTE: all codons are DNA codons and not RNA codons.
     */
    public DefaultCodonTable()
    {
        codons = new HashMap<String, Codon>();
        // 1. Alanine.
        Codon aa = new Codon("Alanine", "Ala", "A", 4, false);
        codons.put("GCA", aa); codons.put("GCC", aa); 
        codons.put("GCG", aa); codons.put("GCT", aa);
        // 2. Arginine.
        aa = new Codon("Arginine", "Arg", "R", 4, false);
        codons.put("CGA", aa); codons.put("CGC", aa); 
        codons.put("CGG", aa); codons.put("CGT", aa);
        aa = new Codon("Arginine", "Arg", "R", 2, false);
        codons.put("AGA", aa); codons.put("AGG", aa); 
        // 3. Asparagine.
        aa = new Codon("Asparagine", "Asn", "N", 2, false);
        codons.put("AAT", aa); codons.put("AAC", aa); 
        // 4. Aspartate.
        aa = new Codon("Aspartate", "Asp", "D", 2, false);
        codons.put("GAT", aa); codons.put("GAC", aa); 
        // 5. Cystein.
        aa = new Codon("Cystein", "Cys", "C", 2, false);
        codons.put("TGT", aa); codons.put("TGC", aa); 
        // 6. Glutamine.
        aa = new Codon("Glutamine", "Gln", "Q", 2, false);
        codons.put("CAA", aa); codons.put("CAG", aa); 
        // 7. Glutamate.
        aa = new Codon("Glutamate", "Glu", "E", 2, false);
        codons.put("GAA", aa); codons.put("GAG", aa); 
        // 8. Glycine.
        aa = new Codon("Glycine", "Gly", "G", 4, false);
        codons.put("GGA", aa); codons.put("GGC", aa); 
        codons.put("GGG", aa); codons.put("GGT", aa);
        // 9. Histidine.
        aa = new Codon("Histidine", "His", "H", 2, false);
        codons.put("CAT", aa); codons.put("CAC", aa);
        // 10. Isoleucine.
        aa = new Codon("Isoleucine", "Ile", "I", 3, false);
        codons.put("ATT", aa); codons.put("ATC", aa); 
        codons.put("ATA", aa); 
        // 11. Leucine.
        aa = new Codon("Leucine", "Leu", "L", 4, false);
        codons.put("CTA", aa); codons.put("CTC", aa); 
        codons.put("CTG", aa); codons.put("CTT", aa);
        aa = new Codon("Leucine", "Leu", "L", 2, false);
        codons.put("TTA", aa); codons.put("TTG", aa);
        // 12. Lysine.
        aa = new Codon("Lysine", "Lys", "K", 2, false);
        codons.put("AAA", aa); codons.put("AAG", aa); 
        // 13. Methionine.
        aa = new Codon("Methionine", "Met", "M", 1, false);
        codons.put("ATG", aa);
        // 14. Phenylalanine.
        aa = new Codon("Phenylalanine", "Phe", "F", 2, false);
        codons.put("TTT", aa); codons.put("TTC", aa);
        // 15. Proline.
        aa = new Codon("Proline", "Pro", "P", 4, false);
        codons.put("CCA", aa); codons.put("CCC", aa); 
        codons.put("CCG", aa); codons.put("CCT", aa);
        // 16. Serine.
        aa = new Codon("Serine", "Ser", "S", 4, false);
        codons.put("TCA", aa); codons.put("TCC", aa); 
        codons.put("TCG", aa); codons.put("TCT", aa);
        aa = new Codon("Serine", "Ser", "S", 2, false);
        codons.put("AGC", aa); codons.put("AGT", aa);
        // 17. Threonine.
        aa = new Codon("Threonine", "Thr", "T", 4, false);
        codons.put("ACA", aa); codons.put("ACC", aa); 
        codons.put("ACG", aa); codons.put("ACT", aa); 
        // 18. Tryptophan.
        aa = new Codon("Tryptophan", "Trp", "W", 1, false);
        codons.put("TGG", aa);
        // 19. Tyrosine.
        aa = new Codon("Tyrosine", "Tyr", "Y", 2, false);
        codons.put("TAC", aa); codons.put("TAT", aa);
        // 20. Valine.
        aa = new Codon("Valine", "Val", "V", 4, false);
        codons.put("GTA", aa); codons.put("GTC", aa); 
        codons.put("GTG", aa); codons.put("GTT", aa); 
        // Terminal codons.
        aa = new Codon("Terminal", "Ter", "*", 2, true);
        codons.put("TAA", aa); codons.put("TAG", aa); 
        aa = new Codon("Terminal", "Ter", "*", 1, true);
        codons.put("TGA", aa);
    }  
    
    public String getName()
    {
        return "Universal genetic code";
    }
        
    /**
     *  Checks whether or not the specified codon is a terminal codon.
     * 
     *  @param codon
     *  @return true, if the codon is terminal and false otherwise.
     * 
     *  Note:
     *      if the codon contains invalid symbols false is returned.
     */
    public boolean isTerminal(String codon)
    {
        Codon aa = codons.get(codon);
        return (aa!=null) ? aa.bTerminal : false;
    }    

    /**
     *  Returns the name/abbreviation of the amino acid.
     *  If the codon is not a valid codon, the method returns null.
     *  
     *  Remarks:
     *  The valid codon can only be built from 4 bases: A,C,G,T.
     *  Complete UIB is not supported.
     *
     *  @param strCodon
     *  @param type
     *  @return
     */
    public String getAminoAcid(String strCodon, TYPE type)
    {
        Codon aa = codons.get(strCodon);
        if(aa==null)
            return null;
        else
        {
            switch(type)
            {
                case FullName:
                    return aa.strName;
                case ThreeLetterCode:
                    return aa.strTLC;
                case OneLetterCode:
                    return aa.strOLC;
                default:
                    return null;
            }
        }
    }
    
    /**
     *  Checks whether or not the codons encode the same amino acid.
     *  If one of the codons is invalid false is returned.
     * 
     *  @param codon1
     *  @param codon2
     *  @return
     */
    public boolean areSynonymous(String codon1, String codon2)
    {
        Codon aa1 = codons.get(codon1.toUpperCase());
        Codon aa2 = codons.get(codon2.toUpperCase());
        return (aa1!=null && aa2!=null) ? aa1.strOLC.equals(aa2.strOLC) : false;
    }    
    
    /**
     *  Returns the fold family of the amino acid encoded by the specified
     *  codon.
     * 
     *  @param strCodon
     *  @return fold of the codon and 0 if the codon is invalid.
     */
    public int getFoldFamily(String strCodon)
    {
        Codon aa = codons.get(strCodon);
        return (aa!=null) ? aa.iFold : 0;
    }

    /**
     *  Returns true if the specified codon is a start codon and false otherwise.
     *
     *  @param strCodon
     *  @return
     */
    public boolean isStartCodon(String strCodon)
    {
        return strCodon.equalsIgnoreCase("ATG");
    }
}
