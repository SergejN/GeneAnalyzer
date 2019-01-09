/*
    File:
        GeneRegion.java
 *   
    Revision:
        2.2.0.2
 * 
    Description:
        Represents one sequence region such as exon, intron etc.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */


package bio.gene;

import algorithms.SequenceBuffer;
import java.util.HashMap;


public class GeneRegion extends DataChunk implements Comparable<GeneRegion>
{
    // Predefined region types.
    public static final String CDS          = "CDS";
    public static final String EXON         = "Exon";
    public static final String INTRON       = "Intron";
    public static final String UTR5         = "5'UTR";
    public static final String UTR3         = "3'UTR";
    public static final String INTERGENIC   = "Intergenic";
    public static final String mRNA         = "mRNA";    
    public static final String UNNAMED      = "Unnamed";
    
    
    /**
     *  Type of the region. This property can be either a
     *  predefined type, such as INTRON or EXON, or be a
     *  user-defined name.
     */
    protected String strType = "";
    
    /**
     *  DNA sequence of the region. The sequence is represented internally as an 
     *  uppercase string.
     *  To save the memory, the internal representation uses SequenceBuffer
     *  instead of a String.
     */
    protected SequenceBuffer sequence = null;
    
    /**
     *  Start and end positions of the region in the complete gene sequence
     *  (StrainEntry instance). The first position is defined to be 1, and not
     *  0.
     */
    protected int iStartPos = 0;
    protected int iEndPos = 0;
    
    
    /**
     *  Constructor. If the strtype parameter is null,
     *  the type of the region is set to "Unnamed region".
     * 
     *  @param strType
     */
    public GeneRegion(String strType)
    {
        properties = new HashMap<String, Object>();
        if(strType==null || strType.isEmpty())
            this.strType = UNNAMED;
        else
            this.strType = strType;
        sequence = new SequenceBuffer();
    }
    
    
    /**
     *  Returns TRUE if the region has the type specified.
     * 
     *  CAUTION: For convenience reasons the method always returns
     *           false if strType is null or empty.
     * 
     *  @param strType
     *  @return true if region has type strType and false otherwise
     */
    public boolean hasType(String strType)
    {
        if((strType==null) || (strType.isEmpty()) )
            return false;
        return this.strType.equalsIgnoreCase(strType);
    }

    
    /**
     *  Assigns the region a DNA sequence with occasional newline symbols
     *  removed. If the sequence is null, the current sequence is removed.
     *
     *  NOTE: the sequence is not checked for allowed characters!
     * 
     *  @param strSequence
     */
    public void setSequence(String strSequence)
    {
        if(strSequence!=null)
            sequence = new SequenceBuffer(strSequence.replaceAll("\n", "").toUpperCase());
        else
            sequence = new SequenceBuffer("");
    }
    
    
    /**
     *  Returns the DNA sequence of the region. If the sequence was not
     *  assigned properly yet, empty String is returned.
     * 
     *  @return DNA sequence
     */
    public String getSequence()
    {
        return (sequence!=null) ? sequence.toString() : "";
    }

    /**
     *  Returns the length of the sequence of the region.
     * 
     *  @return
     */
    public int getSequenceLength()
    {
        return sequence.length();
    }
    
    /**
     *  Sets the starting position of the region. If the position is negative
     *  the start position is set to 0.
     *
     *  NOTE: indexing begins at 1.
     * 
     *  @param i
     */
    public void setStart(int i)
    {
        iStartPos = (i>0) ? i : 0;
    }
    
    
    /**
     *  Returns the starting position of the region in the genomic sequence.
     *  NOTE: The first position of the sequence is indexed by 1. If the region
     *  was not initialized properly 0 is returned.
     * 
     *  @return Starting position of the region. 
     */
    public int getStart()
    {
        return iStartPos;
    }

    
    /**
     *  Sets the end position of the region. If the specified position
     *  is smaller then or equal to the start position, the end position
     *  is set to be 0.
     * 
     *  @param i
     */
    public void setEnd(int i)
    {
        iEndPos = (i>=iStartPos) ? i : 0;
    }
    
    
    /**
     *  Returns the last position of the region in the sequence.
     * 
     *  @return Last position of the region.
     */
    public int getEnd()
    {
        return iEndPos;
    }

    
    /**
     *  Sets the region type. If strType is null, the current type
     *  is removed.
     * 
     *  @param strType
     */
    public void setType(String strType)
    {
        this.strType = (strType==null) ? "" : strType;
    }
    
    /**
     *  Returns the type of the region.
     * 
     *  @return Region type.
     */
    public String getType()
    {
        return strType;
    }
    
    /**
     *  Replaces the base at the specified site. If the site is invalid the method 
     *  does nothing. The base is converted into uppercase.
     * 
     *  @param iSite
     *  @param base
     *  @return true, if the base was replaced and false otherwise
     */
    public boolean replaceBase(int iSite, char base)
    {
        if(iSite>-1 && iSite<sequence.length())
        {
            sequence.setBaseAt(iSite, Character.toUpperCase(base));
            return true;
        }
        return false;
    }

    /**
     *  Removes the specified base and updates the region end position.
     *
     *  Note:
     *  Removing the base changes the end position of the region. Be sure to
     *  update the boundaries of the regions following the current one!
     *
     *  @param iSite
     *  @return true, if the base was removed and false otherwise
     */
    public boolean removeBase(int iSite)
    {
        if(iSite>-1 && iSite<sequence.length())
        {
            sequence.removeBase(iSite);
            iEndPos-=1;
            return true;
        }
        return false;
    }

    /**
     *  Removes the block of bases of the length nCount starting at site iSite.
     *  Returns the number of removed bases or 0 if the start position is invalid.
     *
     *  Note:
     *  If iStart+nCount is greater than the sequence length, nCount is updated
     *  so that iStart+nCount is equal to the sequence length.
     *  Removing the bases changes the end position of the region. Be sure to
     *  update the boundaries of the regions following the current one!
     *
     *  @param iStart
     *  @param nCount
     *  @return
     */
    public int removeBases(int iStart, int nCount)
    {
        int l = sequence.length();
        if(iStart<0 || iStart>=l || nCount<1)
            return 0;
        else
        {
            int nRemoved = sequence.removeBases(iStart, nCount);
            iEndPos-=nRemoved;
            return nRemoved;
        }
    }

    /**
     *  Inserts the base at the specified site. If the site is invalid the method 
     *  does nothing. The base is converted into uppercase.
     *
     *  @param iSite
     *  @param base
     *  @return true, if the base was inserted and false otherwise
     */
    public boolean insertBase(int iSite, char base)
    {
        if(iSite>-1 && iSite<=sequence.length())
        {
            sequence.insertBase(iSite, base);
            iEndPos+=1;
            return true;
        }
        return false;
    }

    /**
     *  Compares two regions and returns 1 if other region starts earlier in the 
     *  sequence than this does, 0 if both start at the same position, and -1 if
     *  this region starts earlier in the sequence than the other one. 
     * 
     * @param other
     * @return -1,0,1 as in other comparison methods.
     */
    public int compareTo(GeneRegion other)
    {
        return (new Integer(this.iStartPos)).compareTo(new Integer(other.getStart()));
    }
}
