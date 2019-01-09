/*
    File:
        SequenceFragment.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Extends the standard GeneRegion class. Sequence fragment does not
        contain the sequence but the fragment start and end positions in the
        alignment instead. This approach helps to keep the memory demand as low
        as possible, since the sequence is not duplicated.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */


package dpgp;

import algorithms.SequenceBuffer;
import bio.gene.GeneRegion;
import java.util.HashMap;
import java.util.Vector;


public class SequenceFragment extends GeneRegion
{
    /**
     *  Since the SequenceFragment class does not contain the sequence but rather
     *  uses the alignment, the sequence is extracted from.
     *  Thus, the the methods replaceBase, removeBase and insertBase  do not change
     *  the sequence directly, but are saved instead and carried out, when the
     *  method getSequence is called.
     */
    /**
     *  Operation type.
     */
    private enum OperationType
    {
        Insert,
        Remove,
        Replace
    };

    /**
     *  Sequence edit operation.
     */
    private class Operation
    {
        public OperationType type;
        int iSite;
        char base;
    };

    private DPGPAlignment ali = null;   // Corresponding alignment
    private int isi       = -1;     // Strain index

    private int iGenStartPos = -1;
    private int iGenEndPos   = -1;

    private Vector<Operation> op_cache = null;  // Operations cache.
    private int iLength                = 0;

    /**
     *  Constructs a new sequence fragment.
     *
     *  Remarks:
     *      The sequence fragment is firmly associated with the alignment.
     *      I.e. if the alignment is changed, the change has a direct
     *      impact on the sequence fragment. Thus, DO NOT delete or change
     *      the alignment while using the sequence fragment instance.
     *
     *  @param strType
     *  @param alignment
     *  @param iIndex
     *  @param iStart
     *  @param iEnd
     */
    public SequenceFragment(String strType, DPGPAlignment alignment, int iIndex, int iStart, int iEnd)
    {
        super(strType);
        this.ali = alignment;
        this.isi = iIndex;
        this.iGenStartPos = iStart;
        this.iGenEndPos = iEnd;
        this.properties = new HashMap<String, Object>();
        this.op_cache = new Vector<Operation>();
        this.iLength = getSequenceLength();
    }

    @Override
    public void setSequence(String strSequence)
    {
        this.sequence = new SequenceBuffer(strSequence);
    }

    /**
     *  Returns the DNA sequence of the region. If the sequence was not
     *  assigned properly yet, the method returns null.
     *
     *  @return DNA sequence
     */
    @Override
    public String getSequence()
    {
        if(!this.sequence.isEmpty())
            return this.sequence.toString();
        else if(ali!=null)
        {
            SequenceBuffer tmp = ali.extractRegion(iGenStartPos, iGenEndPos, isi);
            if(tmp!=null)
                return getEditedSequence(tmp).toString();
            else
                return "";
        }
        else
            return "";
    }

    @Override
    public int getSequenceLength()
    {
        if(ali!=null && iLength==0)
        {
            SequenceBuffer tmp = ali.extractRegion(iGenStartPos, iGenEndPos, isi);
            if(tmp!=null)
                iLength = getEditedSequence(tmp).length();
        }
        return iLength;
    }

    @Override
    public boolean replaceBase(int iSite, char base)
    {
        Operation op = new Operation();
        op.type = OperationType.Replace;
        op.base = base;
        op.iSite = iSite;
        return addCachedOperation(op);
    }

    @Override
    public boolean removeBase(int iSite)
    {
        Operation op = new Operation();
        op.type = OperationType.Remove;
        op.base = 0;
        op.iSite = iSite;
        return addCachedOperation(op);
    }

    @Override
    public int removeBases(int iStart, int nCount)
    {
        if(iStart<0 || iStart>=iLength || nCount<1)
            return 0;
        int nRemoved = 0;
        for(int i=0;i<nCount;i++)
        {
            if(removeBase(iStart+i-nRemoved))
                nRemoved++;
        }
        return nRemoved;
    }

    @Override
    public boolean insertBase(int iSite, char base)
    {
        Operation op = new Operation();
        op.type = OperationType.Insert;
        op.base = base;
        op.iSite = iSite;
        return addCachedOperation(op);
    }

    /**
     *  Returns the genomic start position.
     *
     *  @return
     */
    public int getGenomicStart()
    {
        return iGenStartPos;
    }

    /**
     *  Returns the genomic end position.
     *
     *  @return
     */
    public int getGenomicEnd()
    {
        return iGenEndPos;
    }  
    
    /**
     *  Returns the number of gaps between the specified genomic positions.
     * 
     *  @param iGenomicFrom
     *  @param iGenomicTo
     *  @return
     */
    public int getGapsCount(int iGenomicFrom, int iGenomicTo)
    {
        SequenceBuffer seq = ali.extractRegion(iGenomicFrom, iGenomicTo, 0);
        return seq.getCharacterCount(DPGPAlignment.GAP);
    }

    /**
     *  Returns the edited sequence, i.e. the sequence with all cached operations
     *  integrated.
     *
     *  @param seq
     *  @return
     */
    private SequenceBuffer getEditedSequence(SequenceBuffer seq)
    {
        for(Operation op:op_cache)
        {
            switch(op.type)
            {
                case Insert:
                {
                    seq.insertBase(op.iSite, op.base);
                    break;
                }
                case Remove:
                {
                    seq.removeBase(op.iSite);
                    break;
                }
                case Replace:
                    seq.setBaseAt(op.iSite, op.base);
                    break;
            }
        }
        return seq;
    }

    /**
     *  Checks the operation and, if it is okay, adds it to the list
     *  and updates another operations.
     *
     *  @param op
     *  @return true, if the operation was cached and false otherwise.
     */
    private boolean addCachedOperation(Operation op)
    {
        // The operation must specify a valid site.
        if(op.iSite>-1 && (op.type==OperationType.Insert) ? op.iSite<=iLength : op.iSite<iLength)
        {
            // Update the sequence length and the end position.
            switch(op.type)
            {
                case Insert:
                    iLength++;
                    iEndPos++;
                    break;
                case Remove:
                    iLength--;
                    iEndPos--;
                    break;
                case Replace:
                    break;
            }
            op_cache.add(op);
            return true;
        }
        return false;
    }
}
