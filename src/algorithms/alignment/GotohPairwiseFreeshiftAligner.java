/*
    File:
        GotohPairwiseFreeshiftAligner.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Encapsulates the Gotoh pairwise alignment.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package algorithms.alignment;

import plugin.classes.IAligner;
import algorithms.Matrix2D;
import algorithms.SequenceBuffer;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import gui.IWaitDialog;
import java.awt.Point;
import java.util.Collections;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import kernel.ErrorCode;
import plugin.AInitData;


public class GotohPairwiseFreeshiftAligner implements IAligner
{    
    private AInitData initData = null;
    private String strLastErr  = null;
    private Options options    = null;

    public ErrorCode initialize(AInitData initData)
    {
        this.initData = initData;
        return ErrorCode.Ok;
    }

    public String getName()
    {
        return "Smith-Waterman Gotoh free-shift alignment";
    }

    public String getDescription()
    {
        return "Aligns the sequences using the global Gotoh pairwise alignment";
    }

    public int getSupportedSequenceNumber()
    {
        return 2;
    }    

    public String getParamString()
    {
        return "match='<NUMBER>' mismatch='<NUMBER>' gap_open='<NUMBER>' gap_extend='<NUMBER>'";
    }

    public String getLastError()
    {
        return strLastErr;
    }

    public Alignment alignStrains(StrainEntry[] entries, String strParams)
    {
        if(entries==null || entries.length==0)
        {
            strLastErr = "Entries array null or empty";
            return null;
        }
        if(entries.length!=2)
        {
            strLastErr = "Unsupported number of strains";
            return null;
        }
        // Options.
        Options opt = null;
        if(options==null || strParams==null || (options!=null && strParams!=null && !strParams.isEmpty()))
        {
            opt = getOptions(strParams);
            if(opt==null)
            {
                strLastErr = "Cancelled by user";
                return null;
            }
            options = opt;
        }
        else
            opt = options;
        initData.wd.show(IWaitDialog.TYPE.Aligner);
        // Create alignment matrices.
        String strSeq1 = entries[0].getCompleteSequence();
        String strSeq2 = entries[1].getCompleteSequence();
        int nRows = strSeq1.length()+1;
        int nCols = strSeq2.length()+1;
        Matrix2D g = createGMatrix(nRows, nCols);
        Matrix2D e = createEMatrix(nRows, nCols, opt.gap_extend, opt.gap_open);
        Matrix2D f = createFMatrix(nRows, nCols, opt.gap_extend, opt.gap_open);
        Matrix2D d = createDMatrix(nRows, nCols);
        // Traceback matrix.
        Matrix2D tb = new Matrix2D(nRows, nCols);
        // Align.
        try
        {
            for(int i=1;i<nRows;i++)
            {
                for(int j=1;j<nCols;j++)
                {
                    // E matrix.
                    float val_e = Math.max(g.getValue(i, j-1)+opt.gap_open+opt.gap_extend,
                                           Math.max(e.getValue(i, j-1)+opt.gap_extend, f.getValue(i, j-1)+opt.gap_open+opt.gap_extend));
                    e.setValue(i, j, val_e);
                    // F matrix.
                    float val_f = Math.max(g.getValue(i-1, j)+opt.gap_open+opt.gap_extend,
                                           Math.max(f.getValue(i-1, j)+opt.gap_extend, e.getValue(i-1, j)+opt.gap_open+opt.gap_extend));
                    f.setValue(i, j, val_f);
                    // G matrix.
                    char c1 = strSeq1.charAt(i-1);
                    char c2 = strSeq2.charAt(j-1);
                    float val = (c1==c2) ? opt.match : opt.mismatch;
                    if(c1=='X' || c2=='X' || c1=='N' || c2=='N')
                        val = (opt.match+opt.mismatch)/2.0f;
                    float val_g = Math.max(g.getValue(i-1, j-1)+val,
                                           Math.max(e.getValue(i-1, j-1)+val, f.getValue(i-1, j-1)+val));
                    g.setValue(i, j, val_g);
                    // D matrix.
                    float max = Math.max(val_e, Math.max(val_f, val_g));
                    d.setValue(i, j, max);
                    // Traceback.
                    if(val_g>=val_e)
                        tb.setValue(i, j, 'd');
                    else
                        tb.setValue(i, j, 'l');
                    if(val_f>Math.max(val_g, val_e))
                        tb.setValue(i, j, 't');
                }
            }
        }
        catch(Matrix2D.EIllegalSubscriptException exception) {}
        // Find the maximum.
        float fScore = 0.0f;
        Point start = new Point(0, 0);
        try
        {
            Matrix2D lastRow = d.submatrix(nRows-1, nRows-1, 0, nCols-1);
            float max1 = lastRow.getMaximum();
            Matrix2D lastCol = d.submatrix(0, nRows-1, nCols-1, nCols-1);
            float max2 = lastCol.getMaximum();
            if(max1>max2)
            {
                fScore = max1;
                Point p = lastRow.getMaximumIndices();
                start.x = nRows-1;
                start.y = p.y;
            }
            else
            {
                fScore = max2;
                Point p = lastCol.getMaximumIndices();
                start.x = p.x;
                start.y = nCols-1;
            }
        }
        catch(Matrix2D.EIllegalSubscriptException exception){}
        // Traceback.
        entries = traceback(start, tb, strSeq1, strSeq2, entries);
        String str = String.format("Match: %.3f, mismatch: %.3f, gap open: %.3f, gap extend: %.3f",
                opt.match, opt.mismatch, opt.gap_open, opt.gap_extend);
        initData.wd.close();
        return new Alignment(fScore, entries, "Gotoh freeshift pairwise alignment", str);
    }

    private Matrix2D createGMatrix(int nRows, int nCols)
    {
        Matrix2D m = new Matrix2D(nRows, nCols);
        try
        {
            m.setValue(0, 0, 0.0f);
            for(int i=1;i<nRows;i++)
                m.setValue(i, 0, Float.NEGATIVE_INFINITY);
            for(int i=1;i<nCols;i++)
                m.setValue(0, i, Float.NEGATIVE_INFINITY);
        }
        catch(Matrix2D.EIllegalSubscriptException e){}
        return m;
    }

    private Matrix2D createEMatrix(int nRows, int nCols, float ge, float go)
    {
        Matrix2D m = new Matrix2D(nRows, nCols);
        try
        {
            m.setValue(0, 0, Float.NEGATIVE_INFINITY);
            for(int i=1;i<nRows;i++)
                m.setValue(i, 0, Float.NEGATIVE_INFINITY);
            for(int i=1;i<nCols;i++)
                m.setValue(0, i, -i*ge-go);
        }
        catch(Matrix2D.EIllegalSubscriptException e){}
        return m;
    }

    private Matrix2D createFMatrix(int nRows, int nCols, float ge, float go)
    {
        Matrix2D m = new Matrix2D(nRows, nCols);
        try
        {
            m.setValue(0, 0, Float.NEGATIVE_INFINITY);
            for(int i=1;i<nRows;i++)
                m.setValue(i, 0, -i*ge-go);
            for(int i=1;i<nCols;i++)
                m.setValue(0, i, Float.NEGATIVE_INFINITY);
        }
        catch(Matrix2D.EIllegalSubscriptException e){}
        return m;
    }

    private Matrix2D createDMatrix(int nRows, int nCols)
    {
        Matrix2D m = new Matrix2D(nRows, nCols);
        try
        {
            m.setValue(0, 0, 0.0f);
            for(int i=1;i<nRows;i++)
                m.setValue(i, 0, 0);
            for(int i=1;i<nCols;i++)
                m.setValue(0, i, 0);
        }
        catch(Matrix2D.EIllegalSubscriptException e){}
        return m;
    }

    private StrainEntry[] traceback(Point start, Matrix2D tb, String strSeq1, String strSeq2, StrainEntry[] entries)
    {
        int nRows = strSeq1.length()+1;
        int nCols = strSeq2.length()+1;
        SequenceBuffer seq1 = new SequenceBuffer(nRows+nCols);
        SequenceBuffer seq2 = new SequenceBuffer(nRows+nCols);
        Vector<GeneRegion> regs = new Vector<GeneRegion>(nRows+nCols);
        int i=start.x;
        int j=start.y;
        // Append the overhang: since the maximum is in either last row or last column, one of the sequences is longer,
        // append the trailing gaps to the other sequence.
        // Copy the rest of the first sequence.
        if(i<nRows-1)
        {
            int _i = nRows-1;
            while(_i>i)
            {
                seq1.appendBase(strSeq1.charAt(_i-1));
                seq2.appendBase('-');
                regs.add(entries[0].getRegionAtSite(_i));
                _i--;
            }
        }
        // Copy the rest of the second sequence.
        else if(j<nCols-1)
        {
            int _j = nCols-1;
            while(_j>j)
            {
                seq1.appendBase('-');
                seq2.appendBase(strSeq2.charAt(_j-1));
                regs.add(entries[1].getRegionAtSite(_j));
                _j--;
            }
        }
        // Append the aligned part of the sequences.
        while(i>0 && j>0)
        {
            float tmp = 0.0f;
            try
            {
                tmp = tb.getValue(i, j);
            }
            catch(Matrix2D.EIllegalSubscriptException exception){}
            // Match/mismatch: align sequence1 and sequence2.
            if(tmp=='d')
            {
                seq1.appendBase(strSeq1.charAt(i-1));
                seq2.appendBase(strSeq2.charAt(j-1));
                GeneRegion r1 = entries[0].getRegionAtSite(i);
                GeneRegion r2 = entries[1].getRegionAtSite(j);
                if(r1.getType().equalsIgnoreCase(r2.getType()))
                    regs.add(r1);
                else
                    regs.add(null);
                i--;
                j--;
            }
            else if(tmp=='l') // Insertion in the second sequence: introduce a gap into the first sequence.
            {
                seq1.appendBase('-');
                seq2.appendBase(strSeq2.charAt(j-1));
                regs.add(entries[1].getRegionAtSite(j));
                j--;
            }
            else            // Insertion in the first sequence: introduce a gap into the second sequence.
            {
                seq1.appendBase(strSeq1.charAt(i-1));
                seq2.appendBase('-');
                regs.add(entries[0].getRegionAtSite(i));
                i--;
            }
        }
        // Append the overhang of the sequences.
        if(i>0)
        {
            while(i>0)
            {
                seq1.appendBase(strSeq1.charAt(i-1));
                seq2.appendBase('-');
                regs.add(entries[0].getRegionAtSite(i));
                i--;
            }
        }
        else if(j>0)
        {
            while(j>0)
            {
                seq1.appendBase('-');
                seq2.appendBase(strSeq2.charAt(j-1));
                regs.add(entries[1].getRegionAtSite(j));
                j--;
            }
        }
        seq1.reverse();
        seq2.reverse();
        // Create the alignment.
        Collections.reverse(regs);
        GeneRegion[] regions = restoreRegions(regs);
        StrainEntry[] tmp = entries;        
        entries = new StrainEntry[2];
        for(int n=0;n<2;n++)
        {
            entries[n] = new StrainEntry(tmp[n].getSpeciesName(), tmp[n].getStrainName());
            entries[n].setChromosome(tmp[n].getChromosome());
            entries[n].addPopulations(tmp[n].listPopulations());
            for(GeneRegion r:regions)
            {
                GeneRegion reg = new GeneRegion(r.getType());
                reg.setStart(r.getStart());
                reg.setEnd(r.getEnd());
                reg.setSequence((n==0)? seq1.substring(r.getStart()-1, r.getEnd())
                                      : seq2.substring(r.getStart()-1, r.getEnd()));
                entries[n].addRegion(reg);
            }
        }
        return entries;
    }

    private GeneRegion[] restoreRegions(Vector<GeneRegion> regs)
    {
        Vector<GeneRegion> regions = new Vector<GeneRegion>();
        int nRegs = regs.size();
        int i = 0;
        while(i<nRegs)
        {
            GeneRegion prev = regs.get(i);
            GeneRegion reg = new GeneRegion((prev!=null) ? prev.getType() : AMBIGUOUS_REGION);
            reg.setStart(i+1);
            reg.setEnd(i+1);
            i++;
            // Extend current region.
            while(i<nRegs)
            {
                GeneRegion next = regs.get(i);
                // If the next region is not null, and the type of the next region is the same
                // as the one of the current region, extend the current region.
                if(next!=null && next.getType().equalsIgnoreCase(reg.getType()))
                    reg.setEnd(i+1);
                // If the next region is null, the current region should be ambiguous in order
                // to be extended.
                else if(next==null && reg.hasType(AMBIGUOUS_REGION))
                    reg.setEnd(i+1);
                // If none of the above is true, exit the cycle and create a new region.
                else
                    break;
                i++;
            }
            regions.add(reg);
        }
        return regions.toArray(new GeneRegion[regions.size()]);
    }

    private Options getOptions(String strParams)
    {
        // If the parameters are not specified or empty, show the dialog.
        if(strParams==null || strParams.isEmpty())
        {
            return new OptionsDlg().getOptions();
        }
        Pattern p = Pattern.compile("match='(-{0,1}[0-9]+[,\\.]{0,1}[0-9]*)'\\s+"+        // 1
                                    "mismatch='(-{0,1}[0-9]+[,\\.]{0,1}[0-9]*)'\\s+"+     // 2
                                    "gap_open='(-{0,1}[0-9]+[,\\.]{0,1}[0-9]*)'\\s+"+     // 3
                                    "gap_extend='(-{0,1}[0-9]+[,\\.]{0,1}[0-9]*)'",       // 4
                                    Pattern.CASE_INSENSITIVE);
        Matcher m = p.matcher(strParams);
        if(m.find())
        {
            Options opt = new Options();
            opt.match = Float.parseFloat(m.group(1).contains(",") ? m.group(1).replaceAll(",", ".")
                                                                  : m.group(1));
            opt.mismatch = Float.parseFloat(m.group(2).contains(",") ? m.group(2).replaceAll(",", ".")
                                                                  : m.group(2));
            opt.gap_open = Float.parseFloat(m.group(3).contains(",") ? m.group(3).replaceAll(",", ".")
                                                                  : m.group(3));
            opt.gap_extend = Float.parseFloat(m.group(4).contains(",") ? m.group(4).replaceAll(",", ".")
                                                                  : m.group(4));
            return opt;
        }
        else
        {
            return new OptionsDlg().getOptions();
        }
    }

    /********************************************************************************************
    *                               PRIVATE CLASSES                                             *
    ********************************************************************************************/
    private class Options
    {
        float gap_open = -10.0f;
        float gap_extend = -0.5f;
        float match = 5.0f;
        float mismatch = -4.0f;
    };

    private class OptionsDlg extends javax.swing.JDialog
    {
        private Options options;

        public OptionsDlg()
        {
            super((java.awt.Frame)null, true);
            initComponents();
            setLocationRelativeTo(null);
            options = null;
            jeMatch.setText("5");
            jeMismatch.setText("-4");
            jeGapOpen.setText("-10");
            jeGapExt.setText("-0.5");
        }

        public Options getOptions()
        {
            setVisible(true);
            return options;
        }

        private char getKey(char key, int iTextLength)
        {
            if(key=='-' && iTextLength>0)
                return '\0';
            if(!Character.isDigit(key) && key!='.' && key!=',' && key!='-')
                return '\0';
            else
                return key;
        }

        private void initComponents()
        {
            jLabel1 = new javax.swing.JLabel();
            jButton1 = new javax.swing.JButton();
            jButton2 = new javax.swing.JButton();
            jLabel2 = new javax.swing.JLabel();
            jeMatch = new javax.swing.JTextField();
            jeMismatch = new javax.swing.JTextField();
            jLabel3 = new javax.swing.JLabel();
            jLabel4 = new javax.swing.JLabel();
            jeGapExt = new javax.swing.JTextField();
            jeGapOpen = new javax.swing.JTextField();

            setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
            setTitle("Aligner options");
            setModal(true);
            setResizable(false);

            jLabel1.setText("Match:");

            jButton1.setText("OK");
            jButton1.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    jButton1ActionPerformed(evt);
                }
            });

            jButton2.setText("Cancel");
            jButton2.addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    jButton2ActionPerformed(evt);
                }
            });

            jLabel2.setText("Mismatch:");

            jeMatch.addKeyListener(new java.awt.event.KeyAdapter() {
                    @Override
                public void keyTyped(java.awt.event.KeyEvent evt) {
                    jeMatchKeyTyped(evt);
                }
            });

            jeMismatch.addKeyListener(new java.awt.event.KeyAdapter() {
                    @Override
                public void keyTyped(java.awt.event.KeyEvent evt) {
                    jeMismatchKeyTyped(evt);
                }
            });

            jLabel3.setText("Gap open penalty:");

            jLabel4.setText("Gap extention penalty:");

            jeGapExt.addKeyListener(new java.awt.event.KeyAdapter() {
                    @Override
                public void keyTyped(java.awt.event.KeyEvent evt) {
                    jeGapExtKeyTyped(evt);
                }
            });

            jeGapOpen.addKeyListener(new java.awt.event.KeyAdapter() {
                    @Override
                public void keyTyped(java.awt.event.KeyEvent evt) {
                    jeGapOpenKeyTyped(evt);
                }
            });

            javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
            getContentPane().setLayout(layout);
            layout.setHorizontalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addContainerGap()
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addGroup(layout.createSequentialGroup()
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(jLabel1)
                                .addComponent(jLabel2))
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 16, Short.MAX_VALUE)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                .addComponent(jeMismatch)
                                .addComponent(jeMatch, javax.swing.GroupLayout.DEFAULT_SIZE, 71, Short.MAX_VALUE))
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addGroup(layout.createSequentialGroup()
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 8, Short.MAX_VALUE)
                                    .addComponent(jLabel3)
                                    .addGap(26, 26, 26))
                                .addGroup(layout.createSequentialGroup()
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                    .addComponent(jLabel4)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)))
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                .addComponent(jeGapExt)
                                .addComponent(jeGapOpen, javax.swing.GroupLayout.DEFAULT_SIZE, 75, Short.MAX_VALUE)))
                        .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                            .addComponent(jButton2)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(jButton1)))
                    .addContainerGap())
            );
            layout.setVerticalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addContainerGap()
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                        .addGroup(layout.createSequentialGroup()
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(jLabel1)
                                .addComponent(jLabel3)
                                .addComponent(jeMatch, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGap(18, 18, 18)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(jLabel2)
                                .addComponent(jeMismatch, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addGroup(layout.createSequentialGroup()
                            .addComponent(jeGapOpen, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGap(18, 18, 18)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(jeGapExt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(jLabel4))))
                    .addGap(18, 18, 18)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(jButton1)
                        .addComponent(jButton2))
                    .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            );

            pack();
        }

        private void jeMatchKeyTyped(java.awt.event.KeyEvent evt) {
            evt.setKeyChar(getKey(evt.getKeyChar(), jeMatch.getText().length()));
        }

        private void jeMismatchKeyTyped(java.awt.event.KeyEvent evt) {
            evt.setKeyChar(getKey(evt.getKeyChar(), jeMismatch.getText().length()));
        }

        private void jeGapOpenKeyTyped(java.awt.event.KeyEvent evt) {
            evt.setKeyChar(getKey(evt.getKeyChar(), jeGapOpen.getText().length()));
        }

        private void jeGapExtKeyTyped(java.awt.event.KeyEvent evt) {
            evt.setKeyChar(getKey(evt.getKeyChar(), jeGapExt.getText().length()));
        }

        private void jButton2ActionPerformed(java.awt.event.ActionEvent evt) {
            options = null;
            setVisible(false);
        }

        private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {
            options = new Options();
            options.match = Float.parseFloat(jeMatch.getText().contains(",") ? jeMatch.getText().replaceAll(",", ".")
                                                                             : jeMatch.getText());
            options.mismatch = Float.parseFloat(jeMismatch.getText().contains(",") ? jeMismatch.getText().replaceAll(",", ".")
                                                                                   : jeMismatch.getText());
            options.gap_open = Float.parseFloat(jeGapOpen.getText().contains(",") ? jeGapOpen.getText().replaceAll(",", ".")
                                                                                  : jeGapOpen.getText());
            options.gap_extend = Float.parseFloat(jeGapExt.getText().contains(",") ? jeGapExt.getText().replaceAll(",", ".")
                                                                                   : jeGapExt.getText());
            setVisible(false);
        }

        private javax.swing.JButton jButton1;
        private javax.swing.JButton jButton2;
        private javax.swing.JLabel jLabel1;
        private javax.swing.JLabel jLabel2;
        private javax.swing.JLabel jLabel3;
        private javax.swing.JLabel jLabel4;
        private javax.swing.JTextField jeGapExt;
        private javax.swing.JTextField jeGapOpen;
        private javax.swing.JTextField jeMatch;
        private javax.swing.JTextField jeMismatch;
    }
}
