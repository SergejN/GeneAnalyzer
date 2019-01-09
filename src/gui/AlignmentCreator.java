/*
    File:
        AlignmentCreator.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Shows the dialog allowing the user to align sequences.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui;

import algorithms.alignment.Alignment;
import bio.gene.Dataset;
import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import gui.editor.EditorPanel;
import gui.editor.HeaderRenderer;
import java.awt.Dimension;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Vector;
import javax.swing.DefaultListModel;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreePath;
import kernel.Kernel;
import plugin.classes.IAligner;

public class AlignmentCreator extends javax.swing.JFrame
{
    private class Header extends JComponent
    {
        private final int OFFSET_TOP = 10;

        public Header(HeaderRenderer hr)
        {
            hr.setBounds(0, OFFSET_TOP, hr.getWidth(), hr.getHeight());
            this.add(hr);
            this.setPreferredSize(new Dimension(hr.getWidth(), hr.getHeight()));
        }
    };

    private class DoubleClickListener extends MouseAdapter
    {
        private JTree tree = null;
        
        public DoubleClickListener(JTree tree)
        {
            this.tree = tree;
        }
        
        @Override
        public void mousePressed(MouseEvent e) 
        {
            int selRow = tree.getRowForLocation(e.getX(), e.getY());
            TreePath selPath = tree.getPathForLocation(e.getX(), e.getY());
            if( (selRow != -1) && (e.getClickCount()==2) ) 
            {
                if(selPath.getPathCount()==3) // Leaf.
                {
                    Object[] path = selPath.getPath();
                    // Root.
                    Object root = path[0];
                    // Parent -> gene
                    Object gene = path[1];
                    TreeModel model = tree.getModel();
                    int iGeneIndex = model.getIndexOfChild(root, gene);
                    int iStrainIndex = model.getIndexOfChild(gene, path[2]);
                    Entry entry = new Entry();
                    entry.se = dataset.getGeneEntry(iGeneIndex).getStrainEntry(iStrainIndex);
                    entry.strGene = dataset.getGeneEntry(iGeneIndex).getCommonName();
                    entries.add(entry);
                    updateStrainsList();
                }
            }
        }
    };

    private class Entry
    {
        StrainEntry se = null;
        String strGene = null;
    };

    private IAligner[] aligners = null;
    private Vector<Entry> entries = null;
    private Dataset dataset = null;
    private EditorPanel ep = null;
    private Kernel kernel  = null;
    private IWaitDialog old = null;
    

    public AlignmentCreator(Kernel kernel, IWaitDialog old)
    {
        super("Alignment Creator");
        initComponents();
        setLocationRelativeTo(null);
        entries = new Vector<Entry>();
        // Set mouse listener for the JTree.
        jtGenes.addMouseListener(new DoubleClickListener(jtGenes));
        jlStrains.setModel(new DefaultListModel());
        this.kernel = kernel;
        this.old = old;
        this.setGlassPane(new WaitDialog());
    }

    public void showAlignmentDialog(IAligner[] aligners, Dataset ds)
    {
        this.dataset = ds;
        this.aligners = aligners;
        for(IAligner al:aligners)
            jcbAligners.addItem(al.getName());
        // Add strains.
        DefaultMutableTreeNode root = new DefaultMutableTreeNode();
        int nGenes = ds.getGenesCount();
        for(int i=0;i<nGenes;i++)
        {
            GeneEntry ge = ds.getGeneEntry(i);
            DefaultMutableTreeNode gene = new DefaultMutableTreeNode(ge.getCommonName());
            int nStrains = ge.getStrainsCount();
            for(int n=0;n<nStrains;n++)
            {
                StrainEntry se = ge.getStrainEntry(n);
                String strEntry = String.format("%s (%s)", se.getStrainName(), se.getSpeciesName());
                DefaultMutableTreeNode strain = new DefaultMutableTreeNode(strEntry);
                gene.add(strain);
            }
            root.add(gene);
        }
        DefaultTreeModel model = new DefaultTreeModel(root);
        jtGenes.setModel(model);
        setVisible(true);
    }

    private void updateStrainsList()
    {
        DefaultListModel model = (DefaultListModel)jlStrains.getModel();
        model.removeAllElements();
        int nStrains = entries.size();
        for(int i=0;i<nStrains;i++)
        {
            Entry e = entries.get(i);
            model.addElement(String.format("%s: %s (%s)", e.strGene, e.se.getStrainName(), e.se.getSpeciesName()));
        }
        jlStrains.setModel(model);
        jbRemove.setEnabled(nStrains>0);
        jbAlign.setEnabled(nStrains>0);
    }

    /**
     *  Iterates through the strain entry, removes all gaps and updates the regions
     *  boundaries. The method returns a new StrainEntry instance.
     *
     *  @param se
     *  @return
     */
    private StrainEntry removeGaps(StrainEntry se)
    {
        StrainEntry entry = new StrainEntry(se.getSpeciesName(), se.getStrainName());
        entry.setChromosome(se.getChromosome());
        entry.addPopulations(se.listPopulations());
        int nRegs = se.getRegionsCount();
        int offset = 0;
        for(int i=0;i<nRegs;i++)
        {
            GeneRegion or = se.getRegion(i);
            String orig = or.getSequence();
            String seq = orig.replaceAll("-", "");
            GeneRegion reg = new GeneRegion(or.getType());
            reg.setStart(or.getStart()-offset);
            offset += orig.length()-seq.length();
            reg.setEnd(or.getEnd()-offset);
            reg.setSequence(seq);
            entry.addRegion(reg);
        }
        return entry;
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        buttonGroup1 = new javax.swing.ButtonGroup();
        jPanel2 = new javax.swing.JPanel();
        jPanel5 = new javax.swing.JPanel();
        jToolBar1 = new javax.swing.JToolBar();
        jtbShowAll = new javax.swing.JToggleButton();
        jtbMask = new javax.swing.JToggleButton();
        jSeparator1 = new javax.swing.JToolBar.Separator();
        jtbHighlightRegs = new javax.swing.JToggleButton();
        jspAlignment = new javax.swing.JScrollPane();
        jPanel6 = new javax.swing.JPanel();
        jScrollPane3 = new javax.swing.JScrollPane();
        jtaAlignment = new javax.swing.JTextArea();
        jPanel1 = new javax.swing.JPanel();
        jPanel3 = new javax.swing.JPanel();
        jScrollPane1 = new javax.swing.JScrollPane();
        jtGenes = new javax.swing.JTree();
        jPanel4 = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        jcbAligners = new javax.swing.JComboBox();
        jLabel2 = new javax.swing.JLabel();
        jScrollPane2 = new javax.swing.JScrollPane();
        jlStrains = new javax.swing.JList();
        jbRemove = new javax.swing.JButton();
        jbAlign = new javax.swing.JButton();
        jcbRemoveGaps = new javax.swing.JCheckBox();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Alignment Creator");

        jPanel2.setBorder(javax.swing.BorderFactory.createTitledBorder("Alignment"));

        jToolBar1.setFloatable(false);
        jToolBar1.setRollover(true);

        buttonGroup1.add(jtbShowAll);
        jtbShowAll.setSelected(true);
        jtbShowAll.setText("Show all bases");
        jtbShowAll.setFocusable(false);
        jtbShowAll.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jtbShowAll.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jtbShowAll.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jtbShowAllItemStateChanged(evt);
            }
        });
        jToolBar1.add(jtbShowAll);

        buttonGroup1.add(jtbMask);
        jtbMask.setText("Mask identical bases");
        jtbMask.setFocusable(false);
        jtbMask.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jtbMask.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jToolBar1.add(jtbMask);
        jToolBar1.add(jSeparator1);

        jtbHighlightRegs.setSelected(true);
        jtbHighlightRegs.setText("Highlight regions");
        jtbHighlightRegs.setFocusable(false);
        jtbHighlightRegs.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jtbHighlightRegs.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jtbHighlightRegs.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jtbHighlightRegsItemStateChanged(evt);
            }
        });
        jToolBar1.add(jtbHighlightRegs);

        javax.swing.GroupLayout jPanel5Layout = new javax.swing.GroupLayout(jPanel5);
        jPanel5.setLayout(jPanel5Layout);
        jPanel5Layout.setHorizontalGroup(
            jPanel5Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jToolBar1, javax.swing.GroupLayout.DEFAULT_SIZE, 546, Short.MAX_VALUE)
            .addComponent(jspAlignment, javax.swing.GroupLayout.DEFAULT_SIZE, 546, Short.MAX_VALUE)
        );
        jPanel5Layout.setVerticalGroup(
            jPanel5Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel5Layout.createSequentialGroup()
                .addComponent(jToolBar1, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jspAlignment, javax.swing.GroupLayout.DEFAULT_SIZE, 234, Short.MAX_VALUE))
        );

        jtaAlignment.setColumns(20);
        jtaAlignment.setRows(5);
        jScrollPane3.setViewportView(jtaAlignment);

        javax.swing.GroupLayout jPanel6Layout = new javax.swing.GroupLayout(jPanel6);
        jPanel6.setLayout(jPanel6Layout);
        jPanel6Layout.setHorizontalGroup(
            jPanel6Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane3, javax.swing.GroupLayout.DEFAULT_SIZE, 546, Short.MAX_VALUE)
        );
        jPanel6Layout.setVerticalGroup(
            jPanel6Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane3, javax.swing.GroupLayout.DEFAULT_SIZE, 241, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel5, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jPanel6, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addComponent(jPanel5, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel6, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        jtGenes.setDragEnabled(true);
        jtGenes.setRootVisible(false);
        jtGenes.setShowsRootHandles(true);
        jScrollPane1.setViewportView(jtGenes);

        javax.swing.GroupLayout jPanel3Layout = new javax.swing.GroupLayout(jPanel3);
        jPanel3.setLayout(jPanel3Layout);
        jPanel3Layout.setHorizontalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 241, Short.MAX_VALUE)
                .addContainerGap())
        );
        jPanel3Layout.setVerticalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel3Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 265, Short.MAX_VALUE))
        );

        jLabel1.setText("Alignment algorithm:");

        jLabel2.setText("Strains:");

        jScrollPane2.setViewportView(jlStrains);

        jbRemove.setText("Remove");
        jbRemove.setEnabled(false);
        jbRemove.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbRemoveActionPerformed(evt);
            }
        });

        jbAlign.setText("Align");
        jbAlign.setEnabled(false);
        jbAlign.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbAlignActionPerformed(evt);
            }
        });

        jcbRemoveGaps.setText("Remove original gaps");

        javax.swing.GroupLayout jPanel4Layout = new javax.swing.GroupLayout(jPanel4);
        jPanel4.setLayout(jPanel4Layout);
        jPanel4Layout.setHorizontalGroup(
            jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel4Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 253, Short.MAX_VALUE)
                    .addGroup(jPanel4Layout.createSequentialGroup()
                        .addComponent(jcbRemoveGaps)
                        .addContainerGap())
                    .addGroup(jPanel4Layout.createSequentialGroup()
                        .addGroup(jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel1)
                            .addComponent(jLabel2))
                        .addGap(155, 155, 155))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel4Layout.createSequentialGroup()
                        .addComponent(jcbAligners, 0, 243, Short.MAX_VALUE)
                        .addContainerGap())
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel4Layout.createSequentialGroup()
                        .addComponent(jbRemove)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 117, Short.MAX_VALUE)
                        .addComponent(jbAlign)
                        .addContainerGap())))
        );
        jPanel4Layout.setVerticalGroup(
            jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel4Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jcbAligners, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 117, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jcbRemoveGaps)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jbRemove)
                    .addComponent(jbAlign))
                .addContainerGap())
        );

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jPanel4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jPanel3, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                .addComponent(jPanel3, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jPanel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jbAlignActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbAlignActionPerformed
        final JFrame frame = this;
        Runnable r = new Runnable()
        {
            public void run()
            {
                kernel.setWaitDialog((IWaitDialog)frame.getGlassPane());
                IAligner aligner = aligners[jcbAligners.getSelectedIndex()];
                int nStrains = entries.size();
                if(nStrains>aligner.getSupportedSequenceNumber())
                {
                    String strMsg = String.format("<html><b>%s</b> cannot align %d sequences</html>", aligner.getName(), nStrains);
                    JOptionPane.showMessageDialog(frame, strMsg, "Alignment Creator", JOptionPane.ERROR_MESSAGE);
                    return;
                }
                StrainEntry[] tmp = new StrainEntry[nStrains];
                for(int i=0;i<nStrains;i++)
                {
                    tmp[i] = (jcbRemoveGaps.isSelected()) ? removeGaps(entries.get(i).se) : entries.get(i).se;
                }
                Alignment ali = null;
                try
                {
                    ali = aligner.alignStrains(tmp, null);
                }
                catch(Exception e)
                {
                    String strMsg = String.format("<html><b>%s</b> encountered an error while aligning the sequences<br>" +
                            "This can happen when the strains you want to align have annotation errors.</html>", aligner.getName());
                    JOptionPane.showMessageDialog(frame, strMsg, "Alignment Creator", JOptionPane.ERROR_MESSAGE);
                    return;
                }
                if(ali!=null)
                {
                    StringBuffer sb = new StringBuffer();
                    int length = ali.getLength();
                    int nMatches = ali.getMatchesCount();
                    int nGaps = ali.getGapsCount();
                    sb.append(String.format("Alignment type:    %s\n", ali.getType()));
                    sb.append(String.format("Alignment options: %s\n\n", ali.getOptions()));
                    sb.append(String.format("Length:            %d\n", length));
                    sb.append(String.format("Identity:          %d/%d (%.2f%%)\n", nMatches, length, 100*(float)nMatches/(float)length));
                    sb.append(String.format("Gaps:              %d/%d (%.2f%%)\n", nGaps, length, 100*(float)nGaps/(float)length));
                    sb.append(String.format("Alignment score:   %.2f\n\n\n", ali.getScore()));
                    sb.append(ali.toString(100));
                    jtaAlignment.setText(sb.toString());
                    // Draw the regions.
                    GeneEntry ge = new GeneEntry("Alignment", "Alignment");
                    nStrains = ali.getStrainsCount();
                    for(int i=0;i<nStrains;i++)
                        ge.addStrain(ali.getStrainEntry(i));
                    ep = new EditorPanel();
                    ep.displayAnnotation(ge);
                    jspAlignment.setViewportView(ep);
                    jspAlignment.setRowHeaderView(new Header(new HeaderRenderer(ge)));
                }
                kernel.setWaitDialog(old);
            }
        };
        new Thread(r).start();
}//GEN-LAST:event_jbAlignActionPerformed

    private void jbRemoveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbRemoveActionPerformed
        int iIndex = jlStrains.getSelectedIndex();
        if(iIndex<0)
            return;
        ((DefaultListModel)jlStrains.getModel()).remove(iIndex);
        entries.remove(iIndex);
        updateStrainsList();
}//GEN-LAST:event_jbRemoveActionPerformed

    private void jtbShowAllItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jtbShowAllItemStateChanged
        if(ep!=null)
        {
            if(jtbShowAll.isSelected())
                ep.setDisplayMode(EditorPanel.SM_ALL);
            else
                ep.setDisplayMode(EditorPanel.SM_NO_IDENT);
        }
}//GEN-LAST:event_jtbShowAllItemStateChanged

    private void jtbHighlightRegsItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jtbHighlightRegsItemStateChanged
        if(ep!=null)
        {
            int mask = 0;
            if(jtbHighlightRegs.isSelected())
                mask = EditorPanel.HM_REGIONS;
            ep.setHighlightMask(mask);
        }
}//GEN-LAST:event_jtbHighlightRegsItemStateChanged

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup buttonGroup1;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JPanel jPanel4;
    private javax.swing.JPanel jPanel5;
    private javax.swing.JPanel jPanel6;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JScrollPane jScrollPane3;
    private javax.swing.JToolBar.Separator jSeparator1;
    private javax.swing.JToolBar jToolBar1;
    private javax.swing.JButton jbAlign;
    private javax.swing.JButton jbRemove;
    private javax.swing.JComboBox jcbAligners;
    private javax.swing.JCheckBox jcbRemoveGaps;
    private javax.swing.JList jlStrains;
    private javax.swing.JScrollPane jspAlignment;
    private javax.swing.JTree jtGenes;
    private javax.swing.JTextArea jtaAlignment;
    private javax.swing.JToggleButton jtbHighlightRegs;
    private javax.swing.JToggleButton jtbMask;
    private javax.swing.JToggleButton jtbShowAll;
    // End of variables declaration//GEN-END:variables

}
