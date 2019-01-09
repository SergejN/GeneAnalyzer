/*
    File:
        EditorForm.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        GeneEntry Editor main form. 

 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui.editor;

import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import javax.swing.JComponent;
import javax.swing.JFrame;

public class EditorForm extends JFrame 
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
    
    private EditorPanel ep = null;
    private GeneEntry ge   = null;
    

    public EditorForm() 
    {
        initComponents();
        setLocationRelativeTo(null);
        setExtendedState(JFrame.MAXIMIZED_BOTH);
        setGlassPane(new EditorGlassPanel());
    }
    
    public void editEntry(GeneEntry ge)
    {
        setTitle(ge.getCommonName());
        ep = new EditorPanel();
        ep.displayAnnotation(ge);
        jspAlignment.setViewportView(ep);
        jspAlignment.setRowHeaderView(new Header(new HeaderRenderer(ge)));
        this.ge = ge;
        setVisible(true);
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
        jspAlignment = new javax.swing.JScrollPane();
        jToolBar1 = new javax.swing.JToolBar();
        jtbShowAll = new javax.swing.JToggleButton();
        jtbMask = new javax.swing.JToggleButton();
        jSeparator1 = new javax.swing.JToolBar.Separator();
        jtbHighlightRegs = new javax.swing.JToggleButton();
        jSeparator2 = new javax.swing.JToolBar.Separator();
        jtbEdit = new javax.swing.JButton();

        setTitle("Editor");
        addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyPressed(java.awt.event.KeyEvent evt) {
                formKeyPressed(evt);
            }
        });

        jspAlignment.setBackground(new java.awt.Color(255, 255, 255));
        jspAlignment.setBorder(null);
        jspAlignment.setForeground(new java.awt.Color(255, 255, 255));
        jspAlignment.setDoubleBuffered(true);

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
        jToolBar1.add(jSeparator2);

        jtbEdit.setText("Edit regions");
        jtbEdit.setFocusable(false);
        jtbEdit.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jtbEdit.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jtbEdit.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jtbEditActionPerformed(evt);
            }
        });
        jToolBar1.add(jtbEdit);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jToolBar1, javax.swing.GroupLayout.DEFAULT_SIZE, 540, Short.MAX_VALUE)
            .addComponent(jspAlignment, javax.swing.GroupLayout.DEFAULT_SIZE, 540, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(jToolBar1, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jspAlignment, javax.swing.GroupLayout.DEFAULT_SIZE, 392, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

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

private void jtbEditActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jtbEditActionPerformed
    if(ge.getStrainsCount()<1)
        return;
    StrainEntry se = ge.getStrainEntry(0);
    int nRegs = se.getRegionsCount();
    if(nRegs==0)
        return;
    GeneRegion[] regs = new GeneRegion[nRegs];
    for(int i=0;i<regs.length;i++)
        regs[i] = se.getRegion(i);
    ((EditorGlassPanel)getGlassPane()).displayRegions(ge, ep);
}//GEN-LAST:event_jtbEditActionPerformed

private void formKeyPressed(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_formKeyPressed
    if(evt.getKeyCode()==java.awt.event.KeyEvent.VK_ESCAPE)
        getGlassPane().setVisible(false);
}//GEN-LAST:event_formKeyPressed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup buttonGroup1;
    private javax.swing.JToolBar.Separator jSeparator1;
    private javax.swing.JToolBar.Separator jSeparator2;
    private javax.swing.JToolBar jToolBar1;
    private javax.swing.JScrollPane jspAlignment;
    private javax.swing.JButton jtbEdit;
    private javax.swing.JToggleButton jtbHighlightRegs;
    private javax.swing.JToggleButton jtbMask;
    private javax.swing.JToggleButton jtbShowAll;
    // End of variables declaration//GEN-END:variables


/*****************************************************************************************
 *        PRIVATE CLASSES                                                                *
 ****************************************************************************************/
    public class EditorGlassPanel extends JComponent
    {
        private final int TABLE_WIDTH   = 389;
        private final int TABLE_HEIGHT  = 331;

        private RegionsPanel regs    = null;


        public EditorGlassPanel()
        {
            setLayout(null);
            regs = new RegionsPanel();
            add(regs);
        }

        /**
         *  Displays the panel allows to edit the regions.
         *
         *  @param ge
         *  @param ep
         *  @return
         */
        public void displayRegions(GeneEntry ge, EditorPanel ep)
        {
            regs.editRegions(ge, ep);
            setVisible(true);
        }

        @Override
        public void paint(Graphics g)
        {
            Graphics2D g2d = (Graphics2D)g;
            g2d.setColor(new Color(225,225,225,170));
            g2d.fillRect(0, 0, getWidth(), getHeight());
            // Move the components.
            int w = getWidth();
            int h = getHeight();
            int tx = (w-TABLE_WIDTH)/2;
            int ty = (h-TABLE_HEIGHT)/2;
            regs.setBounds(tx, ty, TABLE_WIDTH, TABLE_HEIGHT);
            super.paint(g);
        }
    }

}