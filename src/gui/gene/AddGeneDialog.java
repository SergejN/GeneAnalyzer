/*
    File:
        AddGeneDialog.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Shows the dialog allowing the user to add a gene to the dataset.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui.gene;

import bio.gene.GeneEntry;
import bio.gene.GeneRegion;
import bio.gene.StrainEntry;
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;
import javax.swing.JOptionPane;
import javax.swing.table.DefaultTableModel;
import kernel.Main;


public class AddGeneDialog extends javax.swing.JDialog
{
    // Gene entry template.
    private GeneEntry ge = null;
    private Vector<GeneRegion> regs = null;

    public AddGeneDialog()
    {
        super((java.awt.Frame)null, true);
        initComponents();
        setLocationRelativeTo(null);
        ge = new GeneEntry("", "");
        regs = new Vector<GeneRegion>();
    }

    public GeneEntry createGeneEntry()
    {
        setVisible(true);
        return ge;
    }

    private void updateStrainsList()
    {
        DefaultTableModel model = (DefaultTableModel)lvStrains.getModel();
        // Remove old lines.
        for(int i=model.getRowCount()-1;i>=0;i--)
            model.removeRow(i);
        // Add the strains.
        int nStrains = ge.getStrainsCount();
        for(int i=0;i<nStrains;i++)
        {
            StrainEntry se = ge.getStrainEntry(i);
            StringBuffer pops = new StringBuffer();
            String[] tmp = se.listPopulations();
            for(String s:tmp)
                pops.append(s+";");
            model.addRow(new Object[]{se.getStrainName(), se.getSpeciesName(), se.getChromosome(), pops.toString()});
        }
        jbRemoveStrain.setEnabled(nStrains>0);
        jbViewAlignment.setEnabled((nStrains>0) && regs.size()>0);
    }

    private void updateRegionsList()
    {
        // Sort the regions.
        Collections.sort(regs, new Comparator<GeneRegion>()
        {
            public int compare(GeneRegion r1, GeneRegion r2)
            {
                int s1 = r1.getStart();
                int s2 = r2.getStart();
                if(s1==s2)
                    return 0;
                else if(s1<s2)
                    return -1;
                else
                    return 1;
            }
        });
        // Add regions.
        DefaultTableModel model = (DefaultTableModel)lvRegs.getModel();
        // Remove old lines.
        for(int i=model.getRowCount()-1;i>=0;i--)
            model.removeRow(i);
        // Add the strains.
        int nRegions = regs.size();
        for(int i=0;i<nRegions;i++)
        {
            GeneRegion r = regs.get(i);
            model.addRow(new Object[]{r.getType(), Integer.toString(r.getStart()), Integer.toString(r.getEnd())});
        }
        jbRemoveRegion.setEnabled(nRegions>0);
        jbViewAlignment.setEnabled((nRegions>0) && ge.getStrainsCount()>0);
    }

    /**
     *  Checks the data and returns true, if the gene entry can be created and false
     *  if there are some problems.
     *
     *  @return
     */
    private boolean checkData()
    {
        if(ge.getStrainsCount()==0)
        {
            JOptionPane.showMessageDialog(this, "No strains added.", Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            return false;
        }
        if(regs.size()==0)
        {
            JOptionPane.showMessageDialog(this, "No regions specified.", Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            return false;
        }
        // Check if all strains have the same sequence length.
        int nStrains = ge.getStrainsCount();
        int l = ge.getStrainEntry(0).getCompleteSequence().length();
        for(int i=1;i<nStrains;i++)
        {
            if(ge.getStrainEntry(i).getCompleteSequence().length()!=l)
            {
                String strMsg = String.format("Strains 1 (%s) and %d (%s) do not have equal sequence length.",
                        ge.getStrainEntry(0).getStrainName(), i+1, ge.getStrainEntry(i).getStrainName());
                JOptionPane.showMessageDialog(this, strMsg, Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
                return false;
            }
        }
        // Check for regions gaps and overlaps.
        String strMsg = "";
        if(regs.get(0).getStart()>1)
            strMsg = "Unnamed sequence region before the first region";
        int nRegs = regs.size();
        for(int i=0;i<nRegs-1;i++)
        {
            GeneRegion r1 = regs.get(i);
            GeneRegion r2 = regs.get(i+1);
            if(r2.getStart()>r1.getEnd()+1)
            {
                if(!strMsg.isEmpty())
                    strMsg += "\n";
                strMsg += String.format("Gap between regions %d and %d", i+1, i+2);
            }
            if(r1.getEnd()>=r2.getStart())
            {
                if(!strMsg.isEmpty())
                    strMsg += "\n";
                strMsg += String.format("Overlap between regions %d and %d", i+1, i+2);
            }
        }
        if(regs.get(nRegs-1).getEnd()!=l)
        {
            if(!strMsg.isEmpty())
                strMsg += "\n";
            strMsg += "Unnamed sequence region after the last region";
        }
        if(!strMsg.isEmpty())
        {
            JOptionPane.showMessageDialog(this, strMsg, Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            return false;
        }
        return true;
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel1 = new javax.swing.JLabel();
        jeCommonName = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        jeAlias = new javax.swing.JTextField();
        jPanel1 = new javax.swing.JPanel();
        jScrollPane1 = new javax.swing.JScrollPane();
        lvStrains = new javax.swing.JTable();
        jbAddStrain = new javax.swing.JButton();
        jbRemoveStrain = new javax.swing.JButton();
        jPanel2 = new javax.swing.JPanel();
        jScrollPane2 = new javax.swing.JScrollPane();
        lvRegs = new javax.swing.JTable();
        jbAddRegion = new javax.swing.JButton();
        jbRemoveRegion = new javax.swing.JButton();
        jButton1 = new javax.swing.JButton();
        jButton2 = new javax.swing.JButton();
        jbViewAlignment = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Add new gene");
        setResizable(false);

        jLabel1.setText("Common name:");

        jLabel2.setText("Alias:");

        jPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Strains"));

        lvStrains.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {

            },
            new String [] {
                "Name", "Species", "Chromosome", "Populations"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class
            };
            boolean[] canEdit = new boolean [] {
                false, false, false, false
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
        });
        jScrollPane1.setViewportView(lvStrains);

        jbAddStrain.setText("Add");
        jbAddStrain.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbAddStrainActionPerformed(evt);
            }
        });

        jbRemoveStrain.setText("Remove");
        jbRemoveStrain.setEnabled(false);
        jbRemoveStrain.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbRemoveStrainActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 460, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(jbAddStrain, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jbRemoveStrain, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addComponent(jbAddStrain)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jbRemoveStrain)
                .addContainerGap(96, Short.MAX_VALUE))
            .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 148, Short.MAX_VALUE)
        );

        jPanel2.setBorder(javax.swing.BorderFactory.createTitledBorder("Regions"));

        lvRegs.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {

            },
            new String [] {
                "Type", "Start", "End"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.String.class, java.lang.String.class, java.lang.String.class
            };
            boolean[] canEdit = new boolean [] {
                false, false, false
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
        });
        jScrollPane2.setViewportView(lvRegs);

        jbAddRegion.setText("Add");
        jbAddRegion.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbAddRegionActionPerformed(evt);
            }
        });

        jbRemoveRegion.setText("Remove");
        jbRemoveRegion.setEnabled(false);
        jbRemoveRegion.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbRemoveRegionActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel2Layout.createSequentialGroup()
                .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 460, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jbRemoveRegion, javax.swing.GroupLayout.DEFAULT_SIZE, 71, Short.MAX_VALUE)
                    .addComponent(jbAddRegion, javax.swing.GroupLayout.DEFAULT_SIZE, 71, Short.MAX_VALUE)))
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addComponent(jbAddRegion)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jbRemoveRegion)
                .addContainerGap(77, Short.MAX_VALUE))
            .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 129, Short.MAX_VALUE)
        );

        jButton1.setText("Cancel");
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });

        jButton2.setText("Add");
        jButton2.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton2ActionPerformed(evt);
            }
        });

        jbViewAlignment.setText("View alignment");
        jbViewAlignment.setEnabled(false);
        jbViewAlignment.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbViewAlignmentActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jeCommonName, javax.swing.GroupLayout.PREFERRED_SIZE, 209, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(jLabel2)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jeAlias, javax.swing.GroupLayout.PREFERRED_SIZE, 196, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addComponent(jbViewAlignment)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 328, Short.MAX_VALUE)
                        .addComponent(jButton2)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jButton1))
                    .addComponent(jPanel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(jeCommonName, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel2)
                    .addComponent(jeAlias, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(1, 1, 1)
                .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jButton1)
                    .addComponent(jButton2)
                    .addComponent(jbViewAlignment))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jbAddStrainActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbAddStrainActionPerformed
        StrainEntry tmp = (new AddStrainDialog()).createStrainEntry();
        if(tmp==null)
            return;
        if(!ge.hasStrain(tmp.getStrainName()))
        {
            ge.addStrain(tmp);
            updateStrainsList();
        }
        else
        {
            String strMsg = String.format("<html>The strain <b>%s</b> already exists.</html>", tmp.getStrainName());
            JOptionPane.showMessageDialog(this, strMsg, Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
        }
}//GEN-LAST:event_jbAddStrainActionPerformed

    private void jbRemoveStrainActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbRemoveStrainActionPerformed
        int[] indices = lvStrains.getSelectedRows();
        int n = indices.length;
        if(n>0)
        {
            for(int i=n-1;i>=0;i--)
                ge.removeStrain(indices[i]);
        }
        updateStrainsList();
    }//GEN-LAST:event_jbRemoveStrainActionPerformed

    private void jbRemoveRegionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbRemoveRegionActionPerformed
        int[] indices = lvRegs.getSelectedRows();
        int n = indices.length;
        if(n>0)
        {
            for(int i=n-1;i>=0;i--)
                regs.remove(indices[i]);
        }
        updateRegionsList();
    }//GEN-LAST:event_jbRemoveRegionActionPerformed

    private void jbAddRegionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbAddRegionActionPerformed
        GeneRegion reg = (new AddGeneRegionDialog()).createGeneRegion();
        if(reg!=null)
        {
            regs.add(reg);
            updateRegionsList();
        }
    }//GEN-LAST:event_jbAddRegionActionPerformed

    private void jbViewAlignmentActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbViewAlignmentActionPerformed
        if(checkData())
        {
            String strName = jeCommonName.getText();
            GeneEntry tmp = new GeneEntry((strName.isEmpty() ? "New gene" : strName), jeAlias.getText());
            int nStrains = ge.getStrainsCount();
            for(int i=0;i<nStrains;i++)
            {
                StrainEntry se = ge.getStrainEntry(i);
                StrainEntry e = new StrainEntry(se.getStrainName(), se.getSpeciesName());
                e.setChromosome(se.getChromosome());
                e.addPopulations(se.listPopulations());
                String strSeq = se.getCompleteSequence();
                int nRegs = regs.size();
                for(int n=0;n<nRegs;n++)
                {
                    GeneRegion reg = regs.get(n);
                    GeneRegion r = new GeneRegion(reg.getType());
                    r.setStart(reg.getStart());
                    r.setEnd(reg.getEnd());
                    r.setSequence(strSeq.substring(r.getStart()-1, r.getEnd()));
                    e.addRegion(r);
                }
                tmp.addStrain(e);
            }
            (new ViewerDialog()).showGeneEntry(tmp);
        }
    }//GEN-LAST:event_jbViewAlignmentActionPerformed

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        ge = null;
        setVisible(false);
    }//GEN-LAST:event_jButton1ActionPerformed

    private void jButton2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton2ActionPerformed
        if(jeCommonName.getText().isEmpty())
        {
            JOptionPane.showMessageDialog(this, "Specify the gene name.", Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            return;
        }
        if(checkData())
        {
            GeneEntry tmp = new GeneEntry((jeCommonName.getText()), jeAlias.getText());
            int nStrains = ge.getStrainsCount();
            for(int i=0;i<nStrains;i++)
            {
                StrainEntry se = ge.getStrainEntry(i);
                StrainEntry e = new StrainEntry(se.getStrainName(), se.getSpeciesName());
                e.setChromosome(se.getChromosome());
                e.addPopulations(se.listPopulations());
                String strSeq = se.getCompleteSequence();
                int nRegs = regs.size();
                for(int n=0;n<nRegs;n++)
                {
                    GeneRegion reg = regs.get(n);
                    GeneRegion r = new GeneRegion(reg.getType());
                    r.setStart(reg.getStart());
                    r.setEnd(reg.getEnd());
                    r.setSequence(strSeq.substring(r.getStart()-1, r.getEnd()));
                    e.addRegion(r);
                }
                tmp.addStrain(e);
            }
            ge = tmp;
            setVisible(false);
        }
    }//GEN-LAST:event_jButton2ActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton jButton1;
    private javax.swing.JButton jButton2;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JButton jbAddRegion;
    private javax.swing.JButton jbAddStrain;
    private javax.swing.JButton jbRemoveRegion;
    private javax.swing.JButton jbRemoveStrain;
    private javax.swing.JButton jbViewAlignment;
    private javax.swing.JTextField jeAlias;
    private javax.swing.JTextField jeCommonName;
    private javax.swing.JTable lvRegs;
    private javax.swing.JTable lvStrains;
    // End of variables declaration//GEN-END:variables

}
