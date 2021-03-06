/*
    File:
        CodonTableCreator.java
 *   
    Revision:
        1.0.0.1
 * 
    Description:
        Shows the dialog allowing the user to create a custom codon table.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui;

import bio.gene.dna.Codon;
import bio.gene.dna.CustomCodonTable;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;


public class CodonTableCreator extends javax.swing.JDialog
{
    private class AminoAcid
    {
        private String  strName;
        private String  strTLC;
        private String  strOLC;
        
        public AminoAcid(String strName, String strTLC, String strOLC)
        {
            this.strName = strName;
            this.strTLC = strTLC;
            this.strOLC = strOLC;
        }
    };
    
    private class MenuItemActionListener implements ActionListener
    {
        private int iIndex = 0;
        
        public MenuItemActionListener(int iIndex)
        {
            this.iIndex = iIndex;
        }
        
        public void actionPerformed(ActionEvent e)
        {
            DefaultTableModel model = (DefaultTableModel)lvCodons.getModel();  
            int[] idx = lvCodons.getSelectedRows();            
            for(int i=0;i<idx.length;i++)
            {
                model.setValueAt(aas[iIndex].strName, idx[i], 3);
                model.setValueAt(aas[iIndex].strTLC, idx[i], 4);
                model.setValueAt(aas[iIndex].strOLC, idx[i], 5);
            }
        }
    };

    
    private AminoAcid[] aas = null;    
    private CustomCodonTable ct  = null;
    

    public CodonTableCreator() 
    {
        initComponents();
        setLocationRelativeTo(null);
        initializeCodons();
        initializeAminoAcids();
        // Menues.
        for(int i=0;i<aas.length;i++)
        {
            JMenuItem mi = new JMenuItem(aas[i].strName);
            mi.addActionListener(new MenuItemActionListener(i));
            jPopupMenu1.add(mi);
        }
        // Columns.
        lvCodons.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        lvCodons.getColumnModel().getColumn(0).setPreferredWidth(40);   // Start
        lvCodons.getColumnModel().getColumn(1).setPreferredWidth(60);   // Terminal
        lvCodons.getColumnModel().getColumn(2).setPreferredWidth(60);   // Codon
        lvCodons.getColumnModel().getColumn(3).setPreferredWidth(70);   // AA
        lvCodons.getColumnModel().getColumn(4).setPreferredWidth(60);   // TLC
        lvCodons.getColumnModel().getColumn(5).setPreferredWidth(60);   // OLC
    }
        
    public CustomCodonTable createCodonTable()
    {
        setVisible(true);
        return ct;
    }
    
    private void initializeCodons()
    {
        DefaultTableModel model = (DefaultTableModel)lvCodons.getModel(); 
        String[] codons = Codon.generateCodonSequences();
        for(int i=0;i<codons.length;i++)
            model.addRow(new Object[]{new Boolean(false), new Boolean(false), codons[i], "", "", ""});        
    }
    
    private void initializeAminoAcids()
    {
        aas = new AminoAcid[20];
        aas[0] = this.new AminoAcid("Alanine", "Ala", "A");
        aas[1] = this.new AminoAcid("Arginine", "Arg", "R");
        aas[2] = this.new AminoAcid("Asparagine", "Asn", "N");
        aas[3] = this.new AminoAcid("Aspartate", "Asp", "D");
        aas[4] = this.new AminoAcid("Cystein", "Cys", "C");
        aas[5] = this.new AminoAcid("Glutamine", "Gln", "Q");
        aas[6] = this.new AminoAcid("Glutamate", "Glu", "E");
        aas[7] = this.new AminoAcid("Glycine", "Gly", "G");
        aas[8] = this.new AminoAcid("Histidine", "His", "H");
        aas[9] = this.new AminoAcid("Isoleucine", "Ile", "I");
        
        aas[10] = this.new AminoAcid("Leucine", "Leu", "L");
        aas[11] = this.new AminoAcid("Lysine", "Lys", "K");
        aas[12] = this.new AminoAcid("Methionine", "Met", "M");
        aas[13] = this.new AminoAcid("Phenylalanine", "Phe", "F");
        aas[14] = this.new AminoAcid("Proline", "Pro", "P");
        aas[15] = this.new AminoAcid("Serine", "Ser", "S");
        aas[16] = this.new AminoAcid("Threonine", "Thr", "T");
        aas[17] = this.new AminoAcid("Tryptophan", "Trp", "W");
        aas[18] = this.new AminoAcid("Tyrosine", "Tyr", "Y");
        aas[19] = this.new AminoAcid("Valine", "Val", "V");
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPopupMenu1 = new javax.swing.JPopupMenu();
        jLabel1 = new javax.swing.JLabel();
        jeName = new javax.swing.JTextField();
        jSeparator1 = new javax.swing.JSeparator();
        jLabel2 = new javax.swing.JLabel();
        jScrollPane1 = new javax.swing.JScrollPane();
        lvCodons = new javax.swing.JTable();
        jbtnSave = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Codon Table Creator");
        setModal(true);

        jLabel1.setText("Name:");

        jLabel2.setText("Codon table");

        lvCodons.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {

            },
            new String [] {
                "Start", "Terminal", "Codon", "Amino acid", "TLC", "OLC"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.Boolean.class, java.lang.Boolean.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class
            };
            boolean[] canEdit = new boolean [] {
                true, true, false, false, false, false
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
        });
        lvCodons.setComponentPopupMenu(jPopupMenu1);
        lvCodons.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                lvCodonsMouseClicked(evt);
            }
        });
        jScrollPane1.setViewportView(lvCodons);

        jbtnSave.setText("Save");
        jbtnSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbtnSaveActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 371, Short.MAX_VALUE)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addComponent(jLabel1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jeName, javax.swing.GroupLayout.DEFAULT_SIZE, 330, Short.MAX_VALUE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel2)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jSeparator1, javax.swing.GroupLayout.DEFAULT_SIZE, 303, Short.MAX_VALUE))
                    .addComponent(jbtnSave, javax.swing.GroupLayout.Alignment.TRAILING))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(jeName, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jLabel2)
                    .addComponent(jSeparator1, javax.swing.GroupLayout.PREFERRED_SIZE, 2, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 280, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jbtnSave)
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jbtnSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbtnSaveActionPerformed
        // First, check whether all codons are specified and terminal and start codons are set.
        boolean bStart = false;
        boolean bStop = false;
        for(int i=0;i<64;i++)
        {
            if(((String)lvCodons.getValueAt(i, 3)).isEmpty())
            {
                JOptionPane.showMessageDialog(null, "Assign the amino acids to all codons", "Custom codon table creator", JOptionPane.ERROR_MESSAGE);
                return;
            }
            if((Boolean)lvCodons.getValueAt(i, 0))
                bStart = true;
            if((Boolean)lvCodons.getValueAt(i, 1))
                bStop = true;
        }
        if(!bStart || !bStop)
        {
            JOptionPane.showMessageDialog(null, "Specify at least one start and one terminal codon", "Custom codon table creator", JOptionPane.ERROR_MESSAGE);
            return;
        }
        // Check table name.
        if(jeName.getText().isEmpty())
        {
            JOptionPane.showMessageDialog(null, "Specify the codon table name", "Custom codon table creator", JOptionPane.ERROR_MESSAGE);
            return;
        }  
        // Create the codon table.
        CustomCodonTable.Codon[] codons = new CustomCodonTable.Codon[64];
        for(int i=0;i<64;i++)
        { 
            bStart = (Boolean)lvCodons.getValueAt(i, 0);
            bStop  = (Boolean)lvCodons.getValueAt(i, 1);
            String strSeq = (String)lvCodons.getValueAt(i, 2);
            String strAa  = (String)lvCodons.getValueAt(i, 3);
            String strTLC = (String)lvCodons.getValueAt(i, 4);
            String strOLC = (String)lvCodons.getValueAt(i, 5);
            codons[i] = new CustomCodonTable.Codon(strAa, strTLC, strOLC, strSeq, bStop, bStart);
        }
        ct = CustomCodonTable.create(codons, jeName.getText());
        setVisible(false);
    }//GEN-LAST:event_jbtnSaveActionPerformed

private void lvCodonsMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_lvCodonsMouseClicked
    if(evt.getClickCount()==1)
    {
        Point p = evt.getPoint();
        int row = lvCodons.rowAtPoint(p);
        int col = lvCodons.columnAtPoint(p);
        if(row>-1 && col>-1)
        {
            boolean bStart = (Boolean)lvCodons.getValueAt(row, 0);
            boolean bTerm  = (Boolean)lvCodons.getValueAt(row, 1);
            // If the column clicked is the first one, check, whether it is
            // selected or not. If the new state is selected, check, whether the second
            // column is selected and unselect it if needed.
            if(col==0 && bStart && bTerm)
                lvCodons.setValueAt(new Boolean(false), row, 1);
            // If the column clicked is the second one, check, whether it is
            // selected or not. If the new state is selected, check, whether the second
            // column is selected and unselect it if needed.
            if(col==1 && bTerm && bStart)
                lvCodons.setValueAt(new Boolean(false), row, 0);
        }
    }
}//GEN-LAST:event_lvCodonsMouseClicked

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JPopupMenu jPopupMenu1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JButton jbtnSave;
    private javax.swing.JTextField jeName;
    private javax.swing.JTable lvCodons;
    // End of variables declaration//GEN-END:variables

}
