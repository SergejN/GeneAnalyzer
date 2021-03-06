/*
    File:
        AddGeneRegionDialog.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Shows the dialog allowing the user to add a gene region to the strains.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package gui.gene;

import bio.gene.GeneRegion;
import javax.swing.JOptionPane;
import kernel.Main;


public class AddGeneRegionDialog extends javax.swing.JDialog
{
    private GeneRegion reg = null;


    public AddGeneRegionDialog()
    {
        super((java.awt.Frame)null, true);
        initComponents();
        setLocationRelativeTo(null);
    }

    /**
     *  Creates a GeneRegion instamce and returns it. If the user clicks CANCEL,
     *  the method returns null.
     *
     *  @return
     */
    public GeneRegion createGeneRegion()
    {
        jcbTypes.addItem(GeneRegion.EXON);
        jcbTypes.addItem(GeneRegion.INTRON);
        jcbTypes.addItem(GeneRegion.UTR5);
        jcbTypes.addItem(GeneRegion.UTR3);
        jcbTypes.addItem(GeneRegion.INTERGENIC);
        setVisible(true);
        return reg;
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jButton1 = new javax.swing.JButton();
        jButton2 = new javax.swing.JButton();
        jLabel1 = new javax.swing.JLabel();
        jcbTypes = new javax.swing.JComboBox();
        jLabel2 = new javax.swing.JLabel();
        jeStart = new javax.swing.JTextField();
        jLabel3 = new javax.swing.JLabel();
        jeEnd = new javax.swing.JTextField();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Add gene region");
        setModal(true);
        setResizable(false);

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

        jLabel1.setText("Type:");

        jcbTypes.setEditable(true);

        jLabel2.setText("Start:");

        jeStart.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent evt) {
                jeStartKeyTyped(evt);
            }
        });

        jLabel3.setText("End:");

        jeEnd.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent evt) {
                jeEndKeyTyped(evt);
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
                        .addComponent(jLabel1)
                        .addGap(18, 18, 18)
                        .addComponent(jcbTypes, 0, 210, Short.MAX_VALUE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel2)
                        .addGap(18, 18, 18)
                        .addComponent(jeStart, javax.swing.GroupLayout.DEFAULT_SIZE, 82, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jLabel3)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jeEnd, javax.swing.GroupLayout.PREFERRED_SIZE, 86, javax.swing.GroupLayout.PREFERRED_SIZE))
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
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(jcbTypes, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(jeStart, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jeEnd, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel3))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jButton1)
                    .addComponent(jButton2))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        reg = null;
        setVisible(false);
    }//GEN-LAST:event_jButton1ActionPerformed

    private void jButton2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton2ActionPerformed
        String type = (String)jcbTypes.getSelectedItem();
        if(type.isEmpty())
        {
            JOptionPane.showMessageDialog(this, "Specify the region type.", Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            return;
        }
        if(jeStart.getText().isEmpty())
        {
            JOptionPane.showMessageDialog(this, "Specify the region start.", Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            return;
        }
        if(jeEnd.getText().isEmpty())
        {
            JOptionPane.showMessageDialog(this, "Specify the region end.", Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            return;
        }
        int iStart = Integer.parseInt(jeStart.getText());
        int iEnd   = Integer.parseInt(jeEnd.getText());
        if(iStart==0)
        {
            JOptionPane.showMessageDialog(this, "Region start must be greater than 0.", Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            return;
        }
        if(iEnd<iStart)
        {
            JOptionPane.showMessageDialog(this, "Region end must be greater than or equal to the start.", Main.APPTITLE, JOptionPane.ERROR_MESSAGE);
            return;
        }
        reg = new GeneRegion(type);
        reg.setStart(iStart);
        reg.setEnd(iEnd);
        setVisible(false);
    }//GEN-LAST:event_jButton2ActionPerformed

    private void jeStartKeyTyped(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_jeStartKeyTyped
        if(!Character.isDigit(evt.getKeyChar()))
            evt.setKeyChar('\0');
    }//GEN-LAST:event_jeStartKeyTyped

    private void jeEndKeyTyped(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_jeEndKeyTyped
        if(!Character.isDigit(evt.getKeyChar()))
            evt.setKeyChar('\0');
    }//GEN-LAST:event_jeEndKeyTyped

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton jButton1;
    private javax.swing.JButton jButton2;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JComboBox jcbTypes;
    private javax.swing.JTextField jeEnd;
    private javax.swing.JTextField jeStart;
    // End of variables declaration//GEN-END:variables

}
