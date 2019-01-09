/*
    File:
        OptionsForm.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Shows the options dialog.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.writers.nfa2;

import javax.swing.DefaultListModel;
import javax.swing.JOptionPane;


public class OptionsForm extends javax.swing.JDialog
{
    private ExportOptions eo = null;


    public OptionsForm()
    {
        super((java.awt.Frame)null, true);
        initComponents();
        setLocationRelativeTo(null);
    }

    public ExportOptions getOptions(String[] pops)
    {
        DefaultListModel model = new DefaultListModel();
        for(String s:pops)
            model.addElement(s);
        jlNames.setModel(model);
        setVisible(true);
        return eo;
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
        jScrollPane1 = new javax.swing.JScrollPane();
        jlNames = new javax.swing.JList();
        jcbLimit = new javax.swing.JCheckBox();
        jeLimit = new javax.swing.JTextField();
        jButton1 = new javax.swing.JButton();
        jButton2 = new javax.swing.JButton();
        jLabel1 = new javax.swing.JLabel();
        jcbNoIntName = new javax.swing.JCheckBox();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Exporter options");
        setResizable(false);

        jScrollPane1.setViewportView(jlNames);

        jcbLimit.setText("Limit strains count");
        jcbLimit.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jcbLimitStateChanged(evt);
            }
        });

        jeLimit.setEnabled(false);

        jButton1.setText("Cancel");
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });

        jButton2.setText("OK");
        jButton2.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton2ActionPerformed(evt);
            }
        });

        jLabel1.setText("Populations to export");

        jcbNoIntName.setSelected(true);
        jcbNoIntName.setText("Do not include internal name");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(58, 58, 58)
                        .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 237, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jLabel1))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(18, 18, 18)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jcbNoIntName)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jcbLimit)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jeLimit, javax.swing.GroupLayout.PREFERRED_SIZE, 73, javax.swing.GroupLayout.PREFERRED_SIZE)))))
                .addContainerGap(14, Short.MAX_VALUE))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap(181, Short.MAX_VALUE)
                .addComponent(jButton2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jButton1)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel1)
                .addGap(7, 7, 7)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jcbLimit)
                    .addComponent(jeLimit, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jcbNoIntName)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jButton1)
                    .addComponent(jButton2))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        eo = null;
        setVisible(false);
    }//GEN-LAST:event_jButton1ActionPerformed

    private void jButton2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton2ActionPerformed
        eo = null;
        if(jlNames.getSelectedIndices().length==0)
        {
            JOptionPane.showMessageDialog(null,
                                      "No populations specified",
                                      "Native FASTA exporter",
                                      JOptionPane.ERROR_MESSAGE);
            return;
        }
        if(jcbLimit.isSelected() && jeLimit.getText().isEmpty())
        {
            JOptionPane.showMessageDialog(null,
                                      "Specify the maximal number of strains",
                                      "Native FASTA exporter",
                                      JOptionPane.ERROR_MESSAGE);
            return;
        }
        eo = new ExportOptions();
        eo.bNoIntNames = jcbNoIntName.isSelected();
        eo.nLimit = (jcbLimit.isSelected()) ? Integer.parseInt(jeLimit.getText()) : Integer.MAX_VALUE;
        int[] indices = jlNames.getSelectedIndices();
        eo.pops = new String[indices.length];
        for(int i=0;i<indices.length;i++)
        {
            eo.pops[i] = (String)jlNames.getModel().getElementAt(indices[i]);
        }
        setVisible(false);
    }//GEN-LAST:event_jButton2ActionPerformed

    private void jcbLimitStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jcbLimitStateChanged
        jeLimit.setEnabled(jcbLimit.isSelected());
    }//GEN-LAST:event_jcbLimitStateChanged

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup buttonGroup1;
    private javax.swing.JButton jButton1;
    private javax.swing.JButton jButton2;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JCheckBox jcbLimit;
    private javax.swing.JCheckBox jcbNoIntName;
    private javax.swing.JTextField jeLimit;
    private javax.swing.JList jlNames;
    // End of variables declaration//GEN-END:variables

}
