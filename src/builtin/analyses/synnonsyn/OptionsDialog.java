/*
    File:
        OptionsDialog.java
 *   
    Revision:
        1.1.0.1
 * 
    Description:
        Allows the user to select options for performing the analysis.
 * 
    Project:
        GeneAnalyzer 2.2
 * 
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.synnonsyn;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileNameExtensionFilter;


public class OptionsDialog extends javax.swing.JDialog 
{
    private AnalysisOptions ao = null; 
    
    /**
     *  Creates the dialog.
     */
    public OptionsDialog(String[] pops) 
    {      
        super((java.awt.Frame)null, true);
        initComponents();
        setLocationRelativeTo(null);
        for(String s:pops)
        {
            jcPop.addItem(s);
            jcOut.addItem(s);
        }        
    }    
    
    /**
     *  Displays the dialog.
     */
    public AnalysisOptions getOptions()
    {
        setVisible(true);
        return ao;
    }  
    

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPanel1 = new javax.swing.JPanel();
        jLabel2 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        jcPop = new javax.swing.JComboBox();
        jcOut = new javax.swing.JComboBox();
        jButton1 = new javax.swing.JButton();
        jButton2 = new javax.swing.JButton();
        jPanel2 = new javax.swing.JPanel();
        jLabel4 = new javax.swing.JLabel();
        jeFilename = new javax.swing.JTextField();
        jButton3 = new javax.swing.JButton();
        jcbShow = new javax.swing.JCheckBox();
        jPanel3 = new javax.swing.JPanel();
        jcbLimit = new javax.swing.JCheckBox();
        jeLimit = new javax.swing.JTextField();
        jcbFreq = new javax.swing.JCheckBox();
        jeFreq = new javax.swing.JTextField();
        jcbUseTerm = new javax.swing.JCheckBox();
        jcbExclTerm = new javax.swing.JCheckBox();
        jPanel4 = new javax.swing.JPanel();
        jcbJC_pi = new javax.swing.JCheckBox();
        jcbJC_theta = new javax.swing.JCheckBox();
        jcbJC_k = new javax.swing.JCheckBox();
        jLabel1 = new javax.swing.JLabel();
        jPanel6 = new javax.swing.JPanel();
        jLabel6 = new javax.swing.JLabel();
        jrbExcludeAll = new javax.swing.JRadioButton();
        jrbExcludeTajD = new javax.swing.JRadioButton();

        setTitle("Synonymous and nonsynonymous analysis options");
        setResizable(false);

        jPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Populations"));

        jLabel2.setText("Population of interest:");

        jLabel3.setText("Outgroup:");

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jLabel2)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jcPop, 0, 211, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                        .addComponent(jLabel3)
                        .addGap(61, 61, 61)
                        .addComponent(jcOut, 0, 211, Short.MAX_VALUE)))
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(jcPop, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(jcOut, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

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

        jPanel2.setBorder(javax.swing.BorderFactory.createTitledBorder("Output"));

        jLabel4.setText("Output file:");

        jButton3.setText("...");
        jButton3.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton3ActionPerformed(evt);
            }
        });

        jcbShow.setSelected(true);
        jcbShow.setText("Show results in a new window");

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel2Layout.createSequentialGroup()
                        .addGap(10, 10, 10)
                        .addComponent(jcbShow)
                        .addContainerGap())
                    .addGroup(jPanel2Layout.createSequentialGroup()
                        .addComponent(jLabel4)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jeFilename, javax.swing.GroupLayout.DEFAULT_SIZE, 232, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jButton3, javax.swing.GroupLayout.PREFERRED_SIZE, 35, javax.swing.GroupLayout.PREFERRED_SIZE))))
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel4)
                    .addComponent(jeFilename, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jButton3))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 9, Short.MAX_VALUE)
                .addComponent(jcbShow))
        );

        jPanel3.setBorder(javax.swing.BorderFactory.createTitledBorder("Analysis options"));

        jcbLimit.setText("Limit the number of strains to");
        jcbLimit.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jcbLimitItemStateChanged(evt);
            }
        });

        jeLimit.setEnabled(false);
        jeLimit.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent evt) {
                jeLimitKeyTyped(evt);
            }
        });

        jcbFreq.setText("Use custom singletons cut-off frequency:");
        jcbFreq.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jcbFreqItemStateChanged(evt);
            }
        });

        jeFreq.setEnabled(false);
        jeFreq.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent evt) {
                jeFreqKeyTyped(evt);
            }
        });

        jcbUseTerm.setText("Use terminal codons when generating the path");
        jcbUseTerm.setToolTipText("<html>Check this option in order not to consider the terminal codons<br>when calculating the number of synonymous sites within a codon");

        jcbExclTerm.setSelected(true);
        jcbExclTerm.setText("Exclude terminal codons from analysis");

        javax.swing.GroupLayout jPanel3Layout = new javax.swing.GroupLayout(jPanel3);
        jPanel3.setLayout(jPanel3Layout);
        jPanel3Layout.setHorizontalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel3Layout.createSequentialGroup()
                        .addComponent(jcbLimit)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jeLimit, javax.swing.GroupLayout.PREFERRED_SIZE, 69, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(jPanel3Layout.createSequentialGroup()
                        .addComponent(jcbFreq)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jeFreq, javax.swing.GroupLayout.PREFERRED_SIZE, 82, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(jcbExclTerm)
                    .addComponent(jcbUseTerm))
                .addContainerGap(25, Short.MAX_VALUE))
        );
        jPanel3Layout.setVerticalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jcbLimit)
                    .addComponent(jeLimit, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jcbFreq)
                    .addComponent(jeFreq, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jcbExclTerm)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jcbUseTerm)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        jPanel4.setBorder(javax.swing.BorderFactory.createTitledBorder("Jukes-Cantor correction"));

        jcbJC_pi.setSelected(true);
        jcbJC_pi.setText("Pi");

        jcbJC_theta.setText("Theta");

        jcbJC_k.setSelected(true);
        jcbJC_k.setText("K");

        javax.swing.GroupLayout jPanel4Layout = new javax.swing.GroupLayout(jPanel4);
        jPanel4.setLayout(jPanel4Layout);
        jPanel4Layout.setHorizontalGroup(
            jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel4Layout.createSequentialGroup()
                .addComponent(jcbJC_pi)
                .addGap(18, 18, 18)
                .addComponent(jcbJC_theta)
                .addGap(18, 18, 18)
                .addComponent(jcbJC_k)
                .addContainerGap(189, Short.MAX_VALUE))
        );
        jPanel4Layout.setVerticalGroup(
            jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel4Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                .addComponent(jcbJC_pi)
                .addComponent(jcbJC_theta)
                .addComponent(jcbJC_k))
        );

        jLabel1.setText("<html>For synonymous and nonsynonymous sites analysis you need to specify<br>the population of interest and the outgroup as well as whether or not<br>to correct the results using Jukes-Cantor correction formula.");

        jPanel6.setBorder(javax.swing.BorderFactory.createTitledBorder("Missing data"));

        jLabel6.setText("<html>Tajima's D cannot be calculated for the sample size less than 4. How<br>should the blocks with the sample size of <4 be treated?");

        jrbExcludeAll.setText("Exclude from the analysis completely");

        jrbExcludeTajD.setSelected(true);
        jrbExcludeTajD.setText("Exclude only when calculating Tajima's D");

        javax.swing.GroupLayout jPanel6Layout = new javax.swing.GroupLayout(jPanel6);
        jPanel6.setLayout(jPanel6Layout);
        jPanel6Layout.setHorizontalGroup(
            jPanel6Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel6Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel6))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel6Layout.createSequentialGroup()
                .addContainerGap(10, Short.MAX_VALUE)
                .addGroup(jPanel6Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jrbExcludeAll)
                    .addComponent(jrbExcludeTajD))
                .addGap(107, 107, 107))
        );
        jPanel6Layout.setVerticalGroup(
            jPanel6Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel6Layout.createSequentialGroup()
                .addComponent(jLabel6)
                .addGap(18, 18, 18)
                .addComponent(jrbExcludeAll)
                .addGap(7, 7, 7)
                .addComponent(jrbExcludeTajD)
                .addContainerGap(12, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jPanel4, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPanel1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jLabel1)
                            .addComponent(jPanel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jPanel3, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 15, Short.MAX_VALUE))
                            .addComponent(jPanel6, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addComponent(jButton2)
                        .addGap(5, 5, 5)
                        .addComponent(jButton1)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jLabel1)
                        .addGap(7, 7, 7)
                        .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jPanel4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jPanel3, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jPanel6, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jButton2)
                    .addComponent(jButton1))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

private void jButton2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton2ActionPerformed
    ao = null;
    setVisible(false);
}//GEN-LAST:event_jButton2ActionPerformed

private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
    ao = new AnalysisOptions();
    // Output file.
    if(!jeFilename.getText().isEmpty())
        ao.strOutput = jeFilename.getText();
    else if(!jcbShow.isSelected())
    {
        JOptionPane.showMessageDialog(null,
                                      "You must either specify an output file or check the box under\n" +
                                      "the file name field to display the result in a new window",
                                      "Synonymous and nonsynonymous sites analysis",
                                      JOptionPane.ERROR_MESSAGE);
        return;
    }
    ao.bShowRes = jcbShow.isSelected();
    // Cut-off frequency.
    if(jcbFreq.isSelected())
    {
        if(jeFreq.getText().isEmpty())
        {
            JOptionPane.showMessageDialog(null, 
                                      "No cut-off frequency specified", 
                                      "Synonymous and nonsynonymous sites analysis",
                                      JOptionPane.ERROR_MESSAGE);
            return;
        }
        else
            ao.cof = Float.parseFloat(jeFreq.getText());
    }
    else
        ao.cof = 1.0f;
    ao.strPop = (String)jcPop.getSelectedItem();
    ao.strOut = (String)jcOut.getSelectedItem();
    ao.bJC_pi = jcbJC_pi.isSelected();
    ao.bJC_t  = jcbJC_theta.isSelected();
    ao.bJC_K  = jcbJC_k.isSelected();
    ao.bExclTer = jcbExclTerm.isSelected();
    ao.bUseTer = jcbUseTerm.isSelected();
    // Strains count limit.
    if(jcbLimit.isSelected() && !jeLimit.getText().isEmpty())
        ao.maxstr = Integer.parseInt(jeLimit.getText());
    else
        ao.maxstr = Integer.MAX_VALUE;
    ao.bExclAll = jrbExcludeAll.isSelected();
    setVisible(false);
}//GEN-LAST:event_jButton1ActionPerformed

private void jButton3ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton3ActionPerformed
    JFileChooser fc = new JFileChooser();
    fc.addChoosableFileFilter(
                    new FileNameExtensionFilter("Tab-separated text", "txt"));
    if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
    {
        try
        {
            jeFilename.setText(fc.getSelectedFile().getCanonicalPath());
        }
        catch (Exception e)
        {
            jeFilename.setText("");
        }
    }
}//GEN-LAST:event_jButton3ActionPerformed

private void jcbLimitItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jcbLimitItemStateChanged
    jeLimit.setEnabled(jcbLimit.isSelected());
}//GEN-LAST:event_jcbLimitItemStateChanged

private void jeLimitKeyTyped(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_jeLimitKeyTyped
    char c = evt.getKeyChar();
    if(!Character.isDigit(c))
        evt.setKeyChar('\0');
}//GEN-LAST:event_jeLimitKeyTyped

private void jcbFreqItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jcbFreqItemStateChanged
    jeFreq.setEnabled(jcbFreq.isSelected());
}//GEN-LAST:event_jcbFreqItemStateChanged

private void jeFreqKeyTyped(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_jeFreqKeyTyped
    char c = evt.getKeyChar();
    // If the key is '.' check, whether the text
    // already contains this character.
    if(c=='.')
    {
        if(jeFreq.getText().contains("."))
            evt.setKeyChar('\0');
        return;
    }
    if(!Character.isDigit(c))
        evt.setKeyChar('\0');
}//GEN-LAST:event_jeFreqKeyTyped

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton jButton1;
    private javax.swing.JButton jButton2;
    private javax.swing.JButton jButton3;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JPanel jPanel4;
    private javax.swing.JPanel jPanel6;
    private javax.swing.JComboBox jcOut;
    private javax.swing.JComboBox jcPop;
    private javax.swing.JCheckBox jcbExclTerm;
    private javax.swing.JCheckBox jcbFreq;
    private javax.swing.JCheckBox jcbJC_k;
    private javax.swing.JCheckBox jcbJC_pi;
    private javax.swing.JCheckBox jcbJC_theta;
    private javax.swing.JCheckBox jcbLimit;
    private javax.swing.JCheckBox jcbShow;
    private javax.swing.JCheckBox jcbUseTerm;
    private javax.swing.JTextField jeFilename;
    private javax.swing.JTextField jeFreq;
    private javax.swing.JTextField jeLimit;
    private javax.swing.JRadioButton jrbExcludeAll;
    private javax.swing.JRadioButton jrbExcludeTajD;
    // End of variables declaration//GEN-END:variables

}
