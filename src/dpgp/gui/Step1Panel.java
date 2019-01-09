/*
    File:
        Step1Panel.java
 *
    Revision:
        1.0.0.0
 *
    Description:
        Displays the options of the first import step.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package dpgp.gui;

import dpgp.DPGPImporter;
import java.io.File;
import java.util.Vector;
import javax.swing.JFileChooser;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.DefaultTableModel;
import kernel.SettingsManager;

public class Step1Panel extends javax.swing.JPanel
{
    private class TableListener implements TableModelListener
    {
        public void tableChanged(TableModelEvent e)
        {
            if(e.getType()==TableModelEvent.UPDATE)
            {
                int n = e.getFirstRow();
                DefaultTableModel model = (DefaultTableModel)lvEntries.getModel();
                updateEntry(n, (String)model.getValueAt(n, 0), (String)model.getValueAt(n, 1));
            }
        }
    };

    private class VMAFileFilter extends FileFilter
    {
        @Override
        public boolean accept(File f)
        {
            // Accept directories.
            if (f.isDirectory())
                return true;
            // Files.
            String fileName = f.getName();
            int i = fileName.lastIndexOf('.');
            if (i>0 && i<fileName.length()-1)
                return fileName.substring(i+1).equalsIgnoreCase("vma");
            return false;
        }

        @Override
        public String getDescription()
        {
            return "Vertical Multiple Alignment file (*.vma)";
        }

    };

    public class Entry
    {
        public String strSpecies = "Drosophila";
        public String strStrain  = "Strain";
        public File   file       = null;
    };

    private SettingsManager sm      = null;
    private Vector<Entry> entries   = null;


    public Step1Panel(SettingsManager sm)
    {
        initComponents();
        this.sm = sm;
        entries = new Vector<Entry>();
        setVisible(false);
        lvEntries.getModel().addTableModelListener(new TableListener());
    }

    public Entry[] getEntries()
    {
        if(entries.size()==0)
            return null;
        else
            return entries.toArray(new Entry[1]);
    }

    public String getChromosomeName()
    {
        return jeChr.getText();
    }

    public int getCutoff()
    {
        if(jeThreshold.getText().isEmpty())
            return 0;
        else
            return Integer.parseInt(jeThreshold.getText());
    }

    private void updateList()
    {
        DefaultTableModel model = (DefaultTableModel)lvEntries.getModel();
        // Remove old lines.
        for(int i=model.getRowCount()-1;i>=0;i--)
            model.removeRow(i);
        // Add lines.
        for(Entry e:entries)
            model.addRow(new Object[]{e.strSpecies, e.strStrain, e.file.getAbsolutePath()});
    }

    private void updateEntry(int i, String strSpecies, String strStrain)
    {
        entries.get(i).strSpecies = strSpecies;
        entries.get(i).strStrain = strStrain;
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
        jLabel3 = new javax.swing.JLabel();
        jeChr = new javax.swing.JTextField();
        jSeparator1 = new javax.swing.JSeparator();
        jPanel1 = new javax.swing.JPanel();
        jScrollPane1 = new javax.swing.JScrollPane();
        lvEntries = new javax.swing.JTable();
        jbUp = new javax.swing.JButton();
        jbDown = new javax.swing.JButton();
        jbDelete = new javax.swing.JButton();
        jbAdd = new javax.swing.JButton();
        jLabel2 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        jeThreshold = new javax.swing.JTextField();

        setPreferredSize(new java.awt.Dimension(681, 416));

        jLabel1.setFont(new java.awt.Font("Georgia", 0, 18));
        jLabel1.setText("STEP 1. Select the chromosome file(s) to import.");

        jLabel3.setText("Chromosome:");

        jPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Strains"));

        lvEntries.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {

            },
            new String [] {
                "Species", "Strain", "Filename"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.String.class, java.lang.String.class, java.lang.String.class
            };
            boolean[] canEdit = new boolean [] {
                true, true, false
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
        });
        lvEntries.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                lvEntriesMouseClicked(evt);
            }
        });
        jScrollPane1.setViewportView(lvEntries);

        jbUp.setText("Up");
        jbUp.setEnabled(false);
        jbUp.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbUpActionPerformed(evt);
            }
        });

        jbDown.setText("Down");
        jbDown.setEnabled(false);
        jbDown.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbDownActionPerformed(evt);
            }
        });

        jbDelete.setText("Delete");
        jbDelete.setEnabled(false);
        jbDelete.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbDeleteActionPerformed(evt);
            }
        });

        jbAdd.setText("Add");
        jbAdd.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbAddActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 562, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jbAdd, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jbDown, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jbDelete, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jbUp, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jbUp)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jbDown)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jbDelete)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jbAdd))
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 256, Short.MAX_VALUE))
                .addContainerGap())
        );

        jLabel2.setText("<html>Specify the chromosome name and then add the files with the sequence. Only *.VMA files can be imported. Also, make sure all<br>sequences are from the same chromosome, otherwise the alignment cannot be done.");

        jLabel4.setText("Threshold (bases with the lower value will be replaced by 'N'):");

        jeThreshold.setText("15");
        jeThreshold.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent evt) {
                jeThresholdKeyTyped(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jPanel1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jLabel1)
                    .addComponent(jSeparator1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 661, Short.MAX_VALUE)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addGroup(javax.swing.GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                            .addComponent(jLabel3)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                            .addComponent(jeChr, javax.swing.GroupLayout.PREFERRED_SIZE, 141, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGap(18, 18, 18)
                            .addComponent(jLabel4)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                            .addComponent(jeThreshold))
                        .addComponent(jLabel2, javax.swing.GroupLayout.Alignment.LEADING)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jSeparator1, javax.swing.GroupLayout.PREFERRED_SIZE, 10, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(jeChr, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel4)
                    .addComponent(jeThreshold, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void jbAddActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbAddActionPerformed
        JFileChooser fc = new JFileChooser();
        String strLastDir = sm.getSetting("", "LastDPGPDir");
        if(strLastDir!=null)
            fc.setCurrentDirectory(new File(strLastDir));
        fc.setMultiSelectionEnabled(true);
        fc.setAcceptAllFileFilterUsed(false);
        fc.addChoosableFileFilter(new VMAFileFilter());
        // Show OpenDialog.
        File[] files = null;
        if(fc.showOpenDialog(null)==JFileChooser.APPROVE_OPTION)
            files = fc.getSelectedFiles();
        // If no files selected, return.
        if(files==null)
            return;
        // Load files.
        for(File f:files)
        {
            if(f.exists())
            {
                Entry en = new Entry();
                en.file = f;
                entries.add(en);
            }
        }
        updateList();
        // Save the directory.
        strLastDir = files[0].getParent();
        sm.addSetting("", "LastDPGPDir", strLastDir);
    }//GEN-LAST:event_jbAddActionPerformed

    private void lvEntriesMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_lvEntriesMouseClicked
        if(evt==null || evt.getClickCount()==1)
        {
            int iIndex = lvEntries.getSelectedRow();
            if(iIndex>-1)
            {
                // First enable all buttons and then disable unnecessary ones.
                jbDelete.setEnabled(true);
                jbUp.setEnabled(true);
                jbDown.setEnabled(true);
                if(iIndex==0)
                    jbUp.setEnabled(false);
                if(iIndex==lvEntries.getRowCount()-1)
                    jbDown.setEnabled(false);
            }
            else    // Disable buttons.
            {
                jbDelete.setEnabled(false);
                jbUp.setEnabled(false);
                jbDown.setEnabled(false);
            }
        }
    }//GEN-LAST:event_lvEntriesMouseClicked

    private void jbDeleteActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbDeleteActionPerformed
        int iIndex = lvEntries.getSelectedRow();
        if(iIndex>-1)
        {
            entries.remove(iIndex);
            updateList();
        }
        lvEntriesMouseClicked(null);
    }//GEN-LAST:event_jbDeleteActionPerformed

    private void jbUpActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbUpActionPerformed
        int iIndex = lvEntries.getSelectedRow();
        if(iIndex>-1)
        {
            Entry tmp = entries.get(iIndex);
            entries.remove(iIndex);
            entries.insertElementAt(tmp, iIndex-1);
        }
        updateList();
        lvEntriesMouseClicked(null);
    }//GEN-LAST:event_jbUpActionPerformed

    private void jbDownActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbDownActionPerformed
        int iIndex = lvEntries.getSelectedRow();
        if(iIndex>-1)
        {
            Entry tmp = entries.get(iIndex);
            entries.remove(iIndex);
            entries.insertElementAt(tmp, iIndex+1);
        }
        updateList();
        lvEntriesMouseClicked(null);
    }//GEN-LAST:event_jbDownActionPerformed

    private void jeThresholdKeyTyped(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_jeThresholdKeyTyped
        char c = evt.getKeyChar();
        if(!Character.isDigit(c))
            evt.setKeyChar('\0');
    }//GEN-LAST:event_jeThresholdKeyTyped


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JButton jbAdd;
    private javax.swing.JButton jbDelete;
    private javax.swing.JButton jbDown;
    private javax.swing.JButton jbUp;
    private javax.swing.JTextField jeChr;
    private javax.swing.JTextField jeThreshold;
    private javax.swing.JTable lvEntries;
    // End of variables declaration//GEN-END:variables

}
