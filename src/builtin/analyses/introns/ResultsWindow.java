/*
    File:
        ResultsWindow.java
 *
    Revision:
        1.0.0.1
 *
    Description:
        Displays the analysis results.
 *
    Project:
        GeneAnalyzer 2.2
 *
    Copyright:
        (c) 2008. Sergej Nowoshilow, Biozentrum, Martinsried, Germany.
 */

package builtin.analyses.introns;

import algorithms.SequenceRoutines;
import java.awt.Toolkit;
import java.awt.datatransfer.StringSelection;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import javax.swing.JFrame;
import javax.swing.table.DefaultTableModel;


public class ResultsWindow extends javax.swing.JFrame
{
    public ResultsWindow()
    {
        initComponents();
        setExtendedState(JFrame.MAXIMIZED_BOTH);
        lvResults.getTableHeader().addMouseListener(new MouseAdapter()
            {
                @Override
                public void mouseClicked(MouseEvent e)
                {
                    int iCol = lvResults.getColumnModel().getColumnIndexAtX(e.getX());
                    lvResults.setColumnSelectionInterval(iCol, iCol);
                    lvResults.setRowSelectionInterval(0, lvResults.getRowCount()-1);
                }
            });
    }

    /**
     *  Sets the results, prepares the table etc. Since depending on the dataset
     *  these operations can take a while, do not display the window yet.
     *
     *  @param strContent
     */
    public void setResults(String strContent)
    {
        DefaultTableModel model = (DefaultTableModel)lvResults.getModel();
        // Remove old lines.
        for(int i=model.getRowCount()-1;i>=0;i--)
            model.removeRow(i);
        String[] lines = strContent.split("\n");
        // Don't display the lines 1-15, since they only contain parameters information.
        for(int i=15;i<lines.length;i++)
            model.addRow(lines[i].split("\t"));
    }

    /**
     *  Display the window, when everything is ready.
     */
    public void displayResults()
    {
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

        jPopupMenu1 = new javax.swing.JPopupMenu();
        miCopy = new javax.swing.JMenuItem();
        miCopyAsIs = new javax.swing.JMenuItem();
        jScrollPane1 = new javax.swing.JScrollPane();
        lvResults = new javax.swing.JTable();

        miCopy.setText("Copy selection");
        miCopy.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miCopyActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miCopy);

        miCopyAsIs.setText("Copy selection and keep the spaces");
        miCopyAsIs.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                miCopyAsIsActionPerformed(evt);
            }
        });
        jPopupMenu1.add(miCopyAsIs);

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Short introns analysis");

        lvResults.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {

            },
            new String [] {
                "Gene", "", "Sample size", "Sites", "pi", "theta", "Tajima's D", "Tajima's D'", "P", "Singletons", "Transitions", "Transversions", "", "Sites", "pi", "theta", "Tajima's D", "Tajima's D'", "P", "Singletons", "Transitions", "Transversions", "", "K", "P", "TSp", "TVp", "D", "TSd", "TVd"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.String.class
            };
            boolean[] canEdit = new boolean [] {
                false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
        });
        lvResults.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        lvResults.setCellSelectionEnabled(true);
        lvResults.setComponentPopupMenu(jPopupMenu1);
        jScrollPane1.setViewportView(lvResults);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 604, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane1, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 435, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void miCopyAsIsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miCopyAsIsActionPerformed
        StringBuffer sb = new StringBuffer();
        int[] rows = lvResults.getSelectedRows();
        int[] cols = lvResults.getSelectedColumns();
        for(int r:rows)
        {
            int ilc = -1;
            for(int c:cols)
            {
                if(ilc>-1)
                    sb.append(SequenceRoutines.generateRepetetiveSequence("\t", c-ilc));
                sb.append(lvResults.getValueAt(r, c));
                ilc = c;
            }
            sb.append("\n");
        }
        StringSelection strsel = new StringSelection(sb.toString());
        Toolkit.getDefaultToolkit().getSystemClipboard().setContents(strsel, strsel);
}//GEN-LAST:event_miCopyAsIsActionPerformed

    private void miCopyActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_miCopyActionPerformed
        StringBuffer sb = new StringBuffer();
        int[] rows = lvResults.getSelectedRows();
        int[] cols = lvResults.getSelectedColumns();
        int iLastIndex = cols[cols.length-1];
        for(int r:rows)
        {
            for(int c:cols)
            {
                sb.append(lvResults.getValueAt(r, c));
                if(c<iLastIndex)
                    sb.append("\t");
            }
            sb.append("\n");
        }
        StringSelection strsel = new StringSelection(sb.toString());
        Toolkit.getDefaultToolkit().getSystemClipboard().setContents(strsel, strsel);
    }//GEN-LAST:event_miCopyActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPopupMenu jPopupMenu1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JTable lvResults;
    private javax.swing.JMenuItem miCopy;
    private javax.swing.JMenuItem miCopyAsIs;
    // End of variables declaration//GEN-END:variables

}
